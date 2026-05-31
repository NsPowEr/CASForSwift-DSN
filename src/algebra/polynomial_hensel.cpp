// polynomial_hensel.cpp — Hensel lifting (quadratic and multi-factor).
// A2 (F2 Block A): quadratic lifting replaces linear lifting.
//
// Algorithm: Geddes-Czapor-Labahn "Algorithms for Computer Algebra" §6.5
//            Algorithm 6.2 (quadratic Hensel lift).
//            von zur Gathen & Gerhard "Modern Computer Algebra" §15.4.
//
// Quadratic lifting: given f ≡ G·H (mod m) and s·G + t·H ≡ 1 (mod m),
// one step maps m → m² while solving the system:
//   e = (f - G·H) / m   (exact integer division, over Z)
//   solve G·ΔH + H·ΔG = e  mod p  using Bézout cofactors s, t
//   G' = G + ΔG·m,  H' = H + ΔH·m
//   then update s', t' so that s'·G' + t'·H' ≡ 1 (mod m²)
// After one step f ≡ G'·H' (mod m²), so the modulus squares each iteration.
//
// Key invariant maintained throughout: all polynomial coefficients are stored
// as symmetric representatives in (−m, m] so that exact division by m is valid.
//
// Correctness certificate: after lifting to modulus M = p^k,
//   G·H ≡ f (mod M).

#include "polynomial_internal.hpp"
#include "cas/numtheory.hpp"
#include <vector>
#include <algorithm>
#include <tuple>

namespace cas::algebra {

// -------------------------------------------------------------------------
// Fp[x] helpers (internal to this file)
// -------------------------------------------------------------------------

static BigInt fp_mod(const BigInt& a, const BigInt& p) {
    BigInt r = a % p;
    if (r.is_negative()) r += p;
    return r;
}

// Balanced (symmetric) representative of a mod m: result in (-m/2, m/2].
// This ensures exact integer division by m after subtraction.
static BigInt balanced_mod(const BigInt& a, const BigInt& m) {
    BigInt r = a % m;
    if (r.is_negative()) r += m;
    // If r > m/2, shift to negative side
    BigInt half = m / BigInt(2);
    if (r > half) r -= m;
    return r;
}

static void poly_balanced_mod(IntPoly& p, const BigInt& m) {
    for (auto& c : p.coefficients()) c = balanced_mod(c, m);
    p.normalize([](const BigInt& v) { return v.is_zero(); });
}

static IntPoly poly_sub_fp(const IntPoly& a, const IntPoly& b, const BigInt& p) {
    std::size_t n = std::max(a.size(), b.size());
    IntPoly r;
    r.resize(n, BigInt(0));
    for (std::size_t i = 0; i < n; ++i) {
        BigInt diff(0);
        if (i < a.size()) diff += a[i];
        if (i < b.size()) diff -= b[i];
        r[i] = fp_mod(diff, p);
    }
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

static IntPoly poly_mul_fp(const IntPoly& a, const IntPoly& b, const BigInt& p) {
    if (a.is_zero() || b.is_zero()) return IntPoly{};
    IntPoly r;
    r.resize(a.size() + b.size() - 1U, BigInt(0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].is_zero()) continue;
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] = fp_mod(r[i + j] + a[i] * b[j], p);
        }
    }
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

// Division with remainder in Fp[x] (p prime, b monic or invertible leading coeff)
static std::pair<IntPoly, IntPoly> poly_divrem_fp(
    const IntPoly& a_in, const IntPoly& b, const BigInt& p) {
    if (b.is_zero()) return {IntPoly{}, a_in};
    IntPoly r = a_in;
    for (auto& c : r.coefficients()) c = fp_mod(c, p);
    r.normalize([](const BigInt& v) { return v.is_zero(); });

    auto inv_lc_res = numtheory::modular_inverse(fp_mod(b.leading_coeff(), p), p);
    if (inv_lc_res.is_error()) return {IntPoly{}, r};
    const BigInt inv_lc = inv_lc_res.value();

    IntPoly q;
    while (!r.is_zero() && r.degree() >= b.degree()) {
        std::size_t dd = r.degree() - b.degree();
        BigInt factor = fp_mod(r.leading_coeff() * inv_lc, p);
        if (q.size() <= dd) q.resize(dd + 1U, BigInt(0));
        q[dd] = factor;
        for (std::size_t i = 0; i < b.size(); ++i) {
            r[i + dd] = fp_mod(r[i + dd] - factor * b[i], p);
        }
        r.normalize([](const BigInt& v) { return v.is_zero(); });
    }
    q.normalize([](const BigInt& v) { return v.is_zero(); });
    return {q, r};
}

// Division with remainder in Z/mZ[x], where m is a prime power.
// Requires lc(b) to be invertible mod m (true for any unit mod m, including monic).
// Coefficients of inputs reduced to [0, m) throughout.
static std::pair<IntPoly, IntPoly> poly_divrem_mod(
    const IntPoly& a_in, const IntPoly& b, const BigInt& m) {
    if (b.is_zero()) return {IntPoly{}, a_in};
    IntPoly r = a_in;
    for (auto& c : r.coefficients()) c = fp_mod(c, m);
    r.normalize([](const BigInt& v) { return v.is_zero(); });

    // Need inv of lc(b) mod m. Since m = p^k and b arises from a Hensel lift
    // where lc(b) is invertible mod p (coprime to p), it is also invertible mod m = p^k.
    // Use extended Euclidean to find inv_lc.
    // For monic b: lc=1, inv_lc=1 (fast path).
    BigInt lc_b = fp_mod(b.is_zero() ? BigInt(0) : b.leading_coeff(), m);
    BigInt inv_lc;
    if (lc_b == BigInt(1)) {
        inv_lc = BigInt(1);
    } else {
        auto inv_res = numtheory::modular_inverse(lc_b, m);
        if (inv_res.is_error()) return {IntPoly{}, r};  // lc not invertible
        inv_lc = inv_res.value();
    }

    IntPoly q;
    while (!r.is_zero() && r.degree() >= b.degree()) {
        std::size_t dd = r.degree() - b.degree();
        BigInt factor = fp_mod(r.leading_coeff() * inv_lc, m);
        if (q.size() <= dd) q.resize(dd + 1U, BigInt(0));
        q[dd] = factor;
        for (std::size_t i = 0; i < b.size(); ++i) {
            r[i + dd] = fp_mod(r[i + dd] - factor * b[i], m);
        }
        r.normalize([](const BigInt& v) { return v.is_zero(); });
    }
    q.normalize([](const BigInt& v) { return v.is_zero(); });
    return {q, r};
}

// Polynomial multiply mod m
static IntPoly poly_mul_mod(const IntPoly& a, const IntPoly& b, const BigInt& m) {
    if (a.is_zero() || b.is_zero()) return IntPoly{};
    IntPoly r;
    r.resize(a.size() + b.size() - 1U, BigInt(0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].is_zero()) continue;
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] = fp_mod(r[i + j] + a[i] * b[j], m);
        }
    }
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

// Polynomial add mod m
static IntPoly poly_add_mod(const IntPoly& a, const IntPoly& b, const BigInt& m) {
    std::size_t n = std::max(a.size(), b.size());
    IntPoly r;
    r.resize(n, BigInt(0));
    for (std::size_t i = 0; i < n; ++i) {
        BigInt sum(0);
        if (i < a.size()) sum += a[i];
        if (i < b.size()) sum += b[i];
        r[i] = fp_mod(sum, m);
    }
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

// Extended GCD in Fp[x]: returns (g, s, t) with s*a + t*b = g, g monic gcd
static std::tuple<IntPoly, IntPoly, IntPoly> poly_xgcd_fp(
    const IntPoly& a_in, const IntPoly& b_in, const BigInt& p) {

    IntPoly r0 = a_in, r1 = b_in;
    for (auto& c : r0.coefficients()) c = fp_mod(c, p);
    for (auto& c : r1.coefficients()) c = fp_mod(c, p);
    r0.normalize([](const BigInt& v) { return v.is_zero(); });
    r1.normalize([](const BigInt& v) { return v.is_zero(); });

    IntPoly s0{std::vector<BigInt>{BigInt(1)}}, s1{std::vector<BigInt>{BigInt(0)}};
    IntPoly t0{std::vector<BigInt>{BigInt(0)}}, t1{std::vector<BigInt>{BigInt(1)}};

    while (!r1.is_zero()) {
        auto [q, rem] = poly_divrem_fp(r0, r1, p);
        r0 = r1; r1 = rem;
        auto tmp_s = poly_sub_fp(s0, poly_mul_fp(q, s1, p), p);
        s0 = s1; s1 = tmp_s;
        auto tmp_t = poly_sub_fp(t0, poly_mul_fp(q, t1, p), p);
        t0 = t1; t1 = tmp_t;
    }

    // Normalize to monic gcd
    if (!r0.is_zero()) {
        auto inv_res = numtheory::modular_inverse(fp_mod(r0.leading_coeff(), p), p);
        if (inv_res.is_ok()) {
            const BigInt inv = inv_res.value();
            for (auto& c : r0.coefficients()) c = fp_mod(c * inv, p);
            for (auto& c : s0.coefficients()) c = fp_mod(c * inv, p);
            for (auto& c : t0.coefficients()) c = fp_mod(c * inv, p);
        }
    }
    return {r0, s0, t0};
}

// Exact integer polynomial multiplication (no modular reduction)
static IntPoly poly_mul_exact(const IntPoly& a, const IntPoly& b) {
    if (a.is_zero() || b.is_zero()) return IntPoly{};
    IntPoly r;
    r.resize(a.size() + b.size() - 1U, BigInt(0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].is_zero()) continue;
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] += a[i] * b[j];
        }
    }
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

// Compute error vector e = (f - G*H) / m as exact integers, then reduce mod p.
// Pre-condition: f ≡ G*H (mod m), so each coefficient of (f - G*H) is
// divisible by m in exact arithmetic.
// G and H must use balanced representatives (values in (-m, m]).
// One quadratic Hensel step (GCL §6.5 Algorithm 6.2).
// Lifts: f ≡ G·H (mod m) to f ≡ G'·H' (mod m²).
// Also lifts Bézout cofactors: s·G + t·H ≡ 1 (mod m) to s'·G' + t'·H' ≡ 1 (mod m²).
//
// Crucially, the Bézout solve is performed modulo m (the current modulus), NOT
// modulo the original prime p.  Only at step 0 (m=p) are they the same.
// Using mod-p arithmetic for higher steps gives only linear convergence.
//
// All polynomials stored with balanced representatives (coefficients in (-m, m]).
static void quadratic_step(
    const IntPoly& f,
    IntPoly& G, IntPoly& H,
    IntPoly& s, IntPoly& t,
    const BigInt& m, const BigInt& m2) {

    // --- Step 1: compute error e_exact = (f - G*H)/m as exact integers, then mod m ---
    // G*H using balanced representatives
    IntPoly GH = poly_mul_exact(G, H);
    std::size_t esz = std::max(f.size(), GH.size());
    IntPoly e_mod;      // e mod m, in [0,m)
    e_mod.resize(esz, BigInt(0));
    for (std::size_t i = 0; i < esz; ++i) {
        BigInt fi = (i < f.size())  ? f[i]  : BigInt(0);
        BigInt gi = (i < GH.size()) ? GH[i] : BigInt(0);
        BigInt diff = fi - gi;
        // diff must be divisible by m
        e_mod[i] = fp_mod(diff / m, m);
    }
    e_mod.normalize([](const BigInt& v) { return v.is_zero(); });

    if (!e_mod.is_zero()) {
        // --- Step 2: solve G_m·ΔG + H_m·ΔH ≡ e_mod (mod m) ---
        // From G'·H' = G·H + (G·ΔH + H·ΔG)·m + O(m²), the lift condition is:
        //   G_m·ΔH + H_m·ΔG ≡ e_mod  (mod m)
        // Apply GCL §6.3 Lemma 6.1 with equation g·σ + h·τ = c (mod m),
        //   where g = G_m, h = H_m, σ = ΔH, τ = ΔG, deg(τ) < deg(g):
        //   τ = ΔG = t_m·e mod G_m   (poly_divrem_mod, monic G_m)
        //   σ = ΔH = s_m·e + q·H_m   (mod m)
        // G and H are balanced mod m — get positive [0,m) form for division
        IntPoly G_m = G;
        for (auto& c : G_m.coefficients()) c = fp_mod(c, m);
        G_m.normalize([](const BigInt& v) { return v.is_zero(); });
        IntPoly H_m = H;
        for (auto& c : H_m.coefficients()) c = fp_mod(c, m);
        H_m.normalize([](const BigInt& v) { return v.is_zero(); });
        IntPoly s_m = s;
        for (auto& c : s_m.coefficients()) c = fp_mod(c, m);
        s_m.normalize([](const BigInt& v) { return v.is_zero(); });
        IntPoly t_m = t;
        for (auto& c : t_m.coefficients()) c = fp_mod(c, m);
        t_m.normalize([](const BigInt& v) { return v.is_zero(); });

        auto [q_g, dg] = poly_divrem_mod(poly_mul_mod(t_m, e_mod, m), G_m, m);
        IntPoly dh = poly_add_mod(poly_mul_mod(s_m, e_mod, m), poly_mul_mod(q_g, H_m, m), m);

        // --- Step 3: update G and H ---
        // G' = G + ΔG·m,  H' = H + ΔH·m
        if (G.size() < dg.size()) G.resize(dg.size(), BigInt(0));
        for (std::size_t i = 0; i < dg.size(); ++i) G[i] += dg[i] * m;
        G.normalize([](const BigInt& v) { return v.is_zero(); });

        if (H.size() < dh.size()) H.resize(dh.size(), BigInt(0));
        for (std::size_t i = 0; i < dh.size(); ++i) H[i] += dh[i] * m;
        H.normalize([](const BigInt& v) { return v.is_zero(); });
    }

    // Reduce G and H to balanced representatives mod m²
    poly_balanced_mod(G, m2);
    poly_balanced_mod(H, m2);

    // --- Step 4: update Bézout cofactors s, t to mod m² ---
    // We need s'·G' + t'·H' ≡ 1 (mod m²).
    // Error: b_err = 1 - s·G' - t·H' (exact integers, divisible by m).
    // δ = b_err / m (exact) mod m.
    // Solve G'_m·Δs + H'_m·Δt ≡ δ (mod m), same pattern:
    //   Δt = t_m·δ mod G'_m,   Δs = s_m·δ + q·H'_m
    {
        IntPoly sG = poly_mul_exact(s, G);
        IntPoly tH_prod = poly_mul_exact(t, H);

        std::size_t bsz = std::max({std::size_t(1), sG.size(), tH_prod.size()});
        IntPoly b_err;
        b_err.resize(bsz, BigInt(0));
        b_err[0] += BigInt(1);
        for (std::size_t i = 0; i < sG.size(); ++i) b_err[i] -= sG[i];
        for (std::size_t i = 0; i < tH_prod.size(); ++i) b_err[i] -= tH_prod[i];

        IntPoly delta;
        delta.resize(b_err.size(), BigInt(0));
        for (std::size_t i = 0; i < b_err.size(); ++i) {
            delta[i] = fp_mod(b_err[i] / m, m);
        }
        delta.normalize([](const BigInt& v) { return v.is_zero(); });

        if (!delta.is_zero()) {
            // G', H' are already balanced mod m² — get [0,m) form for mod-m solve
            IntPoly G_m2 = G;
            for (auto& c : G_m2.coefficients()) c = fp_mod(c, m);
            G_m2.normalize([](const BigInt& v) { return v.is_zero(); });
            IntPoly H_m2 = H;
            for (auto& c : H_m2.coefficients()) c = fp_mod(c, m);
            H_m2.normalize([](const BigInt& v) { return v.is_zero(); });
            IntPoly s_m2 = s;
            for (auto& c : s_m2.coefficients()) c = fp_mod(c, m);
            s_m2.normalize([](const BigInt& v) { return v.is_zero(); });
            IntPoly t_m2 = t;
            for (auto& c : t_m2.coefficients()) c = fp_mod(c, m);
            t_m2.normalize([](const BigInt& v) { return v.is_zero(); });

            // Solve G'_m·Δs + H'_m·Δt ≡ δ (mod m) with deg(Δt) < deg(G'_m):
            auto [q_s, dt] = poly_divrem_mod(poly_mul_mod(t_m2, delta, m), G_m2, m);
            IntPoly ds = poly_add_mod(poly_mul_mod(s_m2, delta, m), poly_mul_mod(q_s, H_m2, m), m);

            // s' = s + Δs·m
            if (s.size() < ds.size()) s.resize(ds.size(), BigInt(0));
            for (std::size_t i = 0; i < ds.size(); ++i) s[i] += ds[i] * m;
            s.normalize([](const BigInt& v) { return v.is_zero(); });

            // t' = t + Δt·m
            if (t.size() < dt.size()) t.resize(dt.size(), BigInt(0));
            for (std::size_t i = 0; i < dt.size(); ++i) t[i] += dt[i] * m;
            t.normalize([](const BigInt& v) { return v.is_zero(); });

            poly_balanced_mod(s, m2);
            poly_balanced_mod(t, m2);
        }
    }
}

// Univariate Quadratic Hensel Lift.
// Given f ∈ Z[x], factors g, h ∈ Z[x] with g*h ≡ f (mod p) and gcd(g mod p, h mod p) = 1,
// lifts to G, H with G*H ≡ f (mod p^k), using quadratic steps (modulus squares each step).
// Number of steps: ceil(log2(k)), each squaring the modulus.
Result<std::pair<IntPoly, IntPoly>> hensel_lift(
    const IntPoly& f,
    const IntPoly& g,
    const IntPoly& h,
    const BigInt& p,
    std::size_t k) {

    if (k == 0) return ok(std::make_pair(g, h));

    // Initialize with balanced representatives mod p
    IntPoly G = g;
    poly_balanced_mod(G, p);
    IntPoly H = h;
    poly_balanced_mod(H, p);

    // Find Bézout cofactors s, t: s*G + t*H ≡ 1 (mod p)
    // (use fp versions for xgcd)
    IntPoly G_fp = G, H_fp = H;
    for (auto& c : G_fp.coefficients()) c = fp_mod(c, p);
    G_fp.normalize([](const BigInt& v) { return v.is_zero(); });
    for (auto& c : H_fp.coefficients()) c = fp_mod(c, p);
    H_fp.normalize([](const BigInt& v) { return v.is_zero(); });

    auto [gcd_p, s, t] = poly_xgcd_fp(G_fp, H_fp, p);
    // s, t are in [0, p) form; convert to balanced
    poly_balanced_mod(s, p);
    poly_balanced_mod(t, p);

    // Compute target modulus: p^k
    BigInt target_modulus(1);
    for (std::size_t i = 0; i < k; ++i) target_modulus *= p;

    // k=1: already at target (balanced mod p = g mod p)
    if (k == 1) {
        return ok(std::make_pair(G, H));
    }

    // Quadratic lifting: modulus doubles in exponent each step
    BigInt modulus = p;
    while (modulus < target_modulus) {
        BigInt next_modulus = modulus * modulus;
        if (next_modulus > target_modulus) next_modulus = target_modulus;

        quadratic_step(f, G, H, s, t, modulus, next_modulus);
        modulus = next_modulus;
    }

    return ok(std::make_pair(G, H));
}

// Multi-factor Hensel lifting (divide-and-conquer on factor groups)
Result<std::vector<IntPoly>> hensel_lift_multi(
    const IntPoly& f,
    const std::vector<IntPoly>& factors,
    const BigInt& p,
    std::size_t k) {

    // Base case: a single factor.  The caller has already lifted the *group
    // product* `f` to modulus p^k (it equals this factor's lift), so the lifted
    // representative of this singleton is `f` itself — NOT the original mod-p
    // input `factors[0]`, which is only correct mod p.  Returning the mod-p
    // factor here breaks the invariant ∏ gᵢ ≡ f (mod p^k) for every recombination.
    if (factors.empty()) return ok(factors);
    if (factors.size() == 1) {
        IntPoly lifted = f;
        poly_balanced_mod(lifted, [&] {
            BigInt m(1);
            for (std::size_t i = 0; i < k; ++i) m *= p;
            return m;
        }());
        return ok(std::vector<IntPoly>{std::move(lifted)});
    }

    // Split factors into two groups
    std::size_t mid = factors.size() / 2;
    std::vector<IntPoly> left_factors(factors.begin(), factors.begin() + mid);
    std::vector<IntPoly> right_factors(factors.begin() + mid, factors.end());

    // Compute product of each group mod p (fp form)
    IntPoly g{std::vector<BigInt>{BigInt(1)}}, h{std::vector<BigInt>{BigInt(1)}};
    for (const auto& poly : left_factors) {
        IntPoly fp_poly = poly;
        for (auto& c : fp_poly.coefficients()) c = fp_mod(c, p);
        fp_poly.normalize([](const BigInt& v) { return v.is_zero(); });
        g = poly_mul_fp(g, fp_poly, p);
    }
    for (const auto& poly : right_factors) {
        IntPoly fp_poly = poly;
        for (auto& c : fp_poly.coefficients()) c = fp_mod(c, p);
        fp_poly.normalize([](const BigInt& v) { return v.is_zero(); });
        h = poly_mul_fp(h, fp_poly, p);
    }

    // Lift f = g*h mod p to f = G*H mod p^k using quadratic lifting
    auto lift_res = hensel_lift(f, g, h, p, k);
    if (lift_res.is_error()) return fail<std::vector<IntPoly>>(lift_res.error());

    auto [G, H] = lift_res.value();

    // Recursively lift groups
    auto left_lifted = hensel_lift_multi(G, left_factors, p, k);
    if (left_lifted.is_error()) return fail<std::vector<IntPoly>>(left_lifted.error());

    auto right_lifted = hensel_lift_multi(H, right_factors, p, k);
    if (right_lifted.is_error()) return fail<std::vector<IntPoly>>(right_lifted.error());

    std::vector<IntPoly> result = left_lifted.value();
    result.insert(result.end(), right_lifted.value().begin(), right_lifted.value().end());
    return ok(result);
}

} // namespace cas::algebra
