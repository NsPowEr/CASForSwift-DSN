// polynomial_gcd_crt.cpp — Modular GCD via multi-prime CRT (B2.1, PLAN F2.1).
//
// Algorithm: von zur Gathen & Gerhard §6.7; Geddes-Czapor-Labahn §7.4.
//
// Invariants tracked:
//   1. deg(gcd_true) = min{ deg(gcd_p) : p lucky }  (degree-stability criterion)
//   2. Bad prime: skip if p | lc(f) OR p | lc(g) OR deg(gcd_p) > min_deg_seen
//   3. CRT accumulates until M = ∏primes > 2 * Mignotte_bound(f,g)
//      Mignotte: |coeff_k(gcd)| ≤ 2^{min(deg f, deg g)} * min(||f||_∞, ||g||_∞)
//   4. Divisibility certificate: gcd | f and gcd | g in Z[x]; otherwise add primes.
//   5. Centered representation: if 2*r > M then r ← r - M.
//
// Wiring: gcd_integer_poly_dispatch() calls this when both polynomials have
// max coefficient bit-length ≥ ctx.modular_gcd_coeff_bits() (default 48).

#include "polynomial_internal.hpp"
#include "cas/numtheory.hpp"
#include "cas/symbolic.hpp"
#include "cas/error_helpers.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace cas::algebra {

// ── helpers ────────────────────────────────────────────────────────────────

static BigInt pos_mod(const BigInt& a, const BigInt& m) {
    BigInt r = a % m;
    if (r.is_negative()) r += m.abs();
    return r;
}

// Reduce polynomial coefficients to [0, p).
static IntPoly reduce_mod_p(const IntPoly& f, const BigInt& p) {
    IntPoly r;
    r.resize(f.size(), BigInt(0));
    for (std::size_t i = 0; i < f.size(); ++i) r[i] = pos_mod(f[i], p);
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

// Infinity norm: max |coeff|.
static BigInt inf_norm(const IntPoly& f) {
    BigInt best(0);
    for (const BigInt& c : f.coefficients()) {
        BigInt ac = c.abs();
        if (ac > best) best = ac;
    }
    return best;
}

// Mignotte GCD bound: |coeff_k(gcd)| ≤ 2^{min_deg} * min(||f||_∞, ||g||_∞).
// Reference: Geddes-Czapor-Labahn §7.4 Lemma 7.3.
static BigInt mignotte_bound(const IntPoly& f, const IntPoly& g) {
    std::size_t min_deg = std::min(f.degree(), g.degree());
    BigInt nf = inf_norm(f);
    BigInt ng = inf_norm(g);
    BigInt mn = (nf < ng) ? nf : ng;
    if (mn.is_zero()) mn = BigInt(1);
    return BigInt(1).shift_left_bits(min_deg) * mn;
}

// GCD in Fp[x], returns monic result; empty poly on degenerate input.
static IntPoly fp_gcd(IntPoly a, IntPoly b, const BigInt& p) {
    a = reduce_mod_p(a, p);
    b = reduce_mod_p(b, p);
    while (!b.is_zero()) {
        // Exact division in Fp: compute a mod b.
        auto inv_lc = numtheory::modular_inverse(pos_mod(b.leading_coeff(), p), p);
        if (inv_lc.is_error()) return IntPoly{};
        BigInt ilc = inv_lc.value();
        IntPoly rem = a;
        while (!rem.is_zero() && rem.degree() >= b.degree()) {
            std::size_t dd = rem.degree() - b.degree();
            BigInt factor = pos_mod(rem.leading_coeff() * ilc, p);
            for (std::size_t i = 0; i < b.size(); ++i) {
                rem[i + dd] = pos_mod(rem[i + dd] - factor * b[i], p);
            }
            rem.normalize([](const BigInt& v) { return v.is_zero(); });
        }
        a = std::move(b);
        b = std::move(rem);
    }
    // Make monic.
    if (!a.is_zero()) {
        auto inv_lc = numtheory::modular_inverse(pos_mod(a.leading_coeff(), p), p);
        if (inv_lc.is_error()) return IntPoly{};
        BigInt ilc = inv_lc.value();
        for (auto& c : a.coefficients()) c = pos_mod(c * ilc, p);
    }
    return a;
}

// Centered representative of r ∈ [0, M): if 2r > M return r - M else r.
static BigInt centered_repr(const BigInt& r, const BigInt& M) {
    BigInt two_r = r + r;
    return (two_r > M) ? (r - M) : r;
}

// Verify exact divisibility in Z[x]: candidate | dividend iff prem = 0.
static bool divides_z(const IntPoly& dividend, const IntPoly& candidate) {
    if (candidate.is_zero()) return false;
    if (dividend.is_zero()) return true;
    if (candidate.degree() == 0) return true;
    if (candidate.degree() > dividend.degree()) return false;
    IntPoly prim = primitive_integer_poly(candidate);
    IntPoly rem = pseudo_remainder_integer_poly(dividend, prim);
    normalize_integer_poly(rem);
    return rem.is_zero();
}

// ── Public API ─────────────────────────────────────────────────────────────

Result<IntPoly> gcd_integer_poly_crt(
    const IntPoly& f, const IntPoly& g, const symbolic::CASContext& ctx)
{
    // Trivial zero cases.
    if (f.is_zero()) {
        if (g.is_zero()) return ok(IntPoly{});
        IntPoly r = primitive_integer_poly(g);
        if (!r.is_zero() && r.leading_coeff().is_negative())
            multiply_integer_coefficients_by_scalar(r, BigInt(-1));
        BigInt cg = integer_content(g);
        multiply_integer_coefficients_by_scalar(r, cg);
        return ok(std::move(r));
    }
    if (g.is_zero()) {
        IntPoly r = primitive_integer_poly(f);
        if (!r.is_zero() && r.leading_coeff().is_negative())
            multiply_integer_coefficients_by_scalar(r, BigInt(-1));
        BigInt cf = integer_content(f);
        multiply_integer_coefficients_by_scalar(r, cf);
        return ok(std::move(r));
    }

    // Content of GCD = gcd of contents.
    BigInt cont_f = integer_content(f);
    BigInt cont_g = integer_content(g);
    BigInt cont_gcd = gcd(cont_f, cont_g);

    // Work with primitive parts.
    IntPoly pf = primitive_integer_poly(f);
    IntPoly pg = primitive_integer_poly(g);

    // Leading coefficient of true GCD divides gcd(lc(pf), lc(pg)).
    BigInt lc_f = pf.leading_coeff();
    BigInt lc_g = pg.leading_coeff();
    BigInt lc_bound = gcd(lc_f.abs(), lc_g.abs());

    // Mignotte bound on primitive part coefficients.
    const BigInt M_target = mignotte_bound(pf, pg);
    // Need M > 2 * lc_bound * M_target to safely recover scaled GCD.
    BigInt M_need = (M_target + M_target) * lc_bound + BigInt(1);

    const std::size_t max_gcd_deg = std::min(pf.degree(), pg.degree());

    // CRT state: one solution per coefficient position 0..max_gcd_deg.
    const std::size_t n_slots = max_gcd_deg + 1U;
    std::vector<BigInt> solutions(n_slots, BigInt(0));
    BigInt M_acc(1);  // accumulated modulus

    std::size_t min_deg_seen = std::numeric_limits<std::size_t>::max();
    bool have_lucky = false;

    // Starting prime: near 2^30 + hash-offset to vary across poly inputs.
    std::size_t h = lc_f.bit_length() * 2654435761ULL;
    h ^= pf.size() * 40503ULL ^ pg.size() * 12347ULL;
    long long start_val = 1073741827LL + static_cast<long long>(h % 65537ULL);
    auto np0 = numtheory::next_prime(BigInt(start_val - 1));
    BigInt current_prime = np0.is_ok() ? np0.value() : BigInt(1073741827LL);

    const std::size_t max_primes = ctx.max_gcd_total_calls();
    std::size_t primes_used = 0;

    while (primes_used < max_primes) {
        // Advance current_prime for next iteration.
        auto np_next = numtheory::next_prime(current_prime);
        const BigInt p = current_prime;
        if (np_next.is_ok()) current_prime = np_next.value();
        else break;

        // Bad prime check: p | lc(pf) or p | lc(pg).
        if ((lc_f % p).is_zero() || (lc_g % p).is_zero()) continue;

        // Compute modular GCD.
        IntPoly gp = fp_gcd(pf, pg, p);
        if (gp.is_zero()) continue;

        const std::size_t deg_p = gp.degree();

        // Bad prime check: degree higher than best seen → unlucky.
        if (have_lucky && deg_p > min_deg_seen) continue;

        // Degree dropped → reset (better degree estimate).
        if (!have_lucky || deg_p < min_deg_seen) {
            min_deg_seen = deg_p;
            have_lucky = true;
            std::fill(solutions.begin(), solutions.end(), BigInt(0));
            M_acc = BigInt(1);
            primes_used = 0;
        }

        // Scale monic gp by lc_bound mod p → matches the true GCD lc scaling.
        BigInt scale = pos_mod(lc_bound, p);
        if (!scale.is_zero()) {
            for (auto& c : gp.coefficients()) c = pos_mod(c * scale, p);
        }

        // CRT-merge each coefficient independently using the *same* prime p.
        // Invariant: all solution[i] are in [0, M_acc) before this prime.
        // After this prime they will be in [0, M_acc * p).
        // We track M_acc only once: do the merge for coeff 0 to get new M_acc,
        // then for i > 0 use the *old* M_acc (before update) consistently.
        const BigInt M_before = M_acc;
        const std::size_t n_acc = min_deg_seen + 1U;
        if (solutions.size() < n_acc) {
            solutions.resize(n_acc, BigInt(0));
        }

        for (std::size_t i = 0; i < n_acc; ++i) {
            BigInt ri = (i < gp.size()) ? gp[i] : BigInt(0);
            // Each slot must merge individually but with same p and SAME M_before.
            // We cannot reuse crt_merge directly because it updates M_acc.
            // Inline the step with M_before:
            BigInt cur_mod_p = pos_mod(solutions[i], p);
            BigInt delta_i = pos_mod(ri - cur_mod_p, p);
            auto inv_m = numtheory::modular_inverse(pos_mod(M_before, p), p);
            if (inv_m.is_error()) goto next_prime_lbl;
            BigInt t_i = pos_mod(delta_i * inv_m.value(), p);
            solutions[i] = solutions[i] + M_before * t_i;
        }
        M_acc = M_before * p;

        ++primes_used;

        // Check sufficiency.
        if (M_acc > M_need) {
            // Reconstruct candidate via centered representation.
            IntPoly candidate;
            candidate.resize(n_acc, BigInt(0));
            for (std::size_t i = 0; i < n_acc; ++i) {
                candidate[i] = centered_repr(solutions[i], M_acc);
            }
            candidate.normalize([](const BigInt& v) { return v.is_zero(); });

            // Remove lc_bound scaling: divide out content.
            IntPoly prim_cand = primitive_integer_poly(candidate);

            // Divisibility certificate (invariant: MAI return senza verifica).
            if (divides_z(pf, prim_cand) && divides_z(pg, prim_cand)) {
                // Re-attach content.
                multiply_integer_coefficients_by_scalar(prim_cand, cont_gcd);
                if (!prim_cand.is_zero() && prim_cand.leading_coeff().is_negative())
                    multiply_integer_coefficients_by_scalar(prim_cand, BigInt(-1));
                return ok(std::move(prim_cand));
            }
            // Certificate failed: keep accumulating, double the target.
            M_need = M_need + M_need;
        }

        continue;
next_prime_lbl:
        continue;
    }

    return make_unimplemented<IntPoly>(
        "algebra", "gcd_integer_poly_crt",
        "prime_budget=" + std::to_string(primes_used),
        "MODULAR_GCD_PRIME_BUDGET_EXHAUSTED",
        "Increase ctx.max_gcd_total_calls() or use subresultant path",
        "B2.1");
}

}  // namespace cas::algebra
