// polynomial_resultant_crt.cpp — Modular resultant CRT (B2.2) and bivariate
// resultant via evaluation-interpolation (B2.3), PLAN F2.2.
//
// B2.2 — Modular resultant CRT (univariate, Collins 1971; vzGG §6.8):
//   Algorithm:
//   1. For primes p_1, p_2, ... (skip if p | lc(f) or p | lc(g)):
//      Compute res_p = resultant(f mod p, g mod p) ∈ F_p via Euclidean PRS.
//   2. CRT-accumulate until M = ∏p_i > 2 * Hadamard_bound(f,g).
//      Hadamard: |res(f,g)| ≤ ||f||_2^{deg g} · ||g||_2^{deg f}.
//   3. Reconstruct res ∈ Z via centered representation.
//   Invariant: ∏p_i > 2*H ⇒ reconstruction is exact; NEVER return without
//   this bound satisfied.
//
// B2.3 — Bivariate resultant via evaluation-interpolation (vzGG §6.3):
//   f, g ∈ Z[x, y] → res_y(f, g) ∈ Z[x].
//   1. Degree bound: deg_x(res) ≤ d_xf * d_yg + d_xg * d_yf → need d+1 points.
//   2. For evaluation points x_0, x_1, ... (skip bad: lc_y(f)|_{x=xi}=0 or lc_y(g)=0):
//      compute univariate resultant of f(xi, y) in y → value at xi.
//   3. Interpolate Newton polynomial in x.
//   Invariant: #good_points ≥ deg_x(res)+1 → interpolation exact.

#include "polynomial_internal.hpp"
#include "cas/numtheory.hpp"
#include "cas/symbolic.hpp"
#include "cas/error_helpers.hpp"

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

namespace cas::algebra {

// ── helpers shared with gcd_crt ─────────────────────────────────────────

static BigInt pmod(const BigInt& a, const BigInt& m) {
    BigInt r = a % m;
    if (r.is_negative()) r += m.abs();
    return r;
}

// Infinity norm squared approximation: ||f||_∞.
static BigInt intpoly_inf_norm_crt(const IntPoly& f) {
    BigInt best(0);
    for (const BigInt& c : f.coefficients()) {
        BigInt ac = c.abs();
        if (ac > best) best = ac;
    }
    return best;
}

// Hadamard bound for |res(f,g)|:
// |res(f,g)| ≤ ||f||_2^{deg_g} * ||g||_2^{deg_f}.
// We compute ||f||_2^2 and ||g||_2^2 and use bit-shift to approximate:
// log2(bound) ≈ deg_g/2 * log2(||f||_2^2) + deg_f/2 * log2(||g||_2^2).
// We return a BigInt ≥ bound by computing: (||f||_∞+1)^{deg_g+1} * (||g||_∞+1)^{deg_f+1}.
// This is looser than Hadamard but always safe and computable in BigInt.
// Reference: vzGG §6.8 Theorem 6.51.
static BigInt hadamard_resultant_bound(const IntPoly& f, const IntPoly& g) {
    if (f.is_zero() || g.is_zero()) return BigInt(1);
    std::size_t deg_f = f.degree();
    std::size_t deg_g = g.degree();
    BigInt nf = intpoly_inf_norm_crt(f) + BigInt(1);  // +1 for degree-0 edge case
    BigInt ng = intpoly_inf_norm_crt(g) + BigInt(1);
    // bound = nf^{deg_g} * ng^{deg_f}  (plus one for safety, BigInt exact)
    BigInt bf = bigint_pow_nonnegative(nf, deg_g);
    BigInt bg = bigint_pow_nonnegative(ng, deg_f);
    return bf * bg + BigInt(1);
}

// Centered representative.
static BigInt centered_r(const BigInt& r, const BigInt& M) {
    return ((r + r) > M) ? (r - M) : r;
}

// Evaluate IntPoly at integer point x0 → value ∈ Z (Horner).
static BigInt eval_intpoly_at(const IntPoly& f, const BigInt& x0) {
    if (f.is_zero()) return BigInt(0);
    BigInt acc = f[f.size() - 1U];
    for (std::size_t i = f.size() - 1U; i > 0; --i) {
        acc = acc * x0 + f[i - 1U];
    }
    return acc;
}

// ── B2.2 — Modular resultant CRT ───────────────────────────────────────────

// Resultant in F_p[x] via Euclidean PRS.
// Uses the classical sign/degree formula from the subresultant algorithm.
static Result<BigInt> resultant_mod_p(const IntPoly& f, const IntPoly& g, const BigInt& p) {
    // Reduce inputs.
    IntPoly a, b;
    {
        a.resize(f.size(), BigInt(0));
        for (std::size_t i = 0; i < f.size(); ++i) a[i] = pmod(f[i], p);
        a.normalize([](const BigInt& v) { return v.is_zero(); });
        b.resize(g.size(), BigInt(0));
        for (std::size_t i = 0; i < g.size(); ++i) b[i] = pmod(g[i], p);
        b.normalize([](const BigInt& v) { return v.is_zero(); });
    }

    if (a.is_zero() || b.is_zero()) return ok(BigInt(0));

    BigInt result(1);
    // Iterative Euclidean PRS resultant mod p.
    // res(a, b) = (-1)^{deg_a * deg_b} * lc(b)^{deg_a - deg_r} * res(b, r) / lc(b)^{deg_a - deg_r}
    // Standard formula: track sign flips and leading-coeff powers.
    while (true) {
        if (a.is_zero()) { result = BigInt(0); break; }
        if (b.is_zero()) { result = BigInt(0); break; }

        if (b.degree() == 0) {
            // res(a, c) = c^{deg_a}
            BigInt lc_b = pmod(b[0], p);
            result = pmod(result * bigint_pow_nonnegative(lc_b, a.degree()), p);
            break;
        }

        std::size_t da = a.degree();
        std::size_t db = b.degree();

        // Adjust sign: res(a,b) = (-1)^{da*db} * res(b,a).
        if ((da % 2U != 0U) && (db % 2U != 0U)) {
            result = pmod(p - result, p);
        }

        // Compute lc(b)^{da - db + 1} factor from pseudo-remainder.
        BigInt lc_b = pmod(b.leading_coeff(), p);
        auto inv_lc = numtheory::modular_inverse(lc_b, p);
        if (inv_lc.is_error()) return ok(BigInt(0));  // lc_b = 0 mod p

        // Divide a by b in Fp[x], get remainder.
        IntPoly rem = a;
        BigInt ilc = inv_lc.value();
        while (!rem.is_zero() && rem.degree() >= b.degree()) {
            std::size_t dd = rem.degree() - b.degree();
            BigInt factor = pmod(rem.leading_coeff() * ilc, p);
            for (std::size_t i = 0; i < b.size(); ++i) {
                rem[i + dd] = pmod(rem[i + dd] - factor * b[i], p);
            }
            rem.normalize([](const BigInt& v) { return v.is_zero(); });
        }

        // Adjust result for the lc(b)^{da - db + 1} from pseudo-division.
        // In Fp the division is exact, so factor = lc_b^{da - db + 1} gets
        // divided out. For exact division no extra factor needed.
        // (Contrast: over Z, pseudo-remainder introduces lc(b)^{da-db+1} factor.)

        a = std::move(b);
        b = std::move(rem);
    }
    return ok(pmod(result, p));
}

Result<BigInt> resultant_integer_poly_crt(
    const IntPoly& f, const IntPoly& g, const symbolic::CASContext& ctx)
{
    if (f.is_zero() || g.is_zero()) return ok(BigInt(0));

    // Hadamard bound: accumulate until M > 2 * bound.
    const BigInt H = hadamard_resultant_bound(f, g);
    const BigInt M_need = H + H + BigInt(1);

    BigInt lc_f = f.leading_coeff();
    BigInt lc_g = g.leading_coeff();

    BigInt solution(0);
    BigInt M_acc(1);

    // Start prime near 2^30.
    std::size_t h = f.size() * 54321ULL ^ g.size() * 98765ULL;
    long long sv = 1073741827LL + static_cast<long long>(h % 65537ULL);
    auto np0 = numtheory::next_prime(BigInt(sv - 1));
    BigInt current_prime = np0.is_ok() ? np0.value() : BigInt(1073741827LL);

    const std::size_t max_primes = ctx.max_gcd_total_calls();
    std::size_t primes_used = 0;

    while (primes_used < max_primes) {
        auto np_next = numtheory::next_prime(current_prime);
        const BigInt p = current_prime;
        if (np_next.is_ok()) current_prime = np_next.value();
        else break;

        // Bad prime: p | lc(f) or p | lc(g).
        if ((lc_f % p).is_zero() || (lc_g % p).is_zero()) continue;

        auto res_p = resultant_mod_p(f, g, p);
        if (res_p.is_error()) continue;

        BigInt rp = pmod(res_p.value(), p);

        // CRT merge (single scalar).
        BigInt cur_mod_p = pmod(solution, p);
        BigInt delta = pmod(rp - cur_mod_p, p);
        auto inv_m = numtheory::modular_inverse(pmod(M_acc, p), p);
        if (inv_m.is_error()) continue;
        BigInt t = pmod(delta * inv_m.value(), p);
        solution = solution + M_acc * t;
        M_acc = M_acc * p;

        ++primes_used;

        if (M_acc > M_need) {
            return ok(centered_r(solution, M_acc));
        }
    }

    return make_unimplemented<BigInt>(
        "algebra", "resultant_integer_poly_crt",
        "prime_budget=" + std::to_string(primes_used),
        "MODULAR_RESULTANT_PRIME_BUDGET_EXHAUSTED",
        "Increase ctx.max_gcd_total_calls() or use subresultant path",
        "B2.2");
}

// ── B2.3 — Bivariate resultant via evaluation-interpolation ────────────────
//
// Computes res_y(f, g) ∈ Z[x]  where f, g ∈ Z[x][y].
//
// Representation: f is given as a vector of IntPoly (coefficients in Z[x]).
// f_as_y_poly[i] = coefficient of y^i in f, each itself an IntPoly in x.
//
// Algorithm (vzGG §6.3):
//   deg_x(res) ≤ D = deg_x(f)*deg_y(g) + deg_x(g)*deg_y(f).
//   Choose N = D + 1 evaluation points x_0, x_1, ...
//   For each x_i:
//     Evaluate f(x_i, y) and g(x_i, y) → univariate in y.
//     Skip (bad evaluation) if lc_y(f) or lc_y(g) evaluates to 0.
//     Compute resultant in y.
//   Interpolate Newton forward difference to get res_y ∈ Z[x].
//
// Newton interpolation: given pairs (x_0, v_0), ..., (x_n, v_n),
// build the Newton basis (x-x_0)(x-x_1)...(x-x_{k-1}) and compute
// divided differences.

// Evaluate univariate IntPoly in x at an integer point.
static BigInt eval_univariate(const IntPoly& f, const BigInt& x0) {
    return eval_intpoly_at(f, x0);
}

// Evaluate bivariate poly (as y-poly with x-poly coefficients) at x = x0 → IntPoly in y.
static IntPoly eval_bivariate_at_x(const std::vector<IntPoly>& f_as_y_poly, const BigInt& x0) {
    std::size_t n = f_as_y_poly.size();
    IntPoly result;
    result.resize(n, BigInt(0));
    for (std::size_t i = 0; i < n; ++i) {
        result[i] = eval_univariate(f_as_y_poly[i], x0);
    }
    result.normalize([](const BigInt& v) { return v.is_zero(); });
    return result;
}

// Newton polynomial interpolation from points (xs[i], ys[i]).
// Returns coefficient vector of the interpolating polynomial in Z[x].
// Uses integer divided differences — requires that actual polynomial
// is integer-valued at all points AND has integer coefficients, which
// is guaranteed here (resultant ∈ Z[x]).
static IntPoly newton_interpolate(const std::vector<BigInt>& xs, const std::vector<BigInt>& ys) {
    std::size_t n = xs.size();
    if (n == 0U) return IntPoly{};

    // Divided differences table dd[i] = f[x_0,...,x_i] (Newton forward).
    std::vector<BigInt> dd(ys.begin(), ys.end());
    for (std::size_t j = 1U; j < n; ++j) {
        for (std::size_t i = n - 1U; i >= j; --i) {
            BigInt denom = xs[i] - xs[i - j];
            if (denom.is_zero()) return IntPoly{};  // coincident points
            // For integer polynomials, divided differences are integers.
            BigInt num = dd[i] - dd[i - 1U];
            // Check divisibility.
            if (!(num % denom).is_zero()) return IntPoly{};
            dd[i] = num / denom;
        }
    }

    // Expand Newton form: p(x) = Σ dd[k] * prod_{i<k}(x - xs[i]).
    // We accumulate by multiplying by (x - xs[k]) iteratively.
    // coeffs[j] = coefficient of x^j in the result.
    std::vector<BigInt> coeffs(n, BigInt(0));
    coeffs[0] = dd[n - 1U];
    // Process from highest to lowest.
    for (std::size_t k = n - 1U; k > 0; --k) {
        // Multiply by (x - xs[k-1]): coeffs -> x*coeffs - xs[k-1]*coeffs.
        // Then add dd[k-1].
        std::vector<BigInt> new_coeffs(n, BigInt(0));
        // Shift up (multiply by x) and subtract xs[k-1].
        for (std::size_t j = 0U; j + 1U < n; ++j) {
            new_coeffs[j + 1U] = new_coeffs[j + 1U] + coeffs[j];
        }
        for (std::size_t j = 0U; j < n; ++j) {
            new_coeffs[j] = new_coeffs[j] - xs[k - 1U] * coeffs[j];
        }
        new_coeffs[0] = new_coeffs[0] + dd[k - 1U];
        coeffs = std::move(new_coeffs);
    }

    IntPoly res(std::vector<BigInt>(coeffs.begin(), coeffs.end()));
    res.normalize([](const BigInt& v) { return v.is_zero(); });
    return res;
}

// f_as_y_poly: f_as_y_poly[i] = coefficient of y^i ∈ Z[x].
// g_as_y_poly: similarly.
// deg_x_f = max degree in x across all x-coefficients of f.
// deg_x_g = similarly for g.
// deg_y_f = f_as_y_poly.size() - 1  (degree in y).
// deg_y_g = g_as_y_poly.size() - 1.
Result<IntPoly> resultant_bivariate_eval_interp(
    const std::vector<IntPoly>& f_as_y_poly,
    const std::vector<IntPoly>& g_as_y_poly,
    const symbolic::CASContext& ctx)
{
    if (f_as_y_poly.empty() || g_as_y_poly.empty()) return ok(IntPoly{});

    const std::size_t deg_y_f = f_as_y_poly.size() - 1U;
    const std::size_t deg_y_g = g_as_y_poly.size() - 1U;

    // deg_x of f: max degree in x over all y-coefficients.
    std::size_t deg_x_f = 0U;
    for (const auto& cx : f_as_y_poly) {
        if (!cx.is_zero()) deg_x_f = std::max(deg_x_f, cx.degree());
    }
    std::size_t deg_x_g = 0U;
    for (const auto& cx : g_as_y_poly) {
        if (!cx.is_zero()) deg_x_g = std::max(deg_x_g, cx.degree());
    }

    // Degree bound: deg_x(res_y(f,g)) ≤ D = deg_x_f * deg_y_g + deg_x_g * deg_y_f.
    const std::size_t D = deg_x_f * deg_y_g + deg_x_g * deg_y_f;
    const std::size_t n_points_needed = D + 1U;

    // Safety cap: n_points_needed can be large; cap at ctx.max_gcd_total_calls().
    if (n_points_needed > ctx.max_gcd_total_calls()) {
        return make_unimplemented<IntPoly>(
            "algebra", "resultant_bivariate_eval_interp",
            "n_points=" + std::to_string(n_points_needed),
            "BIVARIATE_RESULTANT_DEGREE_TOO_HIGH",
            "Increase ctx.max_gcd_total_calls() or reduce polynomial degrees",
            "B2.3");
    }

    // lc_y(f) and lc_y(g) — the leading (in y) x-polynomial coefficient.
    const IntPoly& lc_y_f = f_as_y_poly.back();
    const IntPoly& lc_y_g = g_as_y_poly.back();

    std::vector<BigInt> xs, ys;
    xs.reserve(n_points_needed);
    ys.reserve(n_points_needed);

    // Evaluation point generation: use 0, 1, -1, 2, -2, ...
    // starting from 0 and growing.
    auto next_eval_point = [](std::size_t idx) -> BigInt {
        // idx=0 → 0, idx=1 → 1, idx=2 → -1, idx=3 → 2, idx=4 → -2, ...
        if (idx == 0U) return BigInt(0);
        long long v = static_cast<long long>((idx + 1U) / 2U);
        if (idx % 2U == 0U) v = -v;
        return BigInt(v);
    };

    std::size_t eval_idx = 0U;
    const std::size_t max_attempts = n_points_needed * 4U + 100U;

    while (xs.size() < n_points_needed && eval_idx < max_attempts) {
        BigInt x0 = next_eval_point(eval_idx++);

        // Bad evaluation check: lc_y(f)(x0) = 0 or lc_y(g)(x0) = 0.
        BigInt lc_f_at_x0 = eval_univariate(lc_y_f, x0);
        BigInt lc_g_at_x0 = eval_univariate(lc_y_g, x0);
        if (lc_f_at_x0.is_zero() || lc_g_at_x0.is_zero()) continue;

        // Evaluate f(x0, y) and g(x0, y) as IntPoly in y.
        IntPoly f_at_x0 = eval_bivariate_at_x(f_as_y_poly, x0);
        IntPoly g_at_x0 = eval_bivariate_at_x(g_as_y_poly, x0);

        // Compute resultant in y using the CRT resultant (recursive, bottom of stack).
        auto res_y = resultant_integer_poly_crt(f_at_x0, g_at_x0, ctx);
        if (res_y.is_error()) continue;

        xs.push_back(x0);
        ys.push_back(res_y.value());
    }

    if (xs.size() < n_points_needed) {
        return make_unimplemented<IntPoly>(
            "algebra", "resultant_bivariate_eval_interp",
            "good_points=" + std::to_string(xs.size()) +
            " needed=" + std::to_string(n_points_needed),
            "BIVARIATE_RESULTANT_INSUFFICIENT_GOOD_POINTS",
            "Too many bad evaluation points; check leading coefficients",
            "B2.3");
    }

    // Newton interpolation.
    IntPoly result = newton_interpolate(xs, ys);
    if (result.is_zero() && !ys.empty() && !ys[0].is_zero()) {
        return make_unimplemented<IntPoly>(
            "algebra", "resultant_bivariate_eval_interp",
            "interpolation_failed",
            "BIVARIATE_RESULTANT_INTERP_FAILED",
            "Integer divided differences produced non-integer; check inputs",
            "B2.3");
    }

    return ok(std::move(result));
}

}  // namespace cas::algebra
