// F3.2 — Multivariate (ideal-adic) Hensel lifting for Wang's EEZ algorithm.
// Reference: GCL "Algorithms for Computer Algebra" §6.4 (multivariate Hensel
// construction) and §6.5 (Wang's EEZ).  All arithmetic exact over Z / Q.
//
// Strategy: lift one variable at a time.  After lifting x_2..x_j we hold factors
// U_i in Z[x_1,...,x_j] that are correct modulo the ideal <x_{j+1}-a_{j+1},...>.
// We then lift in x_{j+1} using a Taylor (ideal-adic) Newton iteration:
//   for k = 1 .. deg_{x_{j+1}}(a):
//     error = a - prod(U_i)              (current, in x_1..x_{j+1})
//     c_k   = coefficient of (x_{j+1}-a_{j+1})^k in error
//     solve the multivariate Diophantine  sum_i (a/U_i)|_eval * delta_i = c_k|_eval
//     U_i += delta_i * (x_{j+1}-a_{j+1})^k
//
// The Diophantine system is solved by reducing all coefficients to the ideal of
// already-lifted variables evaluated at the point, giving a UNIVARIATE multi-term
// Diophantine in x_1 solved through the precomputed univariate Bezout data.
//
// Caps: the ideal-adic degree is the EXACT degree of `a` in each variable (a true
// math bound — no truncation of valid math).  Any structural case this routine
// cannot certify returns explicit Unimplemented (never a silent wrong result):
// the driver re-multiplies and checks exact equality regardless.

#include "factor_multivariate_internal.hpp"

#include <algorithm>
#include <numeric>

namespace cas::algebra {

namespace {

// ---- Univariate (in x_1) helpers over Q, using the existing RatPoly toolkit ----

[[nodiscard]] RatPoly intpoly_to_ratpoly(const IntPoly& f) {
    RatPoly r;
    r.coefficients().reserve(f.size());
    for (const auto& c : f.coefficients()) {
        r.push_back(Rational(c));
    }
    normalize_rational_coefficients(r);
    return r;
}

// Reduce an MPoly to a univariate-in-x_1 IntPoly by evaluating x_2..x_n at the
// point.  (Used to project Diophantine coefficients onto the base ideal.)
[[nodiscard]] IntPoly project_to_x1(
    const MPoly& p, const std::vector<BigInt>& point, const WangContext& wc) {
    MPoly cur = p;
    for (std::size_t vi = 1; vi < wc.nvars(); ++vi) {
        cur = mpoly_eval_var(cur, vi, point[vi - 1U]);
    }
    auto ip = mpoly_to_intpoly(cur, 0U);
    return ip.has_value() ? *ip : IntPoly{};
}

// Multiply MPoly by (x_var - a)^k.
[[nodiscard]] MPoly mul_by_shift_power(
    const MPoly& p, std::size_t var, const BigInt& a, unsigned int k,
    std::size_t nvars) {
    // (x_var - a): MPoly with two terms.
    MPoly shift;
    Monomial mx(nvars, 0U);
    mx[var] = 1U;
    shift.terms.emplace(mx, BigInt(1));
    if (!a.is_zero()) {
        shift.terms[Monomial(nvars, 0U)] = -a;
    }
    MPoly acc = p;
    for (unsigned int i = 0; i < k; ++i) {
        acc = mpoly_mul(acc, shift);
    }
    return acc;
}

// Taylor coefficient: coefficient of (x_var - a)^k in p, as an MPoly with the
// x_var exponent removed (i.e. living in the remaining-variable space at that
// order).  Computed by repeated synthetic division by (x_var - a).
[[nodiscard]] MPoly taylor_coeff(
    const MPoly& p, std::size_t var, const BigInt& a, unsigned int k) {
    // Shift p so that x_var -> x_var + a, then take coeff of x_var^k, i.e. the
    // ordinary coefficient_in after substitution.  We perform the substitution
    // x_var := y + a by binomial expansion per term.
    MPoly shifted;
    const std::size_t nvars = p.terms.empty() ? 0U : p.terms.begin()->first.size();
    for (const auto& [mono, coeff] : p.terms) {
        unsigned int e = mono[var];
        if (e == 0U) {
            BigInt& slot = shifted.terms[mono];
            slot += coeff;
            if (slot.is_zero()) shifted.terms.erase(mono);
            continue;
        }
        // (y + a)^e = sum_{j=0}^{e} C(e,j) a^{e-j} y^j
        BigInt binom(1);  // C(e,0)
        for (unsigned int j = 0; j <= e; ++j) {
            Monomial m = mono;
            m[var] = j;
            BigInt term = coeff * binom;
            if (e - j > 0U) {
                term *= bigint_pow_nonnegative(a, static_cast<std::size_t>(e - j));
            }
            BigInt& slot = shifted.terms[m];
            slot += term;
            if (slot.is_zero()) shifted.terms.erase(m);
            // update binom: C(e,j+1) = C(e,j) * (e-j)/(j+1)
            if (j < e) {
                binom *= BigInt(static_cast<long long>(e - j));
                binom /= BigInt(static_cast<long long>(j + 1U));
            }
        }
    }
    (void)nvars;
    return mpoly_coeff_in(shifted, var, k);
}

// Solve the r-term univariate Diophantine over Q:
//   sum_i sigma_i * (prod_{l!=i} u_l) == c     (mod nothing; exact, deg sigma_i < deg u_i)
// using the precomputed cofactors b_i = prod_{l!=i} u_l and a Bezout-style
// partial-fraction decomposition.  Returns sigma_i as RatPolys.
[[nodiscard]] std::optional<std::vector<RatPoly>> univariate_diophantine(
    const std::vector<RatPoly>& u,    // the (monic-ish) univariate factors
    const RatPoly& c) {
    const std::size_t r = u.size();
    if (r == 0U) return std::nullopt;

    // Product of all u.
    RatPoly prod;
    prod.push_back(Rational(BigInt(1)));
    for (const auto& f : u) {
        prod = mul_rational_poly(prod, f);
    }

    // We solve incrementally using the standard pairwise Bezout splitting:
    //   find s,t with s*A + t*B = 1, A = u_i, B = prod_{l!=i} u_l.
    std::vector<RatPoly> sigma(r);
    RatPoly c_work = c;
    // Remaining product for the "rest" factor.
    for (std::size_t i = 0; i < r; ++i) {
        RatPoly A = u[i];
        RatPoly B;
        B.push_back(Rational(BigInt(1)));
        for (std::size_t l = i + 1; l < r; ++l) {
            B = mul_rational_poly(B, u[l]);
        }
        if (i + 1 == r) {
            // last: sigma_i * 1 = c_work  -> reduce mod u_i
            auto [q, rem] = div_rem_rational_poly(c_work, A);
            (void)q;
            sigma[i] = rem;
            break;
        }
        // s*A + t*B = g; require g constant (coprime).
        auto [g, s, t] = extended_gcd_rational_poly(A, B);
        if (g.size() != 1U) {
            return std::nullopt;  // factors not coprime in this image
        }
        Rational ginv = Rational(BigInt(1)) / g[0];
        // c_work = c_work mod (A) gives sigma_i contribution via t? Use:
        // sigma_i ≡ c_work * t * ginv  (mod A);  then peel off and continue with B.
        RatPoly tc = mul_rational_poly(c_work, t);
        for (auto& coeff : tc.coefficients()) coeff *= ginv;
        auto [qi, ri] = div_rem_rational_poly(tc, A);
        sigma[i] = ri;
        // new c_work for the rest: (c_work - sigma_i * B_full_rest)/A ... but our
        // B already excludes u_i, and the rest factors are inside B.  Continue with
        // c' = (c_work - sigma_i * (prod_{l!=i})) / u_i, reduced into B's ring.
        RatPoly siB = mul_rational_poly(sigma[i], B);
        RatPoly diff = sub_rational_poly(c_work, siB);
        auto [cq, cr] = div_rem_rational_poly(diff, A);
        if (cr.size() != 0U) {
            // Not exactly divisible: the lift coefficient isn't representable —
            // signal failure so the driver reports Unimplemented honestly.
            // (Reduce numerically: keep going only if remainder is zero.)
            return std::nullopt;
        }
        c_work = cq;
    }
    return sigma;
}

}  // namespace

Result<std::vector<MPoly>> wang_multivariate_hensel(
    const MPoly& a,
    const std::vector<MPoly>& lifted_lc,
    const std::vector<IntPoly>& univar_factors,
    const std::vector<BigInt>& evaluation_point,
    const WangContext& wc,
    symbolic::CASContext& ctx) {
    (void)ctx;
    const std::size_t r = univar_factors.size();
    const std::size_t n = wc.nvars();
    if (r == 0U) {
        return fail<std::vector<MPoly>>(make_error(
            CASErrorKind::InvalidArgument, "wang hensel: no univariate factors"));
    }

    // Current lifted factors, start as univariate (in x_1) MPolys.
    std::vector<MPoly> U(r);
    for (std::size_t i = 0; i < r; ++i) {
        U[i] = mpoly_from_intpoly(univar_factors[i], 0U, n);
        // Impose the intended leading coefficient: replace the x_1^deg term's
        // integer coeff with lifted_lc[i] (an MPoly in x_2..x_n).
        if (!lifted_lc[i].is_zero()) {
            unsigned int d1 = mpoly_degree_in(U[i], 0U);
            // strip current leading term and add lc * x_1^d1
            Monomial lead(n, 0U);
            lead[0] = d1;
            U[i].terms.erase(lead);
            MPoly lcterm = mul_by_shift_power(lifted_lc[i], 0U, BigInt(0), d1, n);
            U[i] = mpoly_add(U[i], lcterm);
        }
    }

    // Precompute the RatPoly univariate images (for the Diophantine solver).
    std::vector<RatPoly> u_rat(r);
    for (std::size_t i = 0; i < r; ++i) {
        u_rat[i] = intpoly_to_ratpoly(univar_factors[i]);
    }

    // Lift each remaining variable x_{vi} (vi = 1..n-1), one at a time.
    for (std::size_t vi = 1; vi < n; ++vi) {
        const BigInt a_vi = evaluation_point[vi - 1U];
        const unsigned int target_deg = mpoly_degree_in(a, vi);

        for (unsigned int k = 1; k <= target_deg; ++k) {
            // error = a - prod(U)
            MPoly prod = mpoly_constant(BigInt(1), n);
            for (const auto& f : U) prod = mpoly_mul(prod, f);
            MPoly error = mpoly_sub(a, prod);
            if (error.is_zero()) break;

            // c_k = Taylor coeff of (x_vi - a_vi)^k in error, then project the
            // remaining higher variables (vi+1..n-1) to the point so the
            // Diophantine is univariate in x_1.
            MPoly ck = taylor_coeff(error, vi, a_vi, k);
            if (ck.is_zero()) continue;

            IntPoly ck_x1 = project_to_x1(ck, evaluation_point, wc);
            if (ck_x1.is_zero()) continue;
            RatPoly ck_rat = intpoly_to_ratpoly(ck_x1);

            auto sigma = univariate_diophantine(u_rat, ck_rat);
            if (!sigma.has_value()) {
                return fail<std::vector<MPoly>>(make_error(
                    CASErrorKind::Unimplemented,
                    "wang hensel: multivariate Diophantine not solvable in this "
                    "image (non-coprime univariate factors or non-integral lift "
                    "coefficient); leading-coefficient pre-distribution or a "
                    "different evaluation point is required (GCL §6.5)"));
            }

            // delta_i = sigma_i(x_1) * (x_vi - a_vi)^k ; require integer coeffs.
            for (std::size_t i = 0; i < r; ++i) {
                const RatPoly& s = (*sigma)[i];
                IntPoly s_int;
                s_int.coefficients().resize(s.size(), BigInt(0));
                bool integral = true;
                for (std::size_t d = 0; d < s.size(); ++d) {
                    if (!s[d].is_integer()) { integral = false; break; }
                    s_int[d] = s[d].numerator();
                }
                if (!integral) {
                    return fail<std::vector<MPoly>>(make_error(
                        CASErrorKind::Unimplemented,
                        "wang hensel: non-integral Diophantine solution; exact "
                        "integer lift not available with current LC distribution"));
                }
                normalize_integer_poly(s_int);
                if (s_int.is_zero()) continue;
                MPoly delta_x1 = mpoly_from_intpoly(s_int, 0U, n);
                MPoly delta = mul_by_shift_power(delta_x1, vi, a_vi, k, n);
                U[i] = mpoly_add(U[i], delta);
            }
        }
    }

    return ok(std::move(U));
}

}  // namespace cas::algebra
