// F3.2 — Wang leading-coefficient determination (GCL §6.6, Wang 1978 §3).
//
// Given the multivariate poly `a` (primitive & squarefree in main var x_1) and
// the factors of its univariate image at the evaluation point, assign to each
// univariate factor its intended MULTIVARIATE leading coefficient (an MPoly in
// x_2..x_n).  The product of the assigned LCs (evaluated) must reproduce the LC
// of the univariate image, and each must divide lc_{x_1}(a).
//
// Algorithm reference (closing F3.2-WANG-LC-CORRECTION):
//   - Geddes/Czapor/Labahn, "Algorithms for Computer Algebra", §6.6 Alg 6.4.
//   - Paul S. Wang, "An improved multivariate polynomial factoring algorithm",
//     Math. of Computation 32 (1978), pp. 1215–1231, §3.
//   - SymPy reference impl `sympy.polys.factortools.dmp_zz_wang_lead_coeffs`
//     and `dmp_zz_wang_non_divisors`
//     (https://github.com/sympy/sympy/blob/master/sympy/polys/factortools.py).
//     Used here as algorithmic citation only — no code copied; all data
//     structures, types and control flow are native to this project.
//
// Pipeline implemented here:
//   1.  factor lc_x(a) over Z[x_2..x_n] via recursive `factor_multivariate`,
//       collecting irreducible factors T = [t_1..t_m] with multiplicity.
//   2.  evaluate each t_i at the point → integer E[i].
//   3.  non-divisors check (Wang 1978 §3 condition): each |E[i]| must contain
//       a prime not shared with cs·ct and not absorbed by previous E[j], else
//       the assignment is structurally ambiguous → request a different point.
//   4.  for each univariate factor h_j, compute c_j ∈ Z[x_2..x_n] by reducing
//       lc(h_j)·cs against the E[i]'s in reverse and accumulating t_i to the
//       multiplicity observed.  cs = content of the univariate image.
//   5.  surjectivity check (J array): every lc-factor must have been used at
//       least once — otherwise an extraneous factor appeared at this point
//       and we surface Unimplemented for the driver to retry.
//   6.  integer-content correction: rescale each c_j ↔ h_j by the residual
//       integer ratio between lc(h_j) and eval(c_j); if cs ≠ 1 remains,
//       multiply every (c_j, h_j) by cs and rescale the original polynomial
//       by cs^(r-1) via `overall_constant` (the driver scales lift_target).
//
// Every numerical decision is later certified by EXACT reconstruction in the
// Wang driver — an incorrect distribution can never yield a silently-wrong
// factorization; failures surface as explicit Unimplemented diagnostics.

#include "factor_multivariate_internal.hpp"

#include <algorithm>

namespace cas::algebra {

namespace {

[[nodiscard]] BigInt intpoly_leading(const IntPoly& f) {
    return f.is_zero() ? BigInt(0) : f.leading_coeff();
}

// Evaluate an MPoly (in remaining vars; main-var exponent assumed 0) at the
// evaluation point a_2..a_n → a single BigInt.
[[nodiscard]] BigInt eval_mpoly_at_point(
    const MPoly& p,
    const std::vector<BigInt>& point,
    const WangContext& wc) {
    MPoly cur = p;
    for (std::size_t vi = 1; vi < wc.nvars(); ++vi) {
        cur = mpoly_eval_var(cur, vi, point[vi - 1U]);
    }
    if (cur.is_zero()) {
        return BigInt(0);
    }
    // Now cur is a constant (all exponents 0).
    return cur.terms.begin()->second;
}

// Exact integer power of an MPoly.  Used to lift t_i^k where k is the
// divisibility multiplicity observed in step (4).
[[nodiscard]] MPoly mpoly_pow_nonneg(const MPoly& p, unsigned int k,
                                     std::size_t nvars) {
    MPoly r = mpoly_constant(BigInt(1), nvars);
    if (k == 0U) return r;
    MPoly base = p;
    unsigned int e = k;
    while (e > 0U) {
        if (e & 1U) r = mpoly_mul(r, base);
        e >>= 1U;
        if (e > 0U) base = mpoly_mul(base, base);
    }
    return r;
}

// Wang non-divisors condition (Wang 1978 §3 / GCL §6.6, eq. 6.27):
// returns true iff each |E[i]| holds a prime that, after greedy reduction
// against (cs·ct) and previous entries, leaves a non-unit residue.  False
// means the point cannot distinguish the lc-factor images and the driver
// must try a different evaluation point.
//
// Closely mirrors `dmp_zz_wang_non_divisors` (SymPy) which is a faithful
// transcription of Wang's published predicate.
[[nodiscard]] bool wang_non_divisors_ok(
    const std::vector<BigInt>& E_imgs,
    const BigInt& cs,
    const BigInt& ct) {
    std::vector<BigInt> running;
    running.push_back((cs * ct).abs());
    for (const BigInt& e : E_imgs) {
        BigInt q = e.abs();
        if (q.is_zero()) return false;
        for (auto it = running.rbegin(); it != running.rend(); ++it) {
            BigInt r = *it;
            while (!(r == BigInt(1))) {
                BigInt g = gcd(r, q);
                if (g == BigInt(1)) break;
                r = g;
                q = q / r;
            }
            if (q == BigInt(1)) return false;
        }
        running.push_back(q);
    }
    return true;
}

// Distribute the multiplicities of one lc-factor `e` (integer image) over the
// running per-factor divisors `d_j = lc(h_j)*cs` modified in-place.
// Returns the vector k_j = number of times e | d_j (greedy from j = 0..r-1
// is fine because we will sum the multiplicities via reverse pass anyway).
[[nodiscard]] std::vector<unsigned int> peel_multiplicities(
    const BigInt& e, std::vector<BigInt>& d) {
    std::vector<unsigned int> ks(d.size(), 0U);
    if (e.abs() == BigInt(1)) return ks;
    for (std::size_t j = 0; j < d.size(); ++j) {
        while (!d[j].is_zero() && (d[j] % e).is_zero()) {
            d[j] = d[j] / e;
            ++ks[j];
        }
    }
    return ks;
}

}  // namespace

Result<LcDistribution> wang_distribute_leading_coeff(
    const MPoly& a,
    const std::vector<IntPoly>& univar_factors,
    const std::vector<BigInt>& evaluation_point,
    const WangContext& wc,
    symbolic::CASContext& ctx) {
    const std::size_t r = univar_factors.size();
    const std::size_t n = wc.nvars();
    LcDistribution dist;
    dist.lc.assign(r, mpoly_zero());
    dist.adjusted = univar_factors;

    // Leading coefficient of a w.r.t. main variable (MPoly in x_2..x_n).
    const MPoly lc_a = mpoly_leading_coeff_in(a, 0U);

    // ----- Fast path: lc_a constant -----
    const bool lc_constant = lc_a.terms.size() == 1U &&
        lc_a.terms.begin()->first == Monomial(n, 0U);
    if (lc_constant) {
        const BigInt lc_const = lc_a.terms.begin()->second;
        BigInt prod(1);
        for (std::size_t i = 0; i < r; ++i) {
            dist.lc[i] = mpoly_constant(intpoly_leading(univar_factors[i]), n);
            prod *= intpoly_leading(univar_factors[i]);
        }
        if (prod.is_zero()) {
            return fail<LcDistribution>(make_error(
                CASErrorKind::InvalidArgument,
                "wang LC: zero leading coefficient in univariate image"));
        }
        if ((lc_const % prod).is_zero()) {
            dist.overall_constant = lc_const / prod;
            return ok(std::move(dist));
        }
        dist.lc.assign(r, mpoly_zero());  // reset, general path below
    }

    // ----- Single factor with non-constant lc_a -----
    if (r == 1U) {
        dist.lc[0] = lc_a;
        const BigInt lc_eval = eval_mpoly_at_point(lc_a, evaluation_point, wc);
        const BigInt u_lc = intpoly_leading(univar_factors[0]);
        if (u_lc.is_zero() || lc_eval.is_zero()) {
            return fail<LcDistribution>(make_error(
                CASErrorKind::InvalidArgument,
                "wang LC: zero leading coefficient encountered"));
        }
        if (!(lc_eval % u_lc).is_zero()) {
            return fail<LcDistribution>(make_error(
                CASErrorKind::Unimplemented,
                "wang LC: single-factor LC mismatch under evaluation"));
        }
        dist.overall_constant = lc_eval / u_lc;
        return ok(std::move(dist));
    }

    // ----- General path: factor lc_a, distribute by divisibility (Alg 6.4) -----
    auto lc_expr = mpoly_to_expr(lc_a, wc, ctx);
    if (lc_expr.is_error()) return fail<LcDistribution>(lc_expr.error());
    auto lc_fac = factor_multivariate(lc_expr.value(), ctx);
    if (lc_fac.is_error()) {
        return fail<LcDistribution>(make_error(
            CASErrorKind::Unimplemented,
            "wang LC: could not factor the multivariate leading coefficient: " +
                lc_fac.error().message));
    }

    // ct = integer constant in the factorization of lc_a (Wang's "ct").
    BigInt ct(1);
    if (const auto* il = expr_cast<IntegerLit>(lc_fac.value().content)) {
        ct = il->value;
    }

    // Collect non-constant irreducible lc-factors WITHOUT expanding multiplicity.
    // Multiplicity is recovered later by repeated divisibility tests
    // (Wang/EEZ: peel `while d % e == 0`).  Constant lc-factors fold into ct.
    struct LcFactor { MPoly poly; BigInt eval; };
    std::vector<LcFactor> T;
    for (const auto& pf : lc_fac.value().factors) {
        auto mp = mpoly_from_expr(pf.factor, wc, ctx);
        if (mp.is_error()) return fail<LcDistribution>(mp.error());
        const bool is_constant = mp.value().terms.size() == 1U &&
            mp.value().terms.begin()->first == Monomial(n, 0U);
        if (is_constant) {
            BigInt v = mp.value().terms.begin()->second;
            for (unsigned int m = 0; m < pf.multiplicity; ++m) ct *= v;
            continue;
        }
        BigInt ev = eval_mpoly_at_point(mp.value(), evaluation_point, wc);
        if (ev.is_zero()) {
            return fail<LcDistribution>(make_error(
                CASErrorKind::Unimplemented,
                "wang LC: leading-coefficient factor vanishes at evaluation "
                "point; a different good point is required"));
        }
        // (pf.multiplicity carried implicitly — peeling reaches up to its
        // value across the univariate factors since Σ k_j over j equals the
        // total multiplicity of t_i in lc_a.)
        (void)pf.multiplicity;
        T.push_back({mp.value(), ev});
    }

    // cs = integer content of the univariate image (= lc_a_eval / Π lc(h_j))
    // up to the residual; computed exactly from current data.
    const BigInt lc_a_eval = eval_mpoly_at_point(lc_a, evaluation_point, wc);
    BigInt prod_hlc(1);
    for (const auto& h : univar_factors) prod_hlc *= intpoly_leading(h);
    if (prod_hlc.is_zero()) {
        return fail<LcDistribution>(make_error(
            CASErrorKind::InvalidArgument,
            "wang LC: zero leading coefficient in univariate image"));
    }
    if (!(lc_a_eval % prod_hlc).is_zero()) {
        return fail<LcDistribution>(make_error(
            CASErrorKind::Unimplemented,
            "wang LC: lc_a image not divisible by product of univariate LCs"));
    }
    BigInt cs = lc_a_eval / prod_hlc;
    if (cs.is_negative()) cs = -cs;  // sign absorbed by overall_constant below

    // Non-divisors check (Wang 1978 §3 cond.): without it, the distribution
    // cannot be unique on this point.  Surface Unimplemented so the driver
    // may re-run with a different point in a future enhancement.
    std::vector<BigInt> E_imgs;
    E_imgs.reserve(T.size());
    for (const auto& t : T) E_imgs.push_back(t.eval);
    if (!wang_non_divisors_ok(E_imgs, cs, ct)) {
        return fail<LcDistribution>(make_error(
            CASErrorKind::Unimplemented,
            "wang LC: non-divisors condition fails at this evaluation point "
            "(GCL §6.6) — algorithm needs a different evaluation point"));
    }

    // Per-factor working divisor d_j = lc(h_j) * cs.
    std::vector<BigInt> d(r);
    for (std::size_t j = 0; j < r; ++j) d[j] = intpoly_leading(univar_factors[j]) * cs;

    // Per-factor accumulated multivariate LC.
    std::vector<MPoly> C(r, mpoly_constant(BigInt(1), n));
    // Surjectivity tracker J[i]: was t_i used by any factor?
    std::vector<unsigned char> J(T.size(), 0);

    // SymPy/GCL pass: iterate lc-factors in REVERSE order, for each compute
    // multiplicity it divides each d_j (positive guarantee from non-divisors
    // condition), multiply C[j] by t_i^k_j.
    for (std::size_t ii = T.size(); ii-- > 0;) {
        const BigInt& e = T[ii].eval;
        auto ks = peel_multiplicities(e, d);
        for (std::size_t j = 0; j < r; ++j) {
            if (ks[j] == 0U) continue;
            J[ii] = 1;
            C[j] = mpoly_mul(C[j], mpoly_pow_nonneg(T[ii].poly, ks[j], n));
        }
    }

    // Surjectivity check.  If some t_i was never used, an extraneous factor
    // appeared in the univariate factorization at this point (Wang's
    // ExtraneousFactors signal) — the driver should retry with another point.
    for (std::size_t i = 0; i < T.size(); ++i) {
        if (J[i] == 0) {
            return fail<LcDistribution>(make_error(
                CASErrorKind::Unimplemented,
                "wang LC: extraneous factor / unused lc-factor at this "
                "evaluation point; retry with a different point"));
        }
    }

    // Integer-content correction.  Compare eval(C[j]) against lc(h_j) and
    // rebalance.  See GCL Alg 6.4 step 5 / SymPy `dmp_zz_wang_lead_coeffs`.
    std::vector<IntPoly> H = univar_factors;
    for (std::size_t j = 0; j < r; ++j) {
        const BigInt d_eval = eval_mpoly_at_point(C[j], evaluation_point, wc);
        const BigInt lc_h = intpoly_leading(H[j]);
        if (d_eval.is_zero()) {
            return fail<LcDistribution>(make_error(
                CASErrorKind::Unimplemented,
                "wang LC: accumulated lc evaluates to zero at the point"));
        }
        BigInt cc;
        if (cs == BigInt(1)) {
            if (!(lc_h % d_eval).is_zero()) {
                return fail<LcDistribution>(make_error(
                    CASErrorKind::Unimplemented,
                    "wang LC: residual integer correction not exact (cs=1)"));
            }
            cc = lc_h / d_eval;
        } else {
            BigInt g = gcd(lc_h.abs(), d_eval.abs());
            if (g.is_zero()) g = BigInt(1);
            BigInt d_quot = d_eval / g;  // sign on d_quot preserves d's sign
            cc = lc_h / g;
            // Multiply h_j by d_quot, divide cs by d_quot (must be exact).
            if (d_quot.is_zero() || !(cs % d_quot).is_zero()) {
                return fail<LcDistribution>(make_error(
                    CASErrorKind::Unimplemented,
                    "wang LC: cs not divisible by content-correction quotient"));
            }
            multiply_integer_coefficients_by_scalar(H[j], d_quot);
            normalize_integer_poly(H[j]);
            cs = cs / d_quot;
        }
        C[j] = mpoly_scale(C[j], cc);
    }

    // Residual cs distribution: scale each (C_j, H_j) by cs and the lifted
    // input by cs^(r-1) (Wang 1978 §3 final step).  Driver applies that via
    // `overall_constant` on lift_target.
    BigInt overall(1);
    if (!(cs == BigInt(1))) {
        for (std::size_t j = 0; j < r; ++j) {
            C[j] = mpoly_scale(C[j], cs);
            multiply_integer_coefficients_by_scalar(H[j], cs);
            normalize_integer_poly(H[j]);
        }
        // overall = cs^(r-1)
        BigInt acc(1);
        for (std::size_t k = 0; k + 1 < r; ++k) acc *= cs;
        overall = acc;
    }

    // Sign correction: lc_a_eval may be negative; ensure the SIGN of
    // Π eval(C_j) matches lc_a_eval so the lift sees a consistent leading
    // coefficient.  An odd sign mismatch is absorbed into one C_j (and the
    // corresponding H_j) so the certifier reconstructs `a` exactly.
    {
        BigInt prod_eval(1);
        for (std::size_t j = 0; j < r; ++j) {
            prod_eval *= eval_mpoly_at_point(C[j], evaluation_point, wc);
        }
        if (prod_eval.is_zero()) {
            return fail<LcDistribution>(make_error(
                CASErrorKind::Unimplemented,
                "wang LC: zero product of assigned LCs at evaluation"));
        }
        // overall_constant is multiplied into lift_target by the driver, so
        // we need: overall * Π eval(C_j)  ≡  lc_a_eval · (some unit).  When
        // the ratio is not ±1 (e.g. residual integer), surface Unimplemented.
        BigInt lhs = overall * prod_eval;
        if (lhs.is_zero() || !(lc_a_eval % lhs).is_zero()) {
            return fail<LcDistribution>(make_error(
                CASErrorKind::Unimplemented,
                "wang LC: residual leading-coefficient constant not integral"));
        }
        BigInt residual = lc_a_eval / lhs;
        overall = overall * residual;
    }

    dist.lc = std::move(C);
    dist.adjusted = std::move(H);
    dist.overall_constant = overall;
    return ok(std::move(dist));
}

}  // namespace cas::algebra
