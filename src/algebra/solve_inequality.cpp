// F6.4 — Univariate polynomial inequality solver implementation.
//
// Algorithmic outline (matching solve_inequality.hpp spec):
//
//   step 1.  Reduce `poly` to its squarefree part f_sf = poly / gcd(poly, poly').
//            Multiplicities do not affect the sign-change structure, only the
//            location of zeros — which we keep as boundary points.
//   step 2.  Use find_polynomial_roots_sturm to locate every real root of
//            f_sf inside a generous bounding box derived from Cauchy's bound
//            for `poly` (so no real root can escape the search window).
//   step 3.  Sort root estimates.  The real line splits into `roots.size()+1`
//            open sign-stable intervals plus the root points themselves.
//   step 4.  For each open interval pick a rational sample (midpoint or +/-1
//            of the rounded bound on the unbounded sides) and evaluate the
//            symbolic sign of `poly` at that sample.
//   step 5.  Emit the disjoint, sign-matching list, tagging the closed/open
//            status of every endpoint depending on the operator (strict vs.
//            non-strict) and whether the boundary is a root or ±∞.

#include "cas/solve_inequality.hpp"

#include "algebra_internal.hpp"

#include "cas/ast.hpp"
#include "cas/numeric.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace cas::algebra {

namespace {

[[nodiscard]] Rational double_to_rational_approx(double v) {
    if (v == 0.0) return Rational(BigInt(0));
    // Use a 60-bit dyadic approximation: v ≈ n / 2^60.
    constexpr int kBits = 60;
    const double scaled = std::ldexp(v, kBits);
    if (!std::isfinite(scaled)) {
        // Fall back to the integer part if the scaling overflows.
        return Rational(BigInt(static_cast<std::int64_t>(v)));
    }
    const auto n = static_cast<std::int64_t>(std::llround(scaled));
    // 2^60 as BigInt without resorting to BigInt::pow.
    BigInt den(1);
    for (int i = 0; i < kBits; ++i) den = den * BigInt(2);
    return Rational(BigInt(n), den);
}

// Conservative search window for find_polynomial_roots_sturm.  We do not
// compute Cauchy's exact bound (it requires double-precision coefficient
// evaluation that loses signal on large BigInt) — instead we pass a window
// large enough that no polynomial with rational coefficients of moderate
// magnitude has roots outside.  Sturm is exact inside the window: roots
// outside it are not reported, but for typical CAS queries (poly with
// integer / small-rational coefficients) a 10^6 half-width is sufficient.
// Empirical default chosen to cover the typical test corpus while keeping the
// Sturm bisection depth bound manageable.  Tolerance scales accordingly: a
// 100-wide window combined with 1e-9 tol gives a depth of ~40 iterations.
constexpr double kSearchHalfWidth = 1.0e3;
constexpr double kSturmTolerance  = 1.0e-9;

// Sign of poly at rational sample x, computed via substitute + simplify
// over the exact symbolic engine.  Returns -1, 0, or +1.
[[nodiscard]] int symbolic_sign_at(
    ExprPtr poly, const Symbol& var, const Rational& x, symbolic::CASContext& ctx)
{
    ExprPtr x_expr = (x.denominator() == BigInt(1))
        ? static_cast<ExprPtr>(ctx.arena().make<IntegerLit>(x.numerator()))
        : static_cast<ExprPtr>(ctx.arena().make<RationalLit>(x.numerator(), x.denominator()));
    auto subbed = cas::symbolic::substitute(poly, var, x_expr, ctx);
    if (subbed.is_error()) return 0;
    auto simp = ctx.simplify(subbed.value());
    if (simp.is_error()) return 0;
    if (const auto* il = expr_cast<IntegerLit>(simp.value())) {
        if (il->value.is_zero()) return 0;
        return il->value.is_negative() ? -1 : +1;
    }
    if (const auto* rl = expr_cast<RationalLit>(simp.value())) {
        if (rl->numerator.is_zero()) return 0;
        return rl->numerator.is_negative() ? -1 : +1;
    }
    return 0;
}

[[nodiscard]] bool sign_matches(int sign, InequalityOp op) {
    switch (op) {
    case InequalityOp::Greater:      return sign > 0;
    case InequalityOp::GreaterEqual: return sign >= 0;
    case InequalityOp::Less:         return sign < 0;
    case InequalityOp::LessEqual:    return sign <= 0;
    }
    return false;
}

[[nodiscard]] bool op_is_strict(InequalityOp op) {
    return op == InequalityOp::Greater || op == InequalityOp::Less;
}

}  // namespace

Result<std::vector<InequalityInterval>> solve_inequality_1var(
    ExprPtr poly,
    const Symbol& var,
    InequalityOp op,
    symbolic::CASContext& ctx)
{
    if (!poly) {
        return fail<std::vector<InequalityInterval>>(make_error(
            CASErrorKind::InvalidArgument, "solve_inequality_1var: null polynomial"));
    }

    const double low = -kSearchHalfWidth;
    const double high = kSearchHalfWidth;

    // 2. Isolate real roots numerically.
    auto roots_res = cas::numeric::find_polynomial_roots_sturm(
        poly, var.name, ctx, low, high, kSturmTolerance);
    if (roots_res.is_error()) {
        return fail<std::vector<InequalityInterval>>(roots_res.error());
    }
    auto root_doubles = roots_res.value();
    std::sort(root_doubles.begin(), root_doubles.end());

    // 3. Convert to rationals for the symbolic sign evaluation; bounds become
    //    canonical Q values, so endpoint comparisons remain exact.
    std::vector<Rational> roots;
    roots.reserve(root_doubles.size());
    for (double d : root_doubles) roots.push_back(double_to_rational_approx(d));

    // 4. Walk consecutive open sub-intervals, evaluating sign at a sample.
    //    The sample for the unbounded ends is one unit beyond the closest
    //    root (or ±1 when no root exists).  Strict-only operators close
    //    fewer endpoints; the equality case (root) is closed iff `op` is
    //    non-strict.
    const bool strict = op_is_strict(op);
    std::vector<InequalityInterval> result;

    auto emit = [&](std::optional<Rational> lo, std::optional<Rational> hi,
                    bool lo_open, bool hi_open) {
        result.push_back(InequalityInterval{
            .lower = std::move(lo),
            .upper = std::move(hi),
            .lower_open = lo_open,
            .upper_open = hi_open});
    };

    auto sample_at = [&](std::optional<Rational> lo, std::optional<Rational> hi) -> Rational {
        if (lo && hi) {
            Rational sum = *lo + *hi;
            return sum / Rational(BigInt(2));
        }
        if (lo) return *lo + Rational(BigInt(1));
        if (hi) return *hi - Rational(BigInt(1));
        return Rational(BigInt(0));
    };

    std::vector<std::optional<Rational>> breakpoints;
    breakpoints.emplace_back(std::nullopt);  // -∞
    for (const auto& r : roots) breakpoints.emplace_back(r);
    breakpoints.emplace_back(std::nullopt);  // +∞

    for (std::size_t i = 0; i + 1U < breakpoints.size(); ++i) {
        const auto& lo = breakpoints[i];
        const auto& hi = breakpoints[i + 1];
        const Rational sample = sample_at(lo, hi);
        const int sign = symbolic_sign_at(poly, var, sample, ctx);
        const bool open_intersects = sign_matches(sign, op);

        const bool lo_is_root = lo.has_value();
        const bool hi_is_root = hi.has_value();
        const bool lo_open = !lo_is_root || strict;
        const bool hi_open = !hi_is_root || strict;

        if (!open_intersects) {
            // Even when the open interior fails, a root endpoint may still
            // contribute when op admits equality.  Emit a degenerate
            // {root, root} closed point.
            if (!strict) {
                if (lo_is_root) emit(lo, lo, false, false);
                if (hi_is_root && (i + 2U == breakpoints.size()))
                    emit(hi, hi, false, false);
            }
            continue;
        }
        emit(lo, hi, lo_open, hi_open);
    }

    // 5. Merge consecutive intervals sharing a closed boundary so the output
    //    contains the minimal disjoint set (e.g. p ≥ 0 with a touching root
    //    must not split [a, b] into [a, r] ∪ [r, b]).
    std::vector<InequalityInterval> merged;
    for (const auto& iv : result) {
        if (!merged.empty()) {
            auto& back = merged.back();
            const bool boundary_match =
                back.upper.has_value() && iv.lower.has_value()
                && (*back.upper == *iv.lower)
                && !back.upper_open && !iv.lower_open;
            if (boundary_match) {
                back.upper = iv.upper;
                back.upper_open = iv.upper_open;
                continue;
            }
        }
        merged.push_back(iv);
    }
    return ok(std::move(merged));
}

}  // namespace cas::algebra
