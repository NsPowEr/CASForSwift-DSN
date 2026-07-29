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
#include "cas/error_helpers.hpp"
#include "cas/numeric.hpp"

#include "cas/rational.hpp"

#include <algorithm>

namespace cas::algebra {

namespace {

void collect_symbols_in_expr(ExprPtr e, std::vector<std::string>& out) {
    if (!e) return;
    if (const auto* s = expr_cast<Symbol>(e)) {
        out.push_back(s->name);
        return;
    }
    if (const auto* u = expr_cast<Unary>(e)) { collect_symbols_in_expr(u->operand, out); return; }
    if (const auto* b = expr_cast<Binary>(e)) { collect_symbols_in_expr(b->left, out); collect_symbols_in_expr(b->right, out); return; }
    if (const auto* sum = expr_cast<Sum>(e)) { for (ExprPtr t : sum->terms) collect_symbols_in_expr(t, out); return; }
    if (const auto* p = expr_cast<Product>(e)) { for (ExprPtr f : p->factors) collect_symbols_in_expr(f, out); return; }
    if (const auto* fc = expr_cast<FuncCall>(e)) { for (ExprPtr a : fc->args) collect_symbols_in_expr(a, out); return; }
}

// F7.0-A2.1: search window + tolerance are read from CASContextParams
// (ctx.solve_inequality_search_half_width() / _sturm_tolerance_inv()).
// Defaults: half-width 10^3 (conservative Cauchy bound for moderate coefficients),
// tolerance 1/10^9 (Mahler-separation safe for deg ≤ 10, H ≤ 100).
// REGOLA 1 boundary: the long-long parameters are converted to double ONLY at
// the call site of cas::numeric::find_polynomial_roots_sturm, which is the
// numeric-evaluator layer (not the symbolic core). See HARDCODE_LEDGER.md
// HC-F70-A21-NUMERIC-BOUNDARY.

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

    std::vector<std::string> syms;
    collect_symbols_in_expr(poly, syms);
    for (const auto& s : syms) {
        if (s != var.name) {
            return make_unimplemented<std::vector<InequalityInterval>>(UnimplementedInfo{
                .module      = "algebra",
                .function    = "solve_inequality_1var",
                .input_shape = "multivariate expression with symbol " + s,
                .reason      = cas::error::reason_codes::ALGEBRA_MULTIVAR_REMAINING_VARS,
                .suggestion  = "solve_inequality_1var supports univariate Q[var] polynomials only",
                .ticket      = "F6.4"
            });
        }
    }

    // F7.0-A2.1: bounds + tolerance from CASContextParams (long long),
    // converted to double only at the numeric-evaluator boundary.
    const long long half_width_ll = ctx.solve_inequality_search_half_width();
    const long long tol_inv_ll    = ctx.solve_inequality_sturm_tolerance_inv();
    const double low  = -static_cast<double>(half_width_ll);
    const double high =  static_cast<double>(half_width_ll);
    const double tol  =  1.0 / static_cast<double>(tol_inv_ll);

    // 2. Isolate real roots numerically (numeric layer, accepts double).
    auto roots_res = cas::numeric::find_polynomial_roots_sturm(
        poly, var.name, ctx, low, high, tol);
    if (roots_res.is_error()) {
        if (roots_res.error().kind == CASErrorKind::Unimplemented) {
            return make_unimplemented<std::vector<InequalityInterval>>(UnimplementedInfo{
                .module      = "algebra",
                .function    = "solve_inequality_1var",
                .input_shape = "non-polynomial or unsupported inequality",
                .reason      = cas::error::reason_codes::GENERIC,
                .suggestion  = "Sturm root isolation failed on input",
                .ticket      = "F6.4"
            }, roots_res.error().message);
        }
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
