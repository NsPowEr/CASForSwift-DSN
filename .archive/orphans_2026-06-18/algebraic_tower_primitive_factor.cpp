// algebraic_tower_primitive_factor.cpp — Q-factorization and candidate-filter helpers.
// Split from algebraic_tower_primitive_internal.hpp (ANTI-MONOLITH rule).
// Included translation units: this file only.
#include "algebraic_tower_primitive_internal.hpp"
#include "algebra_internal.hpp"

namespace cas {
namespace algebra {
namespace primitive_internal {

std::vector<RatPoly>
collect_irred_factors_over_q(const RatPoly& f, symbolic::CASContext& ctx) {
    // Clear denominators: LCM of all coefficient denominators → Z[y].
    BigInt lcm_den(1);
    for (std::size_t i = 0U; i < f.size(); ++i) {
        const BigInt& d = f[i].denominator();
        BigInt g = gcd(lcm_den, d);
        lcm_den = (lcm_den / g) * d;
    }
    PolyExpr int_poly;
    int_poly.reserve(f.size());
    for (std::size_t i = 0U; i < f.size(); ++i) {
        Rational sc = f[i] * Rational(lcm_den);
        int_poly.push_back(sc.numerator().is_zero()
            ? ExprPtr{}
            : ctx.arena().make<IntegerLit>(sc.numerator()));
    }
    normalize_poly(int_poly);
    // Use a fresh non-user-facing symbol as the polynomial variable.
    Symbol pv = ctx.make_fresh_symbol("pv");
    auto e_res = polynomial_to_expr(int_poly, pv, ctx);
    if (e_res.is_error()) return {f};
    auto fr = factor_over_integers(e_res.value(), pv, ctx);
    if (fr.is_error() || fr.value().factors.size() <= 1U) return {f};
    std::vector<RatPoly> out;
    out.reserve(fr.value().factors.size());
    for (const auto& pf : fr.value().factors) {
        auto pp = parse_polynomial(pf.factor, pv, ctx);
        if (pp.is_error()) return {f};
        auto rp = poly_to_rational_poly(pp.value());
        if (rp.is_error()) return {f};
        out.push_back(ratpoly_make_monic(std::move(rp.value())));
    }
    return out;
}

bool cand_vanishes_at_theta_expr(
    const RatPoly& cand,
    ExprPtr theta_expr,
    symbolic::CASContext& ctx) {
    if (cand.size() == 0U) return false;
    ExprPtr val = ctx.arena().make<RationalLit>(
        cand[cand.size() - 1U].numerator(),
        cand[cand.size() - 1U].denominator());
    for (std::size_t i = cand.size() - 1U; i > 0U; --i) {
        ExprPtr mul = ctx.arena().make<Binary>(BinaryOp::Mul, val, theta_expr);
        ExprPtr ci = ctx.arena().make<RationalLit>(
            cand[i - 1U].numerator(),
            cand[i - 1U].denominator());
        val = ctx.arena().make<Binary>(BinaryOp::Add, mul, ci);
    }
    auto simp = ctx.simplify(val);
    if (simp.is_error()) return false;
    ExprPtr sv = simp.value();
    if (auto* il = expr_cast<IntegerLit>(sv); il != nullptr) {
        return il->value.is_zero();
    }
    if (auto* rl = expr_cast<RationalLit>(sv); rl != nullptr) {
        return rl->numerator.is_zero();
    }
    return false;
}

std::vector<RatPoly> select_candidates_by_theta_expr(
    const std::vector<RatPoly>& all_candidates,
    ExprPtr theta_expr,
    symbolic::CASContext& ctx) {
    std::vector<RatPoly> filtered;
    filtered.reserve(all_candidates.size());
    for (const auto& cq : all_candidates) {
        if (cand_vanishes_at_theta_expr(cq, theta_expr, ctx)) filtered.push_back(cq);
    }
    // No candidate confirmed by symbolic check → fall back to original (any-valid)
    // order. Preserves termination on cases where the simplifier cannot prove
    // vanishing (e.g. deeply nested radicals).
    if (filtered.empty()) return all_candidates;
    return filtered;
}

}  // namespace primitive_internal
}  // namespace algebra
}  // namespace cas
