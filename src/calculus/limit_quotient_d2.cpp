// F7.5.D2 — Pre-MRV quotient/sum transformations.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Gruntz_QuotientCancellation_D2.md
//
// Two pure, architecturally minimal helpers:
//
//   1. try_cancel_product_pow_inverse(expr):
//      If `expr` is Product containing both a factor `f` and a factor
//      Pow(Product(..., f, ...), -1), remove `f` and remove `f` from the
//      inverse Product. Structural (ExprPtr identity) match only.
//      Returns nullopt if no cancellation possible.
//
//   2. try_limit_sum_termwise(expr, var, point, dir, ctx):
//      If `expr` is Sum, compute lim of each addend via top-level `limit()`
//      entry. If all addends yield finite results, return Sum of results.
//      If any addend errors with non-Unimplemented (including Undefined),
//      bail to allow caller to try the original MRV path. ±∞ addends with
//      consistent sign return ±∞.
//
// Both transformations are REGOLA ZERO compliant: no closed pattern table,
// no magic constants, no depth caps beyond the recursion's natural bound
// (Sum AST is finite by construction; recursive limit() call already has
// its own adaptive depth budget).

#include "cas/calculus.hpp"
#include "calculus_internal.hpp"
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "cas/extended_real.hpp"

#include <optional>
#include <vector>

namespace cas::calculus {

namespace {
[[nodiscard]] bool is_unhandled_extended_constant(ExprPtr e) noexcept {
    const auto* c = expr_cast<Constant>(e);
    if (c == nullptr) return false;
    return c->value == MathConstant::ComplexInfinity ||
           c->value == MathConstant::Indeterminate;
}
}  // namespace

namespace {

[[nodiscard]] ExprPtr make_product_or_one(
    AstArena& arena, std::vector<ExprPtr> factors) {
    if (factors.empty()) return limit_make_integer(arena, 1);
    if (factors.size() == 1) return factors[0];
    return arena.make<Product>(std::move(factors));
}

}  // namespace

[[nodiscard]] std::optional<ExprPtr> try_cancel_product_pow_inverse(
    ExprPtr expr, symbolic::CASContext& ctx) {
    if (!expr) return std::nullopt;
    const auto* outer = expr_cast<Product>(expr);
    if (outer == nullptr) return std::nullopt;

    // Find a Pow(inner_product, IntegerLit(-1)) factor.
    std::optional<std::size_t> inv_idx;
    const Product* inner_prod = nullptr;
    for (std::size_t i = 0; i < outer->factors.size(); ++i) {
        const auto* pw = expr_cast<Binary>(outer->factors[i]);
        if (pw == nullptr || pw->op != BinaryOp::Pow) continue;
        const auto* exp_lit = expr_cast<IntegerLit>(pw->right);
        if (exp_lit == nullptr) continue;
        if (exp_lit->value != BigInt(-1)) continue;
        const auto* base_prod = expr_cast<Product>(pw->left);
        if (base_prod == nullptr) continue;
        inv_idx = i;
        inner_prod = base_prod;
        break;
    }
    if (!inv_idx.has_value()) return std::nullopt;

    // Build mutable copies of outer/inner factor lists, then cancel
    // structural matches (ExprPtr identity via canonical_equal).
    std::vector<ExprPtr> outer_keep;
    outer_keep.reserve(outer->factors.size());
    std::vector<ExprPtr> inner_keep(inner_prod->factors.begin(),
                                    inner_prod->factors.end());
    bool any_cancelled = false;

    for (std::size_t i = 0; i < outer->factors.size(); ++i) {
        if (i == *inv_idx) continue;
        ExprPtr f = outer->factors[i];
        bool cancelled = false;
        for (auto it = inner_keep.begin(); it != inner_keep.end(); ++it) {
            if (symbolic::canonical_compare(f, *it) == 0) {
                inner_keep.erase(it);
                cancelled = true;
                any_cancelled = true;
                break;
            }
        }
        if (!cancelled) outer_keep.push_back(f);
    }

    if (!any_cancelled) return std::nullopt;

    AstArena& arena = ctx.arena();
    if (inner_keep.empty()) {
        // Denominator fully cancelled: result is the surviving outer product.
        return make_product_or_one(arena, std::move(outer_keep));
    }

    // Rebuild Pow(new_inner_product, -1) and reassemble outer product.
    ExprPtr new_inner = make_product_or_one(arena, std::move(inner_keep));
    ExprPtr new_inv = arena.make<Binary>(
        BinaryOp::Pow, new_inner, limit_make_integer(arena, -1));
    outer_keep.push_back(new_inv);
    return make_product_or_one(arena, std::move(outer_keep));
}

[[nodiscard]] std::optional<Result<ExprPtr>> try_limit_sum_termwise(
    ExprPtr expr, const Symbol& var, ExprPtr point,
    LimitDirection dir, symbolic::CASContext& ctx) {
    const auto* sum = expr_cast<Sum>(expr);
    if (sum == nullptr) return std::nullopt;
    if (sum->terms.size() < 2) return std::nullopt;

    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> term_limits;
    term_limits.reserve(sum->terms.size());

    int pos_inf_count = 0;
    int neg_inf_count = 0;

    for (ExprPtr term : sum->terms) {
        // Guard against infinite recurse: if term IS the parent Sum (impossible
        // structurally but defensive), bail.
        if (term == expr) return std::nullopt;
        auto r = limit(term, var, point, dir, ctx);
        if (r.is_error()) {
            // Any error (Undefined, Unimplemented, etc.): bail and let caller
            // fall through to MRV path on the original Sum.
            return std::nullopt;
        }
        ExprPtr v = r.value();
        if (is_pos_infinity(v)) {
            ++pos_inf_count;
        } else if (is_neg_infinity(v)) {
            ++neg_inf_count;
        } else if (is_unhandled_extended_constant(v)) {
            // ±∞ mix or undecided — bail.
            return std::nullopt;
        }
        term_limits.push_back(v);
    }

    // Mixed +∞ and -∞ → ∞ - ∞ indeterminate; bail to MRV.
    if (pos_inf_count > 0 && neg_inf_count > 0) {
        return std::nullopt;
    }
    // All +∞ (or some +∞ with finite rest): result is +∞.
    if (pos_inf_count > 0) {
        return ok(arena.make<Constant>(MathConstant::Infinity));
    }
    if (neg_inf_count > 0) {
        return ok(arena.make<Unary>(
            UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
    }
    // All finite: rebuild Sum and simplify.
    ExprPtr result = arena.make<Sum>(std::move(term_limits));
    auto simp = ctx.simplify(result);
    if (simp.is_ok()) return ok(simp.value());
    return ok(result);
}

}  // namespace cas::calculus
