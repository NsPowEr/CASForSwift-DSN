// limit_mrv.cpp — Public entry point and dispatcher for the Gruntz MRV
// limit algorithm (Gruntz 1996 §3.5).
//
// This file contains only:
//   - MRVCompare::operator()            (required by calculus_internal.hpp)
//   - compute_limit_mrv                 (main algorithm dispatcher)
//
// Implementation is split across:
//   limit_mrv_compare.cpp  — compare_growth (Gruntz §3.5 growth ordering)
//   limit_mrv_set.cpp      — mrv_set, rewrite_mrv (MRV set construction)
//   limit_mrv_leading.cpp  — leading_power_w, quotient/leading-power limits
//   limit_mrv_exp.cpp      — exponential-product limit fast path
//   limit_mrv_internal.hpp — shared private types and inline helpers

#include "limit_mrv_internal.hpp"

#include <iostream>

namespace cas::calculus {

bool MRVCompare::operator()(ExprPtr lhs, ExprPtr rhs) const noexcept {
    return symbolic::canonical_compare(lhs, rhs) < 0;
}

Result<ExprPtr> compute_limit_mrv(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    symbolic::CASContext& ctx) {
    if (!limit_is_infinity(point)) {
        return make_unimplemented<ExprPtr>(UnimplementedInfo{
            .module      = "calculus",
            .function    = "compute_limit_mrv",
            .input_shape = "finite limit point",
            .reason      = cas::error::reason_codes::SERIES_GENERAL,
            .suggestion  = "MRV is defined only for infinite limits",
            .ticket      = "F4-MRV"
        });
    }

    if (is_negative_infinity(point)) {
        AstArena& arena = ctx.arena();
        ExprPtr neg_var = arena.make<Unary>(UnaryOp::Neg, arena.make<Symbol>(var.name));
        auto transformed = ctx.substitute(expr, var, neg_var);
        if (transformed.is_error()) {
            return make_unimplemented<ExprPtr>(UnimplementedInfo{
                .module      = "calculus",
                .function    = "compute_limit_mrv",
                .input_shape = "MRV normalization for -infinity",
                .reason      = cas::error::reason_codes::SERIES_GENERAL,
                .suggestion  = "Check variable substitution",
                .ticket      = "F4-MRV"
            });
        }
        auto transformed_simplified = ctx.simplify(transformed.value());
        if (transformed_simplified.is_error()) {
            return make_unimplemented<ExprPtr>(UnimplementedInfo{
                .module      = "calculus",
                .function    = "compute_limit_mrv",
                .input_shape = "MRV normalization simplify for -infinity",
                .reason      = cas::error::reason_codes::SERIES_GENERAL,
                .suggestion  = "Check expression simplification",
                .ticket      = "F4-MRV"
            });
        }
        return compute_limit_mrv(
            transformed_simplified.value(),
            var,
            arena.make<Constant>(MathConstant::Infinity),
            ctx);
    }

    if (auto exp_product_limit = mrv_try_exponential_product_limit(
            expr, var, ctx.arena(), ctx);
        exp_product_limit.has_value()) {
        return *exp_product_limit;
    }

    MRVSet mrv = mrv_set(expr, var, ctx);
    if (mrv.empty()) return ok(expr);

    AstArena& arena = ctx.arena();
    // Fresh MRV variable.  Iterates make_fresh_symbol until the candidate is
    // structurally absent from the input expression — guards against capture
    // even when the user pre-populated the context with similarly-prefixed
    // symbols (Categoria 7 CLAUDE.md).
    Symbol w_var = ctx.make_fresh_symbol("mrv_w");
    while (depends_on(expr, w_var)) {
        w_var = ctx.make_fresh_symbol("mrv_w");
    }
    std::string w_name = w_var.name;
    ExprPtr w_sym = arena.make<Symbol>(w_name);

    ExprPtr replacement = arena.make<Binary>(
        BinaryOp::Div, limit_make_integer(arena, 1), w_sym);
    auto rewritten = rewrite_mrv(expr, mrv, replacement, var, ctx);
    if (rewritten.is_error()) {
        std::cerr << "rewrite_mrv error\n";
        return make_unimplemented<ExprPtr>(UnimplementedInfo{
            .module      = "calculus",
            .function    = "compute_limit_mrv",
            .input_shape = "MRV rewrite expression",
            .reason      = cas::error::reason_codes::SERIES_GENERAL,
            .suggestion  = "Check MRV substitution logic",
            .ticket      = "F4-MRV"
        });
    }

    auto simplified = ctx.simplify(rewritten.value());
    if (simplified.is_error()) {
        return make_unimplemented<ExprPtr>(UnimplementedInfo{
            .module      = "calculus",
            .function    = "compute_limit_mrv",
            .input_shape = "MRV simplified form",
            .reason      = cas::error::reason_codes::SERIES_GENERAL,
            .suggestion  = "MRV simplified form unavailable",
            .ticket      = "F4-MRV"
        });
    }

    if (auto quotient_limit = mrv_try_quotient_valuation_limit(
            simplified.value(), w_var, arena, ctx);
        quotient_limit.has_value()) {
        return *quotient_limit;
    }

    if (auto leading_limit = mrv_try_leading_power_limit(
            simplified.value(), w_var, arena, ctx);
        leading_limit.has_value()) {
        return *leading_limit;
    }

    auto combined_res = algebra::together(simplified.value(), ctx);
    ExprPtr final_rewritten =
        combined_res.is_ok() ? combined_res.value() : simplified.value();

    if (auto quotient_limit = mrv_try_quotient_valuation_limit(
            final_rewritten, w_var, arena, ctx);
        quotient_limit.has_value()) {
        return *quotient_limit;
    }

    auto direct = ctx.substitute(final_rewritten, w_var, limit_make_integer(arena, 0));
    if (direct.is_ok()) {
        auto direct_simplified = ctx.simplify(direct.value());
        if (direct_simplified.is_ok() && !depends_on(direct_simplified.value(), w_var)) {
            return direct_simplified;
        }
    }

    auto leading = mrv_leading_power_w(final_rewritten, w_var, ctx);
    if (!leading.has_value()) {
        return make_unimplemented<ExprPtr>(UnimplementedInfo{
            .module      = "calculus",
            .function    = "compute_limit_mrv",
            .input_shape = "MRV leading power extraction",
            .reason      = cas::error::reason_codes::SERIES_GENERAL,
            .suggestion  = "MRV leading power is not decidable",
            .ticket      = "F4-MRV"
        });
    }
    if (leading->power > 0) {
        return ok(limit_make_integer(arena, 0));
    }
    if (leading->power == 0) {
        auto coeff = ctx.simplify(leading->coefficient);
        if (coeff.is_ok() && !depends_on(coeff.value(), w_var)) {
            return coeff;
        }
        return make_unimplemented<ExprPtr>(UnimplementedInfo{
            .module      = "calculus",
            .function    = "compute_limit_mrv",
            .input_shape = "MRV leading coefficient at w -> 0+",
            .reason      = cas::error::reason_codes::SERIES_GENERAL,
            .suggestion  = "MRV leading coefficient is not decidable at w -> 0+",
            .ticket      = "F4-MRV"
        });
    }

    auto coeff = ctx.simplify(leading->coefficient);
    if (coeff.is_error()) {
        return make_unimplemented<ExprPtr>(UnimplementedInfo{
            .module      = "calculus",
            .function    = "compute_limit_mrv",
            .input_shape = "MRV leading coefficient simplify",
            .reason      = cas::error::reason_codes::SERIES_GENERAL,
            .suggestion  = "MRV leading coefficient simplify failed",
            .ticket      = "F4-MRV"
        });
    }
    auto sign = exact_sign(coeff.value());
    if (!sign.has_value() || *sign == 0) {
        return make_unimplemented<ExprPtr>(UnimplementedInfo{
            .module      = "calculus",
            .function    = "compute_limit_mrv",
            .input_shape = "MRV pole sign determination",
            .reason      = cas::error::reason_codes::SERIES_GENERAL,
            .suggestion  = "MRV pole sign is not exactly decidable",
            .ticket      = "F4-MRV"
        });
    }
    if (*sign > 0) {
        return ok(arena.make<Constant>(MathConstant::Infinity));
    }
    return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
}

} // namespace cas::calculus
