// limit_mrv_leading.cpp — Leading-term extraction for the Gruntz MRV algorithm
// (Gruntz 1996 §3.5).
//
// Provides:
//   - combine_sum_terms         (file-local)
//   - leading_power_w           (file-local; wrapper: mrv_leading_power_w)
//   - try_quotient_valuation_limit (file-local; wrapper: mrv_try_quotient_valuation_limit)
//   - try_leading_power_limit   (file-local; wrapper: mrv_try_leading_power_limit)
//
// Exponential-product helpers live in limit_mrv_exp.cpp.

#include "limit_mrv_internal.hpp"

#include <optional>
#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] std::optional<LeadingPower> leading_power_w(
    ExprPtr expr, const Symbol& w_var, symbolic::CASContext& ctx);

[[nodiscard]] std::optional<LeadingPower> combine_sum_terms(
    const LeadingPower& lhs,
    const LeadingPower& rhs,
    bool is_sub,
    symbolic::CASContext& ctx) {
    if (lhs.power < rhs.power) return lhs;
    if (rhs.power < lhs.power) {
        if (!is_sub) return rhs;
        auto neg_coeff = simplify_unary_neg(rhs.coefficient, ctx);
        if (!neg_coeff.has_value()) return std::nullopt;
        return LeadingPower{.power = rhs.power, .coefficient = *neg_coeff};
    }

    ExprPtr rhs_coefficient = rhs.coefficient;
    if (is_sub) {
        auto neg_coeff = simplify_unary_neg(rhs.coefficient, ctx);
        if (!neg_coeff.has_value()) return std::nullopt;
        rhs_coefficient = *neg_coeff;
    }
    auto combined = simplify_binary(BinaryOp::Add, lhs.coefficient, rhs_coefficient, ctx);
    if (!combined.has_value()) return std::nullopt;
    if (is_exact_zero(*combined)) {
        // A cancelled leading term requires the next series term; do not guess.
        return std::nullopt;
    }
    if (!is_exact_nonzero(*combined)) {
        return std::nullopt;
    }
    return LeadingPower{.power = lhs.power, .coefficient = *combined};
}

[[nodiscard]] std::optional<LeadingPower> leading_power_w(
    ExprPtr expr,
    const Symbol& w_var,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    if (!depends_on(expr, w_var)) {
        return LeadingPower{.power = 0LL, .coefficient = expr};
    }

    if (const auto* symbol = expr_cast<Symbol>(expr)) {
        if (symbol->name == w_var.name) {
            return LeadingPower{
                .power = 1LL,
                .coefficient = limit_make_integer(arena, 1),
            };
        }
        return std::nullopt;
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op != UnaryOp::Neg) return std::nullopt;
        auto inner = leading_power_w(unary->operand, w_var, ctx);
        if (!inner.has_value()) return std::nullopt;
        auto coeff = simplify_unary_neg(inner->coefficient, ctx);
        if (!coeff.has_value()) return std::nullopt;
        return LeadingPower{.power = inner->power, .coefficient = *coeff};
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Add || binary->op == BinaryOp::Sub) {
            auto left = leading_power_w(binary->left, w_var, ctx);
            auto right = leading_power_w(binary->right, w_var, ctx);
            if (!left.has_value() || !right.has_value()) return std::nullopt;
            return combine_sum_terms(*left, *right, binary->op == BinaryOp::Sub, ctx);
        }
        if (binary->op == BinaryOp::Mul || binary->op == BinaryOp::Div) {
            auto left = leading_power_w(binary->left, w_var, ctx);
            auto right = leading_power_w(binary->right, w_var, ctx);
            if (!left.has_value() || !right.has_value()) return std::nullopt;
            auto coeff = simplify_binary(binary->op, left->coefficient, right->coefficient, ctx);
            if (!coeff.has_value()) return std::nullopt;
            const long long power = binary->op == BinaryOp::Mul
                ? left->power + right->power
                : left->power - right->power;
            return LeadingPower{.power = power, .coefficient = *coeff};
        }
        if (binary->op == BinaryOp::Pow) {
            auto exponent = integer_value(binary->right);
            if (!exponent.has_value()) return std::nullopt;
            auto base = leading_power_w(binary->left, w_var, ctx);
            if (!base.has_value()) return std::nullopt;
            auto coeff = simplify_binary(
                BinaryOp::Pow,
                base->coefficient,
                limit_make_integer(arena, *exponent),
                ctx);
            if (!coeff.has_value()) return std::nullopt;
            return LeadingPower{
                .power = base->power * (*exponent),
                .coefficient = *coeff,
            };
        }
        return std::nullopt;
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        LeadingPower result{
            .power = 0LL,
            .coefficient = limit_make_integer(arena, 1),
        };
        for (ExprPtr factor : product->factors) {
            auto term = leading_power_w(factor, w_var, ctx);
            if (!term.has_value()) return std::nullopt;
            result.power += term->power;
            auto coeff = simplify_binary(BinaryOp::Mul, result.coefficient, term->coefficient, ctx);
            if (!coeff.has_value()) return std::nullopt;
            result.coefficient = *coeff;
        }
        return result;
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        if (sum->terms.empty()) {
            return LeadingPower{
                .power = 0LL,
                .coefficient = limit_make_integer(arena, 0),
            };
        }
        auto current = leading_power_w(sum->terms.front(), w_var, ctx);
        if (!current.has_value()) return std::nullopt;
        for (std::size_t index = 1; index < sum->terms.size(); ++index) {
            auto next = leading_power_w(sum->terms[index], w_var, ctx);
            if (!next.has_value()) return std::nullopt;
            auto merged = combine_sum_terms(*current, *next, false, ctx);
            if (!merged.has_value()) {
                // A cancelled dominant term needs the next Laurent term. Factor out
                // the current valuation first, then use Taylor on the regular part.
                const long long scale_power = -current->power;
                ExprPtr scaled_expr = expr;
                if (scale_power != 0LL) {
                    ExprPtr w_expr = arena.make<Symbol>(w_var.name);
                    ExprPtr scale = scale_power == 1LL
                        ? w_expr
                        : arena.make<Binary>(
                            BinaryOp::Pow,
                            w_expr,
                            limit_make_integer(arena, scale_power));
                    auto scaled = ctx.simplify(arena.make<Binary>(BinaryOp::Mul, expr, scale));
                    if (scaled.is_error()) return std::nullopt;
                    scaled_expr = scaled.value();
                }
                auto series = taylor_series(scaled_expr, w_var, limit_make_integer(arena, 0), 6U, ctx);
                if (series.is_error()) return std::nullopt;
                if (structural_equal(series.value().polynomial, scaled_expr)) return std::nullopt;
                auto regular = leading_power_w(series.value().polynomial, w_var, ctx);
                if (!regular.has_value()) return std::nullopt;
                return LeadingPower{
                    .power = current->power + regular->power,
                    .coefficient = regular->coefficient,
                };
            }
            current = merged;
        }
        return current;
    }

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        std::vector<ExprPtr> args_coeffs;
        for (ExprPtr arg : call->args) {
            auto arg_leading = leading_power_w(arg, w_var, ctx);
            if (!arg_leading.has_value()) return std::nullopt;
            if (arg_leading->power < 0) return std::nullopt;
            if (arg_leading->power > 0) {
                args_coeffs.push_back(limit_make_integer(arena, 0));
            } else {
                args_coeffs.push_back(arg_leading->coefficient);
            }
        }
        auto res_coeff = ctx.simplify(arena.make<FuncCall>(call->func_id, std::move(args_coeffs)));
        if (res_coeff.is_error()) return std::nullopt;
        return LeadingPower{.power = 0LL, .coefficient = res_coeff.value()};
    }

    return std::nullopt;
}

[[nodiscard]] std::optional<Result<ExprPtr>> try_quotient_valuation_limit(
    ExprPtr expr,
    const Symbol& w_var,
    AstArena& arena,
    symbolic::CASContext& ctx) {
    auto quotient = extract_quotient_view(expr, arena);
    if (!quotient.has_value()) return std::nullopt;

    auto numerator = leading_power_w(quotient->numerator, w_var, ctx);
    auto denominator = leading_power_w(quotient->denominator, w_var, ctx);
    if (!numerator.has_value() || !denominator.has_value()) return std::nullopt;

    if (numerator->power > denominator->power) {
        return ok(limit_make_integer(arena, 0));
    }

    if (numerator->power == denominator->power) {
        auto ratio = ctx.simplify(arena.make<Binary>(
            BinaryOp::Div,
            numerator->coefficient,
            denominator->coefficient));
        if (ratio.is_error()) return ratio;
        if (depends_on(ratio.value(), w_var)) return std::nullopt;
        return ratio;
    }

    auto ratio = ctx.simplify(arena.make<Binary>(
        BinaryOp::Div,
        numerator->coefficient,
        denominator->coefficient));
    if (ratio.is_error()) return ratio;
    auto sign = exact_sign(ratio.value());
    if (!sign.has_value() || *sign == 0) return std::nullopt;
    if (*sign > 0) return ok(arena.make<Constant>(MathConstant::Infinity));
    return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
}

[[nodiscard]] std::optional<Result<ExprPtr>> try_leading_power_limit(
    ExprPtr expr,
    const Symbol& w_var,
    AstArena& arena,
    symbolic::CASContext& ctx) {
    auto leading = leading_power_w(expr, w_var, ctx);
    if (!leading.has_value()) return std::nullopt;

    if (leading->power > 0) {
        return ok(limit_make_integer(arena, 0));
    }
    if (leading->power == 0) {
        auto coeff = ctx.simplify(leading->coefficient);
        if (coeff.is_ok() && !depends_on(coeff.value(), w_var)) {
            return coeff;
        }
        return std::nullopt;
    }

    auto coeff = ctx.simplify(leading->coefficient);
    if (coeff.is_error()) return coeff;
    auto sign = exact_sign(coeff.value());
    if (!sign.has_value() || *sign == 0) return std::nullopt;
    if (*sign > 0) {
        return ok(arena.make<Constant>(MathConstant::Infinity));
    }
    return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
}

} // namespace

// ---------------------------------------------------------------------------
// Public wrappers — bridge anonymous-namespace implementations to
// compute_limit_mrv in limit_mrv.cpp.
// ---------------------------------------------------------------------------

std::optional<Result<ExprPtr>> mrv_try_quotient_valuation_limit(
    ExprPtr expr, const Symbol& w_var, AstArena& arena, symbolic::CASContext& ctx) {
    return try_quotient_valuation_limit(expr, w_var, arena, ctx);
}

std::optional<Result<ExprPtr>> mrv_try_leading_power_limit(
    ExprPtr expr, const Symbol& w_var, AstArena& arena, symbolic::CASContext& ctx) {
    return try_leading_power_limit(expr, w_var, arena, ctx);
}

std::optional<LeadingPower> mrv_leading_power_w(
    ExprPtr expr, const Symbol& w_var, symbolic::CASContext& ctx) {
    return leading_power_w(expr, w_var, ctx);
}

} // namespace cas::calculus
