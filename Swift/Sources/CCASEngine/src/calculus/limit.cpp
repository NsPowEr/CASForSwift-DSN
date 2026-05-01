#include "cas/calculus.hpp"

#include "calculus_internal.hpp"
#include "cas/error.hpp"

#include <utility>

namespace cas::calculus {
namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

class LimitEngine {
public:
    explicit LimitEngine(symbolic::CASContext& context) noexcept : context_(context), arena_(context.arena()) {}

    [[nodiscard]] Result<ExprPtr> compute(ExprPtr expr, const Symbol& var, ExprPtr point, LimitDirection dir) {
        auto simplified_expr = context_.simplify(expr);
        if (simplified_expr.is_error()) {
            return simplified_expr;
        }

        auto simplified_point = context_.simplify(point);
        if (simplified_point.is_error()) {
            return simplified_point;
        }

        return compute_recursive(simplified_expr.value(), var, simplified_point.value(), dir, 0U);
    }

private:
    [[nodiscard]] Result<ExprPtr> substitute_and_simplify(ExprPtr expr, const Symbol& var, ExprPtr point) {
        auto substituted = context_.substitute(expr, var, point);
        if (substituted.is_error()) {
            return substituted;
        }
        return context_.simplify(substituted.value());
    }

    [[nodiscard]] Result<ExprPtr> compute_recursive(
        ExprPtr expr,
        const Symbol& var,
        ExprPtr point,
        LimitDirection dir,
        unsigned int depth) {
        if (depth > 6U) {
            return fail<ExprPtr>(make_error(
                CASErrorKind::Unimplemented,
                "Il limite richiede piu' iterazioni di quelle supportate"));
        }

        auto direct = substitute_and_simplify(expr, var, point);
        auto quotient = extract_quotient_view(expr, arena_);
        const bool has_quotient = quotient.has_value();

        if (!limit_is_infinity(point)) {
            auto logarithmic_limit = try_logarithmic_root_limit(expr, var, point, dir, arena_);
            if (logarithmic_limit.has_value()) {
                return logarithmic_limit.value();
            }
        }

        if (limit_is_infinity(point)) {
            auto infinite_limit = try_infinite_limit(expr, var, point, arena_);
            if (infinite_limit.is_ok() || infinite_limit.error().kind != CASErrorKind::Unimplemented) {
                return infinite_limit;
            }
        }

        if (direct.is_error() && !has_quotient) {
            return direct;
        }
        if (direct.is_ok() && (!expr_is<Binary>(direct.value()) || direct.value() == expr)) {
            if (!expr_is<Binary>(expr) || (!limit_is_zero(direct.value()) && !limit_is_infinity(direct.value()))) {
                return direct;
            }
        }

        if (limit_is_zero(point) && dir == LimitDirection::Right) {
            if (const auto* power = expr_cast<Binary>(expr);
                power != nullptr &&
                power->op == BinaryOp::Pow &&
                expr_cast<Symbol>(power->left) != nullptr &&
                expr_cast<Symbol>(power->right) != nullptr &&
                expr_ref<Symbol>(power->left).name == var.name &&
                expr_ref<Symbol>(power->right).name == var.name) {
                return ok(limit_make_integer(arena_, 1));
            }
        }

        if (!quotient.has_value()) {
            return direct;
        }

        return compute_quotient_limit(expr, var, point, dir, depth, quotient.value(), direct);
    }

    [[nodiscard]] Result<ExprPtr> compute_quotient_limit(
        ExprPtr,
        const Symbol& var,
        ExprPtr point,
        LimitDirection dir,
        unsigned int depth,
        const QuotientView& quotient,
        const Result<ExprPtr>& direct) {
        auto numerator_at_point = substitute_and_simplify(quotient.numerator, var, point);
        if (numerator_at_point.is_error()) {
            return numerator_at_point;
        }

        auto denominator_at_point = substitute_and_simplify(quotient.denominator, var, point);
        if (denominator_at_point.is_error()) {
            return denominator_at_point;
        }

        if (!limit_is_infinity(point) && limit_is_zero(denominator_at_point.value())) {
            auto pole_limit = try_polynomial_pole_limit(quotient.numerator, quotient.denominator, var, point, dir, arena_);
            if (pole_limit.has_value()) {
                return pole_limit.value();
            }
        }

        const bool zero_over_zero = limit_is_zero(numerator_at_point.value()) && limit_is_zero(denominator_at_point.value());
        const bool infinity_over_infinity =
            limit_is_infinity(numerator_at_point.value()) && limit_is_infinity(denominator_at_point.value());
        if (!zero_over_zero && !infinity_over_infinity) {
            return direct.is_ok()
                ? direct
                : fail<ExprPtr>(make_error(
                    CASErrorKind::Unimplemented,
                    "Il limite non e' una forma indeterminata supportata"));
        }

        if (zero_over_zero) {
            auto canceled = cancel_common_linear_factor(quotient.numerator, quotient.denominator, var, point, arena_);
            if (canceled.has_value()) {
                auto simplified_canceled = context_.simplify(canceled.value());
                if (simplified_canceled.is_ok()) {
                    auto canceled_limit = compute_recursive(simplified_canceled.value(), var, point, dir, depth + 1U);
                    if (canceled_limit.is_ok()) {
                        return canceled_limit;
                    }
                }
            }
        }

        return compute_lhopital_or_taylor(var, point, dir, depth, quotient);
    }

    [[nodiscard]] Result<ExprPtr> compute_lhopital_or_taylor(
        const Symbol& var,
        ExprPtr point,
        LimitDirection dir,
        unsigned int depth,
        const QuotientView& quotient) {
        auto numerator_derivative = diff(quotient.numerator, var, 1U, context_);
        if (numerator_derivative.is_error()) {
            return numerator_derivative;
        }

        auto denominator_derivative = diff(quotient.denominator, var, 1U, context_);
        if (denominator_derivative.is_error()) {
            return denominator_derivative;
        }

        ExprPtr transformed = limit_make_binary(arena_, BinaryOp::Div, numerator_derivative.value(), denominator_derivative.value());
        auto simplified_transformed = context_.simplify(transformed);
        if (simplified_transformed.is_error()) {
            return simplified_transformed;
        }

        auto lhopital = compute_recursive(simplified_transformed.value(), var, point, dir, depth + 1U);
        if (lhopital.is_ok() || limit_is_infinity(point)) {
            return lhopital;
        }

        auto taylor_numerator = taylor_series(quotient.numerator, var, point, 6U, context_);
        if (taylor_numerator.is_error()) {
            return lhopital;
        }

        auto taylor_denominator = taylor_series(quotient.denominator, var, point, 6U, context_);
        if (taylor_denominator.is_error()) {
            return lhopital;
        }

        ExprPtr taylor_quotient = limit_make_binary(
            arena_,
            BinaryOp::Div,
            taylor_numerator.value().polynomial,
            taylor_denominator.value().polynomial);
        auto simplified_quotient = context_.simplify(taylor_quotient);
        if (simplified_quotient.is_error()) {
            return lhopital;
        }

        return compute_recursive(simplified_quotient.value(), var, point, dir, depth + 1U);
    }

    symbolic::CASContext& context_;
    AstArena& arena_;
};

}  // namespace

Result<ExprPtr> limit(ExprPtr expr, const Symbol& var, ExprPtr point, LimitDirection dir, symbolic::CASContext& ctx) {
    auto value = LimitEngine(ctx).compute(expr, var, point, dir);
    if (value.is_error()) {
        return value;
    }
    return symbolic::materialize_expr(value.value(), ctx.arena());
}

}  // namespace cas::calculus
