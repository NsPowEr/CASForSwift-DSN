#include "cas/calculus.hpp"
#include "cas/algebra.hpp"
#include "integrate_definite_patterns.hpp"
#include "integrate_engine.hpp"

#include <optional>

namespace cas::calculus {

namespace {

[[nodiscard]] bool is_pos_infinity(ExprPtr expr) {
    const auto* c = expr_cast<Constant>(expr);
    return c != nullptr && c->value == MathConstant::Infinity;
}

[[nodiscard]] bool is_neg_infinity(ExprPtr expr) {
    const auto* u = expr_cast<Unary>(expr);
    return u != nullptr && u->op == UnaryOp::Neg && is_pos_infinity(u->operand);
}

[[nodiscard]] std::optional<Rational> exact_rational_from_expr(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return Rational(integer->value);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return Rational(rational->numerator, rational->denominator);
    }
    if (const auto* unary = expr_cast<Unary>(expr);
        unary != nullptr && unary->op == UnaryOp::Neg) {
        auto value = exact_rational_from_expr(unary->operand);
        if (value.has_value()) return -value.value();
    }
    return std::nullopt;
}

[[nodiscard]] bool is_between_closed(const Rational& value, const Rational& a, const Rational& b) {
    const Rational& lower = (a <= b) ? a : b;
    const Rational& upper = (a <= b) ? b : a;
    return lower <= value && value <= upper;
}

[[nodiscard]] Result<void> reject_rational_poles_in_closed_interval(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx) {
    auto lower_value = exact_rational_from_expr(lower);
    auto upper_value = exact_rational_from_expr(upper);
    if (!lower_value.has_value() || !upper_value.has_value()) {
        return ok();
    }

    auto together_expr = algebra::together(expr, ctx);
    ExprPtr rational_expr = together_expr.is_ok() ? together_expr.value() : expr;
    auto simplified = ctx.simplify(rational_expr);
    if (simplified.is_ok()) {
        rational_expr = simplified.value();
    }

    auto parts = algebra::apart_num_den(rational_expr, ctx);
    if (parts.is_error()) {
        return ok();
    }

    auto denominator = ctx.simplify(parts.value().denominator);
    if (denominator.is_error()) {
        return ok();
    }
    if (!integrate_detail::depends_on(denominator.value(), var)) {
        return ok();
    }

    auto roots = algebra::solve_polynomial(denominator.value(), var, ctx);
    if (roots.is_error()) {
        return ok();
    }

    for (ExprPtr root : roots.value()) {
        auto root_value = exact_rational_from_expr(root);
        if (!root_value.has_value()) {
            continue;
        }
        if (is_between_closed(root_value.value(), lower_value.value(), upper_value.value())) {
            return fail<void>(integrate_detail::make_error(
                CASErrorKind::Undefined,
                "Definite integral crosses a rational pole; improper/PV handling is not implemented here"));
        }
    }
    return ok();
}

[[nodiscard]] ExprPtr normalize_definite_integrand(ExprPtr expr, symbolic::CASContext& ctx) {
    auto together_expr = algebra::together(expr, ctx);
    ExprPtr normalized = together_expr.is_ok() ? together_expr.value() : expr;
    auto simplified = ctx.simplify(normalized);
    if (simplified.is_ok()) {
        normalized = simplified.value();
    }
    return normalized;
}

} // anonymous namespace
Result<ExprPtr> integrate(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    if (ctx.is_caching_enabled()) {
        auto key = symbolic::CASContext::IntegrateKey{expr, var.name};
        if (auto found = ctx.integrate_cache_.get(key)) {
            return ok(*found);
        }
    }

    auto primitive = integrate_detail::integrate_indefinite_impl(expr, var, ctx);
    if (primitive.is_error()) {
        return primitive;
    }
    auto materialized = symbolic::materialize_expr(primitive.value(), ctx.arena());
    if (materialized.is_error()) {
        return materialized;
    }
    if (ctx.is_caching_enabled()) {
        auto key = symbolic::CASContext::IntegrateKey{expr, var.name};
        ctx.integrate_cache_.put(key, materialized.value());
    }
    return materialized;
}

Result<ExprPtr> definite_integral(ExprPtr expr, const Symbol& var, ExprPtr lower, ExprPtr upper, symbolic::CASContext& ctx) {
    ExprPtr normalized_expr = normalize_definite_integrand(expr, ctx);

    // Extensible pattern table: each matcher returns nullopt to skip, value to commit.
    DefiniteContext dc{
        .integrand = expr,
        .integrand_normalized = normalized_expr,
        .var = var,
        .lower = lower,
        .upper = upper,
        .ctx = ctx,
    };
    for (DefinitePatternFn matcher : definite_patterns()) {
        auto match = matcher(dc);
        if (match.is_error()) return fail<ExprPtr>(match.error());
        if (match.value().has_value()) return ok(match.value().value());
    }

    // Generic infinite-domain fallback: only the Gaussian pattern is currently handled there;
    // anything else over (-inf, +inf) goes to Unimplemented.
    if (is_neg_infinity(lower) && is_pos_infinity(upper)) {
        return fail<ExprPtr>(integrate_detail::make_error(CASErrorKind::Unimplemented,
            "Integrazione su dominio infinito: pattern non riconosciuto."));
    }

    auto pole_check = reject_rational_poles_in_closed_interval(normalized_expr, var, lower, upper, ctx);
    if (pole_check.is_error()) {
        return fail<ExprPtr>(pole_check.error());
    }

    auto primitive = integrate(normalized_expr, var, ctx);
    if (primitive.is_error()) {
        return primitive;
    }

    auto lower_value = ctx.substitute(primitive.value(), var, lower);
    if (lower_value.is_error()) {
        return lower_value;
    }

    auto upper_value = ctx.substitute(primitive.value(), var, upper);
    if (upper_value.is_error()) {
        return upper_value;
    }

    return ctx.simplify(integrate_detail::make_sum(ctx.arena(), {
        upper_value.value(),
        integrate_detail::make_unary(ctx.arena(), UnaryOp::Neg, lower_value.value()),
    }));
}

}  // namespace cas::calculus
