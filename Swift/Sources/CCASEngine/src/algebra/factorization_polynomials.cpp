#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <vector>

namespace cas::algebra {

[[nodiscard]] static Result<RationalParts> normalize_rational_parts(RationalParts parts, symbolic::CASContext& ctx) {
    auto numerator = simplify_expr(parts.numerator, ctx);
    if (numerator.is_error()) {
        return fail<RationalParts>(numerator.error());
    }

    auto denominator = simplify_expr(parts.denominator, ctx);
    if (denominator.is_error()) {
        return fail<RationalParts>(denominator.error());
    }
    if (is_zero_expr(denominator.value())) {
        return fail<RationalParts>(make_error(CASErrorKind::Undefined, "Denominatore nullo in apart_num_den"));
    }

    if (is_zero_expr(numerator.value())) {
        return ok(RationalParts{
            .numerator = make_integer(ctx.arena(), 0),
            .denominator = make_integer(ctx.arena(), 1),
        });
    }

    return ok(RationalParts{
        .numerator = numerator.value(),
        .denominator = denominator.value(),
    });
}

[[nodiscard]] static Result<RationalParts> make_atomic_parts(ExprPtr expr, symbolic::CASContext& ctx) {
    auto cloned = clone_into_context(expr, ctx);
    if (cloned.is_error()) {
        return fail<RationalParts>(cloned.error());
    }
    return normalize_rational_parts(
        RationalParts{
            .numerator = cloned.value(),
            .denominator = make_integer(ctx.arena(), 1),
        },
        ctx);
}

[[nodiscard]] static Result<RationalParts> multiply_parts(const RationalParts& lhs, const RationalParts& rhs, symbolic::CASContext& ctx) {
    auto numerator = multiply_exprs(lhs.numerator, rhs.numerator, ctx);
    if (numerator.is_error()) {
        return fail<RationalParts>(numerator.error());
    }
    auto denominator = multiply_exprs(lhs.denominator, rhs.denominator, ctx);
    if (denominator.is_error()) {
        return fail<RationalParts>(denominator.error());
    }
    return normalize_rational_parts(
        RationalParts{
            .numerator = numerator.value(),
            .denominator = denominator.value(),
        },
        ctx);
}

[[nodiscard]] static Result<RationalParts> divide_parts(const RationalParts& lhs, const RationalParts& rhs, symbolic::CASContext& ctx) {
    if (is_zero_expr(rhs.numerator)) {
        return fail<RationalParts>(make_error(CASErrorKind::Undefined, "Divisione per zero in apart_num_den"));
    }

    auto numerator = multiply_exprs(lhs.numerator, rhs.denominator, ctx);
    if (numerator.is_error()) {
        return fail<RationalParts>(numerator.error());
    }
    auto denominator = multiply_exprs(lhs.denominator, rhs.numerator, ctx);
    if (denominator.is_error()) {
        return fail<RationalParts>(denominator.error());
    }
    return normalize_rational_parts(
        RationalParts{
            .numerator = numerator.value(),
            .denominator = denominator.value(),
        },
        ctx);
}

[[nodiscard]] static Result<RationalParts> add_parts(const RationalParts& lhs, const RationalParts& rhs, symbolic::CASContext& ctx) {
    auto lhs_scaled = multiply_exprs(lhs.numerator, rhs.denominator, ctx);
    if (lhs_scaled.is_error()) {
        return fail<RationalParts>(lhs_scaled.error());
    }
    auto rhs_scaled = multiply_exprs(rhs.numerator, lhs.denominator, ctx);
    if (rhs_scaled.is_error()) {
        return fail<RationalParts>(rhs_scaled.error());
    }
    auto numerator = add_exprs(lhs_scaled.value(), rhs_scaled.value(), ctx);
    if (numerator.is_error()) {
        return fail<RationalParts>(numerator.error());
    }
    auto denominator = multiply_exprs(lhs.denominator, rhs.denominator, ctx);
    if (denominator.is_error()) {
        return fail<RationalParts>(denominator.error());
    }
    return normalize_rational_parts(
        RationalParts{
            .numerator = numerator.value(),
            .denominator = denominator.value(),
        },
        ctx);
}

[[nodiscard]] static Result<RationalParts> subtract_parts(const RationalParts& lhs, const RationalParts& rhs, symbolic::CASContext& ctx) {
    auto lhs_scaled = multiply_exprs(lhs.numerator, rhs.denominator, ctx);
    if (lhs_scaled.is_error()) {
        return fail<RationalParts>(lhs_scaled.error());
    }
    auto rhs_scaled = multiply_exprs(rhs.numerator, lhs.denominator, ctx);
    if (rhs_scaled.is_error()) {
        return fail<RationalParts>(rhs_scaled.error());
    }
    auto numerator = subtract_exprs(lhs_scaled.value(), rhs_scaled.value(), ctx);
    if (numerator.is_error()) {
        return fail<RationalParts>(numerator.error());
    }
    auto denominator = multiply_exprs(lhs.denominator, rhs.denominator, ctx);
    if (denominator.is_error()) {
        return fail<RationalParts>(denominator.error());
    }
    return normalize_rational_parts(
        RationalParts{
            .numerator = numerator.value(),
            .denominator = denominator.value(),
        },
        ctx);
}

[[nodiscard]] static Result<RationalParts> pow_parts(const RationalParts& parts, const IntegerExponent& exponent, symbolic::CASContext& ctx) {
    if (exponent.negative && is_zero_expr(parts.numerator)) {
        return fail<RationalParts>(make_error(
            CASErrorKind::Undefined,
            "Una potenza negativa richiede una base razionale non nulla"));
    }

    auto numerator = pow_expr(parts.numerator, exponent.magnitude, ctx);
    if (numerator.is_error()) {
        return fail<RationalParts>(numerator.error());
    }
    auto denominator = pow_expr(parts.denominator, exponent.magnitude, ctx);
    if (denominator.is_error()) {
        return fail<RationalParts>(denominator.error());
    }

    if (!exponent.negative) {
        return normalize_rational_parts(
            RationalParts{
                .numerator = numerator.value(),
                .denominator = denominator.value(),
            },
            ctx);
    }

    return normalize_rational_parts(
        RationalParts{
            .numerator = denominator.value(),
            .denominator = numerator.value(),
        },
        ctx);
}

Result<RationalParts> split_num_den(ExprPtr expr, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<RationalParts>(make_error(CASErrorKind::InvalidArgument, "Espressione nulla in apart_num_den"));
    }
    if (contains_decimal_literal(expr)) {
        return fail<RationalParts>(make_error(
            CASErrorKind::Unimplemented,
            "I literal decimali non sono supportati in apart_num_den"));
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op == UnaryOp::Neg) {
            auto operand = split_num_den(unary->operand, ctx);
            if (operand.is_error()) {
                return fail<RationalParts>(operand.error());
            }
            auto negated = negate_expr(operand.value().numerator, ctx);
            if (negated.is_error()) {
                return fail<RationalParts>(negated.error());
            }
            return normalize_rational_parts(
                RationalParts{
                    .numerator = negated.value(),
                    .denominator = operand.value().denominator,
                },
                ctx);
        }
        return make_atomic_parts(expr, ctx);
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        auto lhs = split_num_den(binary->left, ctx);
        if (lhs.is_error()) {
            return fail<RationalParts>(lhs.error());
        }

        switch (binary->op) {
        case BinaryOp::Add: {
            auto rhs = split_num_den(binary->right, ctx);
            if (rhs.is_error()) {
                return fail<RationalParts>(rhs.error());
            }
            return add_parts(lhs.value(), rhs.value(), ctx);
        }
        case BinaryOp::Sub: {
            auto rhs = split_num_den(binary->right, ctx);
            if (rhs.is_error()) {
                return fail<RationalParts>(rhs.error());
            }
            return subtract_parts(lhs.value(), rhs.value(), ctx);
        }
        case BinaryOp::Mul: {
            auto rhs = split_num_den(binary->right, ctx);
            if (rhs.is_error()) {
                return fail<RationalParts>(rhs.error());
            }
            return multiply_parts(lhs.value(), rhs.value(), ctx);
        }
        case BinaryOp::Div: {
            auto rhs = split_num_den(binary->right, ctx);
            if (rhs.is_error()) {
                return fail<RationalParts>(rhs.error());
            }
            return divide_parts(lhs.value(), rhs.value(), ctx);
        }
        case BinaryOp::Pow: {
            auto exponent = parse_integer_exponent(binary->right);
            if (exponent.is_error()) {
                return fail<RationalParts>(exponent.error());
            }
            return pow_parts(lhs.value(), exponent.value(), ctx);
        }
        case BinaryOp::Mod:
            return fail<RationalParts>(make_error(
                CASErrorKind::Unimplemented,
                "Il modulo non e' supportato in apart_num_den"));
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        auto result = normalize_rational_parts(
            RationalParts{
                .numerator = make_integer(ctx.arena(), 0),
                .denominator = make_integer(ctx.arena(), 1),
            },
            ctx);
        if (result.is_error()) {
            return fail<RationalParts>(result.error());
        }
        for (ExprPtr term : sum->terms) {
            auto current = split_num_den(term, ctx);
            if (current.is_error()) {
                return fail<RationalParts>(current.error());
            }
            result = add_parts(result.value(), current.value(), ctx);
            if (result.is_error()) {
                return fail<RationalParts>(result.error());
            }
        }
        return result;
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        auto result = normalize_rational_parts(
            RationalParts{
                .numerator = make_integer(ctx.arena(), 1),
                .denominator = make_integer(ctx.arena(), 1),
            },
            ctx);
        if (result.is_error()) {
            return fail<RationalParts>(result.error());
        }
        for (ExprPtr factor : product->factors) {
            auto current = split_num_den(factor, ctx);
            if (current.is_error()) {
                return fail<RationalParts>(current.error());
            }
            result = multiply_parts(result.value(), current.value(), ctx);
            if (result.is_error()) {
                return fail<RationalParts>(result.error());
            }
        }
        return result;
    }

    return make_atomic_parts(expr, ctx);
}

Result<ExprPtr> together(ExprPtr expr, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "together richiede un'espressione non nulla"));
    }

    auto parts = apart_num_den(expr, ctx);
    if (parts.is_error()) {
        return fail<ExprPtr>(parts.error());
    }
    if (is_one_expr(parts.value().denominator)) {
        return ok(parts.value().numerator);
    }
    return divide_exprs(parts.value().numerator, parts.value().denominator, ctx);
}

Result<RationalParts> apart_num_den(ExprPtr expr, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<RationalParts>(make_error(
            CASErrorKind::InvalidArgument,
            "apart_num_den richiede un'espressione non nulla"));
    }
    return split_num_den(expr, ctx);
}

Result<ExprPtr> polynomial_gcd(ExprPtr p, ExprPtr q, const Symbol& var, symbolic::CASContext& ctx) {
    const bool owns_operation = !ctx.operation_active_;
    if (owns_operation) {
        ctx.operation_active_ = true;
        ctx.trace_capture_active_ = ctx.trace_enabled_;
        ctx.trace_.clear();
        ctx.ops_count_ = 0;
        ctx.operation_started_at_ = std::chrono::steady_clock::now();
    }

    auto record_trace = [&](symbolic::RuleId rule_id, ExprPtr result) {
        if (!ctx.trace_enabled_) {
            return;
        }
        ctx.trace_.push_back(symbolic::TraceStep{
            .rule_id = rule_id,
            .depth = 0U,
            .target_before = p,
            .target_after = q,
            .root_after = result,
        });
    };

    auto finalize = [&]() {
        if (owns_operation) {
            ctx.operation_active_ = false;
            ctx.trace_capture_active_ = false;
            ctx.ops_count_ = 0;
        }
    };

    auto result = [&]() -> Result<ExprPtr> {
    if (!p || !q) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::InvalidArgument,
            "polynomial_gcd richiede due espressioni polinomiali non nulle"));
    }

    auto left = parse_polynomial(p, var, ctx);
    if (left.is_error()) {
        return fail<ExprPtr>(left.error());
    }
    auto right = parse_polynomial(q, var, ctx);
    if (right.is_error()) {
        return fail<ExprPtr>(right.error());
    }

    auto left_integer = poly_to_integer_poly(left.value());
    auto right_integer = poly_to_integer_poly(right.value());
    if (left_integer.is_ok() && right_integer.is_ok()) {
        IntegerGcdResult gcd_result =
            gcd_integer_poly_with_subresultant(left_integer.value(), right_integer.value());
        IntPoly gcd_poly = std::move(gcd_result.gcd);
        symbolic::RuleId path_rule = symbolic::RuleId::PolynomialGcdPrimitiveFallback;
        switch (gcd_result.path) {
        case IntegerGcdPath::Subresultant:
            path_rule = symbolic::RuleId::PolynomialGcdSubresultant;
            break;
        case IntegerGcdPath::PrimitiveFallbackPsi:
            path_rule = symbolic::RuleId::PolynomialGcdPrimitiveFallbackPsi;
            break;
        case IntegerGcdPath::PrimitiveFallbackBeta:
            path_rule = symbolic::RuleId::PolynomialGcdPrimitiveFallbackBeta;
            break;
        case IntegerGcdPath::PrimitiveFallback:
            path_rule = symbolic::RuleId::PolynomialGcdPrimitiveFallback;
            break;
        }

        if (is_zero_integer_poly(gcd_poly)) {
            ExprPtr zero = make_integer(ctx.arena(), 0);
            record_trace(path_rule, zero);
            return ok(zero);
        }

        auto gcd_expr = integer_coefficients_to_expr(gcd_poly, var, ctx);
        if (gcd_expr.is_error()) {
            return fail<ExprPtr>(gcd_expr.error());
        }

        auto normalized = parse_polynomial(gcd_expr.value(), var, ctx);
        if (normalized.is_error()) {
            return fail<ExprPtr>(normalized.error());
        }
        auto monic = normalize_poly_monic(normalized.value(), ctx);
        if (monic.is_error()) {
            return fail<ExprPtr>(monic.error());
        }
        auto traced_result = polynomial_to_expr(monic.value(), var, ctx);
        if (traced_result.is_ok()) {
            record_trace(path_rule, traced_result.value());
        }
        return traced_result;
    }

    if (is_zero_poly(left.value()) && is_zero_poly(right.value())) {
        ExprPtr zero = make_integer(ctx.arena(), 0);
        record_trace(symbolic::RuleId::PolynomialGcdSymbolicEuclidean, zero);
        return ok(zero);
    }
    if (is_zero_poly(left.value())) {
        auto normalized = normalize_poly_monic(right.value(), ctx);
        if (normalized.is_error()) {
            return fail<ExprPtr>(normalized.error());
        }
        auto traced_result = polynomial_to_expr(normalized.value(), var, ctx);
        if (traced_result.is_ok()) {
            record_trace(symbolic::RuleId::PolynomialGcdSymbolicEuclidean, traced_result.value());
        }
        return traced_result;
    }
    if (is_zero_poly(right.value())) {
        auto normalized = normalize_poly_monic(left.value(), ctx);
        if (normalized.is_error()) {
            return fail<ExprPtr>(normalized.error());
        }
        auto traced_result = polynomial_to_expr(normalized.value(), var, ctx);
        if (traced_result.is_ok()) {
            record_trace(symbolic::RuleId::PolynomialGcdSymbolicEuclidean, traced_result.value());
        }
        return traced_result;
    }

    PolyExpr a = std::move(left.value());
    PolyExpr b = std::move(right.value());

    while (!is_zero_poly(b)) {
        auto division = divide_poly_with_remainder(a, b, ctx);
        if (division.is_error()) {
            return fail<ExprPtr>(division.error());
        }
        a = std::move(b);
        b = std::move(division.value().remainder);
    }

    auto normalized = normalize_poly_monic(a, ctx);
    if (normalized.is_error()) {
        return fail<ExprPtr>(normalized.error());
    }
    auto traced_result = polynomial_to_expr(normalized.value(), var, ctx);
    if (traced_result.is_ok()) {
        record_trace(symbolic::RuleId::PolynomialGcdSymbolicEuclidean, traced_result.value());
    }
    return traced_result;
    }();

    finalize();
    return result;
}

} // namespace cas::algebra
