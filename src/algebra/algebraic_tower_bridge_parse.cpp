// algebraic_tower_bridge_parse.cpp — Polynomial parsing functions for tower bridge.
#include "algebraic_tower_bridge_internal.hpp"
#include "cas/error.hpp"
#include "cas/ast.hpp"
#include "polynomial_internal.hpp"

namespace cas {
namespace algebra {

Result<ExprPtr> clone_expr_raw(ExprPtr expr, symbolic::CASContext& ctx) {
    return symbolic::materialize_expr(expr, ctx.arena());
}

Result<ExprPtr> canonicalize_root_expr(ExprPtr expr, symbolic::CASContext& ctx) {
    RootOfExplicitDegreeGuard rootof_guard(ctx, 1U);
    return ctx.simplify(expr);
}

Result<PolyExpr> make_constant_poly_raw(ExprPtr coefficient, symbolic::CASContext& ctx) {
    auto simplified = poly_simplify_expr(coefficient, ctx);
    if (simplified.is_error()) return fail<PolyExpr>(simplified.error());
    PolyExpr poly;
    if (!poly_is_zero_expr(simplified.value())) poly.push_back(simplified.value());
    return ok(std::move(poly));
}

Result<PolyExpr> poly_pow_raw(PolyExpr base, std::size_t exponent, symbolic::CASContext& ctx) {
    auto one = make_constant_poly_raw(poly_make_integer(ctx.arena(), 1), ctx);
    if (one.is_error()) return fail<PolyExpr>(one.error());
    PolyExpr result = one.value();
    while (exponent > 0U) {
        if ((exponent & 1U) != 0U) {
            auto multiplied = poly_multiply(result, base, ctx);
            if (multiplied.is_error()) return fail<PolyExpr>(multiplied.error());
            result = std::move(multiplied.value());
        }
        exponent >>= 1U;
        if (exponent == 0U) break;
        auto squared = poly_multiply(base, base, ctx);
        if (squared.is_error()) return fail<PolyExpr>(squared.error());
        base = std::move(squared.value());
    }
    return ok(std::move(result));
}

Result<PolyExpr> parse_polynomial_raw(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<PolyExpr>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "Tower bridge: null polynomial expression",
        });
    }
    if (poly_contains_decimal_literal(expr)) {
        return fail<PolyExpr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Tower bridge: decimal literals are not supported in exact polynomial parsing",
        });
    }
    if (!poly_depends_on(expr, var.name)) {
        auto cloned = clone_expr_raw(expr, ctx);
        if (cloned.is_error()) return fail<PolyExpr>(cloned.error());
        return make_constant_poly_raw(cloned.value(), ctx);
    }
    if (const auto* symbol = expr_cast<Symbol>(expr)) {
        if (symbol->name == var.name) {
            return ok(poly_make_monomial(poly_make_integer(ctx.arena(), 1), 1U));
        }
    }
    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op != UnaryOp::Neg) {
            return fail<PolyExpr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Tower bridge: unsupported unary operator in polynomial parsing",
            });
        }
        auto operand = parse_polynomial_raw(unary->operand, var, ctx);
        if (operand.is_error()) return fail<PolyExpr>(operand.error());
        return poly_negate(operand.value(), ctx);
    }
    if (const auto* binary = expr_cast<Binary>(expr)) {
        switch (binary->op) {
            case BinaryOp::Add: {
                auto lhs = parse_polynomial_raw(binary->left, var, ctx);
                if (lhs.is_error()) return fail<PolyExpr>(lhs.error());
                auto rhs = parse_polynomial_raw(binary->right, var, ctx);
                if (rhs.is_error()) return fail<PolyExpr>(rhs.error());
                return poly_add(lhs.value(), rhs.value(), ctx);
            }
            case BinaryOp::Sub: {
                auto lhs = parse_polynomial_raw(binary->left, var, ctx);
                if (lhs.is_error()) return fail<PolyExpr>(lhs.error());
                auto rhs = parse_polynomial_raw(binary->right, var, ctx);
                if (rhs.is_error()) return fail<PolyExpr>(rhs.error());
                return poly_subtract(lhs.value(), rhs.value(), ctx);
            }
            case BinaryOp::Mul: {
                auto lhs = parse_polynomial_raw(binary->left, var, ctx);
                if (lhs.is_error()) return fail<PolyExpr>(lhs.error());
                auto rhs = parse_polynomial_raw(binary->right, var, ctx);
                if (rhs.is_error()) return fail<PolyExpr>(rhs.error());
                return poly_multiply(lhs.value(), rhs.value(), ctx);
            }
            case BinaryOp::Div: {
                if (poly_depends_on(binary->right, var.name)) {
                    return fail<PolyExpr>(CASError{
                        .kind = CASErrorKind::Unimplemented,
                        .message = "Tower bridge: denominator depending on tower variable is not a polynomial",
                    });
                }
                auto numerator = parse_polynomial_raw(binary->left, var, ctx);
                if (numerator.is_error()) return fail<PolyExpr>(numerator.error());
                auto denominator = clone_expr_raw(binary->right, ctx);
                if (denominator.is_error()) return fail<PolyExpr>(denominator.error());
                return poly_divide_by_scalar(numerator.value(), denominator.value(), ctx);
            }
            case BinaryOp::Pow: {
                auto base = parse_polynomial_raw(binary->left, var, ctx);
                if (base.is_error()) return fail<PolyExpr>(base.error());
                auto exponent = poly_parse_nonnegative_integer_exponent(binary->right);
                if (exponent.is_error()) return fail<PolyExpr>(exponent.error());
                return poly_pow_raw(base.value(), exponent.value(), ctx);
            }
            default:
                return fail<PolyExpr>(CASError{
                    .kind = CASErrorKind::Unimplemented,
                    .message = "Tower bridge: unsupported binary operator in polynomial parsing",
                });
        }
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        PolyExpr result;
        for (ExprPtr term : sum->terms) {
            auto parsed = parse_polynomial_raw(term, var, ctx);
            if (parsed.is_error()) return fail<PolyExpr>(parsed.error());
            auto next = poly_add(result, parsed.value(), ctx);
            if (next.is_error()) return fail<PolyExpr>(next.error());
            result = std::move(next.value());
        }
        return ok(std::move(result));
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        auto one = make_constant_poly_raw(poly_make_integer(ctx.arena(), 1), ctx);
        if (one.is_error()) return fail<PolyExpr>(one.error());
        PolyExpr result = one.value();
        for (ExprPtr factor : product->factors) {
            auto parsed = parse_polynomial_raw(factor, var, ctx);
            if (parsed.is_error()) return fail<PolyExpr>(parsed.error());
            auto next = poly_multiply(result, parsed.value(), ctx);
            if (next.is_error()) return fail<PolyExpr>(next.error());
            result = std::move(next.value());
        }
        return ok(std::move(result));
    }
    return fail<PolyExpr>(CASError{
        .kind = CASErrorKind::Unimplemented,
        .message = "Tower bridge: expression not supported by raw polynomial parser",
    });
}

}  // namespace algebra
}  // namespace cas
