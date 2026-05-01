#include "polynomial_internal.hpp"

#include "cas/symbolic.hpp"

#include <algorithm>
#include <utility>
#include <vector>

namespace cas {
namespace algebra {

namespace {

[[nodiscard]] ExprPtr build_power(AstArena& arena, ExprPtr base, std::size_t exponent) {
    if (exponent == 1U) {
        return base;
    }
    return arena.make<Binary>(BinaryOp::Pow, base, poly_make_integer(arena, static_cast<long long>(exponent)));
}

}  // namespace

Result<PolyExpr> parse_polynomial(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<PolyExpr>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "Espressione polinomiale nulla",
        });
    }
    if (poly_contains_decimal_literal(expr)) {
        return fail<PolyExpr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "I literal decimali non sono supportati nel modello polinomiale esatto",
        });
    }

    if (!poly_depends_on(expr, var.name)) {
        auto cloned = poly_clone_into_context(expr, ctx);
        if (cloned.is_error()) {
            return fail<PolyExpr>(cloned.error());
        }
        return poly_make_constant_poly(cloned.value(), ctx);
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
                .message = "Solo la negazione unaria e' supportata nel parsing polinomiale",
            });
        }
        auto operand = parse_polynomial(unary->operand, var, ctx);
        if (operand.is_error()) {
            return fail<PolyExpr>(operand.error());
        }
        return poly_negate(operand.value(), ctx);
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        switch (binary->op) {
        case BinaryOp::Add: {
            auto lhs = parse_polynomial(binary->left, var, ctx);
            if (lhs.is_error()) {
                return fail<PolyExpr>(lhs.error());
            }
            auto rhs = parse_polynomial(binary->right, var, ctx);
            if (rhs.is_error()) {
                return fail<PolyExpr>(rhs.error());
            }
            return poly_add(lhs.value(), rhs.value(), ctx);
        }
        case BinaryOp::Sub: {
            auto lhs = parse_polynomial(binary->left, var, ctx);
            if (lhs.is_error()) {
                return fail<PolyExpr>(lhs.error());
            }
            auto rhs = parse_polynomial(binary->right, var, ctx);
            if (rhs.is_error()) {
                return fail<PolyExpr>(rhs.error());
            }
            return poly_subtract(lhs.value(), rhs.value(), ctx);
        }
        case BinaryOp::Mul: {
            auto lhs = parse_polynomial(binary->left, var, ctx);
            if (lhs.is_error()) {
                return fail<PolyExpr>(lhs.error());
            }
            auto rhs = parse_polynomial(binary->right, var, ctx);
            if (rhs.is_error()) {
                return fail<PolyExpr>(rhs.error());
            }
            return poly_multiply(lhs.value(), rhs.value(), ctx);
        }
        case BinaryOp::Div: {
            if (poly_depends_on(binary->right, var.name)) {
                return fail<PolyExpr>(CASError{
                    .kind = CASErrorKind::Unimplemented,
                    .message = "I quozienti con denominatore dipendente dalla variabile non sono polinomi",
                });
            }
            auto numerator = parse_polynomial(binary->left, var, ctx);
            if (numerator.is_error()) {
                return fail<PolyExpr>(numerator.error());
            }
            auto denominator = poly_clone_into_context(binary->right, ctx);
            if (denominator.is_error()) {
                return fail<PolyExpr>(denominator.error());
            }
            return poly_divide_by_scalar(numerator.value(), denominator.value(), ctx);
        }
        case BinaryOp::Pow: {
            auto base = parse_polynomial(binary->left, var, ctx);
            if (base.is_error()) {
                return fail<PolyExpr>(base.error());
            }
            auto exponent = poly_parse_nonnegative_integer_exponent(binary->right);
            if (exponent.is_error()) {
                return fail<PolyExpr>(exponent.error());
            }
            return poly_pow(base.value(), exponent.value(), ctx);
        }
        case BinaryOp::Mod:
            return fail<PolyExpr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Il modulo non e' supportato nel modello polinomiale",
            });
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        PolyExpr result;
        for (ExprPtr term : sum->terms) {
            auto parsed_term = parse_polynomial(term, var, ctx);
            if (parsed_term.is_error()) {
                return fail<PolyExpr>(parsed_term.error());
            }
            auto next = poly_add(result, parsed_term.value(), ctx);
            if (next.is_error()) {
                return fail<PolyExpr>(next.error());
            }
            result = std::move(next.value());
        }
        return ok(std::move(result));
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        auto one = poly_make_constant_poly(poly_make_integer(ctx.arena(), 1), ctx);
        if (one.is_error()) {
            return fail<PolyExpr>(one.error());
        }
        PolyExpr result = one.value();
        for (ExprPtr factor : product->factors) {
            auto parsed_factor = parse_polynomial(factor, var, ctx);
            if (parsed_factor.is_error()) {
                return fail<PolyExpr>(parsed_factor.error());
            }
            auto next = poly_multiply(result, parsed_factor.value(), ctx);
            if (next.is_error()) {
                return fail<PolyExpr>(next.error());
            }
            result = std::move(next.value());
        }
        return ok(std::move(result));
    }

    return fail<PolyExpr>(CASError{
        .kind = CASErrorKind::Unimplemented,
        .message = "Espressione non supportata dal parser polinomiale univariato",
        .hint = "Le chiamate di funzione e i nodi dipendenti dalla variabile saranno gestiti in fasi successive.",
    });
}

Result<ExprPtr> polynomial_to_expr(const PolyExpr& poly, const Symbol& var, symbolic::CASContext& ctx) {
    if (poly.empty()) {
        return ok(poly_make_integer(ctx.arena(), 0));
    }

    ExprPtr variable_expr = ctx.arena().make<Symbol>(var.name);
    std::vector<ExprPtr> terms;
    terms.reserve(poly.size());

    for (std::size_t degree = poly.size(); degree > 0U; --degree) {
        const std::size_t index = degree - 1U;
        ExprPtr coefficient = poly[index];
        if (!coefficient || poly_is_zero_expr(coefficient)) {
            continue;
        }

        ExprPtr term = coefficient;
        if (index > 0U) {
            ExprPtr power = build_power(ctx.arena(), variable_expr, index);
            if (poly_is_one_expr(coefficient)) {
                term = power;
            } else if (poly_is_minus_one_expr(coefficient)) {
                term = ctx.arena().make<Unary>(UnaryOp::Neg, power);
            } else {
                term = ctx.arena().make<Product>(std::vector<ExprPtr>{coefficient, power});
            }
        }
        terms.push_back(term);
    }

    if (terms.empty()) {
        return ok(poly_make_integer(ctx.arena(), 0));
    }
    if (terms.size() == 1U) {
        return poly_simplify_expr(terms.front(), ctx);
    }
    return poly_simplify_expr(ctx.arena().make<Sum>(std::move(terms)), ctx);
}

}  // namespace algebra
}  // namespace cas
