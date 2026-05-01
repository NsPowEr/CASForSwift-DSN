#include "polynomial_internal.hpp"

#include "cas/symbolic.hpp"

#include <algorithm>
#include <utility>

namespace cas {
namespace algebra {

void normalize_poly(PolyExpr& poly) {
    poly.normalize([](ExprPtr coefficient) {
        return !coefficient || poly_is_zero_expr(coefficient);
    });
}

PolyExpr poly_make_monomial(ExprPtr coefficient, std::size_t degree) {
    PolyExpr poly;
    poly.resize(degree + 1U, ExprPtr{});
    poly[degree] = coefficient;
    return poly;
}

Result<PolyExpr> poly_make_constant_poly(ExprPtr coefficient, symbolic::CASContext& ctx) {
    auto simplified = poly_simplify_expr(coefficient, ctx);
    if (simplified.is_error()) {
        return fail<PolyExpr>(simplified.error());
    }

    PolyExpr poly;
    if (!poly_is_zero_expr(simplified.value())) {
        poly.push_back(simplified.value());
    }
    return ok(std::move(poly));
}

Result<PolyExpr> poly_add(const PolyExpr& lhs, const PolyExpr& rhs, symbolic::CASContext& ctx) {
    PolyExpr result;
    const std::size_t size = std::max(lhs.size(), rhs.size());
    result.reserve(size);

    for (std::size_t index = 0; index < size; ++index) {
        const ExprPtr lhs_coeff = index < lhs.size() ? lhs[index] : ExprPtr{};
        const ExprPtr rhs_coeff = index < rhs.size() ? rhs[index] : ExprPtr{};
        if (!lhs_coeff && !rhs_coeff) {
            result.push_back(ExprPtr{});
            continue;
        }
        if (!lhs_coeff) {
            result.push_back(rhs_coeff);
            continue;
        }
        if (!rhs_coeff) {
            result.push_back(lhs_coeff);
            continue;
        }

        auto sum = poly_simplify_expr(ctx.arena().make<Sum>(std::vector<ExprPtr>{lhs_coeff, rhs_coeff}), ctx);
        if (sum.is_error()) {
            return fail<PolyExpr>(sum.error());
        }
        result.push_back(sum.value());
    }

    normalize_poly(result);
    return ok(std::move(result));
}

Result<PolyExpr> poly_negate(const PolyExpr& poly, symbolic::CASContext& ctx) {
    PolyExpr result;
    result.reserve(poly.size());
    for (ExprPtr coefficient : poly.coefficients()) {
        if (!coefficient) {
            result.push_back(ExprPtr{});
            continue;
        }
        auto negated = poly_simplify_expr(ctx.arena().make<Unary>(UnaryOp::Neg, coefficient), ctx);
        if (negated.is_error()) {
            return fail<PolyExpr>(negated.error());
        }
        result.push_back(negated.value());
    }
    normalize_poly(result);
    return ok(std::move(result));
}

Result<PolyExpr> poly_subtract(const PolyExpr& lhs, const PolyExpr& rhs, symbolic::CASContext& ctx) {
    auto negated_rhs = poly_negate(rhs, ctx);
    if (negated_rhs.is_error()) {
        return fail<PolyExpr>(negated_rhs.error());
    }
    return poly_add(lhs, negated_rhs.value(), ctx);
}

Result<PolyExpr> poly_multiply(const PolyExpr& lhs, const PolyExpr& rhs, symbolic::CASContext& ctx) {
    if (lhs.empty() || rhs.empty()) {
        return ok(PolyExpr{});
    }

    PolyExpr result;
    result.resize(lhs.size() + rhs.size() - 1U, ExprPtr{});

    for (std::size_t i = 0; i < lhs.size(); ++i) {
        if (!lhs[i] || poly_is_zero_expr(lhs[i])) {
            continue;
        }
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            if (!rhs[j] || poly_is_zero_expr(rhs[j])) {
                continue;
            }

            auto product = poly_simplify_expr(
                ctx.arena().make<Product>(std::vector<ExprPtr>{lhs[i], rhs[j]}),
                ctx);
            if (product.is_error()) {
                return fail<PolyExpr>(product.error());
            }

            ExprPtr& slot = result[i + j];
            if (!slot) {
                slot = product.value();
            } else {
                auto sum = poly_simplify_expr(ctx.arena().make<Sum>(std::vector<ExprPtr>{slot, product.value()}), ctx);
                if (sum.is_error()) {
                    return fail<PolyExpr>(sum.error());
                }
                slot = sum.value();
            }
        }
    }

    normalize_poly(result);
    return ok(std::move(result));
}

Result<PolyExpr> poly_divide_by_scalar(const PolyExpr& poly, ExprPtr scalar, symbolic::CASContext& ctx) {
    auto scalar_value = poly_simplify_expr(scalar, ctx);
    if (scalar_value.is_error()) {
        return fail<PolyExpr>(scalar_value.error());
    }
    if (poly_is_zero_expr(scalar_value.value())) {
        return fail<PolyExpr>(CASError{
            .kind = CASErrorKind::Undefined,
            .message = "Divisione polinomiale per coefficiente nullo",
        });
    }

    PolyExpr result;
    result.reserve(poly.size());
    for (ExprPtr coefficient : poly.coefficients()) {
        if (!coefficient || poly_is_zero_expr(coefficient)) {
            result.push_back(ExprPtr{});
            continue;
        }
        auto quotient = poly_simplify_expr(
            ctx.arena().make<Binary>(BinaryOp::Div, coefficient, scalar_value.value()),
            ctx);
        if (quotient.is_error()) {
            return fail<PolyExpr>(quotient.error());
        }
        result.push_back(quotient.value());
    }

    normalize_poly(result);
    return ok(std::move(result));
}

Result<PolyExpr> poly_pow(PolyExpr base, std::size_t exponent, symbolic::CASContext& ctx) {
    auto one = poly_make_constant_poly(poly_make_integer(ctx.arena(), 1), ctx);
    if (one.is_error()) {
        return fail<PolyExpr>(one.error());
    }
    PolyExpr result = one.value();

    while (exponent > 0U) {
        if ((exponent & 1U) != 0U) {
            auto multiplied = poly_multiply(result, base, ctx);
            if (multiplied.is_error()) {
                return fail<PolyExpr>(multiplied.error());
            }
            result = std::move(multiplied.value());
        }
        exponent >>= 1U;
        if (exponent == 0U) {
            break;
        }
        auto squared = poly_multiply(base, base, ctx);
        if (squared.is_error()) {
            return fail<PolyExpr>(squared.error());
        }
        base = std::move(squared.value());
    }

    return ok(std::move(result));
}

Result<std::size_t> poly_parse_nonnegative_integer_exponent(ExprPtr expr) {
    if (!expr) {
        return fail<std::size_t>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "Null exponent",
        });
    }

    BigInt value;
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        value = integer->value;
    } else if (const auto* rational = expr_cast<RationalLit>(expr)) {
        if (rational->denominator == BigInt(1)) {
            value = rational->numerator;
        } else {
            return fail<std::size_t>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Le potenze polinomiali supportano solo esponenti interi",
            });
        }
    } else {
        return fail<std::size_t>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Le potenze polinomiali supportano solo esponenti numerici esatti",
        });
    }

    if (value.is_negative()) {
        return fail<std::size_t>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Le potenze polinomiali supportano solo esponenti interi non negativi (ricevuto " + value.decimal() + ")",
        });
    }

    const std::string decimal = value.decimal();
    const auto parsed = parse_bounded_unsigned_decimal<std::size_t>(decimal);
    if (!parsed.has_value()) {
        return fail<std::size_t>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Le potenze polinomiali richiedono esponenti rappresentabili nel limite interno",
            .hint = "L'esponente intero non negativo supera la capacita' supportata da size_t.",
        });
    }
    return ok(parsed.value());
}

bool is_zero_poly(const PolyExpr& poly) {
    return poly.is_zero();
}

std::size_t poly_degree(const PolyExpr& poly) {
    return poly.degree();
}

ExprPtr leading_coefficient(const PolyExpr& poly) {
    if (poly.empty()) {
        return ExprPtr{};
    }
    return poly.leading_coeff();
}

Result<PolyDivisionResult> divide_poly_with_remainder(
    const PolyExpr& dividend,
    const PolyExpr& divisor,
    symbolic::CASContext& ctx) {
    if (is_zero_poly(divisor)) {
        return fail<PolyDivisionResult>(CASError{
            .kind = CASErrorKind::Undefined,
            .message = "Divisione polinomiale per divisore nullo",
        });
    }

    PolyExpr quotient;
    PolyExpr remainder = dividend;
    const std::size_t divisor_degree = divisor.degree();
    ExprPtr divisor_leading = divisor.leading_coeff();

    while (!remainder.is_zero() && remainder.degree() >= divisor_degree) {
        ExprPtr remainder_leading = remainder.leading_coeff();
        const std::size_t degree_gap = remainder.degree() - divisor_degree;
        const std::size_t previous_degree = remainder.degree();

        auto factor_coeff = poly_simplify_expr(
            ctx.arena().make<Binary>(BinaryOp::Div, remainder_leading, divisor_leading),
            ctx);
        if (factor_coeff.is_error()) {
            return fail<PolyDivisionResult>(factor_coeff.error());
        }

        PolyExpr factor = poly_make_monomial(factor_coeff.value(), degree_gap);
        auto next_quotient = poly_add(quotient, factor, ctx);
        if (next_quotient.is_error()) {
            return fail<PolyDivisionResult>(next_quotient.error());
        }
        quotient = std::move(next_quotient.value());

        auto subtractor = poly_multiply(divisor, factor, ctx);
        if (subtractor.is_error()) {
            return fail<PolyDivisionResult>(subtractor.error());
        }
        auto next_remainder = poly_subtract(remainder, subtractor.value(), ctx);
        if (next_remainder.is_error()) {
            return fail<PolyDivisionResult>(next_remainder.error());
        }

        if (!next_remainder.value().is_zero() && next_remainder.value().degree() >= previous_degree) {
            return fail<PolyDivisionResult>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "La divisione polinomiale non converge con coefficienti correnti",
                .hint = "Sarà sostituita da PRS subresultante nelle fasi successive.",
            });
        }
        remainder = std::move(next_remainder.value());
    }

    normalize_poly(quotient);
    normalize_poly(remainder);
    return ok(PolyDivisionResult{
        .quotient = std::move(quotient),
        .remainder = std::move(remainder),
    });
}

Result<PolyExpr> normalize_poly_monic(const PolyExpr& poly, symbolic::CASContext& ctx) {
    if (poly.is_zero()) {
        return ok(PolyExpr{});
    }
    return poly_divide_by_scalar(poly, poly.leading_coeff(), ctx);
}

}  // namespace algebra
}  // namespace cas
