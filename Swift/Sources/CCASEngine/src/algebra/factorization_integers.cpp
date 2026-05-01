#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <optional>
#include <vector>

namespace cas::algebra {

[[nodiscard]] static std::optional<std::vector<long long>> try_small_divisors(const BigInt& n) {
    if (n.is_zero()) return std::nullopt;
    const BigInt abs_n = n.abs();
    if (abs_n.bit_length() > 30U) return std::nullopt;
    const long long val = static_cast<long long>(abs_n.to_u64());
    if (val <= 0) return std::nullopt;
    std::vector<long long> divs;
    for (long long d = 1; d * d <= val; ++d) {
        if (val % d == 0) {
            divs.push_back(d);
            if (d != val / d) divs.push_back(val / d);
        }
    }
    return divs;
}

[[nodiscard]] static BigInt eval_int_poly_at_int(const IntPoly& p, long long x) {
    if (p.empty()) return BigInt(0);
    BigInt result(0);
    BigInt x_pow(1);
    const BigInt x_bi(x);
    for (const BigInt& coeff : p.coefficients()) {
        result += coeff * x_pow;
        x_pow *= x_bi;
    }
    return result;
}

[[nodiscard]] static std::optional<IntPoly> try_kronecker_quadratic_factor(
    const IntPoly& f,
    symbolic::CASContext& ctx) {
    if (f.size() < 5U) return std::nullopt;

    const BigInt v0 = eval_int_poly_at_int(f, 0);
    const BigInt v1 = eval_int_poly_at_int(f, 1);
    const BigInt v2 = eval_int_poly_at_int(f, -1);

    auto divs0 = try_small_divisors(v0);
    auto divs1 = try_small_divisors(v1);
    auto divs2 = try_small_divisors(v2);
    if (!divs0 || !divs1 || !divs2) return std::nullopt;

    for (long long d0 : *divs0) {
        for (int s0 : {1, -1}) {
            const long long sd0 = s0 * d0;
            for (long long d1 : *divs1) {
                for (int s1 : {1, -1}) {
                    const long long sd1 = s1 * d1;
                    for (long long d2 : *divs2) {
                        for (int s2 : {1, -1}) {
                            const long long sd2 = s2 * d2;
                            if (((sd1 % 2) + 2) % 2 != ((sd2 % 2) + 2) % 2) continue;
                            long long a_v = (sd1 + sd2) / 2 - sd0;
                            long long b_v = (sd1 - sd2) / 2;
                            long long c_v = sd0;
                            if (a_v == 0) continue;
                            if (a_v < 0) { a_v = -a_v; b_v = -b_v; c_v = -c_v; }
                            IntPoly g(std::vector<BigInt>{BigInt(c_v), BigInt(b_v), BigInt(a_v)});
                            g.normalize([](const BigInt& v) { return v.is_zero(); });
                            if (g.size() != 3U) continue;
                            PolyExpr f_e = integer_coefficients_to_poly(f, ctx.arena());
                            PolyExpr g_e = integer_coefficients_to_poly(g, ctx.arena());
                            auto div = divide_poly_with_remainder(f_e, g_e, ctx);
                            if (div.is_ok() && is_zero_poly(div.value().remainder)) {
                                return g;
                            }
                        }
                    }
                }
            }
        }
    }
    return std::nullopt;
}

[[nodiscard]] Result<BigInt> expr_to_integer_coefficient(ExprPtr expr) {
    if (!expr) {
        return ok(BigInt(0));
    }
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return ok(integer->value);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        if (rational->denominator == BigInt(1)) {
            return ok(rational->numerator);
        }
    }
    return fail<BigInt>(make_error(
        CASErrorKind::Unimplemented,
        "factor_over_integers supporta solo coefficienti interi esatti"));
}

[[nodiscard]] BigInt integer_content(const std::vector<BigInt>& coefficients) {
    BigInt content(0);
    for (const BigInt& coefficient : coefficients) {
        if (coefficient.is_zero()) {
            continue;
        }
        content = content.is_zero() ? coefficient.abs() : gcd(content, coefficient.abs());
    }
    return content.is_zero() ? BigInt(1) : content;
}

[[nodiscard]] BigInt integer_content(const IntPoly& coefficients) {
    return integer_content(coefficients.coefficients());
}

void divide_integer_coefficients_by_scalar(IntPoly& coefficients, const BigInt& scalar) {
    for (BigInt& coefficient : coefficients.coefficients()) {
        coefficient /= scalar;
    }
    coefficients.normalize([](const BigInt& coefficient) {
        return coefficient.is_zero();
    });
}

[[nodiscard]] PolyExpr integer_coefficients_to_poly(const std::vector<BigInt>& coefficients, AstArena& arena) {
    PolyExpr poly;
    poly.reserve(coefficients.size());
    for (const BigInt& coefficient : coefficients) {
        poly.push_back(coefficient.is_zero() ? ExprPtr{} : arena.make<IntegerLit>(coefficient));
    }
    normalize_poly(poly);
    return poly;
}

[[nodiscard]] PolyExpr integer_coefficients_to_poly(const IntPoly& coefficients, AstArena& arena) {
    return integer_coefficients_to_poly(coefficients.coefficients(), arena);
}

[[nodiscard]] Result<ExprPtr> integer_coefficients_to_expr(
    const IntPoly& coefficients,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    return polynomial_to_expr(integer_coefficients_to_poly(coefficients, ctx.arena()), var, ctx);
}

[[nodiscard]] Result<IntPoly> divide_integer_poly_by_linear_factor(
    const IntPoly& coefficients,
    const RationalRootCandidate& root,
    symbolic::CASContext& ctx) {
    PolyExpr dividend = integer_coefficients_to_poly(coefficients, ctx.arena());
    PolyExpr divisor = integer_coefficients_to_poly(
        std::vector<BigInt>{-root.numerator, root.denominator},
        ctx.arena());

    auto division = divide_poly_with_remainder(dividend, divisor, ctx);
    if (division.is_error()) {
        return fail<IntPoly>(division.error());
    }
    if (!is_zero_poly(division.value().remainder)) {
        return fail<IntPoly>(make_error(
            CASErrorKind::Unimplemented,
            "La deflazione lineare non ha prodotto resto nullo"));
    }

    return poly_to_integer_poly(division.value().quotient);
}

void append_factor_with_multiplicity(
    std::vector<PolynomialFactor>& factors,
    ExprPtr factor,
    unsigned int multiplicity) {
    if (!factors.empty() && structural_equal(factors.back().factor, factor)) {
        factors.back().multiplicity += multiplicity;
        return;
    }
    factors.push_back(PolynomialFactor{
        .factor = factor,
        .multiplicity = multiplicity,
    });
}

[[nodiscard]] IntPoly primitive_integer_poly_local(IntPoly coefficients) {
    coefficients.normalize([](const BigInt& coefficient) {
        return coefficient.is_zero();
    });
    if (coefficients.empty()) {
        return coefficients;
    }

    BigInt content = integer_content(coefficients);
    divide_integer_coefficients_by_scalar(coefficients, content);
    if (!coefficients.empty() && coefficients.leading_coeff().is_negative()) {
        for (BigInt& coefficient : coefficients.coefficients()) {
            coefficient = -coefficient;
        }
    }
    coefficients.normalize([](const BigInt& coefficient) {
        return coefficient.is_zero();
    });
    return coefficients;
}

[[nodiscard]] bool is_unit_integer_poly(const IntPoly& coefficients) {
    return coefficients.size() == 1U;
}

[[nodiscard]] IntPoly formal_derivative_integer_poly(const IntPoly& coefficients) {
    if (coefficients.size() <= 1U) {
        return IntPoly{};
    }

    IntPoly derivative;
    derivative.resize(coefficients.size() - 1U, BigInt(0));
    for (std::size_t degree = 1U; degree < coefficients.size(); ++degree) {
        derivative[degree - 1U] = coefficients[degree] * BigInt(static_cast<long long>(degree));
    }
    derivative.normalize([](const BigInt& coefficient) {
        return coefficient.is_zero();
    });
    return derivative;
}

[[nodiscard]] Result<IntPoly> exact_divide_integer_poly(
    const IntPoly& dividend,
    const IntPoly& divisor,
    symbolic::CASContext& ctx) {
    if (divisor.empty()) {
        return fail<IntPoly>(make_error(
            CASErrorKind::Undefined,
            "La divisione esatta polinomiale ha ricevuto un divisore nullo"));
    }

    PolyExpr dividend_expr = integer_coefficients_to_poly(dividend, ctx.arena());
    PolyExpr divisor_expr = integer_coefficients_to_poly(divisor, ctx.arena());
    auto division = divide_poly_with_remainder(dividend_expr, divisor_expr, ctx);
    if (division.is_error()) {
        return fail<IntPoly>(division.error());
    }
    if (!is_zero_poly(division.value().remainder)) {
        return fail<IntPoly>(make_error(
            CASErrorKind::Unimplemented,
            "La divisione polinomiale interna non e' risultata esatta"));
    }
    return poly_to_integer_poly(division.value().quotient);
}

[[nodiscard]] Result<std::vector<IntegerSquareFreeFactor>> square_free_factorize_integer_poly(
    const IntPoly& primitive,
    symbolic::CASContext& ctx) {
    std::vector<IntegerSquareFreeFactor> factors;
    if (primitive.size() <= 1U) {
        return ok(std::move(factors));
    }

    IntPoly derivative = formal_derivative_integer_poly(primitive);
    if (derivative.empty()) {
        factors.push_back(IntegerSquareFreeFactor{
            .factor = primitive,
            .multiplicity = 1U,
        });
        return ok(std::move(factors));
    }

    IntPoly g = gcd_integer_poly_with_subresultant(primitive, derivative).gcd;
    auto w_result = exact_divide_integer_poly(primitive, g, ctx);
    if (w_result.is_error()) {
        return fail<std::vector<IntegerSquareFreeFactor>>(w_result.error());
    }

    IntPoly w = primitive_integer_poly_local(std::move(w_result.value()));
    unsigned int multiplicity = 1U;
    while (!is_unit_integer_poly(w)) {
        IntPoly y = gcd_integer_poly_with_subresultant(w, g).gcd;
        auto z_result = exact_divide_integer_poly(w, y, ctx);
        if (z_result.is_error()) {
            return fail<std::vector<IntegerSquareFreeFactor>>(z_result.error());
        }

        IntPoly z = primitive_integer_poly_local(std::move(z_result.value()));
        if (!is_unit_integer_poly(z)) {
            factors.push_back(IntegerSquareFreeFactor{
                .factor = std::move(z),
                .multiplicity = multiplicity,
            });
        }

        auto next_g = exact_divide_integer_poly(g, y, ctx);
        if (next_g.is_error()) {
            return fail<std::vector<IntegerSquareFreeFactor>>(next_g.error());
        }

        w = primitive_integer_poly_local(std::move(y));
        g = primitive_integer_poly_local(std::move(next_g.value()));
        ++multiplicity;
    }

    return ok(std::move(factors));
}

[[nodiscard]] Result<void> append_integer_factor_component(
    Factorization& factorization,
    const IntPoly& component,
    unsigned int multiplicity,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    if (component.size() <= 1U) {
        return ok();
    }

    IntPoly remaining = component;
    while (remaining.size() > 1U) {
        const auto root = find_rational_root_candidate(remaining);
        if (!root.has_value()) {
            break;
        }

        auto factor_expr = integer_coefficients_to_expr(
            IntPoly(std::vector<BigInt>{-root->numerator, root->denominator}),
            var,
            ctx);
        if (factor_expr.is_error()) {
            return fail<void>(factor_expr.error());
        }
        append_factor_with_multiplicity(factorization.factors, factor_expr.value(), multiplicity);

        auto quotient = divide_integer_poly_by_linear_factor(remaining, *root, ctx);
        if (quotient.is_error()) {
            return fail<void>(quotient.error());
        }
        remaining = primitive_integer_poly_local(std::move(quotient.value()));
    }

    while (remaining.size() >= 5U) {
        auto quad = try_kronecker_quadratic_factor(remaining, ctx);
        if (!quad.has_value()) break;
        auto quad_expr = integer_coefficients_to_expr(*quad, var, ctx);
        if (quad_expr.is_error()) return fail<void>(quad_expr.error());
        append_factor_with_multiplicity(factorization.factors, quad_expr.value(), multiplicity);
        auto quotient = exact_divide_integer_poly(remaining, *quad, ctx);
        if (quotient.is_error()) break;
        remaining = primitive_integer_poly_local(std::move(quotient.value()));
    }

    if (remaining.size() > 1U) {
        auto residual = integer_coefficients_to_expr(remaining, var, ctx);
        if (residual.is_error()) {
            return fail<void>(residual.error());
        }
        append_factor_with_multiplicity(factorization.factors, residual.value(), multiplicity);
    }

    return ok();
}

Result<Factorization> factor_over_integers(ExprPtr poly, const Symbol& var, symbolic::CASContext& ctx) {
    if (!poly) {
        return fail<Factorization>(make_error(
            CASErrorKind::InvalidArgument,
            "factor_over_integers richiede un polinomio non nullo"));
    }

    auto parsed = parse_polynomial(poly, var, ctx);
    if (parsed.is_error()) {
        return fail<Factorization>(parsed.error());
    }

    auto integer_coefficients = poly_to_integer_poly(parsed.value());
    if (integer_coefficients.is_error()) {
        return fail<Factorization>(integer_coefficients.error());
    }
    if (integer_coefficients.value().empty()) {
        return fail<Factorization>(make_error(
            CASErrorKind::InvalidArgument,
            "Il polinomio nullo non ha una fattorizzazione intera canonica"));
    }

    IntPoly primitive = std::move(integer_coefficients.value());
    BigInt content = integer_content(primitive);
    divide_integer_coefficients_by_scalar(primitive, content);

    if (!primitive.empty() && primitive.leading_coeff().is_negative()) {
        content = -content;
        for (BigInt& coefficient : primitive.coefficients()) {
            coefficient = -coefficient;
        }
    }

    Factorization factorization{
        .content = ctx.arena().make<IntegerLit>(content),
        .factors = {},
    };

    if (primitive.size() <= 1U) {
        return ok(std::move(factorization));
    }

    auto square_free = square_free_factorize_integer_poly(primitive, ctx);
    if (square_free.is_error()) {
        return fail<Factorization>(square_free.error());
    }

    for (const IntegerSquareFreeFactor& component : square_free.value()) {
        auto appended = append_integer_factor_component(
            factorization,
            component.factor,
            component.multiplicity,
            var,
            ctx);
        if (appended.is_error()) {
            return fail<Factorization>(appended.error());
        }
    }

    return ok(std::move(factorization));
}

} // namespace cas::algebra
