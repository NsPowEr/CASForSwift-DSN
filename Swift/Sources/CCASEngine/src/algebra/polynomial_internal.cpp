#include "polynomial_internal.hpp"

#include "cas/numtheory.hpp"
#include "cas/symbolic.hpp"

#include <algorithm>
#include <limits>
#include <type_traits>
#include <utility>

namespace cas {
namespace algebra {
namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message, std::optional<std::string> hint = std::nullopt) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::move(hint),
    };
}

[[nodiscard]] PolyExpr make_zero_poly() {
    return PolyExpr{};
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

void normalize_integer_poly(IntPoly& coefficients) {
    coefficients.normalize([](const BigInt& coefficient) {
        return coefficient.is_zero();
    });
}

[[nodiscard]] BigInt integer_content(const IntPoly& coefficients) {
    BigInt content(0);
    for (const BigInt& coefficient : coefficients.coefficients()) {
        if (coefficient.is_zero()) {
            continue;
        }
        content = content.is_zero() ? coefficient.abs() : gcd(content, coefficient.abs());
    }
    return content.is_zero() ? BigInt(1) : content;
}

void divide_integer_coefficients_by_scalar(IntPoly& coefficients, const BigInt& scalar) {
    for (BigInt& coefficient : coefficients.coefficients()) {
        coefficient /= scalar;
    }
    normalize_integer_poly(coefficients);
}

[[nodiscard]] std::size_t integer_poly_degree(const IntPoly& coefficients) {
    return coefficients.degree();
}

[[nodiscard]] const BigInt& integer_poly_leading(const IntPoly& coefficients) {
    return coefficients.leading_coeff();
}

void multiply_integer_coefficients_by_scalar(IntPoly& coefficients, const BigInt& scalar) {
    if (scalar == BigInt(1)) {
        return;
    }
    if (scalar.is_zero()) {
        coefficients = IntPoly{};
        return;
    }
    for (BigInt& coefficient : coefficients.coefficients()) {
        coefficient *= scalar;
    }
    normalize_integer_poly(coefficients);
}

[[nodiscard]] IntPoly primitive_integer_poly(IntPoly coefficients) {
    normalize_integer_poly(coefficients);
    if (coefficients.empty()) {
        return coefficients;
    }

    const BigInt content = integer_content(coefficients);
    divide_integer_coefficients_by_scalar(coefficients, content);
    if (!coefficients.empty() && coefficients.leading_coeff().is_negative()) {
        for (BigInt& coefficient : coefficients.coefficients()) {
            coefficient = -coefficient;
        }
    }
    normalize_integer_poly(coefficients);
    return coefficients;
}

[[nodiscard]] BigInt bigint_pow_nonnegative(BigInt base, std::size_t exponent) {
    BigInt result(1);
    while (exponent > 0U) {
        if ((exponent & 1U) != 0U) {
            result *= base;
        }
        exponent >>= 1U;
        if (exponent > 0U) {
            base *= base;
        }
    }
    return result;
}

[[nodiscard]] IntPoly pseudo_remainder_integer_poly(
    const IntPoly& dividend,
    const IntPoly& divisor) {
    if (divisor.empty()) {
        return IntPoly{};
    }

    IntPoly remainder = dividend;
    normalize_integer_poly(remainder);
    const std::size_t divisor_degree = integer_poly_degree(divisor);
    const BigInt divisor_leading = integer_poly_leading(divisor);

    while (!remainder.empty() && integer_poly_degree(remainder) >= divisor_degree) {
        const std::size_t degree_gap = integer_poly_degree(remainder) - divisor_degree;
        const BigInt remainder_leading = integer_poly_leading(remainder);

        IntPoly scaled_remainder = remainder;
        multiply_integer_coefficients_by_scalar(scaled_remainder, divisor_leading);

        IntPoly subtractor;
        subtractor.resize(divisor.size() + degree_gap, BigInt(0));
        for (std::size_t index = 0; index < divisor.size(); ++index) {
            subtractor[index + degree_gap] = divisor[index] * remainder_leading;
        }

        if (scaled_remainder.size() < subtractor.size()) {
            scaled_remainder.resize(subtractor.size(), BigInt(0));
        }
        for (std::size_t index = 0; index < subtractor.size(); ++index) {
            scaled_remainder[index] -= subtractor[index];
        }

        remainder = std::move(scaled_remainder);
        normalize_integer_poly(remainder);
    }

    return remainder;
}

[[nodiscard]] bool try_divide_integer_coefficients_by_scalar(IntPoly& coefficients, const BigInt& scalar) {
    if (scalar.is_zero()) {
        return coefficients.empty();
    }

    if (scalar == BigInt(1)) {
        return true;
    }
    if (scalar == BigInt(-1)) {
        for (BigInt& coefficient : coefficients.coefficients()) {
            coefficient = -coefficient;
        }
        return true;
    }

    for (const BigInt& coefficient : coefficients.coefficients()) {
        if ((coefficient % scalar) != BigInt(0)) {
            return false;
        }
    }

    for (BigInt& coefficient : coefficients.coefficients()) {
        coefficient /= scalar;
    }
    normalize_integer_poly(coefficients);
    return true;
}

[[nodiscard]] Rational evaluate_integer_polynomial_at_impl(const IntPoly& coefficients, const Rational& value) {
    Rational result(BigInt(0));
    for (std::size_t index = coefficients.size(); index > 0U; --index) {
        result *= value;
        result += Rational(coefficients[index - 1U]);
    }
    return result;
}

[[nodiscard]] Rational evaluate_rational_polynomial_at_impl(const RatPoly& coefficients, const Rational& value) {
    Rational result(BigInt(0));
    for (std::size_t index = coefficients.size(); index > 0U; --index) {
        result *= value;
        result += coefficients[index - 1U];
    }
    return result;
}

[[nodiscard]] std::vector<BigInt> positive_divisors_or_one(const BigInt& value) {
    if (value.is_zero()) {
        return {BigInt(1)};
    }
    auto divisors_result = numtheory::divisors(value.abs());
    if (divisors_result.is_error()) {
        return {BigInt(1)};
    }
    return divisors_result.value();
}

[[nodiscard]] Result<Rational> expr_to_exact_rational_coefficient(ExprPtr expr) {
    if (!expr) {
        return ok(Rational(BigInt(0)));
    }
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return ok(Rational(integer->value));
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return ok(Rational(rational->numerator, rational->denominator));
    }
    return fail<Rational>(make_error(
        CASErrorKind::Unimplemented,
        "partial_fractions richiede coefficienti razionali esatti"));
}

void normalize_rational_coefficients(RatPoly& coefficients) {
    coefficients.normalize([](const Rational& coefficient) {
        return coefficient.numerator().is_zero();
    });
}

struct IntegerSubresultantExecution {
    IntPoly gcd;
    IntegerGcdPath path{IntegerGcdPath::Subresultant};
};

[[nodiscard]] IntegerSubresultantExecution run_integer_subresultant(IntPoly lhs, IntPoly rhs) {
    lhs = primitive_integer_poly(std::move(lhs));
    rhs = primitive_integer_poly(std::move(rhs));

    if (lhs.empty()) {
        return IntegerSubresultantExecution{
            .gcd = std::move(rhs),
            .path = IntegerGcdPath::Subresultant,
        };
    }
    if (rhs.empty()) {
        return IntegerSubresultantExecution{
            .gcd = std::move(lhs),
            .path = IntegerGcdPath::Subresultant,
        };
    }

    if (integer_poly_degree(lhs) < integer_poly_degree(rhs)) {
        std::swap(lhs, rhs);
    }

    IntPoly r_prev = lhs;
    IntPoly r_curr = rhs;
    std::size_t d_prev = integer_poly_degree(r_prev) - integer_poly_degree(r_curr);
    BigInt psi_prev(-1);
    bool first_iteration = true;

    while (!r_curr.empty()) {
        const std::size_t delta = integer_poly_degree(r_prev) - integer_poly_degree(r_curr);
        BigInt beta(0);
        BigInt psi_current(-1);

        if (first_iteration) {
            beta = (delta % 2U == 0U) ? BigInt(-1) : BigInt(1);
        } else {
            const BigInt previous_gamma = integer_poly_leading(r_prev);
            const BigInt numerator = bigint_pow_nonnegative(-previous_gamma, d_prev);
            const BigInt denominator = bigint_pow_nonnegative(psi_prev, d_prev - 1U);
            if (denominator.is_zero() || (numerator % denominator) != BigInt(0)) {
                return IntegerSubresultantExecution{
                    .gcd = IntPoly{},
                    .path = IntegerGcdPath::PrimitiveFallbackPsi,
                };
            }
            psi_current = numerator / denominator;
            beta = -previous_gamma * bigint_pow_nonnegative(psi_current, delta);
        }

        IntPoly remainder = pseudo_remainder_integer_poly(r_prev, r_curr);
        if (!try_divide_integer_coefficients_by_scalar(remainder, beta)) {
            return IntegerSubresultantExecution{
                .gcd = IntPoly{},
                .path = IntegerGcdPath::PrimitiveFallbackBeta,
            };
        }
        normalize_integer_poly(remainder);
        if (remainder.empty()) {
            return IntegerSubresultantExecution{
                .gcd = primitive_integer_poly(std::move(r_curr)),
                .path = IntegerGcdPath::Subresultant,
            };
        }

        r_prev = std::move(r_curr);
        r_curr = std::move(remainder);
        d_prev = delta;
        psi_prev = psi_current;
        first_iteration = false;
    }

    return IntegerSubresultantExecution{
        .gcd = primitive_integer_poly(std::move(r_prev)),
        .path = IntegerGcdPath::Subresultant,
    };
}

}  // namespace

Result<IntPoly> poly_to_integer_poly(const PolyExpr& poly) {
    IntPoly coefficients;
    coefficients.reserve(poly.size());
    for (ExprPtr coefficient : poly.coefficients()) {
        auto integer = expr_to_integer_coefficient(coefficient);
        if (integer.is_error()) {
            return fail<IntPoly>(integer.error());
        }
        coefficients.push_back(integer.value());
    }
    normalize_integer_poly(coefficients);
    return ok(std::move(coefficients));
}

Result<std::vector<BigInt>> poly_to_integer_coefficients(const PolyExpr& poly) {
    auto coefficients = poly_to_integer_poly(poly);
    if (coefficients.is_error()) {
        return fail<std::vector<BigInt>>(coefficients.error());
    }
    return ok(coefficients.value().coefficients());
}

Rational evaluate_integer_polynomial_at(const IntPoly& coefficients, const Rational& value) {
    return evaluate_integer_polynomial_at_impl(coefficients, value);
}

Rational evaluate_rational_polynomial_at(const RatPoly& coefficients, const Rational& value) {
    return evaluate_rational_polynomial_at_impl(coefficients, value);
}

Result<RatPoly> poly_to_rational_poly(const PolyExpr& poly) {
    RatPoly coefficients;
    coefficients.reserve(poly.size());
    for (ExprPtr coefficient : poly.coefficients()) {
        auto rational = expr_to_exact_rational_coefficient(coefficient);
        if (rational.is_error()) {
            return fail<RatPoly>(rational.error());
        }
        coefficients.push_back(rational.value());
    }
    normalize_rational_coefficients(coefficients);
    return ok(std::move(coefficients));
}

ExprPtr poly_make_integer(AstArena& arena, long long value) {
    return arena.make<IntegerLit>(BigInt(value));
}

bool poly_is_zero_expr(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value.is_zero();
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator.is_zero();
    }
    return false;
}

bool poly_is_one_expr(ExprPtr expr) {
    static const BigInt one(1);
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value == one;
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator == one && rational->denominator == one;
    }
    return false;
}

bool poly_is_minus_one_expr(ExprPtr expr) {
    static const BigInt minus_one(-1);
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value == minus_one;
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator == minus_one && rational->denominator == BigInt(1);
    }
    return false;
}

bool poly_depends_on(ExprPtr expr, const std::string& variable_name) {
    if (!expr) {
        return false;
    }

    if (const auto* symbol = expr_cast<Symbol>(expr)) {
        return symbol->name == variable_name;
    }

    switch (expr->kind) {
    case ExprKind::Null:
    case ExprKind::IntegerLit:
    case ExprKind::RationalLit:
    case ExprKind::DecimalLit:
    case ExprKind::Symbol:
    case ExprKind::Constant:
        return false;
    case ExprKind::Unary:
        return poly_depends_on(expr_ref<Unary>(expr).operand, variable_name);
    case ExprKind::Binary: {
        const auto& binary = expr_ref<Binary>(expr);
        return poly_depends_on(binary.left, variable_name) || poly_depends_on(binary.right, variable_name);
    }
    case ExprKind::FuncCall: {
        for (ExprPtr arg : expr_ref<FuncCall>(expr).args) {
            if (poly_depends_on(arg, variable_name)) {
                return true;
            }
        }
        return false;
    }
    case ExprKind::Sum: {
        for (ExprPtr term : expr_ref<Sum>(expr).terms) {
            if (poly_depends_on(term, variable_name)) {
                return true;
            }
        }
        return false;
    }
    case ExprKind::Product: {
        for (ExprPtr factor : expr_ref<Product>(expr).factors) {
            if (poly_depends_on(factor, variable_name)) {
                return true;
            }
        }
        return false;
    }
    case ExprKind::Integral: {
        const auto& node = expr_ref<Integral>(expr);
        return poly_depends_on(node.integrand, variable_name) ||
               (node.lower.has_value() && poly_depends_on(*node.lower, variable_name)) ||
               (node.upper.has_value() && poly_depends_on(*node.upper, variable_name));
    }
    case ExprKind::Derivative:
        return poly_depends_on(expr_ref<Derivative>(expr).expression, variable_name);
    case ExprKind::Limit: {
        const auto& node = expr_ref<Limit>(expr);
        return poly_depends_on(node.expression, variable_name) || poly_depends_on(node.point, variable_name);
    }
    case ExprKind::RootOf:
        return poly_depends_on(expr_ref<RootOf>(expr).polynomial, variable_name);
    case ExprKind::Matrix: {
        for (ExprPtr element : expr_ref<Matrix>(expr).elements) {
            if (poly_depends_on(element, variable_name)) {
                return true;
            }
        }
        return false;
    }
    }

    return false;
}

bool poly_contains_decimal_literal(ExprPtr expr) {
    if (!expr) {
        return false;
    }
    if (expr_is<DecimalLit>(expr)) {
        return true;
    }

    return visit_expr(expr, [](const auto& node) {
        using Node = std::decay_t<decltype(node)>;
        if constexpr (
            std::is_same_v<Node, IntegerLit> ||
            std::is_same_v<Node, RationalLit> ||
            std::is_same_v<Node, Symbol> ||
            std::is_same_v<Node, Constant>) {
            return false;
        } else if constexpr (std::is_same_v<Node, Unary>) {
            return poly_contains_decimal_literal(node.operand);
        } else if constexpr (std::is_same_v<Node, Binary>) {
            return poly_contains_decimal_literal(node.left) || poly_contains_decimal_literal(node.right);
        } else if constexpr (std::is_same_v<Node, FuncCall>) {
            return std::any_of(node.args.begin(), node.args.end(), [](ExprPtr arg) {
                return poly_contains_decimal_literal(arg);
            });
        } else if constexpr (std::is_same_v<Node, Sum>) {
            return std::any_of(node.terms.begin(), node.terms.end(), [](ExprPtr term) {
                return poly_contains_decimal_literal(term);
            });
        } else if constexpr (std::is_same_v<Node, Product>) {
            return std::any_of(node.factors.begin(), node.factors.end(), [](ExprPtr factor) {
                return poly_contains_decimal_literal(factor);
            });
        } else if constexpr (std::is_same_v<Node, Integral>) {
            return poly_contains_decimal_literal(node.integrand) ||
                   (node.lower.has_value() && poly_contains_decimal_literal(*node.lower)) ||
                   (node.upper.has_value() && poly_contains_decimal_literal(*node.upper));
        } else if constexpr (std::is_same_v<Node, Derivative>) {
            return poly_contains_decimal_literal(node.expression);
        } else if constexpr (std::is_same_v<Node, Limit>) {
            return poly_contains_decimal_literal(node.expression) || poly_contains_decimal_literal(node.point);
        } else if constexpr (std::is_same_v<Node, RootOf>) {
            return poly_contains_decimal_literal(node.polynomial);
        } else if constexpr (std::is_same_v<Node, Matrix>) {
            return std::any_of(node.elements.begin(), node.elements.end(), [](ExprPtr element) {
                return poly_contains_decimal_literal(element);
            });
        } else {
            return false;
        }
    });
}

Result<ExprPtr> poly_simplify_expr(ExprPtr expr, symbolic::CASContext& ctx) {
    auto simplified = ctx.simplify(expr);
    if (simplified.is_error()) {
        return fail<ExprPtr>(simplified.error());
    }
    return simplified;
}

Result<ExprPtr> poly_clone_into_context(ExprPtr expr, symbolic::CASContext& ctx) {
    auto cloned = symbolic::materialize_expr(expr, ctx.arena());
    if (cloned.is_error()) {
        return fail<ExprPtr>(cloned.error());
    }
    return poly_simplify_expr(cloned.value(), ctx);
}

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
        return ok(make_zero_poly());
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
        return fail<PolyExpr>(make_error(CASErrorKind::Undefined, "Divisione polinomiale per coefficiente nullo"));
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
        return fail<std::size_t>(make_error(CASErrorKind::InvalidArgument, "Null exponent"));
    }

    BigInt value;
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        value = integer->value;
    } else if (const auto* rational = expr_cast<RationalLit>(expr)) {
        if (rational->denominator == BigInt(1)) {
            value = rational->numerator;
        } else {
            return fail<std::size_t>(make_error(
                CASErrorKind::Unimplemented,
                "Le potenze polinomiali supportano solo esponenti interi"));
        }
    } else {
        return fail<std::size_t>(make_error(
            CASErrorKind::Unimplemented,
            "Le potenze polinomiali supportano solo esponenti numerici esatti"));
    }

    if (value.is_negative()) {
        return fail<std::size_t>(make_error(
            CASErrorKind::Unimplemented,
            "Le potenze polinomiali supportano solo esponenti interi non negativi (ricevuto " + value.decimal() + ")"));
    }

    const std::string decimal = value.decimal();
    const auto parsed = parse_bounded_unsigned_decimal<std::size_t>(decimal);
    if (!parsed.has_value()) {
        return fail<std::size_t>(make_error(
            CASErrorKind::Unimplemented,
            "Le potenze polinomiali richiedono esponenti rappresentabili nel limite interno",
            "L'esponente intero non negativo supera la capacita' supportata da size_t."));
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
        return fail<PolyDivisionResult>(make_error(
            CASErrorKind::Undefined,
            "Divisione polinomiale per divisore nullo"));
    }

    PolyExpr quotient = make_zero_poly();
    PolyExpr remainder = dividend;
    const std::size_t divisor_degree = poly_degree(divisor);
    ExprPtr divisor_leading = leading_coefficient(divisor);

    while (!is_zero_poly(remainder) && poly_degree(remainder) >= divisor_degree) {
        ExprPtr remainder_leading = leading_coefficient(remainder);
        const std::size_t degree_gap = poly_degree(remainder) - divisor_degree;
        const std::size_t previous_degree = poly_degree(remainder);

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

        if (!is_zero_poly(next_remainder.value()) && poly_degree(next_remainder.value()) >= previous_degree) {
            return fail<PolyDivisionResult>(make_error(
                CASErrorKind::Unimplemented,
                "La divisione polinomiale non converge con coefficienti correnti",
                "Sarà sostituita da PRS subresultante nelle fasi successive."));
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
    if (is_zero_poly(poly)) {
        return ok(make_zero_poly());
    }
    return poly_divide_by_scalar(poly, leading_coefficient(poly), ctx);
}

bool is_zero_integer_poly(const IntPoly& coefficients) {
    return coefficients.is_zero();
}

IntPoly gcd_integer_poly_primitive(IntPoly lhs, IntPoly rhs) {
    lhs = primitive_integer_poly(std::move(lhs));
    rhs = primitive_integer_poly(std::move(rhs));

    if (lhs.empty()) {
        return rhs;
    }
    if (rhs.empty()) {
        return lhs;
    }

    while (!rhs.empty()) {
        IntPoly remainder = pseudo_remainder_integer_poly(lhs, rhs);
        remainder = primitive_integer_poly(std::move(remainder));
        lhs = std::move(rhs);
        rhs = std::move(remainder);
    }

    return primitive_integer_poly(std::move(lhs));
}

std::optional<IntPoly> gcd_integer_poly_subresultant(IntPoly lhs, IntPoly rhs) {
    IntegerSubresultantExecution execution = run_integer_subresultant(std::move(lhs), std::move(rhs));
    if (execution.path != IntegerGcdPath::Subresultant) {
        return std::nullopt;
    }
    return std::move(execution.gcd);
}

IntegerGcdResult gcd_integer_poly_with_subresultant(IntPoly lhs, IntPoly rhs) {
    lhs = primitive_integer_poly(std::move(lhs));
    rhs = primitive_integer_poly(std::move(rhs));

    IntegerSubresultantExecution execution = run_integer_subresultant(lhs, rhs);
    if (execution.path == IntegerGcdPath::Subresultant) {
        return IntegerGcdResult{
            .gcd = std::move(execution.gcd),
            .path = execution.path,
        };
    }

    return IntegerGcdResult{
        .gcd = gcd_integer_poly_primitive(std::move(lhs), std::move(rhs)),
        .path = execution.path,
    };
}

std::optional<RationalRootCandidate> find_rational_root_candidate(const IntPoly& coefficients) {
    if (coefficients.size() <= 1U) {
        return std::nullopt;
    }
    if (coefficients.constant_term().is_zero()) {
        return RationalRootCandidate{
            .numerator = BigInt(0),
            .denominator = BigInt(1),
        };
    }

    const std::vector<BigInt> numerator_divisors = positive_divisors_or_one(coefficients.constant_term());
    const std::vector<BigInt> denominator_divisors = positive_divisors_or_one(coefficients.leading_coeff());

    for (const BigInt& denominator_divisor : denominator_divisors) {
        for (const BigInt& numerator_divisor : numerator_divisors) {
            for (int sign : {1, -1}) {
                BigInt numerator = sign > 0 ? numerator_divisor : -numerator_divisor;
                BigInt denominator = denominator_divisor;
                const BigInt common = gcd(numerator.abs(), denominator);
                if (common != BigInt(1)) {
                    numerator /= common;
                    denominator /= common;
                }

                if (evaluate_integer_polynomial_at_impl(coefficients, Rational(numerator, denominator)).numerator().is_zero()) {
                    return RationalRootCandidate{
                        .numerator = std::move(numerator),
                        .denominator = std::move(denominator),
                    };
                }
            }
        }
    }

    return std::nullopt;
}

Result<std::vector<Rational>> poly_to_rational_coefficients(const PolyExpr& poly) {
    auto coefficients = poly_to_rational_poly(poly);
    if (coefficients.is_error()) {
        return fail<std::vector<Rational>>(coefficients.error());
    }
    return ok(coefficients.value().coefficients());
}

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
        return fail<PolyExpr>(make_error(CASErrorKind::InvalidArgument, "Espressione polinomiale nulla"));
    }
    if (poly_contains_decimal_literal(expr)) {
        return fail<PolyExpr>(make_error(
            CASErrorKind::Unimplemented,
            "I literal decimali non sono supportati nel modello polinomiale esatto"));
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
            return fail<PolyExpr>(make_error(
                CASErrorKind::Unimplemented,
                "Solo la negazione unaria e' supportata nel parsing polinomiale"));
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
                return fail<PolyExpr>(make_error(
                    CASErrorKind::Unimplemented,
                    "I quozienti con denominatore dipendente dalla variabile non sono polinomi"));
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
            return fail<PolyExpr>(make_error(
                CASErrorKind::Unimplemented,
                "Il modulo non e' supportato nel modello polinomiale"));
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

    return fail<PolyExpr>(make_error(
        CASErrorKind::Unimplemented,
        "Espressione non supportata dal parser polinomiale univariato",
        "Le chiamate di funzione e i nodi dipendenti dalla variabile saranno gestiti in fasi successive."));
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
