#include "calculus_internal.hpp"

#include "cas/error.hpp"
#include "cas/rational.hpp"

#include <algorithm>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus {
namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

[[nodiscard]] ExprPtr make_integer(AstArena& arena, long long value) {
    return arena.make<IntegerLit>(BigInt(value));
}

[[nodiscard]] ExprPtr make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs) {
    return arena.make<Binary>(op, lhs, rhs);
}

[[nodiscard]] ExprPtr make_sum(AstArena& arena, std::vector<ExprPtr> terms) {
    if (terms.empty()) {
        return make_integer(arena, 0);
    }
    if (terms.size() == 1U) {
        return terms.front();
    }
    return arena.make<Sum>(std::move(terms));
}

[[nodiscard]] ExprPtr make_product(AstArena& arena, std::vector<ExprPtr> factors) {
    if (factors.empty()) {
        return make_integer(arena, 1);
    }
    if (factors.size() == 1U) {
        return factors.front();
    }
    return arena.make<Product>(std::move(factors));
}

[[nodiscard]] ExprPtr make_power(AstArena& arena, ExprPtr base, ExprPtr exponent) {
    return arena.make<Binary>(BinaryOp::Pow, base, exponent);
}

[[nodiscard]] ExprPtr make_rational_expr(AstArena& arena, const Rational& value) {
    if (value.is_integer()) {
        return arena.make<IntegerLit>(value.numerator());
    }
    return arena.make<RationalLit>(value.numerator(), value.denominator());
}

[[nodiscard]] ExprPtr make_signed_infinity(AstArena& arena, int sign) {
    ExprPtr infinity = arena.make<Constant>(MathConstant::Infinity);
    if (sign < 0) {
        return arena.make<Unary>(UnaryOp::Neg, infinity);
    }
    return infinity;
}

[[nodiscard]] std::optional<Rational> rational_from_expr(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return Rational(integer->value);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return Rational(rational->numerator, rational->denominator);
    }
    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op != UnaryOp::Neg) {
            return std::nullopt;
        }
        auto operand = rational_from_expr(unary->operand);
        if (!operand.has_value()) {
            return std::nullopt;
        }
        return -operand.value();
    }
    return std::nullopt;
}

void trim_polynomial(std::vector<Rational>& coefficients) {
    while (coefficients.size() > 1U && coefficients.back() == Rational(BigInt(0))) {
        coefficients.pop_back();
    }
}

[[nodiscard]] std::vector<Rational> add_polynomials(
    const std::vector<Rational>& lhs,
    const std::vector<Rational>& rhs) {
    std::vector<Rational> result(std::max(lhs.size(), rhs.size()), Rational(BigInt(0)));
    for (std::size_t index = 0; index < lhs.size(); ++index) {
        result[index] += lhs[index];
    }
    for (std::size_t index = 0; index < rhs.size(); ++index) {
        result[index] += rhs[index];
    }
    trim_polynomial(result);
    return result;
}

[[nodiscard]] std::vector<Rational> subtract_polynomials(
    const std::vector<Rational>& lhs,
    const std::vector<Rational>& rhs) {
    std::vector<Rational> result(std::max(lhs.size(), rhs.size()), Rational(BigInt(0)));
    for (std::size_t index = 0; index < lhs.size(); ++index) {
        result[index] += lhs[index];
    }
    for (std::size_t index = 0; index < rhs.size(); ++index) {
        result[index] -= rhs[index];
    }
    trim_polynomial(result);
    return result;
}

[[nodiscard]] std::vector<Rational> multiply_polynomials(
    const std::vector<Rational>& lhs,
    const std::vector<Rational>& rhs) {
    std::vector<Rational> result(lhs.size() + rhs.size() - 1U, Rational(BigInt(0)));
    for (std::size_t lhs_index = 0; lhs_index < lhs.size(); ++lhs_index) {
        for (std::size_t rhs_index = 0; rhs_index < rhs.size(); ++rhs_index) {
            result[lhs_index + rhs_index] += lhs[lhs_index] * rhs[rhs_index];
        }
    }
    trim_polynomial(result);
    return result;
}

[[nodiscard]] std::optional<unsigned int> nonnegative_integer_from_expr(ExprPtr expr) {
    const auto* integer = expr_cast<IntegerLit>(expr);
    if (integer == nullptr || integer->value.is_negative()) {
        return std::nullopt;
    }

    unsigned int value = 0U;
    for (char ch : integer->value.decimal()) {
        const unsigned int digit = static_cast<unsigned int>(ch - '0');
        if (value > (std::numeric_limits<unsigned int>::max() - digit) / 10U) {
            return std::nullopt;
        }
        value = value * 10U + digit;
    }
    return value;
}

[[nodiscard]] std::optional<std::vector<Rational>> polynomial_from_expr(ExprPtr expr, const Symbol& var) {
    if (const auto rational = rational_from_expr(expr)) {
        return std::vector<Rational>{rational.value()};
    }

    if (const auto* symbol = expr_cast<Symbol>(expr)) {
        if (symbol->name == var.name) {
            return std::vector<Rational>{Rational(BigInt(0)), Rational(BigInt(1))};
        }
        return std::nullopt;
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op != UnaryOp::Neg) {
            return std::nullopt;
        }
        auto operand = polynomial_from_expr(unary->operand, var);
        if (!operand.has_value()) {
            return std::nullopt;
        }
        for (auto& coefficient : operand.value()) {
            coefficient = -coefficient;
        }
        return operand;
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Add || binary->op == BinaryOp::Sub || binary->op == BinaryOp::Mul) {
            auto lhs = polynomial_from_expr(binary->left, var);
            auto rhs = polynomial_from_expr(binary->right, var);
            if (!lhs.has_value() || !rhs.has_value()) {
                return std::nullopt;
            }
            if (binary->op == BinaryOp::Add) {
                return add_polynomials(lhs.value(), rhs.value());
            }
            if (binary->op == BinaryOp::Sub) {
                return subtract_polynomials(lhs.value(), rhs.value());
            }
            return multiply_polynomials(lhs.value(), rhs.value());
        }

        if (binary->op == BinaryOp::Div) {
            auto numerator = polynomial_from_expr(binary->left, var);
            auto denominator = rational_from_expr(binary->right);
            if (!numerator.has_value() || !denominator.has_value() || denominator->numerator().is_zero()) {
                return std::nullopt;
            }
            for (auto& coefficient : numerator.value()) {
                coefficient /= denominator.value();
            }
            trim_polynomial(numerator.value());
            return numerator;
        }

        if (binary->op == BinaryOp::Pow) {
            auto base = polynomial_from_expr(binary->left, var);
            auto exponent = nonnegative_integer_from_expr(binary->right);
            if (!base.has_value() || !exponent.has_value()) {
                return std::nullopt;
            }
            std::vector<Rational> result{Rational(BigInt(1))};
            for (unsigned int index = 0; index < exponent.value(); ++index) {
                result = multiply_polynomials(result, base.value());
            }
            return result;
        }

        return std::nullopt;
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<Rational> result{Rational(BigInt(0))};
        for (ExprPtr term : sum->terms) {
            auto polynomial = polynomial_from_expr(term, var);
            if (!polynomial.has_value()) {
                return std::nullopt;
            }
            result = add_polynomials(result, polynomial.value());
        }
        return result;
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        std::vector<Rational> result{Rational(BigInt(1))};
        for (ExprPtr factor : product->factors) {
            auto polynomial = polynomial_from_expr(factor, var);
            if (!polynomial.has_value()) {
                return std::nullopt;
            }
            result = multiply_polynomials(result, polynomial.value());
        }
        return result;
    }

    return std::nullopt;
}

[[nodiscard]] ExprPtr polynomial_to_expr(AstArena& arena, const std::vector<Rational>& coefficients, const Symbol& var) {
    std::vector<ExprPtr> terms;
    for (std::size_t degree = 0; degree < coefficients.size(); ++degree) {
        if (coefficients[degree] == Rational(BigInt(0))) {
            continue;
        }
        ExprPtr term = make_rational_expr(arena, coefficients[degree]);
        if (degree >= 1U) {
            ExprPtr power = degree == 1U
                ? arena.make<Symbol>(var)
                : make_power(arena, arena.make<Symbol>(var), make_integer(arena, static_cast<long long>(degree)));
            term = make_product(arena, {term, power});
        }
        terms.push_back(term);
    }
    return make_sum(arena, std::move(terms));
}

[[nodiscard]] int sign_of_rational(const Rational& value) {
    if (value.numerator().is_zero()) {
        return 0;
    }
    return value.numerator().is_negative() ? -1 : 1;
}

[[nodiscard]] Rational evaluate_polynomial(const std::vector<Rational>& coefficients, const Rational& point) {
    Rational result(BigInt(0));
    for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
        result *= point;
        result += *it;
    }
    return result;
}

[[nodiscard]] std::optional<std::vector<Rational>> divide_by_linear_factor(
    const std::vector<Rational>& coefficients,
    const Rational& root) {
    if (coefficients.size() <= 1U) {
        return std::nullopt;
    }

    const std::size_t degree = coefficients.size() - 1U;
    std::vector<Rational> quotient(degree, Rational(BigInt(0)));
    quotient[degree - 1U] = coefficients.back();

    for (std::size_t index = degree - 1U; index > 0U; --index) {
        quotient[index - 1U] = coefficients[index] + quotient[index] * root;
    }

    const Rational remainder = coefficients.front() + quotient.front() * root;
    if (!(remainder == Rational(BigInt(0)))) {
        return std::nullopt;
    }

    trim_polynomial(quotient);
    return quotient;
}

struct RootFactoredPolynomial {
    std::vector<Rational> reduced;
    unsigned int multiplicity{0U};
};

[[nodiscard]] RootFactoredPolynomial factor_out_linear_root(
    std::vector<Rational> coefficients,
    const Rational& root) {
    trim_polynomial(coefficients);

    RootFactoredPolynomial result{
        .reduced = std::move(coefficients),
        .multiplicity = 0U,
    };

    while (result.reduced.size() > 1U) {
        auto quotient = divide_by_linear_factor(result.reduced, root);
        if (!quotient.has_value()) {
            break;
        }
        result.reduced = std::move(quotient.value());
        ++result.multiplicity;
    }

    return result;
}

}  // namespace

std::optional<ExprPtr> cancel_common_linear_factor(
    ExprPtr numerator,
    ExprPtr denominator,
    const Symbol& var,
    ExprPtr point,
    AstArena& arena) {
    auto rational_point = rational_from_expr(point);
    auto numerator_poly = polynomial_from_expr(numerator, var);
    auto denominator_poly = polynomial_from_expr(denominator, var);
    if (!rational_point.has_value() || !numerator_poly.has_value() || !denominator_poly.has_value()) {
        return std::nullopt;
    }

    bool reduced = false;
    while (numerator_poly->size() > 1U && denominator_poly->size() > 1U) {
        auto reduced_numerator = divide_by_linear_factor(numerator_poly.value(), rational_point.value());
        auto reduced_denominator = divide_by_linear_factor(denominator_poly.value(), rational_point.value());
        if (!reduced_numerator.has_value() || !reduced_denominator.has_value()) {
            break;
        }
        numerator_poly = std::move(reduced_numerator);
        denominator_poly = std::move(reduced_denominator);
        reduced = true;
    }

    if (!reduced) {
        return std::nullopt;
    }

    return make_binary(
        arena,
        BinaryOp::Div,
        polynomial_to_expr(arena, numerator_poly.value(), var),
        polynomial_to_expr(arena, denominator_poly.value(), var));
}

std::optional<Result<ExprPtr>> try_polynomial_pole_limit(
    ExprPtr numerator,
    ExprPtr denominator,
    const Symbol& var,
    ExprPtr point,
    LimitDirection dir,
    AstArena& arena) {
    auto rational_point = rational_from_expr(point);
    auto numerator_poly = polynomial_from_expr(numerator, var);
    auto denominator_poly = polynomial_from_expr(denominator, var);
    if (!rational_point.has_value() || !numerator_poly.has_value() || !denominator_poly.has_value()) {
        return std::nullopt;
    }

    auto factored_numerator = factor_out_linear_root(std::move(numerator_poly.value()), rational_point.value());
    auto factored_denominator = factor_out_linear_root(std::move(denominator_poly.value()), rational_point.value());
    if (factored_denominator.multiplicity == 0U ||
        factored_denominator.multiplicity <= factored_numerator.multiplicity) {
        return std::nullopt;
    }

    const Rational numerator_value = evaluate_polynomial(factored_numerator.reduced, rational_point.value());
    const Rational denominator_value = evaluate_polynomial(factored_denominator.reduced, rational_point.value());
    if (denominator_value.numerator().is_zero()) {
        return std::nullopt;
    }
    const int coefficient_sign = sign_of_rational(numerator_value / denominator_value);
    if (coefficient_sign == 0) {
        return std::nullopt;
    }

    const unsigned int residual_power = factored_denominator.multiplicity - factored_numerator.multiplicity;
    auto sign_for_direction = [&](LimitDirection direction) {
        if ((residual_power % 2U) == 0U) {
            return coefficient_sign;
        }
        return direction == LimitDirection::Left ? -coefficient_sign : coefficient_sign;
    };

    if (dir == LimitDirection::Both) {
        const int left_sign = sign_for_direction(LimitDirection::Left);
        const int right_sign = sign_for_direction(LimitDirection::Right);
        if (left_sign != right_sign) {
            return fail<ExprPtr>(make_error(
                CASErrorKind::Undefined,
                "Il limite bilaterale non esiste: i limiti laterali divergono con segno opposto"));
        }
        return ok(make_signed_infinity(arena, left_sign));
    }

    return ok(make_signed_infinity(arena, sign_for_direction(dir)));
}

std::optional<Result<ExprPtr>> try_logarithmic_root_limit(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    LimitDirection dir,
    AstArena& arena) {
    const auto* call = expr_cast<FuncCall>(expr);
    // A41: "log" parses to BuiltinOp::Log, a distinct enum value from Ln, but
    // it is the SAME natural logarithm everywhere else in the engine (see
    // simplify_exp_log.cpp, differentiate_rules.cpp, integrate_elementary.cpp,
    // limit_infinite.cpp, ...) — this was the one call site still checking
    // only Ln, silently excluding log(x) from the one-sided log-pole rule and
    // making limit(x*log(x), x, 0, plus) fail with a raw division-by-zero
    // instead of resolving log(x) -> -infinity.
    if (call == nullptr ||
        (call->func_id != BuiltinOp::Ln && call->func_id != BuiltinOp::Log) ||
        call->args.size() != 1U) {
        return std::nullopt;
    }

    auto rational_point = rational_from_expr(point);
    auto argument_poly = polynomial_from_expr(call->args.front(), var);
    if (!rational_point.has_value() || !argument_poly.has_value()) {
        return std::nullopt;
    }

    auto factored_argument = factor_out_linear_root(std::move(argument_poly.value()), rational_point.value());
    if (factored_argument.multiplicity == 0U) {
        return std::nullopt;
    }

    const Rational coefficient = evaluate_polynomial(factored_argument.reduced, rational_point.value());
    const int coefficient_sign = sign_of_rational(coefficient);
    if (coefficient_sign == 0) {
        return std::nullopt;
    }

    auto argument_sign = [&](LimitDirection direction) {
        if ((factored_argument.multiplicity % 2U) == 0U) {
            return coefficient_sign;
        }
        return direction == LimitDirection::Left ? -coefficient_sign : coefficient_sign;
    };

    auto make_undefined = []() {
        return fail<ExprPtr>(make_error(
            CASErrorKind::Undefined,
            "Il limite del logaritmo non e' definito sul lato richiesto"));
    };

    if (dir == LimitDirection::Both) {
        const int left_sign = argument_sign(LimitDirection::Left);
        const int right_sign = argument_sign(LimitDirection::Right);
        if (left_sign <= 0 || right_sign <= 0) {
            return make_undefined();
        }
        return ok(make_signed_infinity(arena, -1));
    }

    if (argument_sign(dir) <= 0) {
        return make_undefined();
    }

    return ok(make_signed_infinity(arena, -1));
}


}  // namespace cas::calculus
