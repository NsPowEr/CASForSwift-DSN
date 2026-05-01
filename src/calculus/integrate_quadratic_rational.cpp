#include "integrate_engine.hpp"

#include <array>
#include <optional>
#include <vector>

namespace cas::calculus::integrate_detail {
namespace {

struct DegreeTwoPolynomial {
    std::array<Rational, 3> coeffs{
        Rational(BigInt(0)),
        Rational(BigInt(0)),
        Rational(BigInt(0)),
    };
};

[[nodiscard]] bool is_zero_rational(const Rational& value) {
    return value.numerator().is_zero();
}

[[nodiscard]] DegreeTwoPolynomial add_poly(DegreeTwoPolynomial lhs, const DegreeTwoPolynomial& rhs) {
    for (std::size_t i = 0; i < lhs.coeffs.size(); ++i) {
        lhs.coeffs[i] += rhs.coeffs[i];
    }
    return lhs;
}

[[nodiscard]] DegreeTwoPolynomial negate_poly(DegreeTwoPolynomial value) {
    for (Rational& coeff : value.coeffs) {
        coeff = -coeff;
    }
    return value;
}

[[nodiscard]] std::optional<DegreeTwoPolynomial> multiply_poly(
    const DegreeTwoPolynomial& lhs,
    const DegreeTwoPolynomial& rhs) {
    DegreeTwoPolynomial result;
    for (std::size_t i = 0; i < lhs.coeffs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.coeffs.size(); ++j) {
            if (i + j >= result.coeffs.size()) {
                if (!is_zero_rational(lhs.coeffs[i]) && !is_zero_rational(rhs.coeffs[j])) {
                    return std::nullopt;
                }
                continue;
            }
            result.coeffs[i + j] += lhs.coeffs[i] * rhs.coeffs[j];
        }
    }
    return result;
}

[[nodiscard]] std::optional<DegreeTwoPolynomial> parse_degree_two_polynomial(ExprPtr expr, const Symbol& var) {
    if (!expr) {
        return std::nullopt;
    }

    if (auto scalar = exact_scalar_from_expr(expr); scalar.has_value()) {
        DegreeTwoPolynomial result;
        result.coeffs[0] = *scalar;
        return result;
    }

    if (is_same_symbol(expr, var)) {
        DegreeTwoPolynomial result;
        result.coeffs[1] = Rational(BigInt(1));
        return result;
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op != UnaryOp::Neg) {
            return std::nullopt;
        }
        auto inner = parse_degree_two_polynomial(unary->operand, var);
        if (!inner.has_value()) {
            return std::nullopt;
        }
        return negate_poly(*inner);
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Add || binary->op == BinaryOp::Sub) {
            auto left = parse_degree_two_polynomial(binary->left, var);
            auto right = parse_degree_two_polynomial(binary->right, var);
            if (!left.has_value() || !right.has_value()) {
                return std::nullopt;
            }
            if (binary->op == BinaryOp::Sub) {
                right = negate_poly(*right);
            }
            return add_poly(*left, *right);
        }
        if (binary->op == BinaryOp::Mul) {
            auto left = parse_degree_two_polynomial(binary->left, var);
            auto right = parse_degree_two_polynomial(binary->right, var);
            if (!left.has_value() || !right.has_value()) {
                return std::nullopt;
            }
            return multiply_poly(*left, *right);
        }
        if (binary->op == BinaryOp::Div) {
            auto left = parse_degree_two_polynomial(binary->left, var);
            auto right_scalar = exact_scalar_from_expr(binary->right);
            if (!left.has_value() || !right_scalar.has_value() || right_scalar->numerator().is_zero()) {
                return std::nullopt;
            }
            for (Rational& coeff : left->coeffs) {
                coeff /= *right_scalar;
            }
            return left;
        }
        if (binary->op == BinaryOp::Pow && is_same_symbol(binary->left, var)) {
            if (is_rational_value(binary->right, 0, 1)) {
                DegreeTwoPolynomial result;
                result.coeffs[0] = Rational(BigInt(1));
                return result;
            }
            if (is_rational_value(binary->right, 1, 1)) {
                DegreeTwoPolynomial result;
                result.coeffs[1] = Rational(BigInt(1));
                return result;
            }
            if (is_rational_value(binary->right, 2, 1)) {
                DegreeTwoPolynomial result;
                result.coeffs[2] = Rational(BigInt(1));
                return result;
            }
        }
        return std::nullopt;
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        DegreeTwoPolynomial result;
        for (ExprPtr term : sum->terms) {
            auto parsed = parse_degree_two_polynomial(term, var);
            if (!parsed.has_value()) {
                return std::nullopt;
            }
            result = add_poly(result, *parsed);
        }
        return result;
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        DegreeTwoPolynomial result;
        result.coeffs[0] = Rational(BigInt(1));
        for (ExprPtr factor : product->factors) {
            auto parsed = parse_degree_two_polynomial(factor, var);
            if (!parsed.has_value()) {
                return std::nullopt;
            }
            auto multiplied = multiply_poly(result, *parsed);
            if (!multiplied.has_value()) {
                return std::nullopt;
            }
            result = *multiplied;
        }
        return result;
    }

    return std::nullopt;
}

[[nodiscard]] ExprPtr make_linear_polynomial(
    AstArena& arena,
    const Rational& coefficient,
    const Rational& constant,
    const Symbol& var) {
    std::vector<ExprPtr> terms;
    if (!coefficient.numerator().is_zero()) {
        terms.push_back(make_product(arena, {make_rational(arena, coefficient), arena.make<Symbol>(var)}));
    }
    if (!constant.numerator().is_zero()) {
        terms.push_back(make_rational(arena, constant));
    }
    return make_sum(arena, std::move(terms));
}

[[nodiscard]] BigInt bigint_isqrt(const BigInt& n) {
    if (n.is_zero() || n.is_negative()) return BigInt(0);
    const std::size_t bits = n.bit_length();
    BigInt x = n.shift_right_bits(bits / 2U);
    if (x.is_zero()) x = BigInt(1);
    while (true) {
        BigInt x1 = (x + n / x) / BigInt(2);
        if (x1 >= x) return x;
        x = std::move(x1);
    }
}

[[nodiscard]] std::optional<BigInt> perfect_square_root(const BigInt& n) {
    if (n.is_negative()) return std::nullopt;
    if (n.is_zero()) return BigInt(0);
    BigInt root = bigint_isqrt(n);
    if (root * root == n) return root;
    return std::nullopt;
}

[[nodiscard]] ExprPtr make_sqrt(AstArena& arena, const Rational& value) {
    auto numerator_root = perfect_square_root(value.numerator());
    auto denominator_root = perfect_square_root(value.denominator());
    if (numerator_root.has_value() && denominator_root.has_value()) {
        return make_rational(arena, Rational(*numerator_root, *denominator_root));
    }
    return make_function(arena, "sqrt", {make_rational(arena, value)});
}

}  // namespace

std::optional<QuadraticArgument> extract_quadratic_argument(ExprPtr expr, const Symbol& var) {
    auto parsed = parse_degree_two_polynomial(expr, var);
    if (!parsed.has_value()) {
        return std::nullopt;
    }
    return QuadraticArgument{
        .quadratic = parsed->coeffs[2],
        .linear = parsed->coeffs[1],
        .constant = parsed->coeffs[0],
    };
}

Result<ExprPtr> Integrator::integrate_linear_over_quadratic(const Binary& quotient, const Symbol& var) {
    auto numerator = extract_affine_argument(quotient.left, var);
    auto denominator = extract_quadratic_argument(quotient.right, var);
    if (!numerator.has_value() || !denominator.has_value() || denominator->quadratic.numerator().is_zero()) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::Unimplemented,
            "Linear-over-quadratic rational integration requires exact quadratic denominator"));
    }

    const Rational two(BigInt(2));
    const Rational four(BigInt(4));
    const Rational alpha = numerator->coefficient / (two * denominator->quadratic);
    const Rational beta = numerator->constant - alpha * denominator->linear;
    const Rational completed_discriminant =
        four * denominator->quadratic * denominator->constant - denominator->linear * denominator->linear;

    std::vector<ExprPtr> primitive_terms;
    if (!alpha.numerator().is_zero()) {
        primitive_terms.push_back(make_product(arena_, {
            make_rational(arena_, alpha),
            make_function(arena_, "ln", {make_function(arena_, "abs", {quotient.right})}),
        }));
    }

    if (!beta.numerator().is_zero()) {
        if (completed_discriminant.numerator().is_zero()) {
            return fail<ExprPtr>(make_error(
                CASErrorKind::Unimplemented,
                "Repeated quadratic denominator should be handled by partial fractions"));
        }
        if (completed_discriminant.numerator().is_negative()) {
            return fail<ExprPtr>(make_error(
                CASErrorKind::Unimplemented,
                "Real-factorable quadratic denominator should be handled by partial fractions"));
        }

        const Rational two_a = two * denominator->quadratic;
        const ExprPtr affine = make_linear_polynomial(arena_, two_a, denominator->linear, var);

        ExprPtr reciprocal_primitive{};
        const ExprPtr sqrt_delta = make_sqrt(arena_, completed_discriminant);
        reciprocal_primitive = make_product(arena_, {
            make_binary(arena_, BinaryOp::Div, make_integer(arena_, 2), sqrt_delta),
            make_function(arena_, "arctan", {make_binary(arena_, BinaryOp::Div, affine, sqrt_delta)}),
        });

        primitive_terms.push_back(make_product(arena_, {
            make_rational(arena_, beta),
            reciprocal_primitive,
        }));
    }

    if (primitive_terms.empty()) {
        return ok(make_integer(arena_, 0));
    }
    return ok(make_sum(arena_, std::move(primitive_terms)));
}

}  // namespace cas::calculus::integrate_detail
