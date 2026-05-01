#include "polynomial_internal.hpp"

#include "cas/numtheory.hpp"

#include <algorithm>
#include <utility>

namespace cas {
namespace algebra {

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

    if (dividend.empty() || dividend.degree() < divisor.degree()) {
        return dividend;
    }

    const std::size_t n = dividend.degree();
    const std::size_t m = divisor.degree();
    const std::size_t delta = n - m;
    const BigInt divisor_leading = divisor.leading_coeff();

    IntPoly remainder = dividend;
    std::size_t steps_taken = 0;

    while (!remainder.empty() && remainder.degree() >= m) {
        const std::size_t current_degree = remainder.degree();
        const std::size_t degree_gap = current_degree - m;
        const BigInt remainder_leading = remainder.leading_coeff();

        multiply_integer_coefficients_by_scalar(remainder, divisor_leading);
        steps_taken++;

        for (std::size_t index = 0; index < divisor.size(); ++index) {
            remainder[index + degree_gap] -= divisor[index] * remainder_leading;
        }
        normalize_integer_poly(remainder);
    }

    if (steps_taken <= delta) {
        multiply_integer_coefficients_by_scalar(remainder, bigint_pow_nonnegative(divisor_leading, delta + 1 - steps_taken));
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

Rational evaluate_integer_polynomial_at(const IntPoly& coefficients, const Rational& value) {
    return evaluate_integer_polynomial_at_impl(coefficients, value);
}

Rational evaluate_rational_polynomial_at(const RatPoly& coefficients, const Rational& value) {
    return evaluate_rational_polynomial_at_impl(coefficients, value);
}

}  // namespace algebra
}  // namespace cas
