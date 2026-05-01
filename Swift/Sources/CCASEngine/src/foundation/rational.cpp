#include "cas/rational.hpp"

#include "cas/error.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace cas {
namespace {

struct DecimalApprox {
    long double significand{0.0L};
    int exponent10{0};
    bool zero{true};
};

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

[[nodiscard]] DecimalApprox approximate_decimal_magnitude(const BigInt& value) {
    if (value.is_zero()) {
        return DecimalApprox{};
    }

    constexpr std::size_t max_digits = 19U;
    const std::string& digits = value.decimal();
    const std::size_t consumed = std::min(max_digits, digits.size());

    long double significand = 0.0L;
    for (std::size_t index = 0; index < consumed; ++index) {
        significand = significand * 10.0L + static_cast<long double>(digits[index] - '0');
    }

    return DecimalApprox{
        .significand = significand,
        .exponent10 = static_cast<int>(digits.size() - consumed),
        .zero = false,
    };
}

}  // namespace

Rational::Rational() = default;

Rational::Rational(BigInt numerator) : numerator_(std::move(numerator)) {}

Rational::Rational(BigInt numerator, BigInt denominator)
    : numerator_(std::move(numerator)), denominator_(std::move(denominator)) {
    reduce();
}

Result<Rational> Rational::make(BigInt numerator, BigInt denominator) {
    if (denominator.is_zero()) {
        return fail<Rational>(make_error(CASErrorKind::Undefined, "Rational denominator must be non-zero"));
    }
    return ok(Rational(std::move(numerator), std::move(denominator)));
}

const BigInt& Rational::numerator() const noexcept {
    return numerator_;
}

const BigInt& Rational::denominator() const noexcept {
    return denominator_;
}

bool Rational::is_integer() const noexcept {
    return denominator_ == BigInt(1);
}

Result<double> Rational::to_double_checked() const {
    if (denominator_.is_zero()) {
        return fail<double>(make_error(CASErrorKind::Undefined, "Rational denominator must be non-zero"));
    }

    const DecimalApprox numerator = approximate_decimal_magnitude(numerator_.abs());
    if (numerator.zero) {
        return ok(0.0);
    }

    const DecimalApprox denominator = approximate_decimal_magnitude(denominator_.abs());
    const long double ratio = numerator.significand / denominator.significand;
    const int exponent = numerator.exponent10 - denominator.exponent10;
    long double value = ratio * std::pow(10.0L, static_cast<long double>(exponent));

    if (numerator_.is_negative() != denominator_.is_negative()) {
        value = -value;
    }

    if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<double>::max()) {
        return fail<double>(make_error(CASErrorKind::Overflow, "Rational value is outside double range"));
    }

    return ok(static_cast<double>(value));
}

double Rational::to_double() const {
    const auto converted = to_double_checked();
    if (converted.is_ok()) {
        return converted.value();
    }
    return numerator_.is_negative() ? -std::numeric_limits<double>::infinity()
                                    : std::numeric_limits<double>::infinity();
}

Rational Rational::operator-() const {
    return Rational(-numerator_, denominator_);
}

Rational& Rational::operator+=(const Rational& rhs) {
    *this = *this + rhs;
    return *this;
}

Rational& Rational::operator-=(const Rational& rhs) {
    *this = *this - rhs;
    return *this;
}

Rational& Rational::operator*=(const Rational& rhs) {
    *this = *this * rhs;
    return *this;
}

Rational& Rational::operator/=(const Rational& rhs) {
    assert(!rhs.numerator().is_zero());
    *this = Rational(
        numerator_ * rhs.denominator(),
        denominator_ * rhs.numerator());
    return *this;
}

void Rational::reduce() {
    assert(!denominator_.is_zero());

    if (numerator_.is_zero()) {
        denominator_ = BigInt(1);
        return;
    }

    if (denominator_.is_negative()) {
        numerator_ = -numerator_;
        denominator_ = -denominator_;
    }

    const BigInt factor = gcd(numerator_, denominator_);
    numerator_ /= factor;
    denominator_ /= factor;
}

Rational operator+(const Rational& lhs, const Rational& rhs) {
    return Rational(
        lhs.numerator() * rhs.denominator() + rhs.numerator() * lhs.denominator(),
        lhs.denominator() * rhs.denominator());
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
    return Rational(
        lhs.numerator() * rhs.denominator() - rhs.numerator() * lhs.denominator(),
        lhs.denominator() * rhs.denominator());
}

Rational operator*(const Rational& lhs, const Rational& rhs) {
    return Rational(
        lhs.numerator() * rhs.numerator(),
        lhs.denominator() * rhs.denominator());
}

Rational operator/(const Rational& lhs, const Rational& rhs) {
    assert(!rhs.numerator().is_zero());

    return Rational(
        lhs.numerator() * rhs.denominator(),
        lhs.denominator() * rhs.numerator());
}

Result<Rational> checked_divide(const Rational& lhs, const Rational& rhs) {
    if (rhs.numerator().is_zero()) {
        return fail<Rational>(make_error(CASErrorKind::Undefined, "Rational division by zero"));
    }
    return Rational::make(
        lhs.numerator() * rhs.denominator(),
        lhs.denominator() * rhs.numerator());
}

}  // namespace cas
