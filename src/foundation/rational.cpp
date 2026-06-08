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

Rational::Rational() : numerator_(0), denominator_(1) {}

Rational::Rational(BigInt numerator) : numerator_(std::move(numerator)), denominator_(1) {}

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

BigInt Rational::round() const {
    if (is_integer()) return numerator_;
    
    BigInt q = numerator_ / denominator_;
    BigInt r = numerator_ % denominator_;
    
    BigInt abs_r2 = r.abs() * BigInt(2);
    BigInt abs_d = denominator_.abs();
    
    if (abs_r2 > abs_d || (abs_r2 == abs_d && !q.is_zero())) {
        if (numerator_.is_negative() == denominator_.is_negative()) {
            return q + BigInt(1);
        } else {
            return q - BigInt(1);
        }
    }
    return q;
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

bool operator<(const Rational& lhs, const Rational& rhs) {
    return lhs.numerator() * rhs.denominator() < rhs.numerator() * lhs.denominator();
}

bool operator<=(const Rational& lhs, const Rational& rhs) {
    return lhs.numerator() * rhs.denominator() <= rhs.numerator() * lhs.denominator();
}

bool operator>(const Rational& lhs, const Rational& rhs) {
    return lhs.numerator() * rhs.denominator() > rhs.numerator() * lhs.denominator();
}

bool operator>=(const Rational& lhs, const Rational& rhs) {
    return lhs.numerator() * rhs.denominator() >= rhs.numerator() * lhs.denominator();
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

// F1.2-NEW: Euclidean continued-fraction expansion.
// Algorithm (Knuth TAOCP Vol.2 §4.5.3):
//   a0 = floor(p/q), then iterate (p,q) ← (q, p mod q) until q=0 or limit.
// Negative case: standard convention keeps all a_k (k>=1) positive; only
// a0 may be negative (i.e. a0 = floor(rational), which for negative non-
// integers differs from truncation). We implement this via the invariant that
// the remainder at each step is always non-negative.
std::vector<BigInt> Rational::to_continued_fraction(std::size_t n_max) const {
    std::vector<BigInt> quotients;
    quotients.reserve(32U);

    // Work with |numerator_| and |denominator_|, track sign separately.
    // For floor division with negative rationals: floor(-7/3) = -3 (not -2),
    // so we compute floor(p/q) = (p - (q-1)) / q when p is negative and q>0.
    // Simpler: always use the Euclidean floor via BigInt arithmetic.
    BigInt p = numerator_;
    BigInt q = denominator_;   // Always positive after reduce().

    std::size_t count = 0U;
    while (!q.is_zero() && count < n_max) {
        // floor_div: for p>=0 it is p/q; for p<0 it is -(|p|+q-1)/q.
        BigInt a0;
        BigInt r;
        if (!p.is_negative()) {
            a0 = p / q;
            r  = p % q;
        } else {
            // p negative: floor(p/q) = -(|p| + q - 1) / q  (integer ceiling of |p|/q, negated)
            BigInt abs_p = p.abs();
            a0 = -((abs_p + q - BigInt(1)) / q);
            // r = p - a0 * q  (must be in [0, q-1])
            r  = p - a0 * q;
        }
        quotients.push_back(std::move(a0));
        p = q;
        q = std::move(r);
        ++count;
    }

    return quotients;
}

Result<Rational> checked_divide(const Rational& lhs, const Rational& rhs) {
    if (rhs.numerator().is_zero()) {
        return fail<Rational>(make_error(CASErrorKind::Undefined, "Rational division by zero"));
    }
    return Rational::make(
        lhs.numerator() * rhs.denominator(),
        lhs.denominator() * rhs.numerator());
}

Rational double_to_rational_approx(double v) {
    if (v == 0.0) return Rational(BigInt(0));
    // Use a 60-bit dyadic approximation: v ≈ n / 2^60.
    constexpr int kBits = 60;
    const double scaled = std::ldexp(v, kBits);
    if (!std::isfinite(scaled)) {
        // Fall back to the integer part if the scaling overflows.
        return Rational(BigInt(static_cast<std::int64_t>(v)));
    }
    const auto n = static_cast<std::int64_t>(std::llround(scaled));
    // 2^60 as BigInt without resorting to BigInt::pow.
    BigInt den(1);
    for (int i = 0; i < kBits; ++i) den = den * BigInt(2);
    return Rational(BigInt(n), den);
}

}  // namespace cas
