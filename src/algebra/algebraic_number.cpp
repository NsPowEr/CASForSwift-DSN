#include "cas/algebraic_number.hpp"
#include "polynomial_internal.hpp"
#include "algebra_internal.hpp"
#include "cas/error.hpp"

namespace cas {
namespace algebra {

namespace {

[[nodiscard]] RatPoly normalize_monic_min_poly(AlgebraicNumber::CoeffVec coefficients) {
    RatPoly polynomial(std::move(coefficients));
    normalize_rational_coefficients(polynomial);
    if (!polynomial.empty()) {
        const Rational leading = polynomial.leading_coeff();
        if (!leading.numerator().is_zero() && leading != Rational(BigInt(1))) {
            for (auto& coefficient : polynomial.coefficients()) {
                coefficient = coefficient / leading;
            }
            normalize_rational_coefficients(polynomial);
        }
    }
    return polynomial;
}

[[nodiscard]] AlgebraicNumber::CoeffVec reduce_coefficients(
    AlgebraicNumber::CoeffVec value,
    const AlgebraicNumber::CoeffVec& min_poly) {
    RatPoly v(std::move(value));
    normalize_rational_coefficients(v);

    RatPoly mp(min_poly);
    normalize_rational_coefficients(mp);
    if (mp.empty() || mp.is_zero()) {
        return v.coefficients();
    }

    auto [quotient, remainder] = div_rem_rational_poly(v, mp);
    (void)quotient;
    normalize_rational_coefficients(remainder);
    return remainder.coefficients();
}

} // namespace

AlgebraicNumber::AlgebraicNumber(CoeffVec value, CoeffVec min_poly)
    : value_(std::move(value)), min_poly_(std::move(min_poly)) {
    RatPoly mp = normalize_monic_min_poly(std::move(min_poly_));
    min_poly_ = mp.coefficients();
    value_ = reduce_coefficients(std::move(value_), min_poly_);
}

bool AlgebraicNumber::is_zero() const noexcept {
    return value_.empty() || (value_.size() == 1U && value_.front().numerator().is_zero());
}

AlgebraicNumber AlgebraicNumber::operator-() const {
    CoeffVec negated = value_;
    for (auto& coefficient : negated) {
        coefficient = -coefficient;
    }
    return AlgebraicNumber(std::move(negated), min_poly_);
}

AlgebraicNumber AlgebraicNumber::operator+(const AlgebraicNumber& other) const {
    RatPoly a(value_);
    RatPoly b(other.value_);
    RatPoly sum = add_rational_poly(a, b);
    return AlgebraicNumber(sum.coefficients(), min_poly_);
}

AlgebraicNumber AlgebraicNumber::operator-(const AlgebraicNumber& other) const {
    RatPoly a(value_);
    RatPoly b(other.value_);
    RatPoly diff = sub_rational_poly(a, b);
    return AlgebraicNumber(diff.coefficients(), min_poly_);
}

AlgebraicNumber AlgebraicNumber::operator*(const AlgebraicNumber& other) const {
    RatPoly a(value_);
    RatPoly b(other.value_);
    RatPoly prod = mul_rational_poly(a, b);
    RatPoly mp(min_poly_);
    
    auto [q, rem] = div_rem_rational_poly(prod, mp);
    return AlgebraicNumber(rem.coefficients(), min_poly_);
}

bool AlgebraicNumber::operator==(const AlgebraicNumber& other) const {
    return value_ == other.value_ && min_poly_ == other.min_poly_;
}

Result<AlgebraicNumber> AlgebraicNumber::inverse() const {
    RatPoly a(value_);
    RatPoly mp(min_poly_);
    
    if (a.is_zero()) {
        return fail<AlgebraicNumber>(make_error(CASErrorKind::InvalidArgument, "Division by zero in AlgebraicNumber"));
    }

    auto [g, s, t] = extended_gcd_rational_poly(a, mp);
    
    // In Q(alpha), since m(alpha) is irreducible, any p(alpha) != 0 is invertible.
    // gcd(p, m) must be a constant (degree 0).
    if (g.degree() > 0 || g.is_zero()) {
        return fail<AlgebraicNumber>(make_error(CASErrorKind::InvalidArgument, "Element is not invertible (minimal polynomial might not be irreducible)"));
    }

    Rational g_const = g.constant_term();
    CoeffVec inv_coeffs = s.coefficients();
    for (auto& c : inv_coeffs) {
        c = c / g_const;
    }
    
    RatPoly inv_poly(inv_coeffs);
    auto [q, rem] = div_rem_rational_poly(inv_poly, mp);
    return ok(AlgebraicNumber(rem.coefficients(), min_poly_));
}

Result<AlgebraicNumber> AlgebraicNumber::div(const AlgebraicNumber& other) const {
    auto inverse_res = other.inverse();
    if (inverse_res.is_error()) {
        return fail<AlgebraicNumber>(inverse_res.error());
    }
    return ok((*this) * inverse_res.value());
}

Result<AlgebraicNumber> AlgebraicNumber::pow(std::size_t exponent) const {
    AlgebraicNumber result({Rational(BigInt(1))}, min_poly_);
    if (exponent == 0U) {
        return ok(result);
    }

    AlgebraicNumber base = *this;
    std::size_t power = exponent;
    while (power > 0U) {
        if ((power & 1U) != 0U) {
            result = result * base;
        }
        power >>= 1U;
        if (power > 0U) {
            base = base * base;
        }
    }
    return ok(result);
}

} // namespace algebra
} // namespace cas
