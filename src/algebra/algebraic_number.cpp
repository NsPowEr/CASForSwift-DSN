#include "cas/algebraic_number.hpp"
#include "polynomial_internal.hpp"
#include "algebra_internal.hpp"
#include "cas/error.hpp"

namespace cas {
namespace algebra {

AlgebraicNumber::AlgebraicNumber(CoeffVec value, CoeffVec min_poly)
    : value_(std::move(value)), min_poly_(std::move(min_poly)) {
    RatPoly v(value_);
    normalize_rational_coefficients(v);
    value_ = v.coefficients();

    RatPoly mp(min_poly_);
    normalize_rational_coefficients(mp);
    min_poly_ = mp.coefficients();
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

} // namespace algebra
} // namespace cas
