// CAS-F1.6-NEW — ComplexRational arithmetic implementation.
//
// Q[i] is a field (Gaussian rationals). Every non-zero element is
// invertible: (a+bi)^{-1} = (a-bi)/(a²+b²).
//
// Multiplication: (a+bi)(c+di) = (ac-bd) + (ad+bc)i.
// Division:       (a+bi)/(c+di) = (a+bi)(c-di)/(c²+d²).
//
// All operations use Rational arithmetic internally — no double/float.

#include "cas/complex_rational.hpp"
#include "cas/error.hpp"

namespace cas {

namespace {

[[nodiscard]] CASError make_div_zero_error() {
    return CASError{
        .kind    = CASErrorKind::Undefined,
        .message = "ComplexRational: division by zero",
        .hint    = std::nullopt,
    };
}

}  // namespace

// F1.6-NEW: unary negation.
ComplexRational ComplexRational::operator-() const {
    return ComplexRational(-re_, -im_);
}

// F1.6-NEW: Q[i] addition/subtraction.
ComplexRational ComplexRational::operator+(const ComplexRational& o) const {
    return ComplexRational(re_ + o.re_, im_ + o.im_);
}

ComplexRational ComplexRational::operator-(const ComplexRational& o) const {
    return ComplexRational(re_ - o.re_, im_ - o.im_);
}

// F1.6-NEW: Q[i] multiplication: (a+bi)(c+di) = (ac-bd) + (ad+bc)i.
ComplexRational ComplexRational::operator*(const ComplexRational& o) const {
    Rational ac = re_ * o.re_;
    Rational bd = im_ * o.im_;
    Rational ad = re_ * o.im_;
    Rational bc = im_ * o.re_;
    return ComplexRational(ac - bd, ad + bc);
}

// F1.6-NEW: Q[i] division via conjugate denominator.
// (a+bi)/(c+di) = (a+bi)(c-di)/(c²+d²).
Result<ComplexRational> ComplexRational::divide(const ComplexRational& o) const {
    Rational denom = o.norm_sq();
    if (denom == Rational{}) {
        return fail<ComplexRational>(make_div_zero_error());
    }
    ComplexRational num = *this * o.conjugate();
    return ok(ComplexRational(num.re_ / denom, num.im_ / denom));
}

// F1.6-NEW: conjugate a - bi.
ComplexRational ComplexRational::conjugate() const {
    return ComplexRational(re_, -im_);
}

// F1.6-NEW: |z|² = re² + im² ∈ Q.
Rational ComplexRational::norm_sq() const {
    return re_ * re_ + im_ * im_;
}

// F1.6-NEW: unit detection (norm_sq == 1, which for Q[i] means only ±1 and ±i).
bool ComplexRational::is_unit() const noexcept {
    return norm_sq() == Rational(BigInt(1));
}

bool ComplexRational::operator==(const ComplexRational& o) const noexcept {
    return re_ == o.re_ && im_ == o.im_;
}

bool ComplexRational::operator!=(const ComplexRational& o) const noexcept {
    return !(*this == o);
}

// F1.6-NEW: named constants.
ComplexRational ComplexRational::zero() {
    return ComplexRational{};
}

ComplexRational ComplexRational::one() {
    return ComplexRational(Rational(BigInt(1)));
}

ComplexRational ComplexRational::imag_unit() {
    return ComplexRational(Rational{}, Rational(BigInt(1)));
}

}  // namespace cas
