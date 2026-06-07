#include "cas/complex_rational.hpp"

namespace cas {

bool ComplexRational::is_unit() const noexcept {
    // Q[i] is a field, all non-zero elements are units.
    return !is_zero();
}

ComplexRational ComplexRational::operator-() const {
    return {-re_, -im_};
}

ComplexRational ComplexRational::operator+(const ComplexRational& o) const {
    return {re_ + o.re_, im_ + o.im_};
}

ComplexRational ComplexRational::operator-(const ComplexRational& o) const {
    return {re_ - o.re_, im_ - o.im_};
}

ComplexRational ComplexRational::operator*(const ComplexRational& o) const {
    // (a+bi)(c+di) = (ac-bd) + (ad+bc)i
    return {re_ * o.re_ - im_ * o.im_, re_ * o.im_ + im_ * o.re_};
}

Result<ComplexRational> ComplexRational::divide(const ComplexRational& o) const {
    if (o.is_zero()) {
        return fail<ComplexRational>(CASError{CASErrorKind::DivisionByZero, "Complex division by zero", std::nullopt});
    }
    // (a+bi)/(c+di) = (a+bi)(c-di) / (c^2+d^2)
    Rational den = o.norm_sq();
    ComplexRational num = (*this) * o.conjugate();
    return ok(ComplexRational(num.re_ / den, num.im_ / den));
}

ComplexRational ComplexRational::conjugate() const {
    return {re_, -im_};
}

Rational ComplexRational::norm_sq() const {
    return re_ * re_ + im_ * im_;
}

bool ComplexRational::operator==(const ComplexRational& o) const noexcept {
    return re_ == o.re_ && im_ == o.im_;
}

bool ComplexRational::operator!=(const ComplexRational& o) const noexcept {
    return !(*this == o);
}

ComplexRational ComplexRational::zero() { return {}; }
ComplexRational ComplexRational::one() { return {Rational(1), Rational(0)}; }
ComplexRational ComplexRational::imag_unit() { return {Rational(0), Rational(1)}; }

} // namespace cas
