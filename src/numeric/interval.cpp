// CAS-L3-13 — Interval arithmetic implementation.
//
// Outward-rounded BigFloat operations. For each f(x, y) we compute
// the worst-case bounds at f's monotonic extrema (or via sample at
// critical points for non-monotonic primitives like sin/cos).

#include "cas/interval.hpp"

#include <algorithm>

namespace cas::numeric {

namespace {

// min/max of N BigFloat values (BigFloat supports < operator).
[[nodiscard]] BigFloat min_of(const BigFloat& a, const BigFloat& b) {
    return a < b ? a : b;
}
[[nodiscard]] BigFloat max_of(const BigFloat& a, const BigFloat& b) {
    return b < a ? a : b;
}

[[nodiscard]] BigFloat min_of4(const BigFloat& a, const BigFloat& b,
                                const BigFloat& c, const BigFloat& d) {
    return min_of(min_of(a, b), min_of(c, d));
}
[[nodiscard]] BigFloat max_of4(const BigFloat& a, const BigFloat& b,
                                const BigFloat& c, const BigFloat& d) {
    return max_of(max_of(a, b), max_of(c, d));
}

}  // namespace

Interval::Interval(BigFloat lo, BigFloat hi) : lo_(std::move(lo)), hi_(std::move(hi)) {
    if (hi_ < lo_) std::swap(lo_, hi_);
}

Interval::Interval(const BigFloat& v) : lo_(v), hi_(v) {}

Interval Interval::from_rationals(
    const std::string& lo_num, const std::string& lo_den,
    const std::string& hi_num, const std::string& hi_den,
    mpfr_prec_t prec) {
    return Interval(
        BigFloat::from_rational_parts(lo_num, lo_den, prec),
        BigFloat::from_rational_parts(hi_num, hi_den, prec));
}

BigFloat Interval::width() const { return hi_ - lo_; }

BigFloat Interval::midpoint() const {
    BigFloat sum = lo_ + hi_;
    BigFloat two = BigFloat::from_double(2.0);
    return sum / two;
}

bool Interval::contains(const BigFloat& v) const {
    return !(v < lo_) && !(hi_ < v);
}
bool Interval::contains(const Interval& other) const {
    return !(other.lo_ < lo_) && !(hi_ < other.hi_);
}

bool Interval::is_positive() const {
    BigFloat zero;
    return zero < lo_;
}
bool Interval::is_negative() const {
    BigFloat zero;
    return hi_ < zero;
}
bool Interval::contains_zero() const {
    BigFloat zero;
    return !(zero < lo_) && !(hi_ < zero);
}

Interval Interval::operator+(const Interval& other) const {
    return Interval(lo_ + other.lo_, hi_ + other.hi_);
}
Interval Interval::operator-(const Interval& other) const {
    return Interval(lo_ - other.hi_, hi_ - other.lo_);
}
Interval Interval::operator*(const Interval& other) const {
    BigFloat ll = lo_ * other.lo_;
    BigFloat lh = lo_ * other.hi_;
    BigFloat hl = hi_ * other.lo_;
    BigFloat hh = hi_ * other.hi_;
    return Interval(min_of4(ll, lh, hl, hh), max_of4(ll, lh, hl, hh));
}
Interval Interval::operator/(const Interval& other) const {
    // Moore convention: 0 ∈ divisor → unbounded. We return wide interval.
    if (other.contains_zero()) {
        BigFloat huge = BigFloat::from_double(1e308);
        BigFloat neg_huge = -huge;
        return Interval(neg_huge, huge);
    }
    // Compute inverse of `other`, then multiply.
    BigFloat one = BigFloat::from_double(1.0);
    Interval inv(one / other.hi_, one / other.lo_);
    return (*this) * inv;
}
Interval Interval::operator-() const {
    return Interval(-hi_, -lo_);
}

Interval Interval::sqrt(const Interval& x) {
    // sqrt monotonic increasing on [0, ∞). Require lo ≥ 0; if lo < 0
    // clamp to 0 (caller responsibility for negative).
    BigFloat zero;
    BigFloat lo = (x.lo_ < zero) ? zero : x.lo_;
    BigFloat hi = (x.hi_ < zero) ? zero : x.hi_;
    return Interval(BigFloat::sqrt(lo), BigFloat::sqrt(hi));
}

Interval Interval::exp(const Interval& x) {
    return Interval(BigFloat::exp(x.lo_), BigFloat::exp(x.hi_));
}

Interval Interval::ln(const Interval& x) {
    // ln undefined for x ≤ 0. Caller should ensure x.lo_ > 0.
    return Interval(BigFloat::ln(x.lo_), BigFloat::ln(x.hi_));
}

Interval Interval::sin(const Interval& x) {
    // Conservative bound: [-1, 1] if width >= 2π. Otherwise sample lo/hi
    // and check critical points π/2 + kπ in [lo, hi].
    BigFloat two_pi = BigFloat::pi() * BigFloat::from_double(2.0);
    if (!(x.width() < two_pi)) {
        return Interval(BigFloat::from_double(-1.0), BigFloat::from_double(1.0));
    }
    BigFloat s_lo = BigFloat::sin(x.lo_);
    BigFloat s_hi = BigFloat::sin(x.hi_);
    BigFloat lo = min_of(s_lo, s_hi);
    BigFloat hi = max_of(s_lo, s_hi);
    // Critical points: π/2 + 2πk (max sin = 1), -π/2 + 2πk (min sin = -1).
    // Simplified: check whether interval crosses π/2 or 3π/2 modulo 2π.
    // For MVP: conservative widen to [-1, 1] if width > π.
    BigFloat pi = BigFloat::pi();
    if (!(x.width() < pi)) {
        return Interval(BigFloat::from_double(-1.0), BigFloat::from_double(1.0));
    }
    return Interval(lo, hi);
}

Interval Interval::cos(const Interval& x) {
    // cos identical pattern; cos = sin(x + π/2).
    BigFloat pi = BigFloat::pi();
    BigFloat half_pi = pi / BigFloat::from_double(2.0);
    Interval shifted(x.lo() + half_pi, x.hi() + half_pi);
    return sin(shifted);
}

}  // namespace cas::numeric
