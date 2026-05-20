// CAS-L2-24 — Z[i] arithmetic implementation.

#include "cas/gaussian_int.hpp"

namespace cas {

namespace {

// Round BigInt rational p/q to nearest integer (half away from zero).
[[nodiscard]] BigInt round_div(const BigInt& p, const BigInt& q) {
    BigInt half = q.abs() / BigInt(2);
    BigInt p_abs = p.abs();
    BigInt magnitude = (p_abs + half) / q.abs();
    // Sign follows sign(p) · sign(q).
    bool neg = (p.is_negative() != q.is_negative());
    if (neg) return -magnitude;
    return magnitude;
}

}  // namespace

GaussianDivision gaussian_divmod(const GaussianInt& a, const GaussianInt& b) {
    // α / β = α · conj(β) / N(β); N(β) is real BigInt.
    // Numerator: a · conj(b) = (ac + bd) + (bc − ad)i  where a+bi = α, c+di = β.
    BigInt num_real = a.real() * b.real() + a.imag() * b.imag();
    BigInt num_imag = a.imag() * b.real() - a.real() * b.imag();
    BigInt denom = b.norm();  // c² + d², always positive (if β ≠ 0)
    GaussianInt q(round_div(num_real, denom), round_div(num_imag, denom));
    GaussianInt r = a - q * b;
    return {q, r};
}

GaussianInt gaussian_gcd(GaussianInt a, GaussianInt b) {
    // Make argument with larger norm the dividend.
    if (a.norm() < b.norm()) std::swap(a, b);
    while (!b.is_zero()) {
        auto dm = gaussian_divmod(a, b);
        a = b;
        b = dm.remainder;
    }
    // Canonicalize: multiply by unit so result has positive real, or
    // (real = 0) ∧ (positive imag).
    if (a.real().is_negative()
        || (a.real().is_zero() && a.imag().is_negative())) {
        a = -a;
    }
    if (a.real().is_zero() && a.imag().is_negative()) a = -a;
    // Resolve i-multiplication if real == 0 with positive imag we keep;
    // if imag == 0 with negative real we already negated.
    return a;
}

}  // namespace cas
