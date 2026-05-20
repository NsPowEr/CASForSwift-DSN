// CAS-L2-24 — Gaussian integer Z[i] arithmetic.
//
// Implements the Euclidean ring Z[i] = {a + b·i : a, b ∈ Z}.
// Operations: +, -, *, ÷ (with quotient/remainder), norm, gcd.
// Z[i] is a Euclidean domain with N(a+bi) = a² + b², so unique
// factorization holds and Euclidean GCD terminates.
//
// Units of Z[i]: {1, -1, i, -i} — exactly the elements with norm 1.

#pragma once

#include "cas/bigint.hpp"

namespace cas {

class GaussianInt {
public:
    GaussianInt() noexcept = default;
    GaussianInt(BigInt real, BigInt imag) noexcept
        : real_(std::move(real)), imag_(std::move(imag)) {}
    explicit GaussianInt(BigInt real) noexcept
        : real_(std::move(real)), imag_(BigInt(0)) {}

    [[nodiscard]] const BigInt& real() const noexcept { return real_; }
    [[nodiscard]] const BigInt& imag() const noexcept { return imag_; }

    // Norm: N(a + b·i) = a² + b².  Multiplicative: N(αβ) = N(α)·N(β).
    [[nodiscard]] BigInt norm() const noexcept {
        return real_ * real_ + imag_ * imag_;
    }

    // Complex conjugate.
    [[nodiscard]] GaussianInt conjugate() const noexcept {
        return GaussianInt(real_, -imag_);
    }

    [[nodiscard]] bool is_zero() const noexcept {
        return real_.is_zero() && imag_.is_zero();
    }

    // Unit detection: α is a unit iff N(α) = 1, i.e. α ∈ {1, -1, i, -i}.
    [[nodiscard]] bool is_unit() const noexcept {
        return norm() == BigInt(1);
    }

    // Arithmetic
    [[nodiscard]] GaussianInt operator+(const GaussianInt& other) const noexcept {
        return GaussianInt(real_ + other.real_, imag_ + other.imag_);
    }
    [[nodiscard]] GaussianInt operator-(const GaussianInt& other) const noexcept {
        return GaussianInt(real_ - other.real_, imag_ - other.imag_);
    }
    [[nodiscard]] GaussianInt operator-() const noexcept {
        return GaussianInt(-real_, -imag_);
    }
    // (a+bi)(c+di) = (ac − bd) + (ad + bc)i
    [[nodiscard]] GaussianInt operator*(const GaussianInt& other) const noexcept {
        return GaussianInt(
            real_ * other.real_ - imag_ * other.imag_,
            real_ * other.imag_ + imag_ * other.real_);
    }

    [[nodiscard]] bool operator==(const GaussianInt& other) const noexcept {
        return real_ == other.real_ && imag_ == other.imag_;
    }
    [[nodiscard]] bool operator!=(const GaussianInt& other) const noexcept {
        return !(*this == other);
    }

private:
    BigInt real_;
    BigInt imag_;
};

// Euclidean division: returns (q, r) with α = q·β + r and N(r) < N(β).
// Quotient rounds each component to the nearest integer; remainder follows.
struct GaussianDivision {
    GaussianInt quotient;
    GaussianInt remainder;
};

[[nodiscard]] GaussianDivision gaussian_divmod(
    const GaussianInt& a, const GaussianInt& b);

// gcd via Euclidean algorithm. Result canonical up to unit multiplication
// (representative chosen with real > 0, or real = 0 and imag > 0).
[[nodiscard]] GaussianInt gaussian_gcd(GaussianInt a, GaussianInt b);

}  // namespace cas
