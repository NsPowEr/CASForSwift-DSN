// CAS-F1.6-NEW — ComplexRational: exact Q[i] arithmetic.
//
// Represents a + b·i with a, b ∈ Q (Rational). All four arithmetic
// operations are exact over Q[i] — no floating-point, no DecimalLit.
//
// This is the AST-free layer: it operates purely on Rational values and
// is used by the symbolic ComplexLit layer (F1.6 ComplexLit AST node).
//
// Reference: Knuth TAOCP Vol.2 §4.6.4 (complex arithmetic over exact fields).

#pragma once

#include "cas/rational.hpp"
#include "cas/result.hpp"

namespace cas {

class ComplexRational {
public:
    static constexpr std::uint32_t API_VERSION = 1;

    ComplexRational() = default;
    explicit ComplexRational(Rational real)
        : re_(std::move(real)), im_(Rational{}) {}
    ComplexRational(Rational real, Rational imag)
        : re_(std::move(real)), im_(std::move(imag)) {}

    [[nodiscard]] const Rational& real() const noexcept { return re_; }
    [[nodiscard]] const Rational& imag() const noexcept { return im_; }

    [[nodiscard]] bool is_real() const noexcept {
        return im_ == Rational{};
    }
    [[nodiscard]] bool is_zero() const noexcept {
        return re_ == Rational{} && im_ == Rational{};
    }
    [[nodiscard]] bool is_unit() const noexcept;

    // Arithmetic — F1.6-NEW: add/sub/mul/div exact over Q[i].
    [[nodiscard]] ComplexRational operator-() const;
    [[nodiscard]] ComplexRational operator+(const ComplexRational& o) const;
    [[nodiscard]] ComplexRational operator-(const ComplexRational& o) const;
    [[nodiscard]] ComplexRational operator*(const ComplexRational& o) const;

    // Division: multiply by conjugate of denominator.
    // Returns error on division by zero.
    [[nodiscard]] Result<ComplexRational> divide(const ComplexRational& o) const;

    // Conjugate: a - b·i
    [[nodiscard]] ComplexRational conjugate() const;

    // Norm squared: a² + b² ∈ Q (always ≥ 0).
    // abs²(z) = re² + im²  — returned as Rational (exact).
    [[nodiscard]] Rational norm_sq() const;

    [[nodiscard]] bool operator==(const ComplexRational& o) const noexcept;
    [[nodiscard]] bool operator!=(const ComplexRational& o) const noexcept;

    // Unit constants.
    [[nodiscard]] static ComplexRational zero();
    [[nodiscard]] static ComplexRational one();
    [[nodiscard]] static ComplexRational imag_unit();  // i

private:
    Rational re_;
    Rational im_;
};

}  // namespace cas
