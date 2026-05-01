#pragma once

#include "cas/rational.hpp"
#include "cas/result.hpp"
#include <vector>

namespace cas {
namespace algebra {

/**
 * @brief Represents an element in a simple algebraic extension Q(alpha).
 * 
 * The extension is defined by an irreducible minimal polynomial m(x) over Q.
 * Elements are represented as polynomials p(x) with deg(p) < deg(m).
 */
class AlgebraicNumber {
public:
    /// @brief Coefficients vector (from degree 0 up).
    using CoeffVec = std::vector<Rational>;

    /**
     * @brief Construct an algebraic number p(alpha).
     * @param value Coefficients of p(x).
     * @param min_poly Coefficients of m(x).
     */
    AlgebraicNumber(CoeffVec value, CoeffVec min_poly);

    [[nodiscard]] const CoeffVec& value() const noexcept { return value_; }
    [[nodiscard]] const CoeffVec& min_poly() const noexcept { return min_poly_; }

    [[nodiscard]] AlgebraicNumber operator+(const AlgebraicNumber& other) const;
    [[nodiscard]] AlgebraicNumber operator-(const AlgebraicNumber& other) const;
    [[nodiscard]] AlgebraicNumber operator*(const AlgebraicNumber& other) const;

    [[nodiscard]] bool operator==(const AlgebraicNumber& other) const;
    [[nodiscard]] bool operator!=(const AlgebraicNumber& other) const { return !(*this == other); }

    /**
     * @brief Multiplicative inverse using Extended GCD.
     * @return s(x) such that s(x)p(x) + t(x)m(x) = 1.
     */
    [[nodiscard]] Result<AlgebraicNumber> inverse() const;

private:
    CoeffVec value_;
    CoeffVec min_poly_;
};

} // namespace algebra
} // namespace cas
