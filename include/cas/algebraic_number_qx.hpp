#pragma once

#include "cas/rational.hpp"
#include "cas/result.hpp"
#include <vector>
#include <memory>

namespace cas {
namespace algebra {

/**
 * @brief Represents an element in a simple algebraic extension Q(x)[alpha] over the field of rational functions Q(x).
 * 
 * The extension is defined by alpha^2 = r(x), where r(x) is a rational function in Q(x).
 * Elements are represented as a(x) + b(x)*alpha, where a(x), b(x) are rational functions in Q(x).
 */
class AlgebraicNumberQx {
public:
    using CoeffVec = std::vector<Rational>;

    /**
     * @brief Construct an algebraic function a(x) + b(x)*alpha, with alpha^2 = r(x).
     * All rational functions are given by their numerator and denominator coefficients.
     */
    AlgebraicNumberQx(CoeffVec a_num, CoeffVec a_den,
                      CoeffVec b_num, CoeffVec b_den,
                      CoeffVec r_num, CoeffVec r_den);

    // Convenience constructor where denominators default to 1.
    AlgebraicNumberQx(CoeffVec a, CoeffVec b, CoeffVec r);

    [[nodiscard]] CoeffVec a_num() const;
    [[nodiscard]] CoeffVec a_den() const;
    [[nodiscard]] CoeffVec b_num() const;
    [[nodiscard]] CoeffVec b_den() const;
    [[nodiscard]] CoeffVec r_num() const;
    [[nodiscard]] CoeffVec r_den() const;

    [[nodiscard]] bool is_zero() const noexcept;

    [[nodiscard]] AlgebraicNumberQx operator-() const;
    [[nodiscard]] AlgebraicNumberQx operator+(const AlgebraicNumberQx& other) const;
    [[nodiscard]] AlgebraicNumberQx operator-(const AlgebraicNumberQx& other) const;
    [[nodiscard]] AlgebraicNumberQx operator*(const AlgebraicNumberQx& other) const;

    [[nodiscard]] bool operator==(const AlgebraicNumberQx& other) const;
    [[nodiscard]] bool operator!=(const AlgebraicNumberQx& other) const { return !(*this == other); }

    [[nodiscard]] Result<AlgebraicNumberQx> inverse() const;
    [[nodiscard]] Result<AlgebraicNumberQx> div(const AlgebraicNumberQx& other) const;
    [[nodiscard]] Result<AlgebraicNumberQx> pow(std::size_t exponent) const;

private:
    struct Impl;
    std::shared_ptr<Impl> impl_;
};

} // namespace algebra
} // namespace cas
