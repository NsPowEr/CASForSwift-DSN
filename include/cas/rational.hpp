#pragma once

#include "cas/bigint.hpp"
#include "cas/result.hpp"

#include <cstdint>
#include <limits>
#include <vector>

namespace cas {

class Rational {
public:
    static constexpr std::uint32_t API_VERSION = 1;

    Rational();
    explicit Rational(BigInt numerator);
    Rational(BigInt numerator, BigInt denominator);

    [[nodiscard]] static Result<Rational> make(BigInt numerator, BigInt denominator = BigInt(1));

    [[nodiscard]] const BigInt& numerator() const noexcept;
    [[nodiscard]] const BigInt& denominator() const noexcept;

    [[nodiscard]] bool is_integer() const noexcept;
    [[nodiscard]] BigInt round() const;
    [[nodiscard]] Result<double> to_double_checked() const;
    [[nodiscard]] double to_double() const;

    // F1.2-NEW: continued-fraction expansion via the Euclidean algorithm.
    // Returns the partial quotients [a0; a1, a2, ..., ak] with at most
    // n_max terms (default unbounded = terminates when remainder is zero).
    // For a negative rational the sign is absorbed into a0 (a0 may be
    // negative; all subsequent quotients are positive — standard convention).
    // Property: the (n_max)-convergent reconstructed from the returned
    // list satisfies |rational - convergent| <= 1/(q_k * q_{k+1}) where
    // q_k is the k-th convergent denominator.
    [[nodiscard]] std::vector<BigInt> to_continued_fraction(
        std::size_t n_max = std::numeric_limits<std::size_t>::max()) const;

    [[nodiscard]] Rational operator-() const;

    [[nodiscard]] friend bool operator==(const Rational& lhs, const Rational& rhs) noexcept = default;
    friend bool operator<(const Rational& lhs, const Rational& rhs);
    friend bool operator<=(const Rational& lhs, const Rational& rhs);
    friend bool operator>(const Rational& lhs, const Rational& rhs);
    friend bool operator>=(const Rational& lhs, const Rational& rhs);

    friend Rational operator+(const Rational& lhs, const Rational& rhs);
    friend Rational operator-(const Rational& lhs, const Rational& rhs);
    friend Rational operator*(const Rational& lhs, const Rational& rhs);
    friend Rational operator/(const Rational& lhs, const Rational& rhs);

    Rational& operator+=(const Rational& rhs);
    Rational& operator-=(const Rational& rhs);
    Rational& operator*=(const Rational& rhs);
    Rational& operator/=(const Rational& rhs);

private:
    void reduce();

    BigInt numerator_;
    BigInt denominator_{1};
};

[[nodiscard]] Result<Rational> checked_divide(const Rational& lhs, const Rational& rhs);

}  // namespace cas
