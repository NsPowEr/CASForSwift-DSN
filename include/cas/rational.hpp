#pragma once

#include "cas/bigint.hpp"
#include "cas/result.hpp"

#include <cstdint>

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
