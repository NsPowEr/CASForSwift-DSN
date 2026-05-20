#pragma once

#include "cas/result.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace cas {

class BigInt {
public:
    static constexpr std::uint32_t API_VERSION = 1;

    BigInt() = default;
    BigInt(std::int64_t value);
    [[nodiscard]] static BigInt from_u64(std::uint64_t value) noexcept;
    [[nodiscard]] static Result<BigInt> parse(std::string decimal);

    [[nodiscard]] const std::string& decimal() const;
    [[nodiscard]] bool is_negative() const noexcept;
    [[nodiscard]] bool is_zero() const noexcept;
    [[nodiscard]] bool is_positive() const noexcept;
    [[nodiscard]] std::size_t hash() const noexcept;

    [[nodiscard]] BigInt operator-() const;
    [[nodiscard]] BigInt abs() const;

    friend bool operator==(const BigInt& lhs, const BigInt& rhs) noexcept;
    friend bool operator!=(const BigInt& lhs, const BigInt& rhs) noexcept;
    friend bool operator<(const BigInt& lhs, const BigInt& rhs) noexcept;
    friend bool operator<=(const BigInt& lhs, const BigInt& rhs) noexcept;
    friend bool operator>(const BigInt& lhs, const BigInt& rhs) noexcept;
    friend bool operator>=(const BigInt& lhs, const BigInt& rhs) noexcept;
    friend Result<std::pair<BigInt, BigInt>> checked_divide_with_remainder(
        const BigInt& dividend,
        const BigInt& divisor);

    friend BigInt operator+(const BigInt& lhs, const BigInt& rhs);
    friend BigInt operator-(const BigInt& lhs, const BigInt& rhs);
    friend BigInt operator*(const BigInt& lhs, const BigInt& rhs);
    friend BigInt operator/(const BigInt& lhs, const BigInt& rhs);
    friend BigInt operator%(const BigInt& lhs, const BigInt& rhs);

    BigInt& operator+=(const BigInt& rhs);
    BigInt& operator-=(const BigInt& rhs);
    BigInt& operator*=(const BigInt& rhs);
    BigInt& operator/=(const BigInt& rhs);
    BigInt& operator%=(const BigInt& rhs);

    [[nodiscard]] std::size_t bit_length() const noexcept;
    [[nodiscard]] BigInt shift_left_bits(std::size_t bits) const;
    [[nodiscard]] BigInt shift_right_bits(std::size_t bits) const;
    [[nodiscard]] std::uint64_t to_u64() const noexcept;
    [[nodiscard]] double to_double() const noexcept;

private:
    [[nodiscard]] Result<void> assign_decimal_checked(std::string decimal);
    void multiply_by_small(std::uint32_t value) noexcept;
    void add_small(std::uint32_t value) noexcept;
    [[nodiscard]] std::uint32_t divide_by_small(std::uint32_t value) noexcept;
    void normalize() noexcept;
    void invalidate_decimal_cache() noexcept;
    void rebuild_decimal_cache() const;

    [[nodiscard]] static BigInt from_parts(std::vector<std::uint32_t> limbs, bool negative) noexcept;
    [[nodiscard]] static int compare_magnitude(const BigInt& lhs, const BigInt& rhs) noexcept;
    [[nodiscard]] static BigInt add_magnitude(const BigInt& lhs, const BigInt& rhs);
    [[nodiscard]] static BigInt subtract_magnitude(const BigInt& lhs, const BigInt& rhs);
    [[nodiscard]] static BigInt multiply_magnitude(const BigInt& lhs, const BigInt& rhs);
    [[nodiscard]] static std::pair<BigInt, BigInt> divide_magnitude(const BigInt& dividend, const BigInt& divisor);
    [[nodiscard]] static std::pair<BigInt, BigInt> divide_with_remainder(const BigInt& dividend, const BigInt& divisor);

    std::vector<std::uint32_t> limbs_;
    bool negative_{false};
    mutable std::string decimal_cache_{"0"};
    mutable bool decimal_cache_valid_{true};
};

[[nodiscard]] BigInt gcd(BigInt lhs, BigInt rhs);
[[nodiscard]] Result<std::pair<BigInt, BigInt>> checked_divide_with_remainder(
    const BigInt& dividend,
    const BigInt& divisor);
[[nodiscard]] Result<BigInt> checked_divide(const BigInt& lhs, const BigInt& rhs);
[[nodiscard]] Result<BigInt> checked_mod(const BigInt& lhs, const BigInt& rhs);

}  // namespace cas
