#include "cas/bigint.hpp"
#include "cas/error.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas {

namespace {

constexpr std::uint64_t kLimbBase = static_cast<std::uint64_t>(1) << 32U;

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

}  // namespace

BigInt BigInt::operator-() const {
    if (is_zero()) {
        return *this;
    }

    BigInt value = *this;
    value.negative_ = !negative_;
    value.invalidate_decimal_cache();
    return value;
}

BigInt BigInt::abs() const {
    BigInt value = *this;
    if (value.negative_) {
        value.negative_ = false;
        value.invalidate_decimal_cache();
    }
    return value;
}

BigInt& BigInt::operator+=(const BigInt& rhs) {
    *this = *this + rhs;
    return *this;
}

BigInt& BigInt::operator-=(const BigInt& rhs) {
    *this = *this - rhs;
    return *this;
}

BigInt& BigInt::operator*=(const BigInt& rhs) {
    *this = *this * rhs;
    return *this;
}

BigInt& BigInt::operator/=(const BigInt& rhs) {
    assert(!rhs.is_zero());
    *this = divide_with_remainder(*this, rhs).first;
    return *this;
}

BigInt& BigInt::operator%=(const BigInt& rhs) {
    assert(!rhs.is_zero());
    *this = divide_with_remainder(*this, rhs).second;
    return *this;
}

void BigInt::multiply_by_small(std::uint32_t value) noexcept {
    if (value == 0U || is_zero()) {
        if (value == 0U) {
            limbs_.clear();
            negative_ = false;
            invalidate_decimal_cache();
        }
        return;
    }

    std::uint64_t carry = 0U;
    for (std::uint32_t& limb : limbs_) {
        const std::uint64_t product = static_cast<std::uint64_t>(limb) * value + carry;
        limb = static_cast<std::uint32_t>(product);
        carry = product >> 32U;
    }

    if (carry != 0U) {
        limbs_.push_back(static_cast<std::uint32_t>(carry));
    }
    invalidate_decimal_cache();
}

void BigInt::add_small(std::uint32_t value) noexcept {
    std::uint64_t carry = value;
    std::size_t index = 0U;
    while (carry != 0U) {
        if (index == limbs_.size()) {
            limbs_.push_back(0U);
        }

        const std::uint64_t sum = static_cast<std::uint64_t>(limbs_[index]) + carry;
        limbs_[index] = static_cast<std::uint32_t>(sum);
        carry = sum >> 32U;
        ++index;
    }
    invalidate_decimal_cache();
}

std::uint32_t BigInt::divide_by_small(std::uint32_t value) noexcept {
    std::uint64_t remainder = 0U;
    for (std::size_t index = limbs_.size(); index > 0U; --index) {
        const std::uint64_t current =
            (remainder << 32U) | static_cast<std::uint64_t>(limbs_[index - 1U]);
        limbs_[index - 1U] = static_cast<std::uint32_t>(current / value);
        remainder = current % value;
    }
    normalize();
    return static_cast<std::uint32_t>(remainder);
}

BigInt BigInt::from_parts(std::vector<std::uint32_t> limbs, bool negative) noexcept {
    BigInt value;
    value.limbs_ = std::move(limbs);
    value.negative_ = negative;
    value.normalize();
    return value;
}

BigInt BigInt::add_magnitude(const BigInt& lhs, const BigInt& rhs) {
    std::vector<std::uint32_t> result;
    result.reserve(std::max(lhs.limbs_.size(), rhs.limbs_.size()) + 1U);

    const std::size_t count = std::max(lhs.limbs_.size(), rhs.limbs_.size());
    std::uint64_t carry = 0U;
    for (std::size_t index = 0U; index < count; ++index) {
        const std::uint64_t left = index < lhs.limbs_.size() ? lhs.limbs_[index] : 0U;
        const std::uint64_t right = index < rhs.limbs_.size() ? rhs.limbs_[index] : 0U;
        const std::uint64_t sum = left + right + carry;
        result.push_back(static_cast<std::uint32_t>(sum));
        carry = sum >> 32U;
    }

    if (carry != 0U) {
        result.push_back(static_cast<std::uint32_t>(carry));
    }
    return from_parts(std::move(result), false);
}

BigInt BigInt::subtract_magnitude(const BigInt& lhs, const BigInt& rhs) {
    std::vector<std::uint32_t> result;
    result.reserve(lhs.limbs_.size());

    std::int64_t borrow = 0;
    for (std::size_t index = 0U; index < lhs.limbs_.size(); ++index) {
        std::int64_t left = static_cast<std::int64_t>(lhs.limbs_[index]) - borrow;
        const std::int64_t right = index < rhs.limbs_.size()
            ? static_cast<std::int64_t>(rhs.limbs_[index])
            : 0;
        if (left < right) {
            left += static_cast<std::int64_t>(kLimbBase);
            borrow = 1;
        } else {
            borrow = 0;
        }
        result.push_back(static_cast<std::uint32_t>(left - right));
    }

    return from_parts(std::move(result), false);
}

// F7.0-A3.6: BigInt limb-allocation budget.
namespace {
thread_local std::size_t g_bigint_max_limbs = 0U;   // 0 = unlimited
thread_local bool        g_bigint_exhausted = false;
}  // namespace

void BigInt::set_max_limbs(std::size_t n) noexcept {
    g_bigint_max_limbs = n;
    g_bigint_exhausted = false;
}
std::size_t BigInt::max_limbs() noexcept { return g_bigint_max_limbs; }
bool BigInt::budget_exhausted() noexcept { return g_bigint_exhausted; }
void BigInt::clear_budget_exhausted() noexcept { g_bigint_exhausted = false; }

BigInt BigInt::multiply_magnitude(const BigInt& lhs, const BigInt& rhs) {
    if (lhs.is_zero() || rhs.is_zero()) {
        return BigInt(0);
    }

    // F7.0-A3.6: pre-flight budget check. result limbs ≤ lhs.size + rhs.size.
    if (g_bigint_max_limbs > 0U
        && lhs.limbs_.size() + rhs.limbs_.size() > g_bigint_max_limbs) {
        g_bigint_exhausted = true;
        return BigInt(0);
    }

    const std::size_t n = std::max(lhs.limbs_.size(), rhs.limbs_.size());

    // Schoolbook for small n: O(n²) but low overhead.
    // Threshold 32: benchmarked break-even for Karatsuba on x86-64.
    if (n < 32U) {
        std::vector<std::uint32_t> result(lhs.limbs_.size() + rhs.limbs_.size(), 0U);
        for (std::size_t i = 0U; i < lhs.limbs_.size(); ++i) {
            std::uint64_t carry = 0U;
            for (std::size_t j = 0U; j < rhs.limbs_.size(); ++j) {
                const std::size_t index = i + j;
                const std::uint64_t product =
                    static_cast<std::uint64_t>(lhs.limbs_[i]) * rhs.limbs_[j] +
                    result[index] +
                    carry;
                result[index] = static_cast<std::uint32_t>(product);
                carry = product >> 32U;
            }

            std::size_t index = i + rhs.limbs_.size();
            while (carry != 0U) {
                const std::uint64_t sum = static_cast<std::uint64_t>(result[index]) + carry;
                result[index] = static_cast<std::uint32_t>(sum);
                carry = sum >> 32U;
                ++index;
            }
        }
        return from_parts(std::move(result), false);
    }

    // Toom-3 for medium-large n: O(n^1.465).
    // Threshold 64: break-even vs Karatsuba (ref: Brent-Zimmermann §1.3.3).
    // kToom3MaxLimbs = 4096: beyond this, Toom-3 falls back to Karatsuba.
    // HARDCODE-OF-PASSAGE HPP-F1.1-MUL: Schönhage-Strassen for n>=4096 deferred.
    if (n >= 64U && n < 4096U) {
        return multiply_magnitude_toom3(lhs, rhs);
    }

    // Karatsuba for [32,64) and fallback for n>=4096.
    // O(n^1.585). Reference: Knuth TAOCP Vol 2 §4.3.3.A.
    const std::size_t half = (n + 1U) / 2U;
    const std::size_t half_bits = half * 32U;

    BigInt high1 = lhs.shift_right_bits(half_bits);
    BigInt low1 = subtract_magnitude(lhs, high1.shift_left_bits(half_bits));

    BigInt high2 = rhs.shift_right_bits(half_bits);
    BigInt low2 = subtract_magnitude(rhs, high2.shift_left_bits(half_bits));

    BigInt z0 = multiply_magnitude(low1, low2);
    BigInt z2 = multiply_magnitude(high1, high2);
    BigInt z1 = multiply_magnitude(add_magnitude(low1, high1), add_magnitude(low2, high2));

    z1 = subtract_magnitude(z1, add_magnitude(z2, z0));

    BigInt shifted_z2 = z2.shift_left_bits(2U * half_bits);
    BigInt shifted_z1 = z1.shift_left_bits(half_bits);

    return add_magnitude(add_magnitude(shifted_z2, shifted_z1), z0);
}

std::pair<BigInt, BigInt> BigInt::divide_magnitude(const BigInt& dividend, const BigInt& divisor) {
    assert(!divisor.is_zero());

    if (dividend.is_zero() || compare_magnitude(dividend, divisor) < 0) {
        return {BigInt(0), dividend};
    }

    // Single-limb divisor: use divide_by_small (fastest path).
    // divide_by_small modifies the object in-place, so we copy.
    if (divisor.limbs_.size() == 1U) {
        BigInt quotient = dividend;
        const std::uint32_t rem_u32 = quotient.divide_by_small(divisor.limbs_[0]);
        BigInt remainder = from_parts({rem_u32}, false);
        return {std::move(quotient), std::move(remainder)};
    }

    // Multi-limb divisor: Knuth Algorithm D for small/medium divisors,
    // Burnikel-Ziegler for large divisors.
    // - Knuth D (TAOCP Vol 2 §4.3.1) — O(m·n).
    // - Burnikel-Ziegler (1998) — O(M(n)·log n), beats Knuth D once n ≥ 64
    //   limbs (HPP-023 closure). Internal threshold owned by the BZ module.
    if (divisor.limb_count() >= 64U) {
        return divide_burnikel_ziegler(dividend, divisor);
    }
    return divide_knuth_d(dividend, divisor);
}

std::pair<BigInt, BigInt> BigInt::divide_with_remainder(const BigInt& dividend, const BigInt& divisor) {
    assert(!divisor.is_zero());

    if (dividend.is_zero()) {
        return {BigInt(0), BigInt(0)};
    }

    auto [quotient, remainder] = divide_magnitude(dividend.abs(), divisor.abs());
    const bool signs_differ = dividend.is_negative() != divisor.is_negative();

    if (signs_differ && !remainder.is_zero()) {
        quotient += BigInt(1);
        remainder = divisor.abs() - remainder;
    }

    quotient.negative_ = signs_differ && !quotient.is_zero();
    remainder.negative_ = divisor.is_negative() && !remainder.is_zero();
    return {std::move(quotient), std::move(remainder)};
}

BigInt operator+(const BigInt& lhs, const BigInt& rhs) {
    if (lhs.is_negative() == rhs.is_negative()) {
        BigInt result = BigInt::add_magnitude(lhs, rhs);
        result.negative_ = lhs.is_negative() && !result.is_zero();
        return result;
    }

    const int magnitude_cmp = BigInt::compare_magnitude(lhs, rhs);
    if (magnitude_cmp == 0) {
        return BigInt(0);
    }

    BigInt result = magnitude_cmp > 0
        ? BigInt::subtract_magnitude(lhs, rhs)
        : BigInt::subtract_magnitude(rhs, lhs);
    result.negative_ = magnitude_cmp > 0 ? lhs.is_negative() : rhs.is_negative();
    return result;
}

BigInt operator-(const BigInt& lhs, const BigInt& rhs) {
    return lhs + (-rhs);
}

BigInt operator*(const BigInt& lhs, const BigInt& rhs) {
    BigInt result = BigInt::multiply_magnitude(lhs, rhs);
    result.negative_ = lhs.is_negative() != rhs.is_negative() && !result.is_zero();
    return result;
}

BigInt operator/(const BigInt& lhs, const BigInt& rhs) {
    assert(!rhs.is_zero());
    return BigInt::divide_with_remainder(lhs, rhs).first;
}

BigInt operator%(const BigInt& lhs, const BigInt& rhs) {
    assert(!rhs.is_zero());
    return BigInt::divide_with_remainder(lhs, rhs).second;
}

// NOTE: gcd() is now defined in bigint_gcd_lehmer.cpp (F1.1 production algorithms).
// It dispatches: single-limb → Euclidean, [2,16) → Binary GCD, >=16 → Lehmer.
// The old Euclidean implementation here is removed to avoid duplicate symbol.

Result<std::pair<BigInt, BigInt>> checked_divide_with_remainder(
    const BigInt& dividend,
    const BigInt& divisor) {
    if (divisor.is_zero()) {
        return fail<std::pair<BigInt, BigInt>>(
            make_error(CASErrorKind::Undefined, "BigInt division by zero"));
    }
    return ok(BigInt::divide_with_remainder(dividend, divisor));
}

Result<BigInt> checked_divide(const BigInt& lhs, const BigInt& rhs) {
    auto division = checked_divide_with_remainder(lhs, rhs);
    if (division.is_error()) {
        return fail<BigInt>(division.error());
    }
    return ok(std::move(division.value().first));
}

Result<BigInt> checked_mod(const BigInt& lhs, const BigInt& rhs) {
    auto division = checked_divide_with_remainder(lhs, rhs);
    if (division.is_error()) {
        return fail<BigInt>(division.error());
    }
    return ok(std::move(division.value().second));
}

} // namespace cas
