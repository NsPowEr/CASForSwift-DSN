#include "cas/bigint.hpp"

#include "cas/error.hpp"

#include <algorithm>
#include <cassert>
#include <bit>
#include <cctype>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
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

BigInt::BigInt(std::int64_t value) {
    std::uint64_t magnitude = 0;
    if (value < 0) {
        negative_ = true;
        magnitude = static_cast<std::uint64_t>(-(value + 1)) + 1U;
    } else {
        magnitude = static_cast<std::uint64_t>(value);
    }

    while (magnitude != 0U) {
        limbs_.push_back(static_cast<std::uint32_t>(magnitude & 0xFFFFFFFFU));
        magnitude >>= 32U;
    }
    normalize();
}

BigInt BigInt::from_u64(std::uint64_t value) noexcept {
    BigInt result;
    while (value != 0U) {
        result.limbs_.push_back(static_cast<std::uint32_t>(value & 0xFFFFFFFFU));
        value >>= 32U;
    }
    result.normalize();
    return result;
}

Result<BigInt> BigInt::parse(std::string decimal) {
    BigInt value;
    const auto assigned = value.assign_decimal_checked(std::move(decimal));
    if (assigned.is_error()) {
        return fail<BigInt>(assigned.error());
    }
    return ok(std::move(value));
}

const std::string& BigInt::decimal() const {
    if (!decimal_cache_valid_) {
        rebuild_decimal_cache();
    }
    return decimal_cache_;
}

bool BigInt::is_negative() const noexcept {
    return negative_;
}

bool BigInt::is_zero() const noexcept {
    return limbs_.empty();
}

std::size_t BigInt::hash() const noexcept {
    std::size_t seed = 0;
    auto hash_combine = [](std::size_t& s, std::uint32_t v) {
        s ^= static_cast<std::size_t>(v) + 0x9e3779b9 + (s << 6) + (s >> 2);
    };
    hash_combine(seed, static_cast<std::uint32_t>(negative_));
    for (const std::uint32_t limb : limbs_) {
        hash_combine(seed, limb);
    }
    return seed;
}

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

Result<void> BigInt::assign_decimal_checked(std::string decimal) {
    if (decimal.empty()) {
        return fail<void>(make_error(CASErrorKind::InvalidArgument, "BigInt decimal input must not be empty"));
    }

    bool negative = false;
    if (decimal.front() == '+' || decimal.front() == '-') {
        negative = (decimal.front() == '-');
        decimal.erase(decimal.begin());
    }

    if (decimal.empty()) {
        return fail<void>(make_error(CASErrorKind::InvalidArgument, "BigInt decimal input must contain digits"));
    }

    const auto invalid_char = std::find_if(decimal.begin(), decimal.end(), [](unsigned char ch) {
        return !std::isdigit(ch);
    });
    if (invalid_char != decimal.end()) {
        return fail<void>(make_error(
            CASErrorKind::InvalidArgument,
            "BigInt decimal input contains non-digit characters"));
    }

    limbs_.clear();
    negative_ = negative;
    invalidate_decimal_cache();

    for (const char ch : decimal) {
        multiply_by_small(10U);
        add_small(static_cast<std::uint32_t>(ch - '0'));
    }

    normalize();
    return ok();
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

void BigInt::normalize() noexcept {
    while (!limbs_.empty() && limbs_.back() == 0U) {
        limbs_.pop_back();
    }
    if (limbs_.empty()) {
        negative_ = false;
    }
    invalidate_decimal_cache();
}

void BigInt::invalidate_decimal_cache() noexcept {
    decimal_cache_valid_ = false;
}

void BigInt::rebuild_decimal_cache() const {
    if (is_zero()) {
        decimal_cache_ = "0";
        decimal_cache_valid_ = true;
        return;
    }

    constexpr std::uint32_t chunk_base = 1000000000U;
    BigInt value = abs();
    std::vector<std::uint32_t> chunks;
    while (!value.is_zero()) {
        chunks.push_back(value.divide_by_small(chunk_base));
    }

    std::ostringstream out;
    out << chunks.back();
    for (auto it = chunks.rbegin() + 1; it != chunks.rend(); ++it) {
        out << std::setw(9) << std::setfill('0') << *it;
    }

    decimal_cache_ = std::move(out).str();
    decimal_cache_valid_ = true;
}

std::size_t BigInt::bit_length() const noexcept {
    if (is_zero()) {
        return 0U;
    }

    const std::uint32_t top = limbs_.back();
    return (limbs_.size() - 1U) * 32U + (32U - static_cast<std::size_t>(std::countl_zero(top)));
}

BigInt BigInt::shift_left_bits(std::size_t bits) const {
    if (is_zero()) {
        return BigInt(0);
    }

    const std::size_t limb_shift = bits / 32U;
    const std::size_t bit_shift = bits % 32U;

    std::vector<std::uint32_t> shifted(limb_shift, 0U);
    shifted.reserve(limb_shift + limbs_.size() + 1U);

    std::uint64_t carry = 0U;
    for (const std::uint32_t limb : limbs_) {
        const std::uint64_t value = (static_cast<std::uint64_t>(limb) << bit_shift) | carry;
        shifted.push_back(static_cast<std::uint32_t>(value));
        carry = value >> 32U;
    }

    if (carry != 0U) {
        shifted.push_back(static_cast<std::uint32_t>(carry));
    }

    return from_parts(std::move(shifted), negative_);
}

BigInt BigInt::shift_right_bits(std::size_t bits) const {
    if (is_zero() || bits == 0) {
        return *this;
    }

    const std::size_t limb_shift = bits / 32U;
    const std::size_t bit_shift = bits % 32U;

    if (limb_shift >= limbs_.size()) {
        return BigInt(0);
    }

    std::vector<std::uint32_t> shifted;
    shifted.reserve(limbs_.size() - limb_shift);

    for (std::size_t i = limb_shift; i < limbs_.size(); ++i) {
        std::uint32_t current = limbs_[i];
        std::uint32_t next = (i + 1 < limbs_.size()) ? limbs_[i + 1] : 0U;

        std::uint32_t val = current >> bit_shift;
        if (bit_shift > 0 && i + 1 < limbs_.size()) {
            val |= (next << (32U - bit_shift));
        }
        shifted.push_back(val);
    }

    return from_parts(std::move(shifted), negative_);
}

std::uint64_t BigInt::to_u64() const noexcept {
    if (is_zero() || is_negative()) {
        return 0U;
    }

    std::uint64_t result = limbs_[0];
    if (limbs_.size() > 1) {
        result |= (static_cast<std::uint64_t>(limbs_[1]) << 32U);
    }
    return result;
}

BigInt BigInt::from_parts(std::vector<std::uint32_t> limbs, bool negative) noexcept {
    BigInt value;
    value.limbs_ = std::move(limbs);
    value.negative_ = negative;
    value.normalize();
    return value;
}

int BigInt::compare_magnitude(const BigInt& lhs, const BigInt& rhs) noexcept {
    if (lhs.limbs_.size() != rhs.limbs_.size()) {
        return lhs.limbs_.size() < rhs.limbs_.size() ? -1 : 1;
    }

    for (std::size_t index = lhs.limbs_.size(); index > 0U; --index) {
        const std::uint32_t left = lhs.limbs_[index - 1U];
        const std::uint32_t right = rhs.limbs_[index - 1U];
        if (left != right) {
            return left < right ? -1 : 1;
        }
    }
    return 0;
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

BigInt BigInt::multiply_magnitude(const BigInt& lhs, const BigInt& rhs) {
    if (lhs.is_zero() || rhs.is_zero()) {
        return BigInt(0);
    }

    const std::size_t n_max = std::max(lhs.limb_count(), rhs.limb_count());
    if (n_max >= 64U && n_max < 8192U) {
        return multiply_magnitude_toom3(lhs, rhs);
    }

    const std::size_t n = std::max(lhs.limbs_.size(), rhs.limbs_.size());
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
    if (divisor.limb_count() < 2U) {
        BigInt quotient = dividend;
        std::uint32_t rem_val = quotient.divide_by_small(divisor.limb_at(0));
        quotient.invalidate_decimal_cache();
        return {std::move(quotient), BigInt(rem_val)};
    }
    return divide_knuth_d(dividend, divisor);


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

bool operator==(const BigInt& lhs, const BigInt& rhs) noexcept {
    return lhs.negative_ == rhs.negative_ && lhs.limbs_ == rhs.limbs_;
}

bool operator!=(const BigInt& lhs, const BigInt& rhs) noexcept {
    return !(lhs == rhs);
}

bool operator<(const BigInt& lhs, const BigInt& rhs) noexcept {
    if (lhs.is_negative() != rhs.is_negative()) {
        return lhs.is_negative();
    }

    const int magnitude_cmp = BigInt::compare_magnitude(lhs, rhs);
    return lhs.is_negative() ? magnitude_cmp > 0 : magnitude_cmp < 0;
}

bool operator<=(const BigInt& lhs, const BigInt& rhs) noexcept {
    return !(rhs < lhs);
}

bool operator>(const BigInt& lhs, const BigInt& rhs) noexcept {
    return rhs < lhs;
}

bool operator>=(const BigInt& lhs, const BigInt& rhs) noexcept {
    return !(lhs < rhs);
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

BigInt gcd(BigInt lhs, BigInt rhs) {
    lhs = lhs.abs();
    rhs = rhs.abs();
    while (!rhs.is_zero()) {
        BigInt remainder = lhs % rhs;
        lhs = rhs;
        rhs = remainder.abs();
    }
    return lhs;
}

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

double BigInt::to_double() const noexcept {
    if (!decimal_cache_valid_) rebuild_decimal_cache();
    return std::stod(decimal_cache_);
}

}  // namespace cas
