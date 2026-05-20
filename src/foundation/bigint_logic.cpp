#include "cas/bigint.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <vector>

namespace cas {

bool BigInt::is_negative() const noexcept {
    return negative_;
}

bool BigInt::is_zero() const noexcept {
    return limbs_.empty();
}

bool BigInt::is_positive() const noexcept {
    return !negative_ && !limbs_.empty();
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

void BigInt::normalize() noexcept {
    while (!limbs_.empty() && limbs_.back() == 0U) {
        limbs_.pop_back();
    }
    if (limbs_.empty()) {
        negative_ = false;
    }
    invalidate_decimal_cache();
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

} // namespace cas
