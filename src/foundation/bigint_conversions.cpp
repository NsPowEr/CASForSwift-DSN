#include "cas/bigint.hpp"
#include "cas/error.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace cas {

namespace {

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

double BigInt::to_double() const noexcept {
    if (!decimal_cache_valid_) rebuild_decimal_cache();
    return std::stod(decimal_cache_);
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

    std::string result;
    if (negative_) {
        result += '-';
    }
    
    std::ostringstream out;
    out << chunks.back();
    for (auto it = chunks.rbegin() + 1; it != chunks.rend(); ++it) {
        out << std::setw(9) << std::setfill('0') << *it;
    }

    result += std::move(out).str();
    decimal_cache_ = std::move(result);
    decimal_cache_valid_ = true;
}

void BigInt::invalidate_decimal_cache() noexcept {
    decimal_cache_valid_ = false;
}

} // namespace cas
