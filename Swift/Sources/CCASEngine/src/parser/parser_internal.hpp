#pragma once

#include "cas/ast.hpp"
#include "cas/token.hpp"

#include <limits>
#include <optional>
#include <string>
#include <string_view>

namespace cas::parser_internal {

[[nodiscard]] std::optional<MathConstant> parse_math_constant(std::string_view text);
[[nodiscard]] bool is_builtin_function(std::string_view name);
[[nodiscard]] std::string token_text_or_name(const Token& token);
[[nodiscard]] Token eof_token(std::size_t line = 1U, std::size_t column = 1U) noexcept;

template <typename UInt>
[[nodiscard]] std::optional<UInt> parse_unsigned_decimal(const BigInt& value) {
    if (value.is_negative()) {
        return std::nullopt;
    }

    UInt parsed = 0;
    constexpr UInt max_value = std::numeric_limits<UInt>::max();
    for (const char ch : value.decimal()) {
        if (ch < '0' || ch > '9') {
            return std::nullopt;
        }

        const UInt digit = static_cast<UInt>(ch - '0');
        if (parsed > (max_value - digit) / static_cast<UInt>(10)) {
            return std::nullopt;
        }
        parsed = static_cast<UInt>(parsed * static_cast<UInt>(10) + digit);
    }

    return parsed;
}

}  // namespace cas::parser_internal
