#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

namespace cas {

enum class TokenKind : std::uint8_t {
    Integer,
    Rational,
    Float,
    Identifier,

    Plus,
    Minus,
    Star,
    Slash,
    Caret,
    Percent,
    Bang,

    LParen,
    RParen,
    LBracket,
    RBracket,
    LBrace,
    RBrace,
    Comma,
    Semicolon,

    Integral,
    Derivative,
    Limit,
    Sum,
    Product,
    Sqrt,

    EndOfInput,
    Error,
};

struct SourceLocation {
    std::size_t line{1U};
    std::size_t column{1U};
};

struct Token {
    static constexpr std::uint32_t API_VERSION = 1;

    TokenKind kind{TokenKind::Error};
    std::string_view text;
    SourceLocation location;
};

[[nodiscard]] std::string_view token_kind_name(TokenKind kind) noexcept;

}  // namespace cas
