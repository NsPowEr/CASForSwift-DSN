#pragma once

#include "cas/result.hpp"
#include "cas/token.hpp"

#include <string_view>
#include <vector>

namespace cas {

enum class LexerError : std::uint8_t {
    InvalidCharacter,
    InvalidNumber,
};

struct LexerErrorInfo {
    LexerError kind;
    std::size_t line;
    std::size_t column;
    std::string message;
};

class Lexer {
public:
    static constexpr std::uint32_t API_VERSION = 1;

    explicit Lexer(std::string_view input) noexcept;

    /// Token text is a zero-copy view into the original input buffer.
    /// The input string must outlive the returned token vector.
    [[nodiscard]] Result<std::vector<Token>> tokenize() const;

private:
    std::string_view input_;
};

}  // namespace cas
