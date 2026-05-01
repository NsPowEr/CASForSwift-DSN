#include "cas/lexer.hpp"

#include "cas/error.hpp"

#include <cctype>
#include <string>

namespace cas {

namespace {

constexpr std::string_view kUtf8Integral = "\xE2\x88\xAB";
constexpr std::string_view kUtf8Sqrt = "\xE2\x88\x9A";
constexpr std::string_view kUtf8Pi = "\xCF\x80";
constexpr std::string_view kUtf8Sigma = "\xCE\xA3";
constexpr std::string_view kUtf8Product = "\xCE\xA0";
constexpr std::string_view kUtf8Infinity = "\xE2\x88\x9E";

[[nodiscard]] bool is_identifier_start(unsigned char ch) noexcept {
    return std::isalpha(ch) != 0 || ch == '_';
}

[[nodiscard]] bool is_identifier_continue(unsigned char ch) noexcept {
    return std::isalnum(ch) != 0 || ch == '_';
}

[[nodiscard]] CASError make_lexer_error(LexerError kind, std::size_t line, std::size_t column, std::string message) {
    const auto hint = kind == LexerError::InvalidCharacter
        ? std::optional<std::string>{"remove or escape unsupported input"}
        : std::optional<std::string>{"check numeric literal syntax"};

    return CASError{
        .kind = CASErrorKind::ParseError,
        .message = std::move(message) + " at " + std::to_string(line) + ":" + std::to_string(column),
        .hint = hint,
    };
}

}  // namespace

Lexer::Lexer(std::string_view input) noexcept : input_(input) {}

Result<std::vector<Token>> Lexer::tokenize() const {
    std::vector<Token> tokens;
    tokens.reserve(input_.size() + 1U);

    std::size_t position = 0U;
    std::size_t line = 1U;
    std::size_t column = 1U;

    auto push_token = [&](TokenKind kind, std::size_t start, std::size_t length, std::size_t token_line, std::size_t token_column) {
        tokens.push_back(Token{
            .kind = kind,
            .text = input_.substr(start, length),
            .location = SourceLocation{.line = token_line, .column = token_column},
        });
    };

    while (position < input_.size()) {
        const auto token_line = line;
        const auto token_column = column;
        const unsigned char ch = static_cast<unsigned char>(input_[position]);

        if (ch == ' ' || ch == '\t' || ch == '\r') {
            ++position;
            ++column;
            continue;
        }

        if (ch == '\n') {
            ++position;
            ++line;
            column = 1U;
            continue;
        }

        if (input_.substr(position).starts_with(kUtf8Integral)) {
            push_token(TokenKind::Integral, position, kUtf8Integral.size(), token_line, token_column);
            position += kUtf8Integral.size();
            column += 1U;
            continue;
        }

        if (input_.substr(position).starts_with(kUtf8Sqrt)) {
            push_token(TokenKind::Sqrt, position, kUtf8Sqrt.size(), token_line, token_column);
            position += kUtf8Sqrt.size();
            column += 1U;
            continue;
        }

        if (input_.substr(position).starts_with(kUtf8Pi)) {
            push_token(TokenKind::Identifier, position, kUtf8Pi.size(), token_line, token_column);
            position += kUtf8Pi.size();
            column += 1U;
            continue;
        }

        if (input_.substr(position).starts_with(kUtf8Sigma)) {
            push_token(TokenKind::Sum, position, kUtf8Sigma.size(), token_line, token_column);
            position += kUtf8Sigma.size();
            column += 1U;
            continue;
        }

        if (input_.substr(position).starts_with(kUtf8Product)) {
            push_token(TokenKind::Product, position, kUtf8Product.size(), token_line, token_column);
            position += kUtf8Product.size();
            column += 1U;
            continue;
        }

        if (input_.substr(position).starts_with(kUtf8Infinity)) {
            push_token(TokenKind::Identifier, position, kUtf8Infinity.size(), token_line, token_column);
            position += kUtf8Infinity.size();
            column += 1U;
            continue;
        }

        if (input_.substr(position).starts_with("d/d")) {
            push_token(TokenKind::Derivative, position, 3U, token_line, token_column);
            position += 3U;
            column += 3U;
            continue;
        }

        switch (ch) {
        case '+':
            push_token(TokenKind::Plus, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case '-':
            push_token(TokenKind::Minus, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case '*':
            push_token(TokenKind::Star, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case '/':
            push_token(TokenKind::Slash, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case '^':
            push_token(TokenKind::Caret, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case '%':
            push_token(TokenKind::Percent, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case '!':
            push_token(TokenKind::Bang, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case '(':
            push_token(TokenKind::LParen, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case ')':
            push_token(TokenKind::RParen, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case '[':
            push_token(TokenKind::LBracket, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case ']':
            push_token(TokenKind::RBracket, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case '{':
            push_token(TokenKind::LBrace, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case '}':
            push_token(TokenKind::RBrace, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case ',':
            push_token(TokenKind::Comma, position++, 1U, token_line, token_column);
            ++column;
            continue;
        case ';':
            push_token(TokenKind::Semicolon, position++, 1U, token_line, token_column);
            ++column;
            continue;
        default:
            break;
        }

        if (std::isdigit(ch) != 0) {
            const auto start = position;
            while (position < input_.size() && std::isdigit(static_cast<unsigned char>(input_[position])) != 0) {
                ++position;
                ++column;
            }

            bool is_float = false;
            if (position < input_.size() && input_[position] == '.') {
                is_float = true;
                ++position;
                ++column;

                if (position >= input_.size() || std::isdigit(static_cast<unsigned char>(input_[position])) == 0) {
                    return Result<std::vector<Token>>(make_lexer_error(
                        LexerError::InvalidNumber,
                        token_line,
                        token_column,
                        "Invalid decimal literal"));
                }

                while (position < input_.size() && std::isdigit(static_cast<unsigned char>(input_[position])) != 0) {
                    ++position;
                    ++column;
                }
            }

            if (!is_float && position < input_.size() && input_[position] == '/') {
                std::size_t lookahead = position + 1U;
                if (lookahead < input_.size() && std::isdigit(static_cast<unsigned char>(input_[lookahead])) != 0) {
                    while (lookahead < input_.size() && std::isdigit(static_cast<unsigned char>(input_[lookahead])) != 0) {
                        ++lookahead;
                    }

                    const auto denominator = input_.substr(position + 1U, lookahead - position - 1U);
                    if (denominator != "0") {
                        push_token(TokenKind::Rational, start, lookahead - start, token_line, token_column);
                        column += lookahead - position;
                        position = lookahead;
                        continue;
                    }
                }
            }

            push_token(is_float ? TokenKind::Float : TokenKind::Integer, start, position - start, token_line, token_column);
            continue;
        }

        if (is_identifier_start(ch)) {
            const auto start = position;
            while (position < input_.size() && is_identifier_continue(static_cast<unsigned char>(input_[position]))) {
                ++position;
                ++column;
            }

            const auto text = input_.substr(start, position - start);
            TokenKind kind = TokenKind::Identifier;
            if (text == "int") {
                kind = TokenKind::Integral;
            } else if (text == "diff") {
                kind = TokenKind::Derivative;
            } else if (text == "lim") {
                kind = TokenKind::Limit;
            } else if (text == "sum") {
                kind = TokenKind::Sum;
            } else if (text == "prod") {
                kind = TokenKind::Product;
            } else if (text == "sqrt") {
                kind = TokenKind::Sqrt;
            }

            push_token(kind, start, position - start, token_line, token_column);
            continue;
        }

        return Result<std::vector<Token>>(make_lexer_error(
            LexerError::InvalidCharacter,
            token_line,
            token_column,
            "Invalid character in input"));
    }

    tokens.push_back(Token{
        .kind = TokenKind::EndOfInput,
        .text = std::string_view{},
        .location = SourceLocation{.line = line, .column = column},
    });

    return Result<std::vector<Token>>(std::move(tokens));
}

}  // namespace cas
