#pragma once

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/token.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace cas {

enum class ParseError : std::uint8_t {
    EmptyInput,
    UnexpectedToken,
    UnexpectedEof,
    UnmatchedParen,
    InvalidNumber,
    UnknownFunction,
    MissingArgument,
};

struct ParseErrorInfo {
    ParseError kind;
    std::size_t line;
    std::size_t column;
    std::string message;
};

class Parser {
public:
    Parser(const std::vector<Token>& tokens, AstArena& arena) noexcept;

    [[nodiscard]] Result<ExprPtr> parse();

private:
    [[nodiscard]] Result<ExprPtr> parse_expr(int min_precedence = 0);
    [[nodiscard]] Result<ExprPtr> parse_prefix();
    [[nodiscard]] Result<ExprPtr> parse_primary();
    [[nodiscard]] Result<ExprPtr> parse_grouped_expression();
    [[nodiscard]] Result<ExprPtr> parse_function_call(std::string name);
    [[nodiscard]] Result<ExprPtr> parse_integral();
    [[nodiscard]] Result<ExprPtr> parse_derivative();
    [[nodiscard]] Result<ExprPtr> parse_limit();
    [[nodiscard]] Result<ExprPtr> parse_matrix();
    [[nodiscard]] Result<ExprPtr> parse_postfix(ExprPtr left);
    [[nodiscard]] Result<ExprPtr> parse_infix(ExprPtr left, int precedence);
    [[nodiscard]] ExprPtr make_product_expr(ExprPtr left, ExprPtr right);

    [[nodiscard]] bool is_at_end() const noexcept;
    [[nodiscard]] const Token& peek() const noexcept;
    [[nodiscard]] const Token& previous() const noexcept;
    [[nodiscard]] const Token& advance() noexcept;
    [[nodiscard]] bool check(TokenKind kind) const noexcept;
    [[nodiscard]] bool match(TokenKind kind) noexcept;
    [[nodiscard]] bool starts_implicit_multiplication() const noexcept;
    [[nodiscard]] Result<Token> consume(TokenKind kind, ParseError error, std::string message);
    [[nodiscard]] Result<ExprPtr> parse_number_token(const Token& token);
    [[nodiscard]] Result<ExprPtr> parse_identifier_or_call(const Token& token);

    [[nodiscard]] int infix_precedence(TokenKind kind) const noexcept;
    [[nodiscard]] bool is_right_associative(TokenKind kind) const noexcept;
    [[nodiscard]] CASError make_parse_error(ParseError error, const Token& token, std::string message) const;

    const std::vector<Token>& tokens_;
    AstArena& arena_;
    std::size_t current_{0U};
};

[[nodiscard]] Result<std::string> to_round_trip_text(ExprPtr expr);

}  // namespace cas
