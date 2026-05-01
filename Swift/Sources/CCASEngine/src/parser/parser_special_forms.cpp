#include "cas/parser.hpp"

#include "parser_internal.hpp"

#include <optional>

namespace cas {

Result<ExprPtr> Parser::parse_integral() {
    auto open = consume(TokenKind::LParen, ParseError::MissingArgument, "Expected '(' after integral token");
    if (open.is_error()) {
        return Result<ExprPtr>(open.error());
    }

    auto integrand = parse_expr();
    if (integrand.is_error()) {
        return integrand;
    }

    auto comma = consume(TokenKind::Comma, ParseError::MissingArgument, "Expected ',' after integrand");
    if (comma.is_error()) {
        return Result<ExprPtr>(comma.error());
    }

    auto variable_token = consume(TokenKind::Identifier, ParseError::MissingArgument, "Expected integration variable");
    if (variable_token.is_error()) {
        return Result<ExprPtr>(variable_token.error());
    }

    std::optional<ExprPtr> lower;
    std::optional<ExprPtr> upper;
    if (match(TokenKind::Comma)) {
        auto lower_expr = parse_expr();
        if (lower_expr.is_error()) {
            return lower_expr;
        }
        lower = lower_expr.value();

        auto second_comma = consume(TokenKind::Comma, ParseError::MissingArgument, "Expected ',' before upper bound");
        if (second_comma.is_error()) {
            return Result<ExprPtr>(second_comma.error());
        }

        auto upper_expr = parse_expr();
        if (upper_expr.is_error()) {
            return upper_expr;
        }
        upper = upper_expr.value();
    }

    auto close = consume(TokenKind::RParen, ParseError::UnmatchedParen, "Expected ')' after integral");
    if (close.is_error()) {
        return Result<ExprPtr>(close.error());
    }

    return Result<ExprPtr>(arena_.make<Integral>(
        integrand.value(),
        Symbol{std::string(variable_token.value().text)},
        lower,
        upper));
}

Result<ExprPtr> Parser::parse_derivative() {
    if (check(TokenKind::Identifier) && current_ + 1U < tokens_.size() && tokens_[current_ + 1U].kind == TokenKind::LParen) {
        const auto variable = advance();
        auto open = consume(TokenKind::LParen, ParseError::MissingArgument, "Expected '(' after derivative variable");
        if (open.is_error()) {
            return Result<ExprPtr>(open.error());
        }

        auto expression = parse_expr();
        if (expression.is_error()) {
            return expression;
        }

        auto close = consume(TokenKind::RParen, ParseError::UnmatchedParen, "Expected ')' after derivative expression");
        if (close.is_error()) {
            return Result<ExprPtr>(close.error());
        }

        return Result<ExprPtr>(arena_.make<Derivative>(
            expression.value(),
            Symbol{std::string(variable.text)},
            1U));
    }

    auto open = consume(TokenKind::LParen, ParseError::MissingArgument, "Expected '(' after diff");
    if (open.is_error()) {
        return Result<ExprPtr>(open.error());
    }

    auto expression = parse_expr();
    if (expression.is_error()) {
        return expression;
    }

    auto comma = consume(TokenKind::Comma, ParseError::MissingArgument, "Expected ',' after derivative expression");
    if (comma.is_error()) {
        return Result<ExprPtr>(comma.error());
    }

    auto variable_token = consume(TokenKind::Identifier, ParseError::MissingArgument, "Expected derivative variable");
    if (variable_token.is_error()) {
        return Result<ExprPtr>(variable_token.error());
    }

    unsigned int order = 1U;
    if (match(TokenKind::Comma)) {
        auto order_expr = parse_expr();
        if (order_expr.is_error()) {
            return order_expr;
        }

        const auto* order_value = expr_cast<IntegerLit>(order_expr.value());
        if (order_value == nullptr || order_value->value.is_negative()) {
            return Result<ExprPtr>(make_parse_error(ParseError::InvalidNumber, previous(), "Derivative order must be a non-negative integer"));
        }

        const auto parsed_order = parser_internal::parse_unsigned_decimal<unsigned int>(order_value->value);
        if (!parsed_order.has_value()) {
            return Result<ExprPtr>(make_parse_error(ParseError::InvalidNumber, previous(), "Derivative order is out of range"));
        }
        order = *parsed_order;
    }

    auto close = consume(TokenKind::RParen, ParseError::UnmatchedParen, "Expected ')' after derivative");
    if (close.is_error()) {
        return Result<ExprPtr>(close.error());
    }

    return Result<ExprPtr>(arena_.make<Derivative>(
        expression.value(),
        Symbol{std::string(variable_token.value().text)},
        order));
}

Result<ExprPtr> Parser::parse_limit() {
    auto open = consume(TokenKind::LParen, ParseError::MissingArgument, "Expected '(' after limit token");
    if (open.is_error()) {
        return Result<ExprPtr>(open.error());
    }

    auto expression = parse_expr();
    if (expression.is_error()) {
        return expression;
    }

    auto comma = consume(TokenKind::Comma, ParseError::MissingArgument, "Expected ',' after limit expression");
    if (comma.is_error()) {
        return Result<ExprPtr>(comma.error());
    }

    auto variable_token = consume(TokenKind::Identifier, ParseError::MissingArgument, "Expected limit variable");
    if (variable_token.is_error()) {
        return Result<ExprPtr>(variable_token.error());
    }

    auto second_comma = consume(TokenKind::Comma, ParseError::MissingArgument, "Expected ',' after limit variable");
    if (second_comma.is_error()) {
        return Result<ExprPtr>(second_comma.error());
    }

    auto point = parse_expr();
    if (point.is_error()) {
        return point;
    }

    LimitDirection direction = LimitDirection::Both;
    if (match(TokenKind::Comma)) {
        auto dir_token = consume(TokenKind::Identifier, ParseError::MissingArgument, "Expected limit direction");
        if (dir_token.is_error()) {
            return Result<ExprPtr>(dir_token.error());
        }

        if (dir_token.value().text == "left") {
            direction = LimitDirection::Left;
        } else if (dir_token.value().text == "right") {
            direction = LimitDirection::Right;
        } else if (dir_token.value().text == "both") {
            direction = LimitDirection::Both;
        } else {
            return Result<ExprPtr>(make_parse_error(ParseError::UnexpectedToken, dir_token.value(), "Unknown limit direction"));
        }
    }

    auto close = consume(TokenKind::RParen, ParseError::UnmatchedParen, "Expected ')' after limit");
    if (close.is_error()) {
        return Result<ExprPtr>(close.error());
    }

    return Result<ExprPtr>(arena_.make<Limit>(
        expression.value(),
        Symbol{std::string(variable_token.value().text)},
        point.value(),
        direction));
}

Result<ExprPtr> Parser::parse_matrix() {
    std::vector<std::vector<ExprPtr>> rows;
    if (check(TokenKind::RBracket)) {
        const auto& ignored = advance();
        (void)ignored;
        return Result<ExprPtr>(arena_.make<Matrix>(0U, 0U, std::vector<ExprPtr>{}));
    }

    while (true) {
        auto row_open = consume(TokenKind::LBracket, ParseError::UnexpectedToken, "Expected '[' to start a matrix row");
        if (row_open.is_error()) {
            return Result<ExprPtr>(row_open.error());
        }

        std::vector<ExprPtr> row;
        if (!check(TokenKind::RBracket)) {
            while (true) {
                auto value = parse_expr();
                if (value.is_error()) {
                    return value;
                }
                row.push_back(value.value());

                if (!match(TokenKind::Comma)) {
                    break;
                }
            }
        }

        auto row_close = consume(TokenKind::RBracket, ParseError::UnexpectedToken, "Expected ']' to close a matrix row");
        if (row_close.is_error()) {
            return Result<ExprPtr>(row_close.error());
        }

        rows.push_back(std::move(row));
        if (!match(TokenKind::Comma)) {
            break;
        }
    }

    auto matrix_close = consume(TokenKind::RBracket, ParseError::UnexpectedToken, "Expected ']' to close matrix");
    if (matrix_close.is_error()) {
        return Result<ExprPtr>(matrix_close.error());
    }

    const std::size_t cols = rows.front().size();
    for (const auto& row : rows) {
        if (row.size() != cols) {
            return Result<ExprPtr>(make_parse_error(ParseError::UnexpectedToken, matrix_close.value(), "Matrix rows must have equal length"));
        }
    }

    std::vector<ExprPtr> elements;
    for (const auto& row : rows) {
        elements.insert(elements.end(), row.begin(), row.end());
    }

    return Result<ExprPtr>(arena_.make<Matrix>(rows.size(), cols, std::move(elements)));
}

}  // namespace cas
