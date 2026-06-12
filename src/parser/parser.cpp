#include "cas/parser.hpp"

#include "parser_internal.hpp"

namespace cas {

Parser::Parser(const std::vector<Token>& tokens, AstArena& arena) noexcept
    : tokens_(tokens), arena_(arena) {}

Result<ExprPtr> Parser::parse() {
    if (tokens_.empty()) {
        return Result<ExprPtr>(make_parse_error(ParseError::EmptyInput, parser_internal::eof_token(), "Input is empty"));
    }

    if (check(TokenKind::EndOfInput)) {
        return Result<ExprPtr>(make_parse_error(ParseError::EmptyInput, peek(), "Input is empty"));
    }

    auto expression = parse_expr();
    if (expression.is_error()) {
        return expression;
    }

    if (!check(TokenKind::EndOfInput)) {
        return Result<ExprPtr>(make_parse_error(ParseError::UnexpectedToken, peek(), "Unexpected trailing token"));
    }

    return expression;
}

Result<ExprPtr> Parser::parse_expr(int min_precedence) {
    auto left = parse_prefix();
    if (left.is_error()) {
        return left;
    }

    while (true) {
        auto postfix = parse_postfix(left.value());
        if (postfix.is_error()) {
            return postfix;
        }
        left = std::move(postfix);

        constexpr int implicit_multiplication_precedence = 25;
        if (starts_implicit_multiplication()) {
            if (implicit_multiplication_precedence < min_precedence) {
                break;
            }

            auto right = parse_expr(implicit_multiplication_precedence + 1);
            if (right.is_error()) {
                return right;
            }

            left = Result<ExprPtr>(make_product_expr(left.value(), right.value()));
            continue;
        }

        const int precedence = infix_precedence(peek().kind);
        if (precedence < min_precedence) {
            break;
        }

        auto infix = parse_infix(left.value(), precedence);
        if (infix.is_error()) {
            return infix;
        }
        left = std::move(infix);
    }

    return left;
}

Result<ExprPtr> Parser::parse_prefix() {
    if (match(TokenKind::Minus)) {
        auto operand = parse_expr(40);
        if (operand.is_error()) {
            return operand;
        }
        return Result<ExprPtr>(arena_.make<Unary>(UnaryOp::Neg, operand.value()));
    }

    if (match(TokenKind::Integral)) {
        return parse_integral();
    }

    if (match(TokenKind::Derivative)) {
        return parse_derivative();
    }

    if (match(TokenKind::Limit)) {
        return parse_limit();
    }

    return parse_primary();
}

Result<ExprPtr> Parser::parse_primary() {
    const Token token = advance();

    switch (token.kind) {
    case TokenKind::Integer:
    case TokenKind::Rational:
    case TokenKind::Float:
        return parse_number_token(token);
    case TokenKind::Identifier:
        return parse_identifier_or_call(token);
    case TokenKind::Sqrt:
        return parse_function_call("sqrt");
    case TokenKind::LParen:
        return parse_grouped_expression();
    case TokenKind::LBracket:
        return parse_matrix();
    case TokenKind::EndOfInput:
        return Result<ExprPtr>(make_parse_error(ParseError::UnexpectedEof, token, "Unexpected end of input"));
    default:
        return Result<ExprPtr>(make_parse_error(ParseError::UnexpectedToken, token, "Unexpected token"));
    }
}

Result<ExprPtr> Parser::parse_grouped_expression() {
    auto expression = parse_expr();
    if (expression.is_error()) {
        return expression;
    }

    auto closing = consume(TokenKind::RParen, ParseError::UnmatchedParen, "Expected ')'");
    if (closing.is_error()) {
        return Result<ExprPtr>(closing.error());
    }

    return expression;
}

Result<ExprPtr> Parser::parse_function_call(std::string name) {
    auto open = consume(TokenKind::LParen, ParseError::MissingArgument, "Expected '(' after function name");
    if (open.is_error()) {
        return Result<ExprPtr>(open.error());
    }

    std::vector<ExprPtr> args;
    if (!check(TokenKind::RParen)) {
        while (true) {
            auto argument = parse_expr();
            if (argument.is_error()) {
                return argument;
            }
            args.push_back(argument.value());

            if (!match(TokenKind::Comma)) {
                break;
            }
        }
    }

    auto close = consume(TokenKind::RParen, ParseError::UnmatchedParen, "Expected ')' after function arguments");
    if (close.is_error()) {
        return Result<ExprPtr>(close.error());
    }

    if (name == "RootOf") {
        // Accepted arities:
        //   2 args: (polynomial, variable)
        //   3 args: + index
        //   7 args: + index + low_num,low_den,high_num,high_den (bound)
        if (args.size() != 2U && args.size() != 3U && args.size() != 7U) {
            return Result<ExprPtr>(make_parse_error(ParseError::MissingArgument, close.value(),
                "RootOf expects (polynomial, variable [, index [, low_num, low_den, high_num, high_den]])"));
        }

        const auto* variable = expr_cast<Symbol>(args[1]);
        if (variable == nullptr) {
            return Result<ExprPtr>(make_parse_error(ParseError::UnexpectedToken, close.value(), "RootOf variable must be a symbol"));
        }

        std::optional<std::size_t> root_index;
        if (args.size() >= 3U) {
            const auto* index = expr_cast<IntegerLit>(args[2]);
            if (index == nullptr || index->value.is_negative()) {
                return Result<ExprPtr>(make_parse_error(ParseError::InvalidNumber, close.value(), "RootOf index must be a non-negative integer"));
            }
            root_index = parser_internal::parse_unsigned_decimal<std::size_t>(index->value);
            if (!root_index.has_value()) {
                return Result<ExprPtr>(make_parse_error(ParseError::InvalidNumber, close.value(), "RootOf index is out of range"));
            }
        }

        if (args.size() == 7U) {
            // Extract isolating bound from args[3..7]. Each must be a signed
            // integer literal (Unary(Neg, IntegerLit) is accepted to allow
            // negative endpoints).
            auto extract_bigint = [](ExprPtr e) -> std::optional<BigInt> {
                if (const auto* lit = expr_cast<IntegerLit>(e)) return lit->value;
                if (const auto* un = expr_cast<Unary>(e)) {
                    if (un->op == UnaryOp::Neg) {
                        if (const auto* lit = expr_cast<IntegerLit>(un->operand))
                            return -lit->value;
                    }
                }
                return std::nullopt;
            };
            auto low_num  = extract_bigint(args[3]);
            auto low_den  = extract_bigint(args[4]);
            auto high_num = extract_bigint(args[5]);
            auto high_den = extract_bigint(args[6]);
            if (!low_num || !low_den || !high_num || !high_den) {
                return Result<ExprPtr>(make_parse_error(ParseError::InvalidNumber, close.value(),
                    "RootOf bound endpoints must be signed integer literals"));
            }
            if (low_den->is_zero() || high_den->is_zero()) {
                return Result<ExprPtr>(make_parse_error(ParseError::InvalidNumber, close.value(),
                    "RootOf bound denominator must be non-zero"));
            }
            IsolatingBound bound{*low_num, *low_den, *high_num, *high_den};
            return Result<ExprPtr>(arena_.make<RootOf>(
                args[0], *variable, std::move(bound), root_index));
        }

        return Result<ExprPtr>(arena_.make<RootOf>(args[0], *variable, root_index));
    }

    const auto builtin = parser_internal::is_builtin_function(name);
    (void)builtin;
    return Result<ExprPtr>(arena_.make<FuncCall>(std::move(name), std::move(args)));
}

Result<ExprPtr> Parser::parse_postfix(ExprPtr left) {
    while (match(TokenKind::Bang)) {
        left = arena_.make<Unary>(UnaryOp::Factorial, left);
    }
    return Result<ExprPtr>(left);
}

Result<ExprPtr> Parser::parse_infix(ExprPtr left, int precedence) {
    const Token op = advance();
    const int next_precedence = is_right_associative(op.kind) ? precedence : precedence + 1;

    auto right = parse_expr(next_precedence);
    if (right.is_error()) {
        return right;
    }

    switch (op.kind) {
    case TokenKind::Plus: {
        std::vector<ExprPtr> terms;
        if (const auto* sum = expr_cast<Sum>(left)) {
            terms = sum->terms;
        } else {
            terms.push_back(left);
        }

        if (const auto* sum = expr_cast<Sum>(right.value())) {
            terms.insert(terms.end(), sum->terms.begin(), sum->terms.end());
        } else {
            terms.push_back(right.value());
        }
        return Result<ExprPtr>(arena_.make<Sum>(std::move(terms)));
    }
    case TokenKind::Star:
        return Result<ExprPtr>(make_product_expr(left, right.value()));
    case TokenKind::Minus:
        return Result<ExprPtr>(arena_.make<Binary>(BinaryOp::Sub, left, right.value()));
    case TokenKind::Slash:
        return Result<ExprPtr>(arena_.make<Binary>(BinaryOp::Div, left, right.value()));
    case TokenKind::Caret:
        return Result<ExprPtr>(arena_.make<Binary>(BinaryOp::Pow, left, right.value()));
    case TokenKind::Percent:
        return Result<ExprPtr>(arena_.make<Binary>(BinaryOp::Mod, left, right.value()));
    case TokenKind::Less:
        return Result<ExprPtr>(arena_.make<Binary>(BinaryOp::Less, left, right.value()));
    case TokenKind::Greater:
        return Result<ExprPtr>(arena_.make<Binary>(BinaryOp::Greater, left, right.value()));
    case TokenKind::LessEqual:
        return Result<ExprPtr>(arena_.make<Binary>(BinaryOp::LessEqual, left, right.value()));
    case TokenKind::GreaterEqual:
        return Result<ExprPtr>(arena_.make<Binary>(BinaryOp::GreaterEqual, left, right.value()));
    case TokenKind::DoubleEqual:
        return Result<ExprPtr>(arena_.make<Binary>(BinaryOp::Equal, left, right.value()));
    default:
        return Result<ExprPtr>(make_parse_error(ParseError::UnexpectedToken, op, "Unexpected infix operator"));
    }
}

ExprPtr Parser::make_product_expr(ExprPtr left, ExprPtr right) {
    std::vector<ExprPtr> factors;
    if (const auto* product = expr_cast<Product>(left)) {
        factors = product->factors;
    } else {
        factors.push_back(left);
    }

    if (const auto* product = expr_cast<Product>(right)) {
        factors.insert(factors.end(), product->factors.begin(), product->factors.end());
    } else {
        factors.push_back(right);
    }

    return arena_.make<Product>(std::move(factors));
}

}  // namespace cas
