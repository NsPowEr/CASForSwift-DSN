#include "cas/token.hpp"

namespace cas {

std::string_view token_kind_name(TokenKind kind) noexcept {
    switch (kind) {
    case TokenKind::Integer:
        return "Integer";
    case TokenKind::Rational:
        return "Rational";
    case TokenKind::Float:
        return "Float";
    case TokenKind::Identifier:
        return "Identifier";
    case TokenKind::Plus:
        return "Plus";
    case TokenKind::Minus:
        return "Minus";
    case TokenKind::Star:
        return "Star";
    case TokenKind::Slash:
        return "Slash";
    case TokenKind::Caret:
        return "Caret";
    case TokenKind::Percent:
        return "Percent";
    case TokenKind::Bang:
        return "Bang";
    case TokenKind::LParen:
        return "LParen";
    case TokenKind::RParen:
        return "RParen";
    case TokenKind::LBracket:
        return "LBracket";
    case TokenKind::RBracket:
        return "RBracket";
    case TokenKind::LBrace:
        return "LBrace";
    case TokenKind::RBrace:
        return "RBrace";
    case TokenKind::Comma:
        return "Comma";
    case TokenKind::Semicolon:
        return "Semicolon";
    case TokenKind::Integral:
        return "Integral";
    case TokenKind::Derivative:
        return "Derivative";
    case TokenKind::Limit:
        return "Limit";
    case TokenKind::Sum:
        return "Sum";
    case TokenKind::Product:
        return "Product";
    case TokenKind::Sqrt:
        return "Sqrt";
    case TokenKind::EndOfInput:
        return "EndOfInput";
    case TokenKind::Error:
        return "Error";
    }

    return "Unknown";
}

}  // namespace cas
