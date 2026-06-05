#include "cas/parser.hpp"

#include "cas/error.hpp"
#include "parser_internal.hpp"

namespace cas::parser_internal {

std::optional<MathConstant> parse_math_constant(std::string_view text) {
    if (text == "pi" || text == "\xCF\x80") {
        return MathConstant::Pi;
    }
    if (text == "e") {
        return MathConstant::E;
    }
    if (text == "i") {
        return MathConstant::I;
    }
    if (text == "inf" || text == "\xE2\x88\x9E") {
        return MathConstant::Infinity;
    }
    if (text == "nan") {
        return MathConstant::NaN;
    }
    if (text == "EulerGamma" || text == "eulergamma") {
        return MathConstant::EulerGamma;
    }
    return std::nullopt;
}

bool is_builtin_function(std::string_view name) {
    constexpr std::string_view builtins[] = {
        "sin",   "cos",    "tan",    "cot",   "sec",    "csc",    "arcsin", "arccos",
        "arctan","asin",   "acos",   "atan",  "sinh",   "cosh",   "tanh",   "coth",
        "arcsinh","arccosh","arctanh","asinh","acosh",  "atanh",  "exp",    "ln",
        "log",   "log2",   "log10",  "sqrt",  "cbrt",   "floor",  "ceil",   "round",
        "abs",   "sign",   "gamma",  "beta",  "erf",    "erfc",   "zeta",   "Gamma",
        "arg",   "Re",     "Im",     "conj",  "conjugate",
        "bessel_zero", "BesselZero",
        "RootOf",
        "delta", "DiracDelta",
        "theta", "Theta", "heaviside", "Heaviside",
        "Hypergeometric0F1", "hyp0F1",
        "Hypergeometric1F1", "hyp1F1", "Kummer1F1", "KummerM",
        "Hypergeometric2F1", "hyp2F1",
        "EllipticK", "ellipticK",
        "EllipticE", "ellipticE",
        "EllipticPi", "ellipticPi",
        "EllipticF", "ellipticF"
    };

    for (const auto builtin : builtins) {
        if (name == builtin) {
            return true;
        }
    }
    return false;
}

std::string token_text_or_name(const Token& token) {
    if (!token.text.empty()) {
        return std::string(token.text);
    }
    return std::string(token_kind_name(token.kind));
}

Token eof_token(std::size_t line, std::size_t column) noexcept {
    return Token{
        .kind = TokenKind::EndOfInput,
        .text = std::string_view{},
        .location = SourceLocation{.line = line, .column = column},
    };
}

}  // namespace cas::parser_internal

namespace cas {

bool Parser::is_at_end() const noexcept {
    return check(TokenKind::EndOfInput);
}

const Token& Parser::peek() const noexcept {
    return tokens_[current_];
}

const Token& Parser::previous() const noexcept {
    return tokens_[current_ - 1U];
}

const Token& Parser::advance() noexcept {
    const Token& token = current_ < tokens_.size() ? tokens_[current_] : tokens_.back();
    if (!is_at_end()) {
        ++current_;
    }
    return token;
}

bool Parser::check(TokenKind kind) const noexcept {
    if (current_ >= tokens_.size()) {
        return kind == TokenKind::EndOfInput;
    }
    return tokens_[current_].kind == kind;
}

bool Parser::match(TokenKind kind) noexcept {
    if (!check(kind)) {
        return false;
    }

    const auto& ignored = advance();
    (void)ignored;
    return true;
}

bool Parser::starts_implicit_multiplication() const noexcept {
    switch (peek().kind) {
    case TokenKind::Integer:
    case TokenKind::Rational:
    case TokenKind::Float:
    case TokenKind::Identifier:
    case TokenKind::Sqrt:
    case TokenKind::LParen:
    case TokenKind::LBracket:
    case TokenKind::Integral:
    case TokenKind::Derivative:
    case TokenKind::Limit:
        return true;
    default:
        return false;
    }
}

Result<Token> Parser::consume(TokenKind kind, ParseError error, std::string message) {
    if (check(kind)) {
        return Result<Token>(advance());
    }
    return Result<Token>(make_parse_error(error, peek(), std::move(message)));
}

Result<ExprPtr> Parser::parse_number_token(const Token& token) {
    switch (token.kind) {
    case TokenKind::Integer: {
        auto integer = BigInt::parse(std::string(token.text));
        if (integer.is_error()) {
            return Result<ExprPtr>(make_parse_error(ParseError::InvalidNumber, token, integer.error().message));
        }
        return Result<ExprPtr>(arena_.make<IntegerLit>(std::move(integer.value())));
    }
    case TokenKind::Rational: {
        const auto slash = token.text.find('/');
        const auto numerator_text = std::string(token.text.substr(0U, slash));
        const auto denominator_text = std::string(token.text.substr(slash + 1U));
        auto numerator = BigInt::parse(std::move(numerator_text));
        if (numerator.is_error()) {
            return Result<ExprPtr>(make_parse_error(ParseError::InvalidNumber, token, numerator.error().message));
        }

        auto denominator = BigInt::parse(std::move(denominator_text));
        if (denominator.is_error()) {
            return Result<ExprPtr>(make_parse_error(ParseError::InvalidNumber, token, denominator.error().message));
        }
        if (denominator.value().is_zero()) {
            return Result<ExprPtr>(make_parse_error(ParseError::InvalidNumber, token, "Rational denominator must be non-zero"));
        }

        return Result<ExprPtr>(arena_.make<RationalLit>(
            std::move(numerator.value()),
            std::move(denominator.value())));
    }
    case TokenKind::Float: {
        const std::string_view text = token.text;
        const auto dot_pos = text.find('.');
        if (dot_pos == std::string_view::npos) {
            return Result<ExprPtr>(arena_.make<DecimalLit>(std::string(text)));
        }

        const std::string_view int_part = text.substr(0U, dot_pos);
        const std::string_view frac_part = text.substr(dot_pos + 1U);
        const std::size_t k = frac_part.size();

        const std::string num_str = std::string(int_part) + std::string(frac_part);
        auto numerator = BigInt::parse(num_str);
        if (numerator.is_error()) {
            return Result<ExprPtr>(make_parse_error(ParseError::InvalidNumber, token, numerator.error().message));
        }

        BigInt denominator(1);
        const BigInt ten(10);
        for (std::size_t i = 0U; i < k; ++i) {
            denominator *= ten;
        }

        BigInt g = gcd(numerator.value().abs(), denominator);
        BigInt num_reduced = numerator.value() / g;
        BigInt den_reduced = denominator / g;

        return Result<ExprPtr>(arena_.make<RationalLit>(
            std::move(num_reduced),
            std::move(den_reduced)));
    }
    default:
        return Result<ExprPtr>(make_parse_error(ParseError::InvalidNumber, token, "Unexpected numeric token"));
    }
}

Result<ExprPtr> Parser::parse_identifier_or_call(const Token& token) {
    if (check(TokenKind::LParen)) {
        return parse_function_call(std::string(token.text));
    }

    if (const auto constant = parser_internal::parse_math_constant(token.text)) {
        return Result<ExprPtr>(arena_.make<Constant>(*constant));
    }

    return Result<ExprPtr>(arena_.make<Symbol>(std::string(token.text)));
}

int Parser::infix_precedence(TokenKind kind) const noexcept {
    switch (kind) {
    case TokenKind::Less:
    case TokenKind::Greater:
    case TokenKind::LessEqual:
    case TokenKind::GreaterEqual:
    case TokenKind::DoubleEqual:
        return 5;
    case TokenKind::Plus:
    case TokenKind::Minus:
        return 10;
    case TokenKind::Star:
    case TokenKind::Slash:
    case TokenKind::Percent:
        return 20;
    case TokenKind::Caret:
        return 30;
    default:
        return -1;
    }
}

bool Parser::is_right_associative(TokenKind kind) const noexcept {
    return kind == TokenKind::Caret;
}

CASError Parser::make_parse_error(ParseError error, const Token& token, std::string message) const {
    const char* hint = nullptr;
    switch (error) {
    case ParseError::EmptyInput:
        hint = "provide at least one token";
        break;
    case ParseError::UnexpectedToken:
        hint = "check token order and delimiters";
        break;
    case ParseError::UnexpectedEof:
        hint = "expression ended too early";
        break;
    case ParseError::UnmatchedParen:
        hint = "balance parentheses or brackets";
        break;
    case ParseError::InvalidNumber:
        hint = "use a valid integer, rational, or decimal literal";
        break;
    case ParseError::UnknownFunction:
        hint = "register the function downstream if semantics are needed";
        break;
    case ParseError::MissingArgument:
        hint = "provide all required arguments";
        break;
    }

    return CASError{
        .kind = CASErrorKind::ParseError,
        .message = std::move(message) + " at " + std::to_string(token.location.line) + ":" +
            std::to_string(token.location.column) + " near '" + parser_internal::token_text_or_name(token) + "'",
        .hint = hint == nullptr ? std::nullopt : std::optional<std::string>{hint},
    };
}

}  // namespace cas
