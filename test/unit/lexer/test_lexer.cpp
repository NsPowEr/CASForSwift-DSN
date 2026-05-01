#include "cas/lexer.hpp"

#include <gtest/gtest.h>

namespace cas {
namespace {

struct LexerSequenceCase {
    const char* input;
    std::vector<TokenKind> expected;
};

class LexerSequenceTest : public ::testing::TestWithParam<LexerSequenceCase> {};

TEST_P(LexerSequenceTest, ProducesExpectedTokenKinds) {
    const auto& test_case = GetParam();
    Lexer lexer(test_case.input);

    const auto result = lexer.tokenize();

    ASSERT_TRUE(result.is_ok()) << test_case.input;
    const auto& tokens = result.value();
    ASSERT_EQ(tokens.size(), test_case.expected.size()) << test_case.input;
    for (std::size_t i = 0; i < tokens.size(); ++i) {
        EXPECT_EQ(tokens[i].kind, test_case.expected[i]) << test_case.input << " token index " << i;
    }
}

INSTANTIATE_TEST_SUITE_P(
    CoreSequences,
    LexerSequenceTest,
    ::testing::Values(
        LexerSequenceCase{"x", {TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"42", {TokenKind::Integer, TokenKind::EndOfInput}},
        LexerSequenceCase{"3/4", {TokenKind::Rational, TokenKind::EndOfInput}},
        LexerSequenceCase{"1/0", {TokenKind::Integer, TokenKind::Slash, TokenKind::Integer, TokenKind::EndOfInput}},
        LexerSequenceCase{"3.14", {TokenKind::Float, TokenKind::EndOfInput}},
        LexerSequenceCase{"x+y", {TokenKind::Identifier, TokenKind::Plus, TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"x-y", {TokenKind::Identifier, TokenKind::Minus, TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"x*y", {TokenKind::Identifier, TokenKind::Star, TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"x/y", {TokenKind::Identifier, TokenKind::Slash, TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"x^y", {TokenKind::Identifier, TokenKind::Caret, TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"x%y", {TokenKind::Identifier, TokenKind::Percent, TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"x!", {TokenKind::Identifier, TokenKind::Bang, TokenKind::EndOfInput}},
        LexerSequenceCase{"(x)", {TokenKind::LParen, TokenKind::Identifier, TokenKind::RParen, TokenKind::EndOfInput}},
        LexerSequenceCase{"[x]", {TokenKind::LBracket, TokenKind::Identifier, TokenKind::RBracket, TokenKind::EndOfInput}},
        LexerSequenceCase{"{x}", {TokenKind::LBrace, TokenKind::Identifier, TokenKind::RBrace, TokenKind::EndOfInput}},
        LexerSequenceCase{"x,y;z", {TokenKind::Identifier, TokenKind::Comma, TokenKind::Identifier, TokenKind::Semicolon, TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"int(x,y)", {TokenKind::Integral, TokenKind::LParen, TokenKind::Identifier, TokenKind::Comma, TokenKind::Identifier, TokenKind::RParen, TokenKind::EndOfInput}},
        LexerSequenceCase{"diff(x,y)", {TokenKind::Derivative, TokenKind::LParen, TokenKind::Identifier, TokenKind::Comma, TokenKind::Identifier, TokenKind::RParen, TokenKind::EndOfInput}},
        LexerSequenceCase{"lim(x,y,0)", {TokenKind::Limit, TokenKind::LParen, TokenKind::Identifier, TokenKind::Comma, TokenKind::Identifier, TokenKind::Comma, TokenKind::Integer, TokenKind::RParen, TokenKind::EndOfInput}},
        LexerSequenceCase{"sum", {TokenKind::Sum, TokenKind::EndOfInput}},
        LexerSequenceCase{"prod", {TokenKind::Product, TokenKind::EndOfInput}},
        LexerSequenceCase{"sqrt(x)", {TokenKind::Sqrt, TokenKind::LParen, TokenKind::Identifier, TokenKind::RParen, TokenKind::EndOfInput}},
        LexerSequenceCase{"∫(x,x)", {TokenKind::Integral, TokenKind::LParen, TokenKind::Identifier, TokenKind::Comma, TokenKind::Identifier, TokenKind::RParen, TokenKind::EndOfInput}},
        LexerSequenceCase{"√x", {TokenKind::Sqrt, TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"Σ", {TokenKind::Sum, TokenKind::EndOfInput}},
        LexerSequenceCase{"Π", {TokenKind::Product, TokenKind::EndOfInput}},
        LexerSequenceCase{"π", {TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"∞", {TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"d/dx", {TokenKind::Derivative, TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"_x12", {TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"\n\t x", {TokenKind::Identifier, TokenKind::EndOfInput}},
        LexerSequenceCase{"[[1,2],[3,4]]", {TokenKind::LBracket, TokenKind::LBracket, TokenKind::Integer, TokenKind::Comma, TokenKind::Integer, TokenKind::RBracket, TokenKind::Comma, TokenKind::LBracket, TokenKind::Integer, TokenKind::Comma, TokenKind::Integer, TokenKind::RBracket, TokenKind::RBracket, TokenKind::EndOfInput}}));

class LexerInvalidInputTest : public ::testing::TestWithParam<const char*> {};

TEST_P(LexerInvalidInputTest, RejectsMalformedInput) {
    Lexer lexer(GetParam());

    const auto result = lexer.tokenize();

    ASSERT_TRUE(result.is_error()) << GetParam();
    EXPECT_EQ(result.error().kind, CASErrorKind::ParseError);
}

INSTANTIATE_TEST_SUITE_P(
    InvalidInputs,
    LexerInvalidInputTest,
    ::testing::Values(
        "@",
        "x @ y",
        "3.",
        ".5",
        "1..2",
        "x$",
        "\\",
        "\x01",
        "😀"));

TEST(LexerTest, TokenizesOperatorPrecedenceExample) {
    const std::string input = "x + y * z";
    Lexer lexer(input);

    auto result = lexer.tokenize();
    ASSERT_TRUE(result.is_ok());

    const auto& tokens = result.value();
    ASSERT_EQ(tokens.size(), 6U);
    EXPECT_EQ(tokens[0].kind, TokenKind::Identifier);
    EXPECT_EQ(tokens[1].kind, TokenKind::Plus);
    EXPECT_EQ(tokens[2].kind, TokenKind::Identifier);
    EXPECT_EQ(tokens[3].kind, TokenKind::Star);
    EXPECT_EQ(tokens[4].kind, TokenKind::Identifier);
    EXPECT_EQ(tokens[5].kind, TokenKind::EndOfInput);
}

TEST(LexerTest, TokenizesExactAndDecimalNumbersSeparately) {
    const std::string input = "3/4 3.14 1/0";
    Lexer lexer(input);

    auto result = lexer.tokenize();
    ASSERT_TRUE(result.is_ok());

    const auto& tokens = result.value();
    ASSERT_GE(tokens.size(), 6U);
    EXPECT_EQ(tokens[0].kind, TokenKind::Rational);
    EXPECT_EQ(tokens[1].kind, TokenKind::Float);
    EXPECT_EQ(tokens[2].kind, TokenKind::Integer);
    EXPECT_EQ(tokens[3].kind, TokenKind::Slash);
    EXPECT_EQ(tokens[4].kind, TokenKind::Integer);
}

TEST(LexerTest, KeepsMinusAsDedicatedToken) {
    const std::string input = "-3";
    Lexer lexer(input);

    auto result = lexer.tokenize();
    ASSERT_TRUE(result.is_ok());

    const auto& tokens = result.value();
    ASSERT_EQ(tokens.size(), 3U);
    EXPECT_EQ(tokens[0].kind, TokenKind::Minus);
    EXPECT_EQ(tokens[1].kind, TokenKind::Integer);
}

TEST(LexerTest, RecognizesSpecialMathTokens) {
    const std::string input = "∫ diff lim sqrt d/dx Σ Π π ∞";
    Lexer lexer(input);

    auto result = lexer.tokenize();
    ASSERT_TRUE(result.is_ok());

    const auto& tokens = result.value();
    ASSERT_GE(tokens.size(), 10U);
    EXPECT_EQ(tokens[0].kind, TokenKind::Integral);
    EXPECT_EQ(tokens[1].kind, TokenKind::Derivative);
    EXPECT_EQ(tokens[2].kind, TokenKind::Limit);
    EXPECT_EQ(tokens[3].kind, TokenKind::Sqrt);
    EXPECT_EQ(tokens[4].kind, TokenKind::Derivative);
    EXPECT_EQ(tokens[5].kind, TokenKind::Identifier);
    EXPECT_EQ(tokens[6].kind, TokenKind::Sum);
    EXPECT_EQ(tokens[7].kind, TokenKind::Product);
    EXPECT_EQ(tokens[8].kind, TokenKind::Identifier);
    EXPECT_EQ(tokens[9].kind, TokenKind::Identifier);
}

TEST(LexerTest, TracksOneBasedSourceLocationsAcrossLines) {
    const std::string input = "x\n  + y";
    Lexer lexer(input);

    auto result = lexer.tokenize();
    ASSERT_TRUE(result.is_ok());

    const auto& tokens = result.value();
    EXPECT_EQ(tokens[0].location.line, 1U);
    EXPECT_EQ(tokens[0].location.column, 1U);
    EXPECT_EQ(tokens[1].location.line, 2U);
    EXPECT_EQ(tokens[1].location.column, 3U);
}

TEST(LexerTest, RejectsInvalidCharacters) {
    const std::string input = "x @ y";
    Lexer lexer(input);

    auto result = lexer.tokenize();
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::ParseError);
}

TEST(LexerTest, RejectsMalformedDecimalLiteral) {
    const std::string input = "3.";
    Lexer lexer(input);

    auto result = lexer.tokenize();
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::ParseError);
}

}  // namespace
}  // namespace cas
