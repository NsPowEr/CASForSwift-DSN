#include "cas/token.hpp"

#include <gtest/gtest.h>

namespace cas {
namespace {

TEST(TokenModelTest, MinusRemainsDedicatedOperatorToken) {
    Token token{
        .kind = TokenKind::Minus,
        .text = "-",
        .location = SourceLocation{.line = 1U, .column = 1U},
    };

    EXPECT_EQ(token.kind, TokenKind::Minus);
    EXPECT_EQ(token.text, "-");
    EXPECT_EQ(token_kind_name(token.kind), "Minus");
}

TEST(TokenModelTest, FloatTokenExistsOnlyAtLexerFrontEndBoundary) {
    Token token{
        .kind = TokenKind::Float,
        .text = "3.14",
        .location = SourceLocation{.line = 2U, .column = 4U},
    };

    EXPECT_EQ(token.kind, TokenKind::Float);
    EXPECT_EQ(token.text, "3.14");
}

TEST(TokenModelTest, DedicatedSpecialTokensAreAddressable) {
    EXPECT_EQ(token_kind_name(TokenKind::Integral), "Integral");
    EXPECT_EQ(token_kind_name(TokenKind::Derivative), "Derivative");
    EXPECT_EQ(token_kind_name(TokenKind::Limit), "Limit");
    EXPECT_EQ(token_kind_name(TokenKind::Sqrt), "Sqrt");
}

TEST(TokenModelTest, SourceLocationDefaultsToOneBasedCoordinates) {
    Token token;

    EXPECT_EQ(token.location.line, 1U);
    EXPECT_EQ(token.location.column, 1U);
}

}  // namespace
}  // namespace cas
