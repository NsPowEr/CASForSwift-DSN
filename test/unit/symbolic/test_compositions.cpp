#include <gtest/gtest.h>
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::symbolic;

namespace {
[[nodiscard]] Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}
}

class CompositionTest : public ::testing::Test {
protected:
    CASContext ctx;
};

TEST_F(CompositionTest, SinAsin) {
    auto e = parse_expr("sin(asin(x))", ctx.arena());
    auto simplified = ctx.simplify(e.value());
    ASSERT_TRUE(simplified.is_ok());
    auto text = to_round_trip_text(simplified.value());
    ASSERT_TRUE(text.is_ok());
    EXPECT_EQ(text.value(), "x");
}

TEST_F(CompositionTest, AsinSin) {
    auto x = ctx.arena().make<Symbol>("x");
    // Assume -pi/2 <= x <= pi/2
    ExprPtr pi_2 = ctx.arena().make<Binary>(BinaryOp::Div, ctx.arena().make<Constant>(MathConstant::Pi), ctx.arena().make<IntegerLit>(2));
    ExprPtr neg_pi_2 = ctx.arena().make<Unary>(UnaryOp::Neg, pi_2);
    ctx.assumptions().assume_greater_equal(x, neg_pi_2);
    ctx.assumptions().assume_greater_equal(pi_2, x);

    auto e = parse_expr("asin(sin(x))", ctx.arena());
    auto simplified = ctx.simplify(e.value());
    ASSERT_TRUE(simplified.is_ok());
    auto text = to_round_trip_text(simplified.value());
    ASSERT_TRUE(text.is_ok());
    EXPECT_EQ(text.value(), "x");
}

TEST_F(CompositionTest, SqrtSqrt) {
    auto e = parse_expr("sqrt(sqrt(x))", ctx.arena());
    auto simplified = ctx.simplify(e.value());
    ASSERT_TRUE(simplified.is_ok());
    auto text = to_round_trip_text(simplified.value());
    ASSERT_TRUE(text.is_ok());
    EXPECT_EQ(text.value(), "x^(1/4)");
}
