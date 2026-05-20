#include <gtest/gtest.h>
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
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

class AbsSignTest : public ::testing::Test {
protected:
    CASContext ctx;
};

TEST_F(AbsSignTest, AbsPositive) {
    auto x = ctx.arena().make<Symbol>("x");
    ctx.assumptions().assume_positive(expr_ref<Symbol>(x));
    
    auto e = parse_expr("abs(x)", ctx.arena());
    auto simplified = ctx.simplify(e.value());
    ASSERT_TRUE(simplified.is_ok());
    auto text = to_round_trip_text(simplified.value());
    ASSERT_TRUE(text.is_ok());
    EXPECT_EQ(text.value(), "x");
}

TEST_F(AbsSignTest, AbsNegative) {
    auto x = ctx.arena().make<Symbol>("x");
    ctx.assumptions().assume_greater(ctx.arena().make<IntegerLit>(0), x); // x < 0
    
    auto e = parse_expr("abs(x)", ctx.arena());
    auto simplified = ctx.simplify(e.value());
    ASSERT_TRUE(simplified.is_ok());
    auto text = to_round_trip_text(simplified.value());
    ASSERT_TRUE(text.is_ok());
    EXPECT_EQ(text.value(), "-x");
}

TEST_F(AbsSignTest, AbsSquare) {
    auto x = ctx.arena().make<Symbol>("x");
    ctx.assumptions().assume_real(expr_ref<Symbol>(x));
    auto e = parse_expr("abs(x^2)", ctx.arena());
    auto simplified = ctx.simplify(e.value());
    ASSERT_TRUE(simplified.is_ok());
    auto text = to_round_trip_text(simplified.value());
    ASSERT_TRUE(text.is_ok());
    EXPECT_EQ(text.value(), "x^2");
}

TEST_F(AbsSignTest, AbsAbs) {
    auto e = parse_expr("abs(abs(x))", ctx.arena());
    auto simplified = ctx.simplify(e.value());
    ASSERT_TRUE(simplified.is_ok());
    auto text = to_round_trip_text(simplified.value());
    ASSERT_TRUE(text.is_ok());
    EXPECT_EQ(text.value(), "abs(x)");
}

TEST_F(AbsSignTest, SignPositive) {
    auto x = ctx.arena().make<Symbol>("x");
    ctx.assumptions().assume_positive(expr_ref<Symbol>(x));
    
    auto e = parse_expr("sign(x)", ctx.arena());
    auto simplified = ctx.simplify(e.value());
    ASSERT_TRUE(simplified.is_ok());
    auto text = to_round_trip_text(simplified.value());
    ASSERT_TRUE(text.is_ok());
    EXPECT_EQ(text.value(), "1");
}

TEST_F(AbsSignTest, SignNegative) {
    auto x = ctx.arena().make<Symbol>("x");
    ctx.assumptions().assume_greater(ctx.arena().make<IntegerLit>(0), x); // x < 0
    
    auto e = parse_expr("sign(x)", ctx.arena());
    auto simplified = ctx.simplify(e.value());
    ASSERT_TRUE(simplified.is_ok());
    auto text = to_round_trip_text(simplified.value());
    ASSERT_TRUE(text.is_ok());
    EXPECT_EQ(text.value(), "-1");
}
