#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include <gtest/gtest.h>

namespace cas::symbolic {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tok = Lexer(input).tokenize();
    if (tok.is_error()) return fail<ExprPtr>(tok.error());
    return Parser(tok.value(), arena).parse();
}

TEST(AssumptionsRelationalTest, BinaryRelationalAssumptions) {
    CASContext ctx;
    auto cond = parse_expr("x > 0", ctx.arena());
    ASSERT_TRUE(cond.is_ok()) << cond.error().message;
    
    ctx.assumptions().assume(cond.value());
    
    auto x = ctx.arena().make<Symbol>("x");
    EXPECT_TRUE(ctx.assumptions().is_positive(x));
    EXPECT_TRUE(ctx.assumptions().is_nonzero(x));
    EXPECT_FALSE(ctx.assumptions().is_negative(x));
}

TEST(AssumptionsRelationalTest, LessEqualAssumption) {
    CASContext ctx;
    auto cond = parse_expr("x <= 5", ctx.arena());
    ASSERT_TRUE(cond.is_ok());
    
    ctx.assumptions().assume(cond.value());
    
    auto x = ctx.arena().make<Symbol>("x");
    auto five = ctx.arena().make<IntegerLit>(BigInt(5));
    auto six = ctx.arena().make<IntegerLit>(BigInt(6));
    
    EXPECT_TRUE(ctx.assumptions().is_greater_equal(five, x));
    EXPECT_TRUE(ctx.assumptions().is_greater(six, x));
}

TEST(AssumptionsRelationalTest, TransitiveProperty) {
    CASContext ctx;
    ctx.assumptions().assume(parse_expr("x > y", ctx.arena()).value());
    ctx.assumptions().assume(parse_expr("y > z", ctx.arena()).value());
    
    auto x = ctx.arena().make<Symbol>("x");
    auto z = ctx.arena().make<Symbol>("z");
    
    EXPECT_TRUE(ctx.assumptions().is_greater(x, z));
}

// L1-09: x>0, y>0 → x*y>0 must be inferred automatically (no hardcode)
TEST(AssumptionsRelationalTest, L1_09_ProductOfPositivesIsPositive) {
    CASContext ctx;
    ctx.assumptions().assume(parse_expr("x > 0", ctx.arena()).value());
    ctx.assumptions().assume(parse_expr("y > 0", ctx.arena()).value());

    auto x = ctx.arena().make<Symbol>("x");
    auto y = ctx.arena().make<Symbol>("y");
    // Binary(Mul, x, y)
    auto xy = ctx.arena().make<Binary>(BinaryOp::Mul, x, y);
    EXPECT_TRUE(ctx.assumptions().is_positive(xy)) << "x*y must be positive when x>0 and y>0";
    EXPECT_TRUE(ctx.assumptions().is_nonzero(xy)) << "x*y must be nonzero when x>0 and y>0";
}

// L1-09: is_nonzero derives from relational assume(x>0)
TEST(AssumptionsRelationalTest, L1_09_NonzeroDerivesFromGreaterThanZero) {
    CASContext ctx;
    ctx.assumptions().assume(parse_expr("x > 0", ctx.arena()).value());

    auto x = ctx.arena().make<Symbol>("x");
    EXPECT_TRUE(ctx.assumptions().is_positive(x));
    EXPECT_TRUE(ctx.assumptions().is_nonzero(x)) << "is_nonzero must derive from relational graph x>0";
}

// L1-09: transitivity with 3 hops
TEST(AssumptionsRelationalTest, L1_09_TransitiveChain3Hops) {
    CASContext ctx;
    ctx.assumptions().assume(parse_expr("a > b", ctx.arena()).value());
    ctx.assumptions().assume(parse_expr("b > c", ctx.arena()).value());
    ctx.assumptions().assume(parse_expr("c > d", ctx.arena()).value());

    auto a = ctx.arena().make<Symbol>("a");
    auto d = ctx.arena().make<Symbol>("d");
    EXPECT_TRUE(ctx.assumptions().is_greater(a, d)) << "a>b>c>d must imply a>d transitively";
}

} // namespace
} // namespace cas::symbolic
