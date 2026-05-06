#include <gtest/gtest.h>
#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::algebra;
using namespace cas::symbolic;

namespace {
[[nodiscard]] Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = cas::Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}
}

class ResultantTest : public ::testing::Test {
protected:
    CASContext ctx;
};

TEST_F(ResultantTest, QuadraticDiscriminant) {
    auto x = ctx.arena().make<Symbol>("x");
    // p = x^2 + x + 1 => Disc = 1^2 - 4(1)(1) = -3
    auto p = parse_expr("x^2 + x + 1", ctx.arena());
    ASSERT_TRUE(p.is_ok());

    auto disc = polynomial_discriminant(p.value(), expr_ref<Symbol>(x), ctx);
    ASSERT_TRUE(disc.is_ok());
    
    auto simplified = ctx.simplify(disc.value());
    ASSERT_TRUE(simplified.is_ok());
    
    const auto* lit = expr_cast<IntegerLit>(simplified.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(-3));
}

TEST_F(ResultantTest, CubicDiscriminant) {
    auto x = ctx.arena().make<Symbol>("x");
    // p = x^3 + a*x + b => Disc = -4a^3 - 27b^2
    auto p = parse_expr("x^3 + a*x + b", ctx.arena());
    ASSERT_TRUE(p.is_ok());

    auto disc = polynomial_discriminant(p.value(), expr_ref<Symbol>(x), ctx);
    ASSERT_TRUE(disc.is_ok());
    
    // Check if it matches -4a^3 - 27b^2 when expanded
    auto expanded = expand(disc.value(), ctx);
    ASSERT_TRUE(expanded.is_ok());
    
    // We expect -4*a^3 - 27*b^2
    // Let's check by substitution a=0, b=1 => Disc = -27
    auto a = ctx.arena().make<Symbol>("a");
    auto b = ctx.arena().make<Symbol>("b");
    
    auto sub1 = ctx.substitute(expanded.value(), expr_ref<Symbol>(a), ctx.arena().make<IntegerLit>(0)).value();
    auto sub2 = ctx.substitute(sub1, expr_ref<Symbol>(b), ctx.arena().make<IntegerLit>(1)).value();
    auto final_val = ctx.simplify(sub2).value();
    
    const auto* lit = expr_cast<IntegerLit>(final_val);
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(-27));
}

TEST_F(ResultantTest, ResultantSimple) {
    auto x = ctx.arena().make<Symbol>("x");
    // f = x - a, g = x - b => Res = b - a (wait, standard is a - b or b - a?)
    // Sylvester: [ 1 -a ]
    //            [ 1 -b ]
    // Det = -b - (-a) = a - b
    auto f = parse_expr("x - a", ctx.arena());
    auto g = parse_expr("x - b", ctx.arena());

    auto res = polynomial_resultant(f.value(), g.value(), expr_ref<Symbol>(x), ctx);
    ASSERT_TRUE(res.is_ok());

    auto a = ctx.arena().make<Symbol>("a");
    auto b = ctx.arena().make<Symbol>("b");
    auto sub1 = ctx.substitute(res.value(), expr_ref<Symbol>(a), ctx.arena().make<IntegerLit>(1)).value();
    auto sub2 = ctx.substitute(sub1, expr_ref<Symbol>(b), ctx.arena().make<IntegerLit>(0)).value();
    auto final_val = ctx.simplify(sub2).value();

    const auto* lit = expr_cast<IntegerLit>(final_val);
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(1)); // 1 - 0 = 1
}

// L1-15: resultant output must be normalized (IntegerLit, not (-1)*N or unevaluated)
TEST_F(ResultantTest, L1_15_ResultantOutputNormalized) {
    // Res(x+1, x+2) = (2-(-1)) err: Res = lc(g)^(deg f) * prod_{f(xi)=0} g(xi)
    // f = x+1 root -1; g(-1) = -1+2 = 1; Res = 1^1 * 1 = 1
    auto f = parse_expr("x + 1", ctx.arena());
    auto g = parse_expr("x + 2", ctx.arena());
    ASSERT_TRUE(f.is_ok());
    ASSERT_TRUE(g.is_ok());

    Symbol x("x");
    auto res = polynomial_resultant(f.value(), g.value(), x, ctx);
    ASSERT_TRUE(res.is_ok());

    // Must be IntegerLit(1) directly, not a Mul/Neg tree
    const auto* ilit = expr_cast<IntegerLit>(res.value());
    ASSERT_NE(ilit, nullptr) << "L1-15: resultant must normalize to IntegerLit, not unevaluated tree";
    EXPECT_EQ(ilit->value, BigInt(1));
}

// L1-15: discriminant of x^2+2x+1 = (x+1)^2 must be 0 as IntegerLit
TEST_F(ResultantTest, L1_15_DiscriminantNormalized) {
    auto p = parse_expr("x^2 + 2*x + 1", ctx.arena());
    ASSERT_TRUE(p.is_ok());

    Symbol x("x");
    auto disc = polynomial_discriminant(p.value(), x, ctx);
    ASSERT_TRUE(disc.is_ok());

    // Disc of (x+1)^2 = 0
    const auto* ilit = expr_cast<IntegerLit>(disc.value());
    ASSERT_NE(ilit, nullptr) << "L1-15: discriminant must normalize to IntegerLit";
    EXPECT_EQ(ilit->value, BigInt(0));
}
