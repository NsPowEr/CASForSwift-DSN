#include <gtest/gtest.h>
#include "cas/symbolic.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/formatter.hpp"

using namespace cas;
using namespace cas::symbolic;

class AdvancedSimplificationTest : public ::testing::Test {
protected:
    CASContext ctx;
};

TEST_F(AdvancedSimplificationTest, TranscendentalIdentities) {
    auto& arena = ctx.arena();
    
    // ln(e^x) -> x
    auto x = arena.make<Symbol>("x");
    auto exp_x = arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{x});
    auto ln_exp_x = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{exp_x});
    auto res1 = ctx.simplify(ln_exp_x);
    ASSERT_TRUE(res1.is_ok());
    EXPECT_EQ(res1.value(), x);

    // exp(ln(x)) -> x
    auto ln_x = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{x});
    auto exp_ln_x = arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{ln_x});
    auto res2 = ctx.simplify(exp_ln_x);
    ASSERT_TRUE(res2.is_ok());
    EXPECT_EQ(res2.value(), x);
}

TEST_F(AdvancedSimplificationTest, LogExpansion) {
    auto& arena = ctx.arena();
    auto x = arena.make<Symbol>("x");
    auto y = arena.make<Symbol>("y");
    
    ctx.assumptions().assume_positive(Symbol("x"));
    ctx.assumptions().assume_positive(Symbol("y"));
    
    // ln(x*y) -> ln(x) + ln(y)
    auto prod = arena.make<Product>(std::vector<ExprPtr>{x, y});
    auto ln_prod = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{prod});
    auto res = ctx.simplify(ln_prod);
    ASSERT_TRUE(res.is_ok());
    
    // Check if result is Sum of two Ln
    EXPECT_EQ(res.value()->kind, ExprKind::Sum);
}

TEST_F(AdvancedSimplificationTest, ExpExpansion) {
    auto& arena = ctx.arena();
    auto x = arena.make<Symbol>("x");
    auto y = arena.make<Symbol>("y");
    
    // exp(x+y) -> exp(x) * exp(y)
    auto sum = arena.make<Sum>(std::vector<ExprPtr>{x, y});
    auto exp_sum = arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{sum});
    auto res = ctx.simplify(exp_sum);
    ASSERT_TRUE(res.is_ok());
    
    // Check if result is Product of two Exp
    EXPECT_EQ(res.value()->kind, ExprKind::Product);
}

TEST_F(AdvancedSimplificationTest, Denesting) {
    auto& arena = ctx.arena();
    
    // sqrt(3 + 2*sqrt(2)) -> 1 + sqrt(2)
    // 3 + 2*sqrt(2)
    auto sqrt2 = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{arena.make<IntegerLit>(BigInt(2))});
    auto two_sqrt2 = arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(2)), sqrt2);
    auto inner = arena.make<Sum>(std::vector<ExprPtr>{arena.make<IntegerLit>(BigInt(3)), two_sqrt2});
    auto outer = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{inner});
    
    auto res = ctx.simplify(outer);
    ASSERT_TRUE(res.is_ok());
    
    // Expected: 1 + sqrt(2) (canonical order is sqrt(2) + 1)
    auto expected_raw = arena.make<Sum>(std::vector<ExprPtr>{arena.make<IntegerLit>(BigInt(1)), sqrt2});
    auto expected = ctx.simplify(expected_raw);
    ASSERT_TRUE(expected.is_ok());
    
    EXPECT_TRUE(structural_equal(res.value(), expected.value()));
}
