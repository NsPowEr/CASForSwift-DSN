#include <gtest/gtest.h>
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

using namespace cas::symbolic;
using namespace cas;

TEST(AssumptionsAdvancedTest, DeductionChain3Vars) {
    CASContext ctx;
    auto& arena = ctx.arena();
    
    auto x = arena.make<Symbol>("x");
    auto y = arena.make<Symbol>("y");
    auto z = arena.make<Symbol>("z");
    
    ctx.assumptions().assume_greater(x, y); // x > y
    ctx.assumptions().assume_greater(y, z); // y > z
    ctx.assumptions().assume_greater(z, nullptr); // z > 0
    
    EXPECT_TRUE(ctx.assumptions().is_positive(x));
    EXPECT_TRUE(ctx.assumptions().is_positive(y));
    EXPECT_TRUE(ctx.assumptions().is_positive(z));
    EXPECT_TRUE(ctx.assumptions().is_greater(x, z));
    EXPECT_TRUE(ctx.assumptions().is_greater(x, nullptr));
}

TEST(AssumptionsAdvancedTest, LinearDeduction) {
    CASContext ctx;
    auto& arena = ctx.arena();
    
    auto x = arena.make<Symbol>("x");
    auto y = arena.make<Symbol>("y");
    
    ctx.assumptions().assume_greater(x, y); // x > y
    
    // Simplification of |x - y| -> x - y
    auto x_minus_y = arena.make<Binary>(BinaryOp::Sub, x, y);
    auto abs_expr = arena.make<FuncCall>("abs", std::vector<ExprPtr>{x_minus_y});
    
    auto result = ctx.simplify(abs_expr);
    ASSERT_TRUE(result.is_ok());
    // Simplifier canonicalizes x-y to Sum{x,-y}, so use mathematical equality
    auto eq = mathematically_equal(result.value(), x_minus_y, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

TEST(AssumptionsAdvancedTest, MixedChain) {
    CASContext ctx;
    auto& arena = ctx.arena();
    
    auto a = arena.make<Symbol>("a");
    auto b = arena.make<Symbol>("b");
    
    ctx.assumptions().assume_greater_equal(a, b); // a >= b
    ctx.assumptions().assume_greater(b, nullptr); // b > 0
    
    EXPECT_TRUE(ctx.assumptions().is_positive(a));
    EXPECT_TRUE(ctx.assumptions().is_greater_equal(a, nullptr));
}

TEST(AssumptionsAdvancedTest, SumOfPositive) {
    CASContext ctx;
    auto& arena = ctx.arena();
    
    auto x = arena.make<Symbol>("x");
    auto y = arena.make<Symbol>("y");
    
    ctx.assumptions().assume_positive(expr_ref<Symbol>(x));
    ctx.assumptions().assume_positive(expr_ref<Symbol>(y));
    
    auto sum = arena.make<Sum>(std::vector<ExprPtr>{x, y});
    EXPECT_TRUE(ctx.assumptions().is_positive(sum));
}

TEST(AssumptionsAdvancedTest, ProductOfPositive) {
    CASContext ctx;
    auto& arena = ctx.arena();
    
    auto x = arena.make<Symbol>("x");
    auto y = arena.make<Symbol>("y");
    
    ctx.assumptions().assume_positive(expr_ref<Symbol>(x));
    ctx.assumptions().assume_positive(expr_ref<Symbol>(y));
    
    auto prod = arena.make<Product>(std::vector<ExprPtr>{x, y});
    EXPECT_TRUE(ctx.assumptions().is_positive(prod));
}

// P1-005: assumption-driven simplifications

TEST(P1_005_AssumptionsSimplifier, AbsEliminatedWhenPositive) {
    CASContext ctx;
    auto& arena = ctx.arena();
    auto x = arena.make<Symbol>("x");
    ctx.assumptions().assume_positive(expr_ref<Symbol>(x));

    auto abs_x = arena.make<FuncCall>("abs", std::vector<ExprPtr>{x});
    auto result = ctx.simplify(abs_x);
    ASSERT_TRUE(result.is_ok());
    // abs(x) → x when x > 0
    EXPECT_EQ(result.value(), x);
}

TEST(P1_005_AssumptionsSimplifier, SqrtSquareToSelfWhenPositive) {
    CASContext ctx;
    auto& arena = ctx.arena();
    auto x = arena.make<Symbol>("x");
    ctx.assumptions().assume_positive(expr_ref<Symbol>(x));

    // sqrt(x^2) → x when x > 0
    auto x2 = arena.make<Binary>(BinaryOp::Pow, x, arena.make<IntegerLit>(BigInt(2)));
    auto sqrt_x2 = arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{x2});
    auto result = ctx.simplify(sqrt_x2);
    ASSERT_TRUE(result.is_ok());
    EXPECT_EQ(result.value(), x);
}

TEST(P1_005_AssumptionsSimplifier, SqrtSquareToAbsWhenUnknownSign) {
    CASContext ctx;
    auto& arena = ctx.arena();
    auto x = arena.make<Symbol>("x");
    // No assumption on x — sqrt(x^2) = |x|
    auto x2 = arena.make<Binary>(BinaryOp::Pow, x, arena.make<IntegerLit>(BigInt(2)));
    auto sqrt_x2 = arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{x2});
    auto result = ctx.simplify(sqrt_x2);
    ASSERT_TRUE(result.is_ok());
    const auto* func = expr_cast<FuncCall>(result.value());
    ASSERT_NE(func, nullptr);
    EXPECT_EQ(func->func_id, BuiltinOp::Abs);
}

TEST(P1_005_AssumptionsSimplifier, LnAbsCollapsesWhenPositive) {
    CASContext ctx;
    auto& arena = ctx.arena();
    auto x = arena.make<Symbol>("x");
    ctx.assumptions().assume_positive(expr_ref<Symbol>(x));

    // ln(abs(x)) → ln(x): abs(x) simplifies to x, then ln(x) stays
    auto abs_x = arena.make<FuncCall>("abs", std::vector<ExprPtr>{x});
    auto ln_abs_x = arena.make<FuncCall>("ln", std::vector<ExprPtr>{abs_x});
    auto result = ctx.simplify(ln_abs_x);
    ASSERT_TRUE(result.is_ok());
    const auto* ln = expr_cast<FuncCall>(result.value());
    ASSERT_NE(ln, nullptr);
    EXPECT_EQ(ln->func_id, BuiltinOp::Ln);
    EXPECT_EQ(ln->args.size(), 1U);
    EXPECT_EQ(ln->args[0], x);
}

TEST(P1_005_AssumptionsSimplifier, AbsDifferenceEliminatedByRelation) {
    CASContext ctx;
    auto& arena = ctx.arena();
    auto x = arena.make<Symbol>("x");
    auto y = arena.make<Symbol>("y");
    ctx.assumptions().assume_greater(x, y);

    // abs(x - y) → x - y when x > y
    auto x_minus_y = arena.make<Binary>(BinaryOp::Sub, x, y);
    auto abs_diff = arena.make<FuncCall>("abs", std::vector<ExprPtr>{x_minus_y});
    auto result = ctx.simplify(abs_diff);
    ASSERT_TRUE(result.is_ok());
    auto eq = mathematically_equal(result.value(), x_minus_y, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}
