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
    // structural equality check
    EXPECT_EQ(result.value(), x_minus_y);
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
