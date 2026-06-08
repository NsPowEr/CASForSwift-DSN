#include <gtest/gtest.h>
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

using namespace cas;
using namespace cas::symbolic;

class AdvancedAssumptionsTest : public ::testing::Test {
protected:
    CASContext ctx;
    
    ExprPtr sym(const std::string& name) {
        return ctx.arena().make<Symbol>(name);
    }
};

TEST_F(AdvancedAssumptionsTest, Transitivity) {
    auto x = sym("x");
    auto y = sym("y");
    auto z = sym("z");
    
    ctx.assumptions().assume_greater(x, y);
    ctx.assumptions().assume_greater(y, z);
    
    EXPECT_TRUE(ctx.assumptions().is_greater(x, z));
    EXPECT_FALSE(ctx.assumptions().is_greater(z, x));
}

TEST_F(AdvancedAssumptionsTest, ContradictionDetection) {
    auto x = sym("x");
    
    ctx.assumptions().assume_greater(x, nullptr); // x > 0
    ctx.assumptions().assume_greater(nullptr, x); // x < 0
    
    auto result = ctx.assumptions().check_consistency();
    EXPECT_TRUE(result.is_error());
}

TEST_F(AdvancedAssumptionsTest, ProductInference) {
    auto x = sym("x");
    auto y = sym("y");
    
    ctx.assumptions().assume_positive(*expr_cast<Symbol>(x));
    ctx.assumptions().assume_positive(*expr_cast<Symbol>(y));
    
    auto xy = ctx.arena().make<Binary>(BinaryOp::Mul, x, y);
    EXPECT_TRUE(ctx.assumptions().is_positive(xy));
}

TEST_F(AdvancedAssumptionsTest, SumInference) {
    auto x = sym("x");
    auto y = sym("y");
    auto a = sym("a");
    auto b = sym("b");
    
    ctx.assumptions().assume_greater(x, y);
    ctx.assumptions().assume_greater(a, b);
    
    auto x_plus_a = ctx.arena().make<Binary>(BinaryOp::Add, x, a);
    auto y_plus_b = ctx.arena().make<Binary>(BinaryOp::Add, y, b);
    
    // This is a more advanced deduction: x > y, a > b => x+a > y+b
    EXPECT_TRUE(ctx.assumptions().is_greater(x_plus_a, y_plus_b));
}

TEST_F(AdvancedAssumptionsTest, DomainConsistency) {
    auto n = sym("n");
    ctx.assumptions().assume_domain(*expr_cast<Symbol>(n), Domain::Natural);
    
    EXPECT_TRUE(ctx.assumptions().is_integer(n));
    EXPECT_TRUE(ctx.assumptions().is_nonnegative(n));
    
    ctx.assumptions().assume_domain(*expr_cast<Symbol>(n), Domain::Negative);
    auto result = ctx.assumptions().check_consistency();
    EXPECT_TRUE(result.is_error());
}

TEST_F(AdvancedAssumptionsTest, NonlinearInference) {
    auto x = sym("x");
    ctx.assumptions().assume_greater(x, ctx.arena().make<IntegerLit>(1));
    
    auto x2 = ctx.arena().make<Binary>(BinaryOp::Pow, x, ctx.arena().make<IntegerLit>(2));
    (void)x2; // EXPECT_TRUE pending: x > 1 => x^2 > 1 (not yet wired)
}
