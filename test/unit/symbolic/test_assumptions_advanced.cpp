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

// A32 prerequisite: 0^x = 0 when x is PROVABLY positive (assumption), exact and
// unconditional. Uses the proven assumption, not an A31 side-condition.
TEST_F(AdvancedAssumptionsTest, ZeroPower_PositiveExponent_FoldsToZero) {
    auto x = sym("x");
    ctx.assumptions().assume_positive(*expr_cast<Symbol>(x));
    ExprPtr zero = ctx.arena().make<IntegerLit>(BigInt(0));
    ExprPtr pow = ctx.arena().make<Binary>(BinaryOp::Pow, zero, x);
    auto s = ctx.simplify(pow);
    ASSERT_TRUE(s.is_ok());
    const auto* il = expr_cast<IntegerLit>(s.value());
    EXPECT_TRUE(il != nullptr && il->value.is_zero())
        << "0^x with x>0 did not fold to 0";
}

// Literal positive non-integer exponent: 0^(1/2) = 0 (now covered by the
// is_known_positive branch, previously integer-only).
TEST_F(AdvancedAssumptionsTest, ZeroPower_PositiveRationalExponent_FoldsToZero) {
    ExprPtr zero = ctx.arena().make<IntegerLit>(BigInt(0));
    ExprPtr half = ctx.arena().make<RationalLit>(BigInt(1), BigInt(2));
    auto s = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Pow, zero, half));
    ASSERT_TRUE(s.is_ok());
    const auto* il = expr_cast<IntegerLit>(s.value());
    EXPECT_TRUE(il != nullptr && il->value.is_zero());
}

// Soundness guard: with NO sign knowledge, 0^x must NOT fold to 0 (0^0=1,
// 0^negative undefined) — stays symbolic.
TEST_F(AdvancedAssumptionsTest, ZeroPower_UnknownExponent_StaysSymbolic) {
    auto x = sym("x");
    ExprPtr zero = ctx.arena().make<IntegerLit>(BigInt(0));
    auto s = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Pow, zero, x));
    ASSERT_TRUE(s.is_ok());
    EXPECT_EQ(expr_cast<IntegerLit>(s.value()), nullptr)
        << "0^x with unknown-sign x wrongly folded";
}

TEST_F(AdvancedAssumptionsTest, NonlinearInference) {
    auto x = sym("x");
    ctx.assumptions().assume_greater(x, ctx.arena().make<IntegerLit>(1));
    
    auto x2 = ctx.arena().make<Binary>(BinaryOp::Pow, x, ctx.arena().make<IntegerLit>(2));
    (void)x2; // EXPECT_TRUE pending: x > 1 => x^2 > 1 (not yet wired)
}
