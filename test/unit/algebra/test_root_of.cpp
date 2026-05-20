#include "cas/symbolic.hpp"
#include "cas/algebra.hpp"
#include "cas/numeric.hpp"
#include "cas/ast_debug.hpp"
#include <gtest/gtest.h>

namespace cas::test {

class RootOfTest : public ::testing::Test {
protected:
    void SetUp() override {
        ctx = std::make_unique<symbolic::CASContext>();
    }
    std::unique_ptr<symbolic::CASContext> ctx;
};

TEST_F(RootOfTest, LinearSimplification) {
    // RootOf(x - 2, x, 0) -> 2
    auto x = ctx->arena().make<Symbol>("x");
    auto poly = ctx->arena().make<Binary>(BinaryOp::Sub, x, ctx->arena().make<IntegerLit>(BigInt(2)));
    auto root = ctx->arena().make<RootOf>(poly, Symbol("x"), 0U);
    
    auto simplified = ctx->simplify(root);
    ASSERT_TRUE(simplified.is_ok());
    
    auto expected = ctx->arena().make<IntegerLit>(BigInt(2));
    EXPECT_TRUE(structural_equal(simplified.value(), expected));
}

TEST_F(RootOfTest, AlgebraicEquality) {
    // RootOf(x^2 - 2, x, 0)^2 -> 2
    auto x = ctx->arena().make<Symbol>("x");
    auto x2 = ctx->arena().make<Binary>(BinaryOp::Pow, x, ctx->arena().make<IntegerLit>(BigInt(2)));
    auto poly = ctx->arena().make<Binary>(BinaryOp::Sub, x2, ctx->arena().make<IntegerLit>(BigInt(2)));
    auto root = ctx->arena().make<RootOf>(poly, Symbol("x"), 0U);
    auto root2 = ctx->arena().make<Binary>(BinaryOp::Pow, root, ctx->arena().make<IntegerLit>(BigInt(2)));
    
    auto simplified = ctx->simplify(root2);
    ASSERT_TRUE(simplified.is_ok());
    
    auto expected = ctx->arena().make<IntegerLit>(BigInt(2));
    EXPECT_TRUE(structural_equal(simplified.value(), expected)) << "Got: " << debug_print(simplified.value());
}

TEST_F(RootOfTest, AlgebraicEqualityCubic) {
    // RootOf(x^3 - x - 1, x, 0)^3 -> RootOf(x^3 - x - 1, x, 0) + 1
    auto x = ctx->arena().make<Symbol>("x");
    auto x3 = ctx->arena().make<Binary>(BinaryOp::Pow, x, ctx->arena().make<IntegerLit>(BigInt(3)));
    auto poly = ctx->arena().make<Binary>(BinaryOp::Sub, 
                    ctx->arena().make<Binary>(BinaryOp::Sub, x3, x),
                    ctx->arena().make<IntegerLit>(BigInt(1)));
    auto root = ctx->arena().make<RootOf>(poly, Symbol("x"), 0U);
    auto root3 = ctx->arena().make<Binary>(BinaryOp::Pow, root, ctx->arena().make<IntegerLit>(BigInt(3)));
    
    auto simplified = ctx->simplify(root3);
    ASSERT_TRUE(simplified.is_ok());
    
    auto expected = ctx->arena().make<Binary>(BinaryOp::Add, root, ctx->arena().make<IntegerLit>(BigInt(1)));
    auto eq = symbolic::mathematically_equal(simplified.value(), expected, *ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value()) << "Got: " << debug_print(simplified.value());
}

TEST_F(RootOfTest, NumericalEvaluation) {
    // RootOf(x^2 - 2, x, 0) -> sqrt(2) approx 1.41421356
    auto x = ctx->arena().make<Symbol>("x");
    auto x2 = ctx->arena().make<Binary>(BinaryOp::Pow, x, ctx->arena().make<IntegerLit>(BigInt(2)));
    auto poly = ctx->arena().make<Binary>(BinaryOp::Sub, x2, ctx->arena().make<IntegerLit>(BigInt(2)));
    auto root = ctx->arena().make<RootOf>(poly, Symbol("x"), 0U); // index 0 is usually positive root
    
    auto val = numeric::eval(root);
    ASSERT_TRUE(val.is_ok());
    EXPECT_NEAR(val.value(), 1.41421356, 1e-7);
    
    auto root_neg = ctx->arena().make<RootOf>(poly, Symbol("x"), 1U);
    auto val_neg = numeric::eval(root_neg);
    ASSERT_TRUE(val_neg.is_ok());
    EXPECT_NEAR(val_neg.value(), -1.41421356, 1e-7);
}

} // namespace cas::test
