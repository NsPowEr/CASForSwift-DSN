// test_csolve_linear.cpp — Coverage for the linear fast-path in csolve.
//
// Closes HC-F4-GOSPER-CONSTANT-HANG (csolve audit branch):
// before the fix, csolve dispatched every system through F4 (Buchberger).
// On linear underdetermined systems (rank deficit) the Buchberger pipeline
// did not terminate in any reasonable time bound. The Gosper polynomial
// ansatz (one free additive constant) is the canonical trigger of this
// path. The linear fast-path computes a particular solution by setting
// every free variable to 0, returning a Matrix(1, n) of rationals.

#include <gtest/gtest.h>
#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

namespace cas::algebra {

class CsolveLinearTest : public ::testing::Test {
protected:
    void SetUp() override {
        ctx = std::make_unique<symbolic::CASContext>();
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

namespace {

ExprPtr int_lit(symbolic::CASContext& c, long long v) {
    return c.arena().make<IntegerLit>(BigInt(v));
}

}  // namespace

TEST_F(CsolveLinearTest, UnderdeterminedReturnsParametricParticular) {
    // System: u_1 + u_2 = 0, 2*u_2 - 1 = 0, over vars {u_0, u_1, u_2}.
    // Rank 2, free var u_0 → set to 0.
    // Particular solution: u_0 = 0, u_1 = -1/2, u_2 = 1/2.
    Symbol u0("u_0"), u1("u_1"), u2("u_2");

    ExprPtr eq1 = ctx->arena().make<Sum>(std::vector<ExprPtr>{
        ctx->arena().make<Symbol>(u1), ctx->arena().make<Symbol>(u2)});
    ExprPtr eq2 = ctx->arena().make<Sum>(std::vector<ExprPtr>{
        ctx->arena().make<Binary>(BinaryOp::Mul, int_lit(*ctx, 2), ctx->arena().make<Symbol>(u2)),
        int_lit(*ctx, -1)});

    auto eqs_matrix = ctx->arena().make<Matrix>(2U, 1U, std::vector<ExprPtr>{eq1, eq2});
    auto vars_matrix = ctx->arena().make<Matrix>(3U, 1U, std::vector<ExprPtr>{
        ctx->arena().make<Symbol>(u0), ctx->arena().make<Symbol>(u1), ctx->arena().make<Symbol>(u2)});

    auto sol = csolve(eqs_matrix, vars_matrix, *ctx);
    ASSERT_TRUE(sol.is_ok()) << sol.error().message;

    const auto* m = expr_cast<Matrix>(sol.value());
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->rows, 1U);
    EXPECT_EQ(m->cols, 3U);
    ASSERT_EQ(m->elements.size(), 3U);

    // u_0 = 0 (free).
    const auto* u0v = expr_cast<IntegerLit>(m->elements[0]);
    ASSERT_NE(u0v, nullptr);
    EXPECT_EQ(u0v->value, BigInt(0));

    // u_1 = -1/2.
    const auto* u1v = expr_cast<RationalLit>(m->elements[1]);
    ASSERT_NE(u1v, nullptr);
    EXPECT_EQ(u1v->numerator, BigInt(-1));
    EXPECT_EQ(u1v->denominator, BigInt(2));

    // u_2 = 1/2.
    const auto* u2v = expr_cast<RationalLit>(m->elements[2]);
    ASSERT_NE(u2v, nullptr);
    EXPECT_EQ(u2v->numerator, BigInt(1));
    EXPECT_EQ(u2v->denominator, BigInt(2));
}

TEST_F(CsolveLinearTest, InconsistentLinearReturnsEmptyMatrix) {
    // System: x + y = 1, 2*(x+y) = 1 (i.e. 2x + 2y - 1 = 0).
    // Inconsistent (first says x+y=1, second says x+y=1/2).
    // Expected: empty Matrix (no solutions).
    Symbol x("x"), y("y");

    ExprPtr eq1 = ctx->arena().make<Sum>(std::vector<ExprPtr>{
        ctx->arena().make<Symbol>(x), ctx->arena().make<Symbol>(y),
        int_lit(*ctx, -1)});
    ExprPtr eq2 = ctx->arena().make<Sum>(std::vector<ExprPtr>{
        ctx->arena().make<Binary>(BinaryOp::Mul, int_lit(*ctx, 2), ctx->arena().make<Symbol>(x)),
        ctx->arena().make<Binary>(BinaryOp::Mul, int_lit(*ctx, 2), ctx->arena().make<Symbol>(y)),
        int_lit(*ctx, -1)});

    auto eqs_matrix = ctx->arena().make<Matrix>(2U, 1U, std::vector<ExprPtr>{eq1, eq2});
    auto vars_matrix = ctx->arena().make<Matrix>(2U, 1U, std::vector<ExprPtr>{
        ctx->arena().make<Symbol>(x), ctx->arena().make<Symbol>(y)});

    auto sol = csolve(eqs_matrix, vars_matrix, *ctx);
    ASSERT_TRUE(sol.is_ok()) << sol.error().message;

    const auto* m = expr_cast<Matrix>(sol.value());
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->rows, 0U);
}

TEST_F(CsolveLinearTest, SquareDeterminedLinearSolvedExactly) {
    // System: 2x + 3y = 7, x - y = 1.
    // Solution: x = 2, y = 1.
    Symbol x("x"), y("y");

    // 2x + 3y - 7
    ExprPtr eq1 = ctx->arena().make<Sum>(std::vector<ExprPtr>{
        ctx->arena().make<Binary>(BinaryOp::Mul, int_lit(*ctx, 2), ctx->arena().make<Symbol>(x)),
        ctx->arena().make<Binary>(BinaryOp::Mul, int_lit(*ctx, 3), ctx->arena().make<Symbol>(y)),
        int_lit(*ctx, -7)});
    // x - y - 1
    ExprPtr eq2 = ctx->arena().make<Sum>(std::vector<ExprPtr>{
        ctx->arena().make<Symbol>(x),
        ctx->arena().make<Binary>(BinaryOp::Mul, int_lit(*ctx, -1), ctx->arena().make<Symbol>(y)),
        int_lit(*ctx, -1)});

    auto eqs_matrix = ctx->arena().make<Matrix>(2U, 1U, std::vector<ExprPtr>{eq1, eq2});
    auto vars_matrix = ctx->arena().make<Matrix>(2U, 1U, std::vector<ExprPtr>{
        ctx->arena().make<Symbol>(x), ctx->arena().make<Symbol>(y)});

    auto sol = csolve(eqs_matrix, vars_matrix, *ctx);
    ASSERT_TRUE(sol.is_ok()) << sol.error().message;

    const auto* m = expr_cast<Matrix>(sol.value());
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->rows, 1U);
    EXPECT_EQ(m->cols, 2U);

    const auto* xv = expr_cast<IntegerLit>(m->elements[0]);
    ASSERT_NE(xv, nullptr);
    EXPECT_EQ(xv->value, BigInt(2));
    const auto* yv = expr_cast<IntegerLit>(m->elements[1]);
    ASSERT_NE(yv, nullptr);
    EXPECT_EQ(yv->value, BigInt(1));
}

}  // namespace cas::algebra
