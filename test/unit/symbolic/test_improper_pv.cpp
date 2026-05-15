// L2-11: tests for classify_improper_convergence and cauchy_principal_value.

#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include "cas/improper_integral.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <string>

namespace cas::test {

namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

}  // namespace

class ImproperPVTest : public ::testing::Test {
protected:
    void SetUp() override { ctx = std::make_unique<symbolic::CASContext>(); }

    [[nodiscard]] ExprPtr E(const std::string& src) {
        auto r = parse_expr(src, ctx->arena());
        EXPECT_TRUE(r.is_ok()) << "Parse failed: " << src
                               << "  err=" << (r.is_error() ? r.error().message : std::string{});
        return r.is_ok() ? r.value() : ExprPtr{};
    }

    [[nodiscard]] ExprPtr pos_inf() {
        return ctx->arena().make<Constant>(MathConstant::Infinity);
    }
    [[nodiscard]] ExprPtr neg_inf() {
        return ctx->arena().make<Unary>(UnaryOp::Neg, pos_inf());
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

TEST_F(ImproperPVTest, ConvergesAtBothEnds) {
    // 1/(1+x^2) from -oo to +oo.  Order at +-oo is -2.  Convergent.
    auto report = calculus::classify_improper_convergence(
        E("1/(1 + x^2)"), Symbol("x"), neg_inf(), pos_inf(), *ctx);
    ASSERT_TRUE(report.is_ok()) << report.error().message;
    EXPECT_EQ(report.value().status, calculus::ConvergenceStatus::Convergent)
        << report.value().diagnostic;
    EXPECT_EQ(report.value().leading_order_at_lower, -2);
    EXPECT_EQ(report.value().leading_order_at_upper, -2);
}

TEST_F(ImproperPVTest, DivergesAtZero) {
    // 1/x^2 from 0 to 1.  Order -2 at x=0 -> divergent (algebraic).
    auto report = calculus::classify_improper_convergence(
        E("1/x^2"), Symbol("x"), E("0"), E("1"), *ctx);
    ASSERT_TRUE(report.is_ok()) << report.error().message;
    EXPECT_EQ(report.value().status, calculus::ConvergenceStatus::DivergentAtLowerEnd)
        << report.value().diagnostic;
    EXPECT_EQ(report.value().leading_order_at_lower, -2);
}

TEST_F(ImproperPVTest, LogarithmicDivergenceAtZero) {
    // 1/x from 0 to 1.  Order -1 at x=0 -> log divergence.
    auto report = calculus::classify_improper_convergence(
        E("1/x"), Symbol("x"), E("0"), E("1"), *ctx);
    ASSERT_TRUE(report.is_ok()) << report.error().message;
    EXPECT_EQ(report.value().status, calculus::ConvergenceStatus::DivergentAtLowerEnd)
        << report.value().diagnostic;
    EXPECT_EQ(report.value().leading_order_at_lower, -1);
}

TEST_F(ImproperPVTest, CauchyPV_SimplePole_OneOverX) {
    // PV ∫_{-1}^{1} 1/x dx = 0.
    auto pv = calculus::cauchy_principal_value(
        E("1/x"), Symbol("x"), E("-1"), E("1"), E("0"), *ctx);
    ASSERT_TRUE(pv.is_ok()) << pv.error().message;
    auto zero = ctx->simplify(E("0"));
    ASSERT_TRUE(zero.is_ok());
    auto eq = symbolic::mathematically_equal(pv.value(), zero.value(), *ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value())
        << "PV(1/x, -1, 1) expected 0, got " << debug_print(pv.value());
}

TEST_F(ImproperPVTest, CauchyPV_AlgebraicAtPole) {
    // PV ∫_{0}^{2} 1/(x-1) dx = ln|1| - ln|1| = 0.
    auto pv = calculus::cauchy_principal_value(
        E("1/(x-1)"), Symbol("x"), E("0"), E("2"), E("1"), *ctx);
    ASSERT_TRUE(pv.is_ok()) << pv.error().message;
    auto zero = ctx->simplify(E("0"));
    ASSERT_TRUE(zero.is_ok());
    auto eq = symbolic::mathematically_equal(pv.value(), zero.value(), *ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value())
        << "PV(1/(x-1), 0, 2) expected 0, got " << debug_print(pv.value());
}

}  // namespace cas::test
