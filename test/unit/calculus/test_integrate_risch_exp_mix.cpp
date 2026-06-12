// L1-02 STEP 5 — Risch DE dispatch for ∫f(x)·exp(g(x))dx pattern.
// Verifies Bronstein cap. 6 Risch DE solver wired as top-level
// dispatch in integrate_risch.cpp BEFORE Hermite/Rothstein-Trager.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class IntegrateRischExpMixTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};
    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
    [[nodiscard]] bool verify_antider(ExprPtr F, ExprPtr expr) {
        auto D = calculus::diff(F, x, 1U, ctx);
        if (D.is_error()) return false;
        ExprPtr delta = ctx.arena().make<Binary>(BinaryOp::Sub, D.value(), expr);
        auto t = algebra::together(delta, ctx);
        if (t.is_error()) return false;
        auto simp = ctx.simplify(t.value());
        if (simp.is_error()) return false;
        auto* lit = expr_cast<IntegerLit>(simp.value());
        return lit != nullptr && lit->value.is_zero();
    }
};

TEST_F(IntegrateRischExpMixTest, IntegralOfXExp) {
    // ∫ x·exp(x) dx = (x-1)·exp(x)
    auto e = parse("x * exp(x)");
    auto r = calculus::integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok()) << (r.is_error() ? r.error().message : "");
    EXPECT_TRUE(verify_antider(r.value(), e));
}

TEST_F(IntegrateRischExpMixTest, IntegralOfXSquaredExp) {
    // ∫ x²·exp(x) dx = (x²-2x+2)·exp(x)
    auto e = parse("x^2 * exp(x)");
    auto r = calculus::integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(verify_antider(r.value(), e));
}

TEST_F(IntegrateRischExpMixTest, IntegralOfPolyTimesExpQuadratic) {
    // ∫ (2x+1)·exp(x²+x) dx = exp(x²+x)  (since D(x²+x) = 2x+1)
    auto e = parse("(2*x + 1) * exp(x^2 + x)");
    auto r = calculus::integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(verify_antider(r.value(), e));
}

TEST_F(IntegrateRischExpMixTest, IntegralOfXExpXSquared) {
    // ∫ x·exp(x²) dx = exp(x²)/2  (since D(x²) = 2x, y = 1/2)
    auto e = parse("x * exp(x^2)");
    auto r = calculus::integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(verify_antider(r.value(), e));
}

TEST_F(IntegrateRischExpMixTest, AntiHardcodeHigherDegreeExpQuadratic) {
    // ∫ (3x³+x)·exp(x²) dx, y satisfies y' + 2x·y = 3x³+x
    auto e = parse("(3*x^3 + x) * exp(x^2)");
    auto r = calculus::integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(verify_antider(r.value(), e));
}

TEST_F(IntegrateRischExpMixTest, IntegralOfErfTimesExp) {
    // ∫ x·erf(x)·exp(x²) dx = 1/2 · erf(x)·exp(x²) − x/sqrt(pi)
    auto e = parse("x * erf(x) * exp(x^2)");
    auto r = calculus::integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok()) << (r.is_error() ? r.error().message : "");
    EXPECT_TRUE(verify_antider(r.value(), e));
}

}  // namespace
