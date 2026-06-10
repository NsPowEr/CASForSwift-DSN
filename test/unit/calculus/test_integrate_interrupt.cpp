// F7.5.A3 / HC-F75-A3-HARD-TIMEOUT — interrupt cancellation coverage.
//
// The golden runner's per-entry SIGALRM handler calls ctx.interrupt(); if
// the integrator does not poll the interrupt flag it keeps running until
// the simplifier (which does poll) or another downstream check fires.
// For heavy integrands the simplifier may not be entered for many seconds,
// so the entry-point poll-points added in src/calculus/integrate_core.cpp
// (both Integrator::integrate and Integrator::integrate_once) are the
// canonical cancellation surface. These tests pin that contract.

#include <gtest/gtest.h>

#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class IntegrateInterruptTest : public ::testing::Test {
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
};

TEST_F(IntegrateInterruptTest, PreInterruptedIntegrateReturnsTimeout) {
    // Setting the interrupt flag before calling integrate must short-circuit
    // at the first entry-point poll. The exact integrand is irrelevant — the
    // contract is that no algebraic work happens once the flag is set.
    ctx.interrupt();
    auto integrand = parse("x^2 + sin(x) + exp(x)");
    auto r = calculus::integrate(integrand, x, ctx);
    ASSERT_FALSE(r.is_ok());
    EXPECT_EQ(r.error().kind, CASErrorKind::Timeout);
}

TEST_F(IntegrateInterruptTest, ClearInterruptRestoresNormalIntegration) {
    // After interrupt + clear_interrupt the integrator must resume normal
    // operation. Guards against accidentally caching the interrupted state.
    ctx.interrupt();
    ctx.clear_interrupt();
    auto integrand = parse("x");
    auto r = calculus::integrate(integrand, x, ctx);
    ASSERT_TRUE(r.is_ok()) << (r.is_ok() ? "" : r.error().message);
}

TEST_F(IntegrateInterruptTest, InterruptObservedOnEachRecursiveEntry) {
    // The poll lives in Integrator::integrate AND integrate_once, so even if
    // the flag is set between two integrate calls inside the same context
    // the next call must observe it. Validates that the cancellation flag
    // is not consumed/cleared by integrate itself.
    auto a = parse("x^2");
    auto first = calculus::integrate(a, x, ctx);
    ASSERT_TRUE(first.is_ok());
    ctx.interrupt();
    auto b = parse("sin(x)");
    auto second = calculus::integrate(b, x, ctx);
    ASSERT_FALSE(second.is_ok());
    EXPECT_EQ(second.error().kind, CASErrorKind::Timeout);
    ctx.clear_interrupt();
    auto third = calculus::integrate(b, x, ctx);
    EXPECT_TRUE(third.is_ok());
}

}  // namespace
