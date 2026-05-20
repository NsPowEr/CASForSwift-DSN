// DEBT-002 smoke coverage for src/calculus/integrate_risch.cpp.
// Targets Risch entry covering: rational integration, log derivative,
// exp integration, polynomial-in-t with exp ext, and roundtrip.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class IntegrateRischSmokeTest : public ::testing::Test {
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

    [[nodiscard]] bool verify_antiderivative(ExprPtr F, ExprPtr expr) {
        auto D = calculus::diff(F, x, 1U, ctx);
        EXPECT_TRUE(D.is_ok());
        if (!D.is_ok()) return false;
        auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, D.value(), expr);
        auto t = algebra::together(delta, ctx);
        if (t.is_error()) return false;
        auto s = ctx.simplify(t.value());
        if (s.is_error()) return false;
        auto* lit = expr_cast<IntegerLit>(s.value());
        return lit != nullptr && lit->value.is_zero();
    }
};

TEST_F(IntegrateRischSmokeTest, IntegralOfReciprocalX) {
    auto e = parse("1 / x");
    auto r = calculus::integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(verify_antiderivative(r.value(), e));
}

TEST_F(IntegrateRischSmokeTest, IntegralOfExp) {
    auto e = parse("exp(x)");
    auto r = calculus::integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(verify_antiderivative(r.value(), e));
}

TEST_F(IntegrateRischSmokeTest, IntegralOfXOverXSquaredPlusOne) {
    auto e = parse("x / (x^2 + 1)");
    auto r = calculus::integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(verify_antiderivative(r.value(), e));
}

TEST_F(IntegrateRischSmokeTest, IntegralOfXTimesExp) {
    auto e = parse("x * exp(x)");
    auto r = calculus::integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(verify_antiderivative(r.value(), e));
}

TEST_F(IntegrateRischSmokeTest, IntegralOfReciprocalXLnX) {
    // DEBT-004 regression sentinel — ln(ln(x)) via logarithmic-derivative.
    auto e = parse("1 / (x * ln(x))");
    auto r = calculus::integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(verify_antiderivative(r.value(), e));
}

}  // namespace
