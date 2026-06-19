// T-016 / HC-F75-B1-IBP-DOUBLE-APPLY — regression guard for integration by parts.
//
// The bug: integrate(x*log(x)) re-applied by-parts on an already-reduced
// sub-integral, emitting 4 redundant non-cancelling terms instead of the
// correct ½x²·log(x) - ¼x². The fix (simplify the v·du sub-integrand before
// recursing) was verified only against the golden corpus — no unit test pinned
// it. This adds that guard.
//
// Verification is by MATHEMATICAL EQUIVALENCE (CLAUDE.md: never validate via
// toString): for an antiderivative I of f, d/dx(I) - f must simplify to 0.

#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas::calculus {

class IntegrateByPartsTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};

    ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }

    // Assert that ∫f is a valid antiderivative: d/dx(∫f) - f ≡ 0.
    void expect_integral_roundtrips(const std::string& integrand) {
        ExprPtr f = parse(integrand);
        auto integral = integrate(f, x, ctx);
        ASSERT_TRUE(integral.is_ok())
            << "integrate(" << integrand << ") failed: " << integral.error().message;

        auto back = diff(integral.value(), x, 1U, ctx);
        ASSERT_TRUE(back.is_ok())
            << "diff failed: " << back.error().message;

        ExprPtr residual = ctx.arena().make<Binary>(BinaryOp::Sub, back.value(), f);
        auto simplified = ctx.simplify(residual);
        ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

        const auto* lit = expr_cast<IntegerLit>(simplified.value());
        ASSERT_NE(lit, nullptr)
            << "d/dx(∫ " << integrand << ") - (" << integrand
            << ") did not simplify to a literal 0 — IBP likely produced a wrong antiderivative";
        EXPECT_EQ(lit->value, BigInt(0))
            << "round-trip residual for " << integrand << " is non-zero";
    }
};

TEST_F(IntegrateByPartsTest, XLogX_RoundTrips) {
    expect_integral_roundtrips("x*log(x)");
}

TEST_F(IntegrateByPartsTest, XLogXSquared_RoundTrips) {
    expect_integral_roundtrips("x*log(x)^2");
}

TEST_F(IntegrateByPartsTest, LogXCubed_RoundTrips) {
    expect_integral_roundtrips("log(x)^3");
}

}  // namespace cas::calculus
