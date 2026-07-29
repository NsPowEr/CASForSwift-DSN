// DEBT-002 smoke coverage for src/calculus/limit.cpp.
// Happy-path checks for finite and infinite limits using LimitEngine.

#include <gtest/gtest.h>

#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class LimitSmokeTest : public ::testing::Test {
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

    [[nodiscard]] bool limit_equals(const std::string& expr_str, ExprPtr point,
                                     LimitDirection dir,
                                     const std::string& expected_str) {
        auto e = parse(expr_str);
        auto exp = parse(expected_str);
        auto r = calculus::limit(e, x, point, dir, ctx);
        EXPECT_TRUE(r.is_ok()) << expr_str;
        if (!r.is_ok()) return false;
        return structural_equal(ctx.simplify(r.value()).value(), exp);
    }
};

TEST_F(LimitSmokeTest, FiniteLimitDirectSubstitution) {
    auto pt = parse("0");
    EXPECT_TRUE(limit_equals("x^2 + 3*x + 1", pt, LimitDirection::Both, "1"));
}

TEST_F(LimitSmokeTest, LimitSinXOverXAtZero) {
    auto pt = parse("0");
    auto e = parse("sin(x) / x");
    auto r = calculus::limit(e, x, pt, LimitDirection::Both, ctx);
    ASSERT_TRUE(r.is_ok());
    auto* lit = expr_cast<IntegerLit>(ctx.simplify(r.value()).value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(1));
}

TEST_F(LimitSmokeTest, InfiniteLimitPolynomialRatio) {
    auto pt = ctx.arena().make<Constant>(MathConstant::Infinity);
    auto e = parse("(2*x^2 + x) / (x^2 + 5)");
    auto r = calculus::limit(e, x, pt, LimitDirection::Right, ctx);
    ASSERT_TRUE(r.is_ok());
    auto simp = ctx.simplify(r.value());
    ASSERT_TRUE(simp.is_ok());
    if (auto* il = expr_cast<IntegerLit>(simp.value())) {
        EXPECT_EQ(il->value, BigInt(2));
    } else if (auto* rl = expr_cast<RationalLit>(simp.value())) {
        EXPECT_EQ(rl->numerator, BigInt(2));
        EXPECT_EQ(rl->denominator, BigInt(1));
    } else {
        FAIL() << "limit returned non-literal kind=" << (int)simp.value()->kind;
    }
}

TEST_F(LimitSmokeTest, LHopitalIndeterminateForm) {
    auto pt = parse("0");
    auto e = parse("(exp(x) - 1) / x");
    auto r = calculus::limit(e, x, pt, LimitDirection::Both, ctx);
    ASSERT_TRUE(r.is_ok());
    auto* lit = expr_cast<IntegerLit>(ctx.simplify(r.value()).value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(1));
}

TEST_F(LimitSmokeTest, RemovableSingularityCancellation) {
    auto pt = parse("1");
    auto e = parse("(x^2 - 1) / (x - 1)");
    auto r = calculus::limit(e, x, pt, LimitDirection::Both, ctx);
    ASSERT_TRUE(r.is_ok());
    auto* lit = expr_cast<IntegerLit>(ctx.simplify(r.value()).value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(2));
}

// Helper: classify a limit result as -inf regardless of canonical form
// (Unary(Neg, Infinity) or Constant(NegInfinity)).
static bool is_negative_infinity(ExprPtr v) {
    if (const auto* c = expr_cast<Constant>(v))
        return c->value == MathConstant::NegInfinity;
    if (const auto* u = expr_cast<Unary>(v); u && u->op == UnaryOp::Neg)
        if (const auto* c = expr_cast<Constant>(u->operand))
            return c->value == MathConstant::Infinity;
    return false;
}

// A signed-pole regression: tan(x)/(x-pi/2) ~ -1/(x-pi/2)^2 near pi/2, so the
// two-sided limit is -inf (the pole has even order, same sign on both sides).
// The engine previously returned an unsigned +inf for any nonzero/0 pole; the
// general reciprocal-based sign analysis must recover the -inf.
TEST_F(LimitSmokeTest, TanPoleSignedInfinityAtHalfPi) {
    auto pt = parse("pi/2");
    auto e = parse("tan(x) / (x - pi/2)");
    auto r = calculus::limit(e, x, pt, LimitDirection::Both, ctx);
    ASSERT_TRUE(r.is_ok()) << (r.is_ok() ? "" : r.error().message);
    auto s = ctx.simplify(r.value());
    ASSERT_TRUE(s.is_ok());
    EXPECT_TRUE(is_negative_infinity(s.value()))
        << "expected -inf for tan(x)/(x-pi/2) at pi/2";
}

}  // namespace
