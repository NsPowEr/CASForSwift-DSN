// F7.5.D1 — Gruntz §3.5 nested log tower tests.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Gruntz_Nested_Log.md

#include <gtest/gtest.h>

#include "cas/calculus.hpp"
#include "cas/extended_real.hpp"
#include "cas/ast_debug.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class GruntzNestedLogTest : public ::testing::Test {
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

    [[nodiscard]] ExprPtr limit_at_pos_inf(const std::string& expr_str) {
        ExprPtr pt = ctx.arena().make<Constant>(MathConstant::Infinity);
        auto e = parse(expr_str);
        auto r = calculus::limit(e, x, pt, LimitDirection::Right, ctx);
        if (!r.is_ok()) {
            ADD_FAILURE() << expr_str << " — error kind=" << static_cast<int>(r.error().kind)
                          << " msg=" << r.error().message;
            return nullptr;
        }
        return ctx.simplify(r.value()).value();
    }

    [[nodiscard]] bool is_zero(ExprPtr e) {
        if (!e) return false;
        if (const auto* lit = expr_cast<IntegerLit>(e)) return lit->value == BigInt(0);
        return false;
    }

    [[nodiscard]] bool is_one(ExprPtr e) {
        if (!e) return false;
        if (const auto* lit = expr_cast<IntegerLit>(e)) return lit->value == BigInt(1);
        return false;
    }
};

// A11 / F5.2 — fractional-power growth vs logarithm. A positive rational power of
// x dominates every logarithm, so log(x)/sqrt(x) → 0 (not ∞). Regression for the
// integer-only growth predicates + Log-vs-Ln logarithm recognition.
TEST_F(GruntzNestedLogTest, LogOverSqrtX) {
    EXPECT_TRUE(is_zero(limit_at_pos_inf("log(x) / sqrt(x)")));
    EXPECT_TRUE(is_zero(limit_at_pos_inf("log(x) / x^(1/2)")));
    EXPECT_TRUE(is_zero(limit_at_pos_inf("sqrt(x) * log(x) / x")));
    // Same with the natural-log spelling.
    EXPECT_TRUE(is_zero(limit_at_pos_inf("ln(x) / sqrt(x)")));
}

TEST_F(GruntzNestedLogTest, FractionalPowerBeatsLog) {
    // sqrt(x)·log(x) → +∞ (positive fractional power dominates the logarithm).
    EXPECT_TRUE(is_pos_infinity(limit_at_pos_inf("sqrt(x) * log(x)")));
    // NOTE: the harder x^(3/2)/(x·log(x)) form still needs the full Gruntz
    // slowly-varying-coefficient track (a ComplexRational coefficient hits a
    // division-by-zero in the leading-power path) — tracked under A11 as the
    // remaining nested-coefficient work, not covered by this fix.
}

// Case 1: lim x→∞ log(log(x)) / x = 0
TEST_F(GruntzNestedLogTest, LogLogOverX) {
    auto r = limit_at_pos_inf("log(log(x)) / x");
    EXPECT_TRUE(is_zero(r));
}

// Direct simplify check: log(Infinity) should reduce to Infinity.
TEST_F(GruntzNestedLogTest, DirectSimplifyLogInfinity) {
    auto inf = ctx.arena().make<Constant>(MathConstant::Infinity);
    auto log_inf = ctx.arena().make<FuncCall>(
        BuiltinOp::Ln, std::vector<ExprPtr>{inf});
    auto s = ctx.simplify(log_inf);
    ASSERT_TRUE(s.is_ok());
    EXPECT_TRUE(is_pos_infinity(s.value())) << "got: " << debug_print(s.value());
}

TEST_F(GruntzNestedLogTest, DirectSimplifyLogOpInfinityViaLog) {
    auto inf = ctx.arena().make<Constant>(MathConstant::Infinity);
    auto log_inf = ctx.arena().make<FuncCall>(
        BuiltinOp::Log, std::vector<ExprPtr>{inf});
    auto s = ctx.simplify(log_inf);
    ASSERT_TRUE(s.is_ok());
    EXPECT_TRUE(is_pos_infinity(s.value())) << "got: " << debug_print(s.value());
}

TEST_F(GruntzNestedLogTest, DirectSimplifyLogLogInfinity) {
    auto inf = ctx.arena().make<Constant>(MathConstant::Infinity);
    auto log_inf = ctx.arena().make<FuncCall>(
        BuiltinOp::Ln, std::vector<ExprPtr>{inf});
    auto log_log_inf = ctx.arena().make<FuncCall>(
        BuiltinOp::Ln, std::vector<ExprPtr>{log_inf});
    auto s = ctx.simplify(log_log_inf);
    ASSERT_TRUE(s.is_ok());
    EXPECT_TRUE(is_pos_infinity(s.value())) << "got: " << debug_print(s.value());
}

// Case 2: lim x→∞ log(log(log(x))) = +∞
TEST_F(GruntzNestedLogTest, LogLogLogXGoesInfinity) {
    auto r = limit_at_pos_inf("log(log(log(x)))");
    EXPECT_TRUE(is_pos_infinity(r)) << "got: " << debug_print(r);
}

// Case 3: lim x→∞ exp(sqrt(log(x))) / x = 0
TEST_F(GruntzNestedLogTest, ExpSqrtLogOverX) {
    auto r = limit_at_pos_inf("exp(sqrt(log(x))) / x");
    EXPECT_TRUE(is_zero(r));
}

// Case 4: lim x→∞ exp(log(x)^2) / exp(x) = 0
TEST_F(GruntzNestedLogTest, ExpLogSquaredOverExpX) {
    auto r = limit_at_pos_inf("exp(log(x)^2) / exp(x)");
    EXPECT_TRUE(is_zero(r));
}

// Case 5: lim x→∞ log(x) * log(log(x)) / log(x)^2 = 0
TEST_F(GruntzNestedLogTest, LogTimesLogLogOverLogSquared) {
    auto r = limit_at_pos_inf("log(x) * log(log(x)) / log(x)^2");
    EXPECT_TRUE(is_zero(r));
}

// Case 6 — F7.5.D2 closure: sum-termwise pre-MRV dispatch.
TEST_F(GruntzNestedLogTest, LogPlusLogLogOverLogX) {
    auto r = limit_at_pos_inf("(log(x) + log(log(x))) / log(x)");
    EXPECT_TRUE(is_one(r));
}

// Case 7 — F7.5.D2 closure: Product · Pow(Product, -1) cancellation pre-MRV.
TEST_F(GruntzNestedLogTest, XLogLogOverXLogX) {
    auto r = limit_at_pos_inf("x * log(log(x)) / (x * log(x))");
    EXPECT_TRUE(is_zero(r));
}

// Case 8: lim x→∞ log(x + log(x)) / log(x) = 1
TEST_F(GruntzNestedLogTest, LogXPlusLogXOverLogX) {
    auto r = limit_at_pos_inf("log(x + log(x)) / log(x)");
    EXPECT_TRUE(is_one(r));
}

}  // namespace
