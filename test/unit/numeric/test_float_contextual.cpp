// CAS-L3-03 — Float contestuale: precisione default + overload ctx-aware.
//
// Verifica che `eval_mpfr(expr, ctx)` usa `ctx.numeric_precision_digits()`
// e che `set_numeric_precision_digits` applica clamp [6, 10000].

#include <gtest/gtest.h>

#include "cas/lexer.hpp"
#include "cas/numeric.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::numeric;

namespace {

class FloatContextualTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
};

TEST_F(FloatContextualTest, DefaultPrecisionIs15Digits) {
    EXPECT_EQ(ctx.numeric_precision_digits(), 15U);
}

TEST_F(FloatContextualTest, SetPrecisionRoundtrip) {
    ctx.set_numeric_precision_digits(100);
    EXPECT_EQ(ctx.numeric_precision_digits(), 100U);
    ctx.set_numeric_precision_digits(50);
    EXPECT_EQ(ctx.numeric_precision_digits(), 50U);
}

TEST_F(FloatContextualTest, PrecisionClampLow) {
    ctx.set_numeric_precision_digits(1);
    EXPECT_EQ(ctx.numeric_precision_digits(), 6U);  // clamped to min 6
}

TEST_F(FloatContextualTest, PrecisionClampHigh) {
    ctx.set_numeric_precision_digits(20000);
    EXPECT_EQ(ctx.numeric_precision_digits(), 10000U);  // clamped to max
}

TEST_F(FloatContextualTest, ContextAwareEvalUsesContextPrecision) {
    auto e = parse("pi");
    ctx.set_numeric_precision_digits(50);
    auto r = eval_mpfr(e, ctx);
    ASSERT_TRUE(r.is_ok());
    // First 15 digits should match pi.
    EXPECT_TRUE(r.value().starts_with("3.14159265358979"))
        << "got: " << r.value();
    // Roughly 50 digits of significance — value length depends on
    // formatter; verify it's substantially longer than default 15.
    EXPECT_GE(r.value().size(), 45U);
}

TEST_F(FloatContextualTest, ContextAwareSwitchesPrecision) {
    auto e = parse("sqrt(2)");
    ctx.set_numeric_precision_digits(10);
    auto r10 = eval_mpfr(e, ctx);
    ctx.set_numeric_precision_digits(40);
    auto r40 = eval_mpfr(e, ctx);
    ASSERT_TRUE(r10.is_ok() && r40.is_ok());
    EXPECT_GT(r40.value().size(), r10.value().size())
        << "Higher precision should yield longer string";
    // Both must start with sqrt(2) prefix.
    EXPECT_TRUE(r10.value().starts_with("1.41421"));
    EXPECT_TRUE(r40.value().starts_with("1.41421356237"));
}

}  // namespace
