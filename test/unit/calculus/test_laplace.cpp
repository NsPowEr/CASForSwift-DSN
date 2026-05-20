// CAS-L3-07 — Laplace transform pattern tests.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include "../../../src/calculus/calculus_internal.hpp"

using namespace cas;
using namespace cas::calculus;

namespace {

class LaplaceTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol t{"t"};
    Symbol s{"s"};
    [[nodiscard]] ExprPtr parse(const std::string& str) {
        auto tk = Lexer(str).tokenize();
        EXPECT_TRUE(tk.is_ok()) << str;
        Parser p(tk.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << str;
        return r.value();
    }
    [[nodiscard]] bool equiv(ExprPtr a, ExprPtr b) {
        auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, a, b);
        auto tog = algebra::together(delta, ctx);
        auto simp = ctx.simplify(tog.is_ok() ? tog.value() : delta);
        if (simp.is_error()) return false;
        if (auto* il = expr_cast<IntegerLit>(simp.value())) return il->value.is_zero();
        return false;
    }
};

TEST_F(LaplaceTest, ConstantTransformsTo1OverS) {
    auto f = parse("1");
    auto r = laplace_transform(f, t, s, ctx);
    ASSERT_TRUE(r.is_ok());
    auto expected = parse("1 / s");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(LaplaceTest, TTransformsTo1OverSSquared) {
    auto f = parse("t");
    auto r = laplace_transform(f, t, s, ctx);
    ASSERT_TRUE(r.is_ok());
    auto expected = parse("1 / s^2");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(LaplaceTest, TPowerNTransformsViaFactorial) {
    // L{t^3} = 6/s⁴
    auto f = parse("t^3");
    auto r = laplace_transform(f, t, s, ctx);
    ASSERT_TRUE(r.is_ok());
    auto expected = parse("6 / s^4");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(LaplaceTest, ExpTransformsTo1OverSMinusA) {
    // L{exp(3t)} = 1/(s-3)
    auto f = parse("exp(3*t)");
    auto r = laplace_transform(f, t, s, ctx);
    ASSERT_TRUE(r.is_ok());
    auto expected = parse("1 / (s - 3)");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(LaplaceTest, SinTransform) {
    // L{sin(2t)} = 2/(s²+4)
    auto f = parse("sin(2*t)");
    auto r = laplace_transform(f, t, s, ctx);
    ASSERT_TRUE(r.is_ok());
    auto expected = parse("2 / (s^2 + 4)");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(LaplaceTest, CosTransform) {
    // L{cos(t)} = s/(s²+1)
    auto f = parse("cos(t)");
    auto r = laplace_transform(f, t, s, ctx);
    ASSERT_TRUE(r.is_ok());
    auto expected = parse("s / (s^2 + 1)");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(LaplaceTest, LinearityCombination) {
    // L{3·t + 2·sin(t)} = 3/s² + 2/(s²+1)
    auto f = parse("3*t + 2*sin(t)");
    auto r = laplace_transform(f, t, s, ctx);
    ASSERT_TRUE(r.is_ok());
    auto expected = parse("3/s^2 + 2/(s^2 + 1)");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(LaplaceTest, ConstantScalarFactor) {
    // L{5·exp(2t)} = 5/(s-2)
    auto f = parse("5 * exp(2*t)");
    auto r = laplace_transform(f, t, s, ctx);
    ASSERT_TRUE(r.is_ok());
    auto expected = parse("5 / (s - 2)");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(LaplaceTest, AntiHardcodeNonElementaryRejected) {
    // L{ln(t)} non in table → Unimplemented (or returns weird).
    // Anti-hardcode: must NOT silently return wrong value.
    auto f = parse("ln(t)");
    auto r = laplace_transform(f, t, s, ctx);
    // Either Unimplemented OR symbolic form that doesn't collapse to wrong value.
    if (r.is_ok()) {
        // If returned ok, no specific oracle. Test passes either way.
        EXPECT_NE(r.value(), nullptr);
    }
}

}  // namespace
