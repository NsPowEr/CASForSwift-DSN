// F7.5 follow-up — exact-value identities for sinh/cosh/tanh/coth and
// boundary values of the inverse trig builtins. Closes the limit-area
// regressions where the dispatcher returned an unevaluated FuncCall
// after the canonical-value substitution (e.g. `cosh(0)`, `asin(0)`).

#include <gtest/gtest.h>

#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class HyperbolicZeroTest : public ::testing::Test {
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

    [[nodiscard]] ExprPtr simp(const std::string& s) {
        auto r = ctx.simplify(parse(s));
        EXPECT_TRUE(r.is_ok()) << s;
        return r.is_ok() ? r.value() : nullptr;
    }

    [[nodiscard]] bool is_int(ExprPtr e, long long n) {
        if (!e) return false;
        const auto* lit = expr_cast<IntegerLit>(e);
        return lit != nullptr && lit->value == BigInt(n);
    }
};

TEST_F(HyperbolicZeroTest, SinhZeroEqualsZero) {
    EXPECT_TRUE(is_int(simp("sinh(0)"), 0));
}

TEST_F(HyperbolicZeroTest, CoshZeroEqualsOne) {
    EXPECT_TRUE(is_int(simp("cosh(0)"), 1));
}

TEST_F(HyperbolicZeroTest, TanhZeroEqualsZero) {
    EXPECT_TRUE(is_int(simp("tanh(0)"), 0));
}

TEST_F(HyperbolicZeroTest, CothZeroLeftUnchanged) {
    // coth(0) is a true singularity and must NOT be silently rewritten
    // to a finite value. The simplifier returns the unevaluated FuncCall;
    // a future Extended-Real pass can map it to ComplexInfinity.
    auto r = simp("coth(0)");
    ASSERT_NE(r, nullptr);
    const auto* fc = expr_cast<FuncCall>(r);
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->func_id, BuiltinOp::Coth);
}

TEST_F(HyperbolicZeroTest, CoshIsEvenViaNegation) {
    auto lhs = simp("cosh(-x)");
    auto rhs = simp("cosh(x)");
    auto eq = cas::symbolic::mathematically_equal(lhs, rhs, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

TEST_F(HyperbolicZeroTest, SinhIsOddViaNegation) {
    // sinh(-x) = -sinh(x): rewriting must produce the canonical -sinh(x)
    // form, structurally distinct from sinh(-x).
    auto lhs = simp("sinh(-x)");
    auto neg_rhs = simp("-sinh(x)");
    auto eq = cas::symbolic::mathematically_equal(lhs, neg_rhs, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

TEST_F(HyperbolicZeroTest, AsinZeroAndOneBoundaries) {
    EXPECT_TRUE(is_int(simp("asin(0)"), 0));
    // asin(1) = π/2; structurally Pi/2 after simplify.
    auto r = simp("asin(1)");
    auto pi_2 = ctx.arena().make<Binary>(BinaryOp::Div,
        ctx.arena().make<Constant>(MathConstant::Pi),
        ctx.arena().make<IntegerLit>(BigInt(2)));
    auto eq = cas::symbolic::mathematically_equal(r, pi_2, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

TEST_F(HyperbolicZeroTest, AcosBoundariesZeroOneMinusOne) {
    auto pi = ctx.arena().make<Constant>(MathConstant::Pi);
    auto pi_2 = ctx.arena().make<Binary>(BinaryOp::Div, pi,
        ctx.arena().make<IntegerLit>(BigInt(2)));
    for (auto [input, expected] : std::vector<std::pair<std::string, ExprPtr>>{
            {"acos(0)", pi_2},
            {"acos(1)", ctx.arena().make<IntegerLit>(BigInt(0))},
            {"acos(-1)", pi},
        }) {
        auto eq = cas::symbolic::mathematically_equal(simp(input), expected, ctx);
        ASSERT_TRUE(eq.is_ok()) << input;
        EXPECT_TRUE(eq.value()) << input;
    }
}

TEST_F(HyperbolicZeroTest, AtanAtInfinityIsHalfPi) {
    auto inf = ctx.arena().make<Constant>(MathConstant::Infinity);
    auto atan_inf = ctx.arena().make<FuncCall>(
        BuiltinOp::Atan, std::vector<ExprPtr>{inf});
    auto r = ctx.simplify(atan_inf);
    ASSERT_TRUE(r.is_ok());
    auto pi_2 = ctx.arena().make<Binary>(BinaryOp::Div,
        ctx.arena().make<Constant>(MathConstant::Pi),
        ctx.arena().make<IntegerLit>(BigInt(2)));
    auto eq = cas::symbolic::mathematically_equal(r.value(), pi_2, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

TEST_F(HyperbolicZeroTest, AtanAtNegativeInfinityIsMinusHalfPi) {
    auto neg_inf = ctx.arena().make<Unary>(UnaryOp::Neg,
        ctx.arena().make<Constant>(MathConstant::Infinity));
    auto atan_neg_inf = ctx.arena().make<FuncCall>(
        BuiltinOp::Atan, std::vector<ExprPtr>{neg_inf});
    auto r = ctx.simplify(atan_neg_inf);
    ASSERT_TRUE(r.is_ok());
    auto pi_2 = ctx.arena().make<Binary>(BinaryOp::Div,
        ctx.arena().make<Constant>(MathConstant::Pi),
        ctx.arena().make<IntegerLit>(BigInt(2)));
    auto expected = ctx.arena().make<Unary>(UnaryOp::Neg, pi_2);
    auto eq = cas::symbolic::mathematically_equal(r.value(), expected, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

}  // namespace
