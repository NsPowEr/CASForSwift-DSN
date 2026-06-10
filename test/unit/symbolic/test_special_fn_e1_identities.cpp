// F7.5.E1 — Γ/ζ/erf canonical identity coverage.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Special_Fn_Identities.md

#include <gtest/gtest.h>

#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/extended_real.hpp"
#include "cas/symbolic.hpp"

namespace cas {
namespace {

class SpecialFnE1Test : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] ExprPtr call1(BuiltinOp op, ExprPtr arg) {
        return ctx.arena().make<FuncCall>(op, std::vector<ExprPtr>{arg});
    }
    [[nodiscard]] ExprPtr int_lit(long long v) {
        return ctx.arena().make<IntegerLit>(BigInt(v));
    }
    [[nodiscard]] ExprPtr neg(ExprPtr e) {
        return ctx.arena().make<Unary>(UnaryOp::Neg, e);
    }
    [[nodiscard]] ExprPtr simplify(ExprPtr e) {
        auto r = ctx.simplify(e);
        EXPECT_TRUE(r.is_ok());
        return r.is_ok() ? r.value() : nullptr;
    }
    [[nodiscard]] bool is_integer_value(ExprPtr e, long long v) {
        if (!e) return false;
        if (const auto* lit = expr_cast<IntegerLit>(e)) return lit->value == BigInt(v);
        return false;
    }
};

// Γ(0) → ComplexInfinity (pole)
TEST_F(SpecialFnE1Test, GammaAtZeroIsComplexInfinity) {
    auto r = simplify(call1(BuiltinOp::Gamma, int_lit(0)));
    EXPECT_TRUE(is_complex_infinity(r)) << "got: " << debug_print(r);
}

// Γ(-1) → ComplexInfinity (pole)
TEST_F(SpecialFnE1Test, GammaAtNegOneIsComplexInfinity) {
    auto r = simplify(call1(BuiltinOp::Gamma, int_lit(-1)));
    EXPECT_TRUE(is_complex_infinity(r)) << "got: " << debug_print(r);
}

// Γ(-5) → ComplexInfinity (pole)
TEST_F(SpecialFnE1Test, GammaAtNegFiveIsComplexInfinity) {
    auto r = simplify(call1(BuiltinOp::Gamma, int_lit(-5)));
    EXPECT_TRUE(is_complex_infinity(r)) << "got: " << debug_print(r);
}

// Γ(+∞) → +∞ (Stirling growth)
TEST_F(SpecialFnE1Test, GammaAtPosInfIsPosInf) {
    auto inf = ctx.arena().make<Constant>(MathConstant::Infinity);
    auto r = simplify(call1(BuiltinOp::Gamma, inf));
    EXPECT_TRUE(is_pos_infinity(r)) << "got: " << debug_print(r);
}

// Γ(5) = 24 (regression — should still work)
TEST_F(SpecialFnE1Test, GammaAtFiveIsTwentyFour) {
    auto r = simplify(call1(BuiltinOp::Gamma, int_lit(5)));
    EXPECT_TRUE(is_integer_value(r, 24)) << "got: " << debug_print(r);
}

// erf(+∞) = 1
TEST_F(SpecialFnE1Test, ErfAtPosInfIsOne) {
    auto inf = ctx.arena().make<Constant>(MathConstant::Infinity);
    auto r = simplify(call1(BuiltinOp::Erf, inf));
    EXPECT_TRUE(is_integer_value(r, 1)) << "got: " << debug_print(r);
}

// erf(-∞) = -1
TEST_F(SpecialFnE1Test, ErfAtNegInfIsMinusOne) {
    auto neg_inf = ctx.arena().make<Constant>(MathConstant::NegInfinity);
    auto r = simplify(call1(BuiltinOp::Erf, neg_inf));
    // Expect -1: either canonical IntegerLit(-1) or Unary(Neg, 1).
    bool match = is_integer_value(r, -1);
    if (!match && r) {
        if (const auto* un = expr_cast<Unary>(r)) {
            match = un->op == UnaryOp::Neg && is_integer_value(un->operand, 1);
        }
    }
    EXPECT_TRUE(match) << "got: " << debug_print(r);
}

// erf(0) = 0 (regression)
TEST_F(SpecialFnE1Test, ErfAtZeroIsZero) {
    auto r = simplify(call1(BuiltinOp::Erf, int_lit(0)));
    EXPECT_TRUE(is_integer_value(r, 0)) << "got: " << debug_print(r);
}

// ζ(0) = -1/2 (regression)
TEST_F(SpecialFnE1Test, ZetaAtZeroIsMinusHalf) {
    auto r = simplify(call1(BuiltinOp::Zeta, int_lit(0)));
    ASSERT_NE(r, nullptr);
    const auto* rat = expr_cast<RationalLit>(r);
    ASSERT_NE(rat, nullptr) << "got: " << debug_print(r);
    EXPECT_EQ(rat->numerator, BigInt(-1));
    EXPECT_EQ(rat->denominator, BigInt(2));
}

// ζ(1) → ComplexInfinity (Riemann zeta pole)
TEST_F(SpecialFnE1Test, ZetaAtOneIsComplexInfinity) {
    auto r = simplify(call1(BuiltinOp::Zeta, int_lit(1)));
    EXPECT_TRUE(is_complex_infinity(r)) << "got: " << debug_print(r);
}

// ζ(-1) = -1/12 (regression — Bernoulli)
TEST_F(SpecialFnE1Test, ZetaAtMinusOneIsMinusOneTwelfth) {
    auto r = simplify(call1(BuiltinOp::Zeta, int_lit(-1)));
    ASSERT_NE(r, nullptr);
    const auto* rat = expr_cast<RationalLit>(r);
    ASSERT_NE(rat, nullptr) << "got: " << debug_print(r);
    EXPECT_EQ(rat->numerator, BigInt(-1));
    EXPECT_EQ(rat->denominator, BigInt(12));
}

// ζ(-2) = 0 (trivial zero)
TEST_F(SpecialFnE1Test, ZetaAtMinusTwoIsZero) {
    auto r = simplify(call1(BuiltinOp::Zeta, int_lit(-2)));
    EXPECT_TRUE(is_integer_value(r, 0)) << "got: " << debug_print(r);
}

}  // namespace
}  // namespace cas
