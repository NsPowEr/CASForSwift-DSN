// F7.5.E2 — Bessel canonical identities at x = 0.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Bessel_Identities.md

#include <gtest/gtest.h>

#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/extended_real.hpp"
#include "cas/symbolic.hpp"

namespace cas {
namespace {

class BesselE2Test : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] ExprPtr bessel(BuiltinOp op, long long order, ExprPtr x_arg) {
        return ctx.arena().make<FuncCall>(op, std::vector<ExprPtr>{
            ctx.arena().make<IntegerLit>(BigInt(order)), x_arg});
    }
    [[nodiscard]] ExprPtr zero() { return ctx.arena().make<IntegerLit>(BigInt(0)); }
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

// J_0(0) = 1
TEST_F(BesselE2Test, BesselJOrderZeroAtZeroIsOne) {
    auto r = simplify(bessel(BuiltinOp::BesselJ, 0, zero()));
    EXPECT_TRUE(is_integer_value(r, 1)) << "got: " << debug_print(r);
}

// J_1(0) = 0
TEST_F(BesselE2Test, BesselJOrderOneAtZeroIsZero) {
    auto r = simplify(bessel(BuiltinOp::BesselJ, 1, zero()));
    EXPECT_TRUE(is_integer_value(r, 0)) << "got: " << debug_print(r);
}

// J_5(0) = 0
TEST_F(BesselE2Test, BesselJOrderFiveAtZeroIsZero) {
    auto r = simplify(bessel(BuiltinOp::BesselJ, 5, zero()));
    EXPECT_TRUE(is_integer_value(r, 0)) << "got: " << debug_print(r);
}

// I_0(0) = 1
TEST_F(BesselE2Test, BesselIOrderZeroAtZeroIsOne) {
    auto r = simplify(bessel(BuiltinOp::BesselI, 0, zero()));
    EXPECT_TRUE(is_integer_value(r, 1)) << "got: " << debug_print(r);
}

// I_2(0) = 0
TEST_F(BesselE2Test, BesselIOrderTwoAtZeroIsZero) {
    auto r = simplify(bessel(BuiltinOp::BesselI, 2, zero()));
    EXPECT_TRUE(is_integer_value(r, 0)) << "got: " << debug_print(r);
}

// Y_0(0) = -∞
TEST_F(BesselE2Test, BesselYOrderZeroAtZeroIsNegInfinity) {
    auto r = simplify(bessel(BuiltinOp::BesselY, 0, zero()));
    EXPECT_TRUE(is_neg_infinity(r)) << "got: " << debug_print(r);
}

// Y_3(0) = -∞
TEST_F(BesselE2Test, BesselYOrderThreeAtZeroIsNegInfinity) {
    auto r = simplify(bessel(BuiltinOp::BesselY, 3, zero()));
    EXPECT_TRUE(is_neg_infinity(r)) << "got: " << debug_print(r);
}

// K_0(0) = +∞
TEST_F(BesselE2Test, BesselKOrderZeroAtZeroIsPosInfinity) {
    auto r = simplify(bessel(BuiltinOp::BesselK, 0, zero()));
    EXPECT_TRUE(is_pos_infinity(r)) << "got: " << debug_print(r);
}

// K_2(0) = +∞
TEST_F(BesselE2Test, BesselKOrderTwoAtZeroIsPosInfinity) {
    auto r = simplify(bessel(BuiltinOp::BesselK, 2, zero()));
    EXPECT_TRUE(is_pos_infinity(r)) << "got: " << debug_print(r);
}

}  // namespace
}  // namespace cas
