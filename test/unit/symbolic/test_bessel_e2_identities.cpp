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

// ─── F7.5.E2 extension: half-integer reductions via recurrence ───────────────

namespace {

[[nodiscard]] ExprPtr make_symbol(symbolic::CASContext& ctx, const char* name) {
    return ctx.arena().make<Symbol>(name);
}

[[nodiscard]] ExprPtr make_half_order(symbolic::CASContext& ctx, long long num) {
    return ctx.arena().make<RationalLit>(BigInt(num), BigInt(2));
}

[[nodiscard]] bool structurally_equal_after_simplify(
    symbolic::CASContext& ctx, ExprPtr a, ExprPtr b) {
    auto ra = ctx.simplify(a);
    auto rb = ctx.simplify(b);
    if (!ra.is_ok() || !rb.is_ok()) return false;
    // Compare via simplifying (a - b) and checking against zero.
    ExprPtr diff = ctx.arena().make<Sum>(std::vector<ExprPtr>{
        ra.value(),
        ctx.arena().make<Unary>(UnaryOp::Neg, rb.value())});
    auto rd = ctx.simplify(diff);
    if (!rd.is_ok()) return false;
    const auto* lit = expr_cast<IntegerLit>(rd.value());
    return lit != nullptr && lit->value.is_zero();
}

}  // namespace

// J_{1/2}(x) = sqrt(2/(π·x))·sin(x): direct base case
TEST_F(BesselE2Test, BesselJ_OneHalf_IsClosedForm) {
    ExprPtr x = make_symbol(ctx, "x");
    ExprPtr call = ctx.arena().make<FuncCall>(BuiltinOp::BesselJ,
        std::vector<ExprPtr>{make_half_order(ctx, 1), x});
    auto r = simplify(call);
    ASSERT_NE(r, nullptr);
    EXPECT_EQ(expr_cast<FuncCall>(r), nullptr)
        << "J_{1/2} must reduce to elementary form, got: " << debug_print(r);
}

// J_{3/2}(x) - ((1/x)·J_{1/2}(x) - J_{-1/2}(x)) = 0   (Knuth/DLMF 10.6.1)
TEST_F(BesselE2Test, BesselJ_ThreeHalves_MatchesRecurrence) {
    ExprPtr x = make_symbol(ctx, "x");
    ExprPtr lhs = ctx.arena().make<FuncCall>(BuiltinOp::BesselJ,
        std::vector<ExprPtr>{make_half_order(ctx, 3), x});
    ExprPtr j_half = ctx.arena().make<FuncCall>(BuiltinOp::BesselJ,
        std::vector<ExprPtr>{make_half_order(ctx, 1), x});
    ExprPtr j_minus_half = ctx.arena().make<FuncCall>(BuiltinOp::BesselJ,
        std::vector<ExprPtr>{make_half_order(ctx, -1), x});
    ExprPtr inv_x = ctx.arena().make<Binary>(BinaryOp::Div,
        ctx.arena().make<IntegerLit>(BigInt(1)), x);
    ExprPtr rhs = ctx.arena().make<Sum>(std::vector<ExprPtr>{
        ctx.arena().make<Product>(std::vector<ExprPtr>{inv_x, j_half}),
        ctx.arena().make<Unary>(UnaryOp::Neg, j_minus_half)});
    EXPECT_TRUE(structurally_equal_after_simplify(ctx, lhs, rhs))
        << "J_{3/2} != (1/x)·J_{1/2} - J_{-1/2}";
}

// I_{3/2}(x) = I_{-1/2}(x) - (1/x)·I_{1/2}(x)  (DLMF 10.29.1)
TEST_F(BesselE2Test, BesselI_ThreeHalves_MatchesRecurrence) {
    ExprPtr x = make_symbol(ctx, "x");
    ExprPtr lhs = ctx.arena().make<FuncCall>(BuiltinOp::BesselI,
        std::vector<ExprPtr>{make_half_order(ctx, 3), x});
    ExprPtr i_half = ctx.arena().make<FuncCall>(BuiltinOp::BesselI,
        std::vector<ExprPtr>{make_half_order(ctx, 1), x});
    ExprPtr i_minus_half = ctx.arena().make<FuncCall>(BuiltinOp::BesselI,
        std::vector<ExprPtr>{make_half_order(ctx, -1), x});
    ExprPtr inv_x = ctx.arena().make<Binary>(BinaryOp::Div,
        ctx.arena().make<IntegerLit>(BigInt(1)), x);
    ExprPtr rhs = ctx.arena().make<Sum>(std::vector<ExprPtr>{
        i_minus_half,
        ctx.arena().make<Unary>(UnaryOp::Neg,
            ctx.arena().make<Product>(std::vector<ExprPtr>{inv_x, i_half}))});
    EXPECT_TRUE(structurally_equal_after_simplify(ctx, lhs, rhs))
        << "I_{3/2} != I_{-1/2} - (1/x)·I_{1/2}";
}

// J_{-3/2}(x) reduces to elementary closed form (no residual BesselJ FuncCall)
TEST_F(BesselE2Test, BesselJ_NegThreeHalves_ReducesToClosedForm) {
    ExprPtr x = make_symbol(ctx, "x");
    ExprPtr call = ctx.arena().make<FuncCall>(BuiltinOp::BesselJ,
        std::vector<ExprPtr>{make_half_order(ctx, -3), x});
    auto r = simplify(call);
    ASSERT_NE(r, nullptr);
    // Any residual BesselJ FuncCall indicates the half-integer recurrence did not fire.
    const auto* fc = expr_cast<FuncCall>(r);
    if (fc != nullptr) {
        EXPECT_NE(fc->func_id, BuiltinOp::BesselJ)
            << "J_{-3/2} must reduce to elementary form, got: " << debug_print(r);
    }
}

// K_{±1/2} symmetry: K_{-1/2}(x) = K_{1/2}(x)
TEST_F(BesselE2Test, BesselK_HalfIntegerSymmetry) {
    ExprPtr x = make_symbol(ctx, "x");
    ExprPtr k_plus = ctx.arena().make<FuncCall>(BuiltinOp::BesselK,
        std::vector<ExprPtr>{make_half_order(ctx, 1), x});
    ExprPtr k_minus = ctx.arena().make<FuncCall>(BuiltinOp::BesselK,
        std::vector<ExprPtr>{make_half_order(ctx, -1), x});
    EXPECT_TRUE(structurally_equal_after_simplify(ctx, k_plus, k_minus))
        << "K_{-1/2} must equal K_{+1/2}";
}

// Symmetry for negative integer order: J_{-2}(x) - J_2(x) = 0
TEST_F(BesselE2Test, BesselJ_NegativeIntegerSymmetry_Even) {
    ExprPtr x = make_symbol(ctx, "x");
    ExprPtr j_pos = bessel(BuiltinOp::BesselJ, 2, x);
    ExprPtr j_neg = bessel(BuiltinOp::BesselJ, -2, x);
    EXPECT_TRUE(structurally_equal_after_simplify(ctx, j_pos, j_neg))
        << "J_{-2}(x) must equal J_{2}(x) (even order)";
}

// I_{-3}(x) = I_3(x): even/odd-agnostic symmetry for I
TEST_F(BesselE2Test, BesselI_NegativeIntegerSymmetry) {
    ExprPtr x = make_symbol(ctx, "x");
    ExprPtr i_pos = bessel(BuiltinOp::BesselI, 3, x);
    ExprPtr i_neg = bessel(BuiltinOp::BesselI, -3, x);
    EXPECT_TRUE(structurally_equal_after_simplify(ctx, i_pos, i_neg))
        << "I_{-3}(x) must equal I_{3}(x) for all integer n";
}

// K_{-1}(x) = K_1(x)
TEST_F(BesselE2Test, BesselK_NegativeIntegerSymmetry) {
    ExprPtr x = make_symbol(ctx, "x");
    ExprPtr k_pos = bessel(BuiltinOp::BesselK, 1, x);
    ExprPtr k_neg = bessel(BuiltinOp::BesselK, -1, x);
    EXPECT_TRUE(structurally_equal_after_simplify(ctx, k_pos, k_neg))
        << "K_{-1}(x) must equal K_{1}(x)";
}

// I_n integer recurrence opt-in expands a Bessel call into a Sum containing
// lower-order Bessel calls (no residual order-3 Bessel call after one expansion).
TEST_F(BesselE2Test, BesselI_IntegerRecurrence_OptIn_Expands) {
    ctx.set_expand_bessel_recurrence(true);
    ExprPtr x = make_symbol(ctx, "x");
    ExprPtr call = bessel(BuiltinOp::BesselI, 3, x);
    auto r = simplify(call);
    ASSERT_NE(r, nullptr);
    // After recurrence the top-level node must no longer be a direct BesselI FuncCall.
    const auto* top_fc = expr_cast<FuncCall>(r);
    EXPECT_TRUE(top_fc == nullptr || top_fc->func_id != BuiltinOp::BesselI)
        << "I_3 recurrence did not expand, got: " << debug_print(r);
    ctx.set_expand_bessel_recurrence(false);
}

// K_n integer recurrence opt-in expands.
TEST_F(BesselE2Test, BesselK_IntegerRecurrence_OptIn_Expands) {
    ctx.set_expand_bessel_recurrence(true);
    ExprPtr x = make_symbol(ctx, "x");
    ExprPtr call = bessel(BuiltinOp::BesselK, 3, x);
    auto r = simplify(call);
    ASSERT_NE(r, nullptr);
    const auto* top_fc = expr_cast<FuncCall>(r);
    EXPECT_TRUE(top_fc == nullptr || top_fc->func_id != BuiltinOp::BesselK)
        << "K_3 recurrence did not expand, got: " << debug_print(r);
    ctx.set_expand_bessel_recurrence(false);
}

}  // namespace
}  // namespace cas
