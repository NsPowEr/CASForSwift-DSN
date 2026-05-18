#include "cas/algebraic_number_bridge.hpp"
#include "cas/ast_debug.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <vector>

namespace cas::test {

class AlgebraicNumberBridgeTest : public ::testing::Test {
protected:
    void SetUp() override { ctx = std::make_unique<symbolic::CASContext>(); }

    [[nodiscard]] ExprPtr make_integer(long long n) {
        return ctx->arena().make<IntegerLit>(BigInt(n));
    }

    [[nodiscard]] ExprPtr make_sqrt2_rootof() {
        ExprPtr x_sym = ctx->arena().make<Symbol>("x");
        ExprPtr x2 = ctx->arena().make<Binary>(BinaryOp::Pow, x_sym, make_integer(2));
        ExprPtr poly = ctx->arena().make<Binary>(BinaryOp::Sub, x2, make_integer(2));
        return ctx->arena().make<RootOf>(poly, Symbol("x"), 0U);
    }

    // RootOf(x^3 - 2, x, 0)  →  generator alpha with alpha^3 = 2.
    [[nodiscard]] ExprPtr make_cuberoot2_rootof() {
        ExprPtr x_sym = ctx->arena().make<Symbol>("x");
        ExprPtr x3 = ctx->arena().make<Binary>(BinaryOp::Pow, x_sym, make_integer(3));
        ExprPtr poly = ctx->arena().make<Binary>(BinaryOp::Sub, x3, make_integer(2));
        return ctx->arena().make<RootOf>(poly, Symbol("x"), 0U);
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

TEST_F(AlgebraicNumberBridgeTest, MinPolyExtractedAndCanonicalized) {
    ExprPtr root_expr = make_sqrt2_rootof();
    const auto& root = expr_ref<RootOf>(root_expr);

    auto mp = cas::algebra::rootof_min_poly(root, *ctx);
    ASSERT_TRUE(mp.is_ok());

    const auto& coeffs = mp.value();
    ASSERT_EQ(coeffs.size(), 3U);
    EXPECT_EQ(coeffs[0], Rational(BigInt(-2)));
    EXPECT_EQ(coeffs[1], Rational(BigInt(0)));
    EXPECT_EQ(coeffs[2], Rational(BigInt(1)));
}

TEST_F(AlgebraicNumberBridgeTest, AlphaFromRootOfIsXItself) {
    ExprPtr root_expr = make_sqrt2_rootof();
    const auto& root = expr_ref<RootOf>(root_expr);

    auto alpha_res = cas::algebra::alpha_from_rootof(root, *ctx);
    ASSERT_TRUE(alpha_res.is_ok());

    const auto& v = alpha_res.value().value();
    ASSERT_EQ(v.size(), 2U);
    EXPECT_EQ(v[0], Rational(BigInt(0)));
    EXPECT_EQ(v[1], Rational(BigInt(1)));
}

TEST_F(AlgebraicNumberBridgeTest, ExpressIntegerAsConstant) {
    ExprPtr root_expr = make_sqrt2_rootof();
    const auto& root = expr_ref<RootOf>(root_expr);
    auto mp = cas::algebra::rootof_min_poly(root, *ctx).value();

    auto seven = make_integer(7);
    auto res = cas::algebra::try_express_in_q_alpha(seven, root_expr, mp, *ctx);
    ASSERT_TRUE(res.is_ok());
    ASSERT_TRUE(res.value().has_value());

    const auto& v = res.value().value().value();
    ASSERT_EQ(v.size(), 1U);
    EXPECT_EQ(v[0], Rational(BigInt(7)));
}

TEST_F(AlgebraicNumberBridgeTest, ExpressAlphaPlusOne) {
    // alpha + 1 in Q(sqrt(2))  →  value [1, 1]
    ExprPtr root_expr = make_sqrt2_rootof();
    const auto& root = expr_ref<RootOf>(root_expr);
    auto mp = cas::algebra::rootof_min_poly(root, *ctx).value();

    ExprPtr e = ctx->arena().make<Binary>(BinaryOp::Add, root_expr, make_integer(1));
    auto res = cas::algebra::try_express_in_q_alpha(e, root_expr, mp, *ctx);
    ASSERT_TRUE(res.is_ok());
    ASSERT_TRUE(res.value().has_value());

    const auto& v = res.value().value().value();
    ASSERT_EQ(v.size(), 2U);
    EXPECT_EQ(v[0], Rational(BigInt(1)));
    EXPECT_EQ(v[1], Rational(BigInt(1)));
}

TEST_F(AlgebraicNumberBridgeTest, ExpressDifferenceOfSquaresReducesToRational) {
    // (alpha+1)*(alpha-1) = alpha^2 - 1 = 2 - 1 = 1 over Q(sqrt(2))
    ExprPtr root_expr = make_sqrt2_rootof();
    const auto& root = expr_ref<RootOf>(root_expr);
    auto mp = cas::algebra::rootof_min_poly(root, *ctx).value();

    ExprPtr alpha_plus_1  = ctx->arena().make<Binary>(BinaryOp::Add, root_expr, make_integer(1));
    ExprPtr alpha_minus_1 = ctx->arena().make<Binary>(BinaryOp::Sub, root_expr, make_integer(1));
    ExprPtr product = ctx->arena().make<Binary>(BinaryOp::Mul, alpha_plus_1, alpha_minus_1);

    auto res = cas::algebra::try_express_in_q_alpha(product, root_expr, mp, *ctx);
    ASSERT_TRUE(res.is_ok());
    ASSERT_TRUE(res.value().has_value());

    const auto& v = res.value().value().value();
    ASSERT_EQ(v.size(), 1U);
    EXPECT_EQ(v[0], Rational(BigInt(1)));
}

TEST_F(AlgebraicNumberBridgeTest, ExpressDivisionByAlpha) {
    // 1 / alpha in Q(sqrt(2)) = sqrt(2)/2  →  value [0, 1/2]
    ExprPtr root_expr = make_sqrt2_rootof();
    const auto& root = expr_ref<RootOf>(root_expr);
    auto mp = cas::algebra::rootof_min_poly(root, *ctx).value();

    ExprPtr inv = ctx->arena().make<Binary>(BinaryOp::Div, make_integer(1), root_expr);
    auto res = cas::algebra::try_express_in_q_alpha(inv, root_expr, mp, *ctx);
    ASSERT_TRUE(res.is_ok());
    ASSERT_TRUE(res.value().has_value());

    const auto& v = res.value().value().value();
    ASSERT_EQ(v.size(), 2U);
    EXPECT_EQ(v[0], Rational(BigInt(0)));
    EXPECT_EQ(v[1], Rational(BigInt(1), BigInt(2)));
}

TEST_F(AlgebraicNumberBridgeTest, ExpressRejectsUnknownSymbol) {
    ExprPtr root_expr = make_sqrt2_rootof();
    const auto& root = expr_ref<RootOf>(root_expr);
    auto mp = cas::algebra::rootof_min_poly(root, *ctx).value();

    ExprPtr y = ctx->arena().make<Symbol>("y");
    ExprPtr e = ctx->arena().make<Binary>(BinaryOp::Add, root_expr, y);
    auto res = cas::algebra::try_express_in_q_alpha(e, root_expr, mp, *ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_FALSE(res.value().has_value());
}

TEST_F(AlgebraicNumberBridgeTest, ExpressRejectsTranscendental) {
    ExprPtr root_expr = make_sqrt2_rootof();
    const auto& root = expr_ref<RootOf>(root_expr);
    auto mp = cas::algebra::rootof_min_poly(root, *ctx).value();

    ExprPtr x_sym = ctx->arena().make<Symbol>("x");
    ExprPtr sinx  = ctx->arena().make<FuncCall>(BuiltinOp::Sin, std::vector<ExprPtr>{x_sym});
    ExprPtr e     = ctx->arena().make<Binary>(BinaryOp::Add, root_expr, sinx);
    auto res = cas::algebra::try_express_in_q_alpha(e, root_expr, mp, *ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_FALSE(res.value().has_value());
}

TEST_F(AlgebraicNumberBridgeTest, CubicAlphaPlusOneCubedReducesViaMinPoly) {
    // (alpha+1)^3 over Q(alpha) with alpha^3 = 2
    //   = alpha^3 + 3*alpha^2 + 3*alpha + 1
    //   = 2 + 3*alpha^2 + 3*alpha + 1
    //   = 3 + 3*alpha + 3*alpha^2     →  value [3, 3, 3]
    ExprPtr root_expr = make_cuberoot2_rootof();
    const auto& root = expr_ref<RootOf>(root_expr);
    auto mp = cas::algebra::rootof_min_poly(root, *ctx).value();

    ExprPtr alpha_plus_1 = ctx->arena().make<Binary>(BinaryOp::Add, root_expr, make_integer(1));
    ExprPtr cube = ctx->arena().make<Binary>(BinaryOp::Pow, alpha_plus_1, make_integer(3));

    auto res = cas::algebra::try_express_in_q_alpha(cube, root_expr, mp, *ctx);
    ASSERT_TRUE(res.is_ok());
    ASSERT_TRUE(res.value().has_value());

    const auto& v = res.value().value().value();
    ASSERT_EQ(v.size(), 3U);
    EXPECT_EQ(v[0], Rational(BigInt(3)));
    EXPECT_EQ(v[1], Rational(BigInt(3)));
    EXPECT_EQ(v[2], Rational(BigInt(3)));
}

TEST_F(AlgebraicNumberBridgeTest, RoundTripAlgebraicNumberToExprAndBack) {
    // Build  e = alpha^2 + alpha + 3   over Q(sqrt(2)) → value [3, 1, 1]
    // alpha^2 = 2 so reduced: value [5, 1].
    ExprPtr root_expr = make_sqrt2_rootof();
    const auto& root = expr_ref<RootOf>(root_expr);
    auto mp = cas::algebra::rootof_min_poly(root, *ctx).value();

    ExprPtr alpha_sq = ctx->arena().make<Binary>(BinaryOp::Pow, root_expr, make_integer(2));
    ExprPtr sum1     = ctx->arena().make<Binary>(BinaryOp::Add, alpha_sq, root_expr);
    ExprPtr e        = ctx->arena().make<Binary>(BinaryOp::Add, sum1, make_integer(3));

    auto an_res = cas::algebra::try_express_in_q_alpha(e, root_expr, mp, *ctx);
    ASSERT_TRUE(an_res.is_ok());
    ASSERT_TRUE(an_res.value().has_value());

    const auto& v = an_res.value().value().value();
    ASSERT_EQ(v.size(), 2U);
    EXPECT_EQ(v[0], Rational(BigInt(5)));
    EXPECT_EQ(v[1], Rational(BigInt(1)));

    auto back = cas::algebra::algebraic_number_to_expr(an_res.value().value(), root_expr, *ctx);
    ASSERT_TRUE(back.is_ok());

    // Verify mathematical equivalence with the original expression.
    auto equal = symbolic::mathematically_equal(back.value(), e, *ctx);
    ASSERT_TRUE(equal.is_ok());
    EXPECT_TRUE(equal.value()) << "Round-trip mismatch. Got: " << debug_print(back.value());
}

TEST_F(AlgebraicNumberBridgeTest, ExpressZeroIsZero) {
    ExprPtr root_expr = make_sqrt2_rootof();
    const auto& root = expr_ref<RootOf>(root_expr);
    auto mp = cas::algebra::rootof_min_poly(root, *ctx).value();

    auto res = cas::algebra::try_express_in_q_alpha(make_integer(0), root_expr, mp, *ctx);
    ASSERT_TRUE(res.is_ok());
    ASSERT_TRUE(res.value().has_value());
    EXPECT_TRUE(res.value().value().is_zero());
}

TEST_F(AlgebraicNumberBridgeTest, ReduceInQAlpha_ProductDiffSqrt2) {
    // (alpha + 1) * (alpha - 1) over Q(sqrt(2)) -> 1
    ExprPtr root_expr = make_sqrt2_rootof();
    ExprPtr alpha_plus_1  = ctx->arena().make<Binary>(BinaryOp::Add, root_expr, make_integer(1));
    ExprPtr alpha_minus_1 = ctx->arena().make<Binary>(BinaryOp::Sub, root_expr, make_integer(1));
    ExprPtr product = ctx->arena().make<Binary>(BinaryOp::Mul, alpha_plus_1, alpha_minus_1);

    auto reduced = cas::algebra::simplify_in_q_alpha(product, *ctx);
    ASSERT_TRUE(reduced.is_ok());

    auto expected = make_integer(1);
    auto eq = symbolic::mathematically_equal(reduced.value(), expected, *ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value()) << "Got: " << debug_print(reduced.value());
}

TEST_F(AlgebraicNumberBridgeTest, ReduceInQAlpha_CubicReductionAlphaCubedPlusAlpha) {
    // alpha^3 + alpha   over Q(alpha) with alpha^3 = 2   -> 2 + alpha
    ExprPtr root_expr = make_cuberoot2_rootof();
    ExprPtr alpha_cubed = ctx->arena().make<Binary>(BinaryOp::Pow, root_expr, make_integer(3));
    ExprPtr sum = ctx->arena().make<Binary>(BinaryOp::Add, alpha_cubed, root_expr);

    auto reduced = cas::algebra::simplify_in_q_alpha(sum, *ctx);
    ASSERT_TRUE(reduced.is_ok());

    ExprPtr expected = ctx->arena().make<Binary>(BinaryOp::Add, make_integer(2), root_expr);
    auto eq = symbolic::mathematically_equal(reduced.value(), expected, *ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value()) << "Got: " << debug_print(reduced.value());
}

TEST_F(AlgebraicNumberBridgeTest, ReduceInQAlpha_LeavesExpressionWithoutRootOfUnchanged) {
    ExprPtr x = ctx->arena().make<Symbol>("x");
    ExprPtr e = ctx->arena().make<Binary>(BinaryOp::Add, x, make_integer(3));

    auto reduced = cas::algebra::simplify_in_q_alpha(e, *ctx);
    ASSERT_TRUE(reduced.is_ok());

    auto eq = symbolic::mathematically_equal(reduced.value(), e, *ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

TEST_F(AlgebraicNumberBridgeTest, ReduceInQAlpha_DoesNothingWithTwoDistinctRootOfs) {
    // sqrt(2) + cuberoot(2)  -> currently bridge handles only single extension
    ExprPtr a = make_sqrt2_rootof();
    ExprPtr b = make_cuberoot2_rootof();
    ExprPtr e = ctx->arena().make<Binary>(BinaryOp::Add, a, b);

    auto reduced = cas::algebra::simplify_in_q_alpha(e, *ctx);
    ASSERT_TRUE(reduced.is_ok());

    auto eq = symbolic::mathematically_equal(reduced.value(), e, *ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value()) << "Got: " << debug_print(reduced.value());
}

TEST_F(AlgebraicNumberBridgeTest, ReduceInQAlpha_DivisionByAlphaProducesHalfAlpha) {
    // 1 / alpha   over Q(sqrt(2))   -> alpha / 2
    ExprPtr root_expr = make_sqrt2_rootof();
    ExprPtr inv = ctx->arena().make<Binary>(BinaryOp::Div, make_integer(1), root_expr);

    auto reduced = cas::algebra::simplify_in_q_alpha(inv, *ctx);
    ASSERT_TRUE(reduced.is_ok());

    // Build expected:  alpha / 2  (or equivalently (1/2) * alpha)
    ExprPtr half = ctx->arena().make<RationalLit>(BigInt(1), BigInt(2));
    ExprPtr expected = ctx->arena().make<Binary>(BinaryOp::Mul, half, root_expr);
    auto eq = symbolic::mathematically_equal(reduced.value(), expected, *ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value()) << "Got: " << debug_print(reduced.value());
}

TEST_F(AlgebraicNumberBridgeTest, MinPolyRejectsRootOfWithSymbolicCoefficients) {
    // RootOf(x^2 - y, x)  →  non-rational min_poly → bridge must refuse.
    ExprPtr x_sym = ctx->arena().make<Symbol>("x");
    ExprPtr y_sym = ctx->arena().make<Symbol>("y");
    ExprPtr x2 = ctx->arena().make<Binary>(BinaryOp::Pow, x_sym, make_integer(2));
    ExprPtr poly = ctx->arena().make<Binary>(BinaryOp::Sub, x2, y_sym);
    ExprPtr root_expr = ctx->arena().make<RootOf>(poly, Symbol("x"), 0U);
    const auto& root = expr_ref<RootOf>(root_expr);

    auto mp = cas::algebra::rootof_min_poly(root, *ctx);
    EXPECT_TRUE(mp.is_error());
}

// ── L1-05: register_algebraic_simplify_hook auto-trigger tests ──────────────

class AlgebraicSimplifyHookTest : public ::testing::Test {
protected:
    void SetUp() override {
        ctx = std::make_unique<symbolic::CASContext>();
        cas::algebra::register_algebraic_simplify_hook(*ctx);
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

TEST_F(AlgebraicSimplifyHookTest, HookRegistered) {
    EXPECT_TRUE(ctx->has_post_simplify_hook());
}

TEST_F(AlgebraicSimplifyHookTest, ClearHookDisables) {
    ctx->clear_post_simplify_hook();
    EXPECT_FALSE(ctx->has_post_simplify_hook());
}

TEST_F(AlgebraicSimplifyHookTest, CubeRoot2_AlphaCubed_ReducesToTwo) {
    // alpha = RootOf(x^3 - 2, x, 0)  →  alpha^3 should reduce to 2 via hook.
    // Degree-3 RootOf stays as RootOf (not auto-solved), so the hook fires.
    auto& ar = ctx->arena();
    ExprPtr x_sym = ar.make<Symbol>("x");
    ExprPtr x3 = ar.make<Binary>(BinaryOp::Pow, x_sym, ar.make<IntegerLit>(BigInt(3)));
    ExprPtr poly = ar.make<Binary>(BinaryOp::Sub, x3, ar.make<IntegerLit>(BigInt(2)));
    ExprPtr alpha = ar.make<RootOf>(poly, Symbol("x"), 0U);

    // alpha^3
    ExprPtr expr = ar.make<Binary>(
        BinaryOp::Pow, alpha, ar.make<IntegerLit>(BigInt(3)));

    auto res = ctx->simplify(expr);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    // Expected: IntegerLit(2)
    bool is_two = expr_is<IntegerLit>(res.value())
        && expr_ref<IntegerLit>(res.value()).value == BigInt(2);
    EXPECT_TRUE(is_two) << "Got: " << debug_print(res.value());
}

TEST_F(AlgebraicSimplifyHookTest, CubeRoot2_AlphaFourth_ReducesToTwoAlpha) {
    // alpha^4 = alpha * alpha^3 = 2*alpha.
    auto& ar = ctx->arena();
    ExprPtr x_sym = ar.make<Symbol>("x");
    ExprPtr x3 = ar.make<Binary>(BinaryOp::Pow, x_sym, ar.make<IntegerLit>(BigInt(3)));
    ExprPtr poly = ar.make<Binary>(BinaryOp::Sub, x3, ar.make<IntegerLit>(BigInt(2)));
    ExprPtr alpha = ar.make<RootOf>(poly, Symbol("x"), 0U);

    // alpha^4
    ExprPtr expr = ar.make<Binary>(
        BinaryOp::Pow, alpha, ar.make<IntegerLit>(BigInt(4)));

    auto res = ctx->simplify(expr);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    // Expected: 2*alpha  →  diff from (result - 2*alpha) = 0
    ExprPtr two_alpha = ar.make<Binary>(
        BinaryOp::Mul, ar.make<IntegerLit>(BigInt(2)), alpha);
    auto diff = ar.make<Binary>(BinaryOp::Sub, res.value(), two_alpha);
    auto diff_s = ctx->simplify(diff);
    ASSERT_TRUE(diff_s.is_ok());
    bool is_zero = expr_is<IntegerLit>(diff_s.value())
        && expr_ref<IntegerLit>(diff_s.value()).value == BigInt(0);
    EXPECT_TRUE(is_zero) << "Got: " << debug_print(res.value())
                         << "  diff: " << debug_print(diff_s.value());
}

TEST_F(AlgebraicSimplifyHookTest, NoRootOf_ExprUnchanged) {
    // Without RootOf, hook is a no-op: 2 + 3 = 5
    auto& ar = ctx->arena();
    ExprPtr two = ar.make<IntegerLit>(BigInt(2));
    ExprPtr three = ar.make<IntegerLit>(BigInt(3));
    ExprPtr expr = ar.make<Binary>(BinaryOp::Add, two, three);
    auto res = ctx->simplify(expr);
    ASSERT_TRUE(res.is_ok());
    ASSERT_TRUE(expr_is<IntegerLit>(res.value()));
    EXPECT_EQ(expr_ref<IntegerLit>(res.value()).value, BigInt(5));
}

}  // namespace cas::test
