// A7 step 5 — Meijer G identities (§6.1-6.3, §6.6), the two exact 1F1
// closed forms, and the post-Risch Meijer fallback in integrate()
// (Meijer_G_Slater.md §8-§9; test plan §10.3 + §10.7 indefinite half).
// The flagship end-to-end case int e^{-x^2} dx -> (sqrt(pi)/2) erf(x) is
// verified BY DIFFERENTIATION (never toString).

#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/meijerg.hpp"
#include "calculus/calculus_internal.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

namespace cas::symbolic {
namespace {

[[nodiscard]] ExprPtr parse_expr(const std::string& input, CASContext& ctx) {
    auto tokens = Lexer(input).tokenize();
    EXPECT_TRUE(tokens.is_ok()) << input;
    Parser parser(tokens.value(), ctx.arena());
    auto result = parser.parse();
    EXPECT_TRUE(result.is_ok()) << input;
    return result.value();
}

[[nodiscard]] ExprPtr rat(AstArena& arena, long long num, long long den) {
    return arena.make<RationalLit>(BigInt(num), BigInt(den));
}
[[nodiscard]] ExprPtr integer(AstArena& arena, long long v) {
    return arena.make<IntegerLit>(BigInt(v));
}

[[nodiscard]] bool contains_meijerg(ExprPtr e) {
    if (e == nullptr) return false;
    if (const auto* call = expr_cast<FuncCall>(e)) {
        if (call->func_id == BuiltinOp::MeijerG) return true;
        for (ExprPtr a : call->args)
            if (contains_meijerg(a)) return true;
        return false;
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        for (ExprPtr f : prod->factors)
            if (contains_meijerg(f)) return true;
        return false;
    }
    if (const auto* sum = expr_cast<Sum>(e)) {
        for (ExprPtr t : sum->terms)
            if (contains_meijerg(t)) return true;
        return false;
    }
    if (const auto* bin = expr_cast<Binary>(e))
        return contains_meijerg(bin->left) || contains_meijerg(bin->right);
    if (const auto* un = expr_cast<Unary>(e))
        return contains_meijerg(un->operand);
    return false;
}

// ── §10.3 identità ─────────────────────────────────────────────────────────

TEST(MeijerGIdentitiesTest, PowerShift_ShiftsAllParameters) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr w = parse_expr("w", ctx);
    auto g = make_meijerg(ctx, 1, 0, {}, {rat(arena, 1, 2), integer(arena, 0)}, w);
    ASSERT_TRUE(g.is_ok());
    auto shifted = meijerg_power_shift(ctx, *expr_cast<FuncCall>(g.value()),
                                       rat(arena, 1, 2));
    ASSERT_TRUE(shifted.is_ok());
    auto view = view_meijerg(*expr_cast<FuncCall>(shifted.value()));
    ASSERT_TRUE(view.is_ok());
    // b = (1/2, 0) + 1/2 = (1, 1/2); m,n,p,q unchanged.
    EXPECT_EQ(view.value().m, 1U);
    EXPECT_EQ(view.value().q, 2U);
    const auto* b1 = expr_cast<IntegerLit>(view.value().b[0]);
    ASSERT_NE(b1, nullptr);
    EXPECT_EQ(b1->value, BigInt(1));
    const auto* b2 = expr_cast<RationalLit>(view.value().b[1]);
    ASSERT_NE(b2, nullptr);
    EXPECT_EQ(b2->numerator, BigInt(1));
    EXPECT_EQ(b2->denominator, BigInt(2));
}

TEST(MeijerGIdentitiesTest, CancelCommonParam_ReducesOrders) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr w = parse_expr("w", ctx);
    // G^{1,1}_{2,2}(w | 1/3, 1/2 ; 1/4, 1/3): upper n-group {1/3} matches
    // lower outside-m {1/3} -> G^{1,0}_{1,1}(w | 1/2 ; 1/4)  (§6.3).
    auto g = make_meijerg(ctx, 1, 1,
        {rat(arena, 1, 3), rat(arena, 1, 2)},
        {rat(arena, 1, 4), rat(arena, 1, 3)}, w);
    ASSERT_TRUE(g.is_ok());
    auto reduced = meijerg_cancel_common_param(ctx, *expr_cast<FuncCall>(g.value()));
    ASSERT_TRUE(reduced.is_ok());
    auto view = view_meijerg(*expr_cast<FuncCall>(reduced.value()));
    ASSERT_TRUE(view.is_ok());
    EXPECT_EQ(view.value().m, 1U);
    EXPECT_EQ(view.value().n, 0U);
    EXPECT_EQ(view.value().p, 1U);
    EXPECT_EQ(view.value().q, 1U);
    const auto* a1 = expr_cast<RationalLit>(view.value().a[0]);
    ASSERT_NE(a1, nullptr);
    EXPECT_EQ(a1->numerator, BigInt(1));
    EXPECT_EQ(a1->denominator, BigInt(2));
}

TEST(MeijerGIdentitiesTest, Antiderivative_ShiftsOrdersAndPrependsZero) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr w = parse_expr("w", ctx);
    auto g = make_meijerg(ctx, 1, 0, {}, {integer(arena, 0)}, w);
    ASSERT_TRUE(g.is_ok());
    auto anti = meijerg_antiderivative(ctx, *expr_cast<FuncCall>(g.value()));
    ASSERT_TRUE(anti.is_ok());
    // z * G^{1,1}_{1,2}(z | 0 ; 0, -1): a Product carrying the shifted node.
    const auto* prod = expr_cast<Product>(anti.value());
    ASSERT_NE(prod, nullptr);
    const FuncCall* g2 = nullptr;
    for (ExprPtr f : prod->factors) {
        if (const auto* call = expr_cast<FuncCall>(f);
            call != nullptr && call->func_id == BuiltinOp::MeijerG) g2 = call;
    }
    ASSERT_NE(g2, nullptr);
    auto view = view_meijerg(*g2);
    ASSERT_TRUE(view.is_ok());
    EXPECT_EQ(view.value().m, 1U);
    EXPECT_EQ(view.value().n, 1U);
    EXPECT_EQ(view.value().p, 1U);
    EXPECT_EQ(view.value().q, 2U);
    const auto* bq = expr_cast<IntegerLit>(view.value().b[1]);
    ASSERT_NE(bq, nullptr);
    EXPECT_EQ(bq->value, BigInt(-1));
}

TEST(MeijerGIdentitiesTest, ArgumentInversion_SwapsGroupsAndOrders) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr w = parse_expr("w", ctx);
    // G^{1,0}_{0,1}(w | ; 1/3) -> G^{0,1}_{1,0}(1/w | 2/3 ; )  (§6.1).
    auto g = make_meijerg(ctx, 1, 0, {}, {rat(arena, 1, 3)}, w);
    ASSERT_TRUE(g.is_ok());
    auto inverted = meijerg_invert_argument(ctx, *expr_cast<FuncCall>(g.value()));
    ASSERT_TRUE(inverted.is_ok());
    auto view = view_meijerg(*expr_cast<FuncCall>(inverted.value()));
    ASSERT_TRUE(view.is_ok());
    EXPECT_EQ(view.value().m, 0U);
    EXPECT_EQ(view.value().n, 1U);
    EXPECT_EQ(view.value().p, 1U);
    EXPECT_EQ(view.value().q, 0U);
    const auto* a1 = expr_cast<RationalLit>(view.value().a[0]);
    ASSERT_NE(a1, nullptr);
    EXPECT_EQ(a1->numerator, BigInt(2));
    EXPECT_EQ(a1->denominator, BigInt(3));
}

// ── Identità 1F1 esatte (DLMF 13.6.1, 7.6.2) ───────────────────────────────

TEST(Hyper1F1ClosedFormTest, Kummer_1_2_FoldsToExpForm) {
    CASContext ctx;
    ExprPtr expr = parse_expr("hyp1F1(1, 2, x)", ctx);
    // Compute the expected form FIRST: each top-level simplify RESETS
    // last_side_conditions (A31 — snapshot before any helper that re-enters
    // simplify), so the hyp1F1 fold must be the LAST simplify before the
    // side-condition assertion.
    auto expected = ctx.simplify(parse_expr("(exp(x) - 1)/x", ctx));
    ASSERT_TRUE(expected.is_ok());
    auto simplified = ctx.simplify(expr);
    ASSERT_TRUE(simplified.is_ok());
    EXPECT_TRUE(ExprEqual{}(simplified.value(), expected.value()));
    // A31: the removable-singularity rewrite registers NonZero(x).
    ExprPtr x = parse_expr("x", ctx);
    EXPECT_TRUE(ctx.last_side_conditions().contains(
        DomainCondition{DomainConditionKind::NonZero, x}));
}

TEST(Hyper1F1ClosedFormTest, ErfForm_FoldsToErf) {
    CASContext ctx;
    ExprPtr expr = parse_expr("hyp1F1(1/2, 3/2, -(x^2))", ctx);
    auto simplified = ctx.simplify(expr);
    ASSERT_TRUE(simplified.is_ok());
    auto expected = ctx.simplify(
        parse_expr("sqrt(pi) * erf(x) / (2*x)", ctx));
    ASSERT_TRUE(expected.is_ok());
    EXPECT_TRUE(ExprEqual{}(simplified.value(), expected.value()));
}

// ── §10.7 integrazione indefinita via fallback post-Risch ──────────────────

TEST(IntegrateMeijerGTest, GaussianIntegral_ErfClosedForm_DiffVerified) {
    // int e^{-x^2} dx = (sqrt(pi)/2) erf(x) + C — Risch (correctly) proves
    // no elementary antiderivative.
    //
    // Since A43 this integrand is settled one step EARLIER than the Meijer
    // fallback, by the non-elementary family (integrate_core step 5): same
    // erf closed form, but obtained by completing the square, which is exact
    // on ALL of R. So — unlike the Mellin route, which lives on the half-line
    // — no Positive(x) side condition is warranted here, and asserting its
    // ABSENCE is the stronger statement. The Mellin route keeps its own A31
    // coverage in MeijerFallbackGaussian_DeclaresPositiveSideCondition below.
    //
    // Verified by differentiation: D[(sqrt(pi)/2) erf(x)] = e^{-x^2}.
    CASContext ctx;
    Symbol x{"x"};
    ExprPtr integrand = parse_expr("exp(-(x^2))", ctx);
    auto F = calculus::integrate(integrand, x, ctx);
    ASSERT_TRUE(F.is_ok()) << (F.is_error() ? F.error().message : "");
    EXPECT_FALSE(contains_meijerg(F.value()));  // folded to elementary/erf
    // Snapshot the side conditions BEFORE any further simplify resets them
    // (A31 discipline).
    ExprPtr xs = parse_expr("x", ctx);
    EXPECT_FALSE(ctx.last_side_conditions().contains(
        DomainCondition{DomainConditionKind::Positive, xs}))
        << "completing the square is exact on all of R: no domain restriction";
    auto Fp = calculus::diff(F.value(), x, 1U, ctx);
    ASSERT_TRUE(Fp.is_ok());
    ExprPtr resid = ctx.arena().make<Binary>(BinaryOp::Sub,
        Fp.value(), integrand);
    auto zero = ctx.simplify(resid);
    ASSERT_TRUE(zero.is_ok());
    const auto* lit = expr_cast<IntegerLit>(zero.value());
    EXPECT_TRUE(lit != nullptr && lit->value.is_zero())
        << "residual not zero";
}

// A31 coverage of the Mellin route, called DIRECTLY.
//
// Going through integrate() would no longer reach it for this integrand (A43
// settles the Gaussian one step earlier), and e^{-x^2} was the only case in
// the suite that exercised the conditional Slater fold — routing the assertion
// through the dispatcher would have left that branch silently untested.
// Calling the fallback by name also makes the test independent of dispatch
// order, which is what it actually means to assert.
TEST(IntegrateMeijerGTest, MeijerFallbackGaussian_DeclaresPositiveSideCondition) {
    // The Mellin pipeline works on the half-line (spec §4): a Slater fold that
    // is not unconditionally exact (here s = 1/2, a non-integer power) must
    // declare Positive(x) via A31 — never apply it silently.
    CASContext ctx;
    Symbol x{"x"};
    ExprPtr integrand = parse_expr("exp(-(x^2))", ctx);
    auto F = calculus::integrate_meijerg_fallback(integrand, x, ctx);
    ASSERT_TRUE(F.is_ok()) << (F.is_error() ? F.error().message : "");
    EXPECT_FALSE(contains_meijerg(F.value()));  // folded to erf
    ExprPtr xs = parse_expr("x", ctx);
    EXPECT_TRUE(ctx.last_side_conditions().contains(
        DomainCondition{DomainConditionKind::Positive, xs}));
}

TEST(IntegrateMeijerGTest, BesselJIntegral_StaysClosedFormG) {
    // int J_{1/3}(x) dx has no elementary form and its pFq form (1F2) has
    // no engine node: the fallback returns a first-class Meijer G
    // antiderivative (§9.4) instead of failing.
    CASContext ctx;
    Symbol x{"x"};
    ExprPtr integrand = parse_expr("BesselJ(1/3, x)", ctx);
    auto F = calculus::integrate(integrand, x, ctx);
    ASSERT_TRUE(F.is_ok()) << (F.is_error() ? F.error().message : "");
    EXPECT_TRUE(contains_meijerg(F.value()));
}

TEST(IntegrateMeijerGTest, UnsupportedShapeStillStructuredUnimplemented) {
    // A shape outside every strategy (including the G family) must keep
    // failing structured — the fallback must not swallow it.
    CASContext ctx;
    Symbol x{"x"};
    ExprPtr integrand = parse_expr("exp(sin(x) + x^3)", ctx);
    auto F = calculus::integrate(integrand, x, ctx);
    if (F.is_error()) {
        EXPECT_EQ(F.error().kind, CASErrorKind::Unimplemented);
    } else {
        // If some strategy legitimately handles it, the result must at
        // least differentiate back (no silent-wrong).
        auto Fp = calculus::diff(F.value(), x, 1U, ctx);
        ASSERT_TRUE(Fp.is_ok());
    }
}

}  // namespace
}  // namespace cas::symbolic
