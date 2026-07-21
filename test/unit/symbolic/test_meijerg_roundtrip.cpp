// A7 step 4 — from_meijerg / slater_expand / expand_meijerg_nodes
// (Meijer_G_Slater.md §3.2, §5 inverse, §9.4; test plan §10.1 round-trips +
// §10.6 Slater). Sign convention: round-trip equality holds on the Mellin
// half-line (positive elementary argument) — see meijerg_from.cpp header.

#include "cas/lexer.hpp"
#include "cas/meijerg.hpp"
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

// to_meijerg(f) -> expand_meijerg_nodes -> simplify == simplify(f)?
void expect_roundtrip(const char* input) {
    CASContext ctx;
    ExprPtr f = parse_expr(input, ctx);
    auto g_form = to_meijerg(ctx, f);
    ASSERT_TRUE(g_form.is_ok()) << input;
    auto back = expand_meijerg_nodes(ctx, g_form.value());
    ASSERT_TRUE(back.is_ok()) << input;
    auto lhs = ctx.simplify(back.value());
    auto rhs = ctx.simplify(f);
    ASSERT_TRUE(lhs.is_ok() && rhs.is_ok()) << input;
    EXPECT_TRUE(ExprEqual{}(lhs.value(), rhs.value()))
        << input << "\n  back: " << (lhs.is_ok() ? "ok" : "err");
}

// ── §10.1 round-trip elementare → G → elementare ───────────────────────────

TEST(MeijerGRoundTripTest, Exp)      { expect_roundtrip("exp(x)"); }
TEST(MeijerGRoundTripTest, Sin)      { expect_roundtrip("sin(x)"); }
TEST(MeijerGRoundTripTest, Cos)      { expect_roundtrip("cos(x)"); }
TEST(MeijerGRoundTripTest, Sinh)     { expect_roundtrip("sinh(x)"); }
TEST(MeijerGRoundTripTest, Cosh)     { expect_roundtrip("cosh(x)"); }
TEST(MeijerGRoundTripTest, LnOnePlus){ expect_roundtrip("ln(1 + x)"); }
TEST(MeijerGRoundTripTest, Atan)     { expect_roundtrip("arctan(x)"); }
TEST(MeijerGRoundTripTest, Asin)     { expect_roundtrip("arcsin(x)"); }
TEST(MeijerGRoundTripTest, Erf)      { expect_roundtrip("erf(x)"); }
TEST(MeijerGRoundTripTest, BesselJ)  { expect_roundtrip("BesselJ(1/3, x)"); }

// §10.4 — incomplete gamma (Meijer_G_Slater.md §5.9). Order 1/3 chosen so the
// b-group (0, 1/3) is non-integer-spaced (the clean, non-confluent case).
TEST(MeijerGRoundTripTest, GammaUpper) {
    expect_roundtrip("gamma_incomplete(1/3, x)");
}
TEST(MeijerGRoundTripTest, GammaLower) {
    expect_roundtrip("gamma_incomplete_lower(1/3, x)");
}

// from_meijerg direct: the §5.9 G shapes fold back to Γ/γ (not just via the
// forward table). Γ(a,z) = G^{2,0}_{1,2}(z | 1 ; 0,a).
TEST(MeijerGFromTest, UpperGammaShapeFolds) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr z = parse_expr("z", ctx);
    auto g = make_meijerg(ctx, 2, 0, {integer(arena, 1)},
        {integer(arena, 0), rat(arena, 1, 3)}, z);
    ASSERT_TRUE(g.is_ok());
    const auto* g_call = expr_cast<FuncCall>(g.value());
    ASSERT_NE(g_call, nullptr);
    auto folded = from_meijerg(ctx, *g_call);
    ASSERT_TRUE(folded.is_ok());
    auto expected = ctx.simplify(parse_expr("gamma_incomplete(1/3, z)", ctx));
    ASSERT_TRUE(expected.is_ok());
    EXPECT_TRUE(ExprEqual{}(folded.value(), expected.value()));
}

// γ(a,z) = G^{1,1}_{1,2}(z | 1 ; a,0) with a ≠ 1/2 (distinct from the erf entry
// that shares this shape only for a = 1/2 AND a squared argument).
TEST(MeijerGFromTest, LowerGammaShapeFolds) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr z = parse_expr("z", ctx);
    auto g = make_meijerg(ctx, 1, 1, {integer(arena, 1)},
        {rat(arena, 1, 3), integer(arena, 0)}, z);
    ASSERT_TRUE(g.is_ok());
    const auto* g_call = expr_cast<FuncCall>(g.value());
    ASSERT_NE(g_call, nullptr);
    auto folded = from_meijerg(ctx, *g_call);
    ASSERT_TRUE(folded.is_ok());
    auto expected = ctx.simplify(parse_expr("gamma_incomplete_lower(1/3, z)", ctx));
    ASSERT_TRUE(expected.is_ok());
    EXPECT_TRUE(ExprEqual{}(folded.value(), expected.value()));
}

// The erf entry (§5.7) must NOT be captured by the γ recognizer: erf(x)
// still round-trips to erf, not to γ(1/2, x²).
TEST(MeijerGFromTest, ErfStillWinsOverLowerGamma) {
    expect_roundtrip("erf(x)");
}

TEST(MeijerGRoundTripTest, PowerExpProduct) {
    // z^{1/3} e^{-z}: the (1,0,0,1) node has NO table entry — this
    // round-trip exercises the GENERAL Slater path (0F0 closed form).
    CASContext ctx;
    auto input = ctx.simplify(parse_expr("x^(1/3) * exp(-x)", ctx));
    ASSERT_TRUE(input.is_ok());
    auto g_form = to_meijerg(ctx, input.value());
    ASSERT_TRUE(g_form.is_ok());
    auto back = expand_meijerg_nodes(ctx, g_form.value());
    ASSERT_TRUE(back.is_ok());
    auto lhs = ctx.simplify(back.value());
    ASSERT_TRUE(lhs.is_ok());
    EXPECT_TRUE(ExprEqual{}(lhs.value(), input.value()));
}

TEST(MeijerGRoundTripTest, BinomialPower) {
    // (1+x)^{-1/3}: (1,1,1,1) node, general Slater path (1F0 closed form).
    CASContext ctx;
    ExprPtr f = parse_expr("(1 + x)^(-1/3)", ctx);
    auto g_form = to_meijerg(ctx, f);
    ASSERT_TRUE(g_form.is_ok());
    auto back = expand_meijerg_nodes(ctx, g_form.value());
    ASSERT_TRUE(back.is_ok());
    auto lhs = ctx.simplify(back.value());
    auto rhs = ctx.simplify(f);
    ASSERT_TRUE(lhs.is_ok() && rhs.is_ok());
    EXPECT_TRUE(ExprEqual{}(lhs.value(), rhs.value()));
}

// ── §10.6 Slater inverso ───────────────────────────────────────────────────

TEST(MeijerGSlaterTest, ArctanG_ExpandsTo2F1) {
    // G^{1,2}_{2,2}(w | 1/2,1 ; 1/2 ; 0) = 2 sqrt(w) 2F1(1, 1/2; 3/2; -w)
    // (§3.2 worked example, verified analytically + numerically in the spec).
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr w = parse_expr("w", ctx);
    auto g = make_meijerg(ctx, 1, 2,
        {rat(arena, 1, 2), integer(arena, 1)},
        {rat(arena, 1, 2), integer(arena, 0)}, w);
    ASSERT_TRUE(g.is_ok());
    const auto* g_call = expr_cast<FuncCall>(g.value());
    ASSERT_NE(g_call, nullptr);
    auto expanded = slater_expand(ctx, *g_call);
    ASSERT_TRUE(expanded.is_ok());
    auto expected = ctx.simplify(
        parse_expr("2 * sqrt(w) * hyp2F1(1, 1/2, 3/2, -w)", ctx));
    ASSERT_TRUE(expected.is_ok());
    auto actual = ctx.simplify(expanded.value());
    ASSERT_TRUE(actual.is_ok());
    EXPECT_TRUE(ExprEqual{}(actual.value(), expected.value()));
}

TEST(MeijerGSlaterTest, TwoTermSum_MEquals2) {
    // m = 2 with non-integer-spaced b = (0, 1/3): the expansion is a SUM of
    // two pFq terms (§10.6 second case).
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr w = parse_expr("w", ctx);
    auto g = make_meijerg(ctx, 2, 0, {integer(arena, 1)},
        {integer(arena, 0), rat(arena, 1, 3)}, w);
    ASSERT_TRUE(g.is_ok());
    const auto* g_call = expr_cast<FuncCall>(g.value());
    ASSERT_NE(g_call, nullptr);
    auto expanded = slater_expand(ctx, *g_call);
    ASSERT_TRUE(expanded.is_ok());
    // Structural sanity: a Sum whose terms each contain a 1F1 node
    // (p=1, q-1=1); exact coefficients are Gamma-valued and engine-reduced.
    const auto* sum = expr_cast<Sum>(expanded.value());
    ASSERT_NE(sum, nullptr);
    EXPECT_EQ(sum->terms.size(), 2U);
}

TEST(MeijerGSlaterTest, ConfluentPolesRefused) {
    // b = (0, 1) inside the m-group: integer spacing -> logarithmic terms,
    // the Slater series must be REFUSED (structured), never emitted wrong.
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr w = parse_expr("w", ctx);
    auto g = make_meijerg(ctx, 2, 0, {},
        {integer(arena, 0), integer(arena, 1)}, w);
    ASSERT_TRUE(g.is_ok());
    const auto* g_call = expr_cast<FuncCall>(g.value());
    ASSERT_NE(g_call, nullptr);
    auto expanded = slater_expand(ctx, *g_call);
    ASSERT_TRUE(expanded.is_error());
    EXPECT_EQ(expanded.error().kind, CASErrorKind::Unimplemented);
}

TEST(MeijerGSlaterTest, PGreaterThanQRefused) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr w = parse_expr("w", ctx);
    auto g = make_meijerg(ctx, 1, 0,
        {rat(arena, 1, 3), rat(arena, 1, 4)}, {integer(arena, 0)}, w);
    ASSERT_TRUE(g.is_ok());
    const auto* g_call = expr_cast<FuncCall>(g.value());
    ASSERT_NE(g_call, nullptr);
    auto expanded = slater_expand(ctx, *g_call);
    ASSERT_TRUE(expanded.is_error());
    EXPECT_EQ(expanded.error().kind, CASErrorKind::Unimplemented);
}

// ── §9.4: nodo G senza espansione RESTA (non è un errore) ──────────────────

TEST(MeijerGExpandTest, UnexpandableGStaysIntact) {
    // p=3, q=4 -> pF(q-1) = 3F3: no engine node, Slater refuses, and
    // expand_meijerg_nodes keeps the G node unchanged (same pointer).
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr w = parse_expr("w", ctx);
    auto g = make_meijerg(ctx, 1, 0,
        {rat(arena, 1, 5), rat(arena, 1, 7), rat(arena, 2, 7)},
        {integer(arena, 0), rat(arena, 1, 3), rat(arena, 1, 5),
         rat(arena, 2, 5)},
        w);
    ASSERT_TRUE(g.is_ok());
    auto mapped = expand_meijerg_nodes(ctx, g.value());
    ASSERT_TRUE(mapped.is_ok());
    EXPECT_EQ(mapped.value().get(), g.value().get());
}

TEST(MeijerGExpandTest, MapsInsideCompositeExpressions) {
    // A G node buried in a Sum/Product tree is expanded in place; the
    // surrounding structure survives.
    CASContext ctx;
    ExprPtr f = parse_expr("sin(x)", ctx);
    auto g_form = to_meijerg(ctx, f);
    ASSERT_TRUE(g_form.is_ok());
    ExprPtr wrapped = ctx.arena().make<Sum>(std::vector<ExprPtr>{
        g_form.value(), parse_expr("y^2", ctx)});
    auto back = expand_meijerg_nodes(ctx, wrapped);
    ASSERT_TRUE(back.is_ok());
    auto lhs = ctx.simplify(back.value());
    auto rhs = ctx.simplify(parse_expr("sin(x) + y^2", ctx));
    ASSERT_TRUE(lhs.is_ok() && rhs.is_ok());
    EXPECT_TRUE(ExprEqual{}(lhs.value(), rhs.value()));
}

}  // namespace
}  // namespace cas::symbolic
