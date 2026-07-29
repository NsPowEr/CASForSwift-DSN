// A7 step 3 — to_meijerg / pfq_to_meijerg (Meijer_G_Slater.md §5, §3.1;
// test plan §10.1 structural half + §10.5). Each expected form below is a
// DIRECT TRANSCRIPTION of the spec's §5.x / §3.1 line (mpmath-verified in
// Appendice C) — the test checks the implementation against the spec, not
// against itself. Full round-trips (§10.1 complete) need from_meijerg
// (next step).

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

// Recursively locates the (single) MeijerG node inside a prefactor*G tree.
[[nodiscard]] const FuncCall* find_meijerg(ExprPtr e) {
    if (e == nullptr) return nullptr;
    if (const auto* call = expr_cast<FuncCall>(e)) {
        if (call->func_id == BuiltinOp::MeijerG) return call;
        for (ExprPtr arg : call->args)
            if (const auto* found = find_meijerg(arg)) return found;
        return nullptr;
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        for (ExprPtr f : prod->factors)
            if (const auto* found = find_meijerg(f)) return found;
        return nullptr;
    }
    if (const auto* bin = expr_cast<Binary>(e)) {
        if (const auto* found = find_meijerg(bin->left)) return found;
        return find_meijerg(bin->right);
    }
    if (const auto* un = expr_cast<Unary>(e)) return find_meijerg(un->operand);
    return nullptr;
}

// Asserts the G node inside `actual` has the given indices and parameter
// lists (after the factory's within-group canonicalization) and that the
// WHOLE expression equals prefactor*G structurally after simplification.
void expect_g_structure(CASContext& ctx, ExprPtr actual,
                        std::size_t m, std::size_t n,
                        const std::vector<ExprPtr>& a,
                        const std::vector<ExprPtr>& b,
                        ExprPtr z) {
    const auto* g = find_meijerg(actual);
    ASSERT_NE(g, nullptr);
    auto view = view_meijerg(*g);
    ASSERT_TRUE(view.is_ok());
    EXPECT_EQ(view.value().m, m);
    EXPECT_EQ(view.value().n, n);
    ASSERT_EQ(view.value().p, a.size());
    ASSERT_EQ(view.value().q, b.size());
    // Compare against a factory-built twin so both sides carry the same
    // canonical within-group ordering (§7.2).
    auto expected_g = make_meijerg(ctx, m, n, a, b, z);
    ASSERT_TRUE(expected_g.is_ok());
    const auto* eg = expr_cast<FuncCall>(expected_g.value());
    ASSERT_NE(eg, nullptr);
    auto eview = view_meijerg(*eg);
    ASSERT_TRUE(eview.is_ok());
    for (std::size_t i = 0; i < a.size(); ++i)
        EXPECT_TRUE(ExprEqual{}(view.value().a[i], eview.value().a[i])) << "a[" << i << "]";
    for (std::size_t i = 0; i < b.size(); ++i)
        EXPECT_TRUE(ExprEqual{}(view.value().b[i], eview.value().b[i])) << "b[" << i << "]";
    EXPECT_TRUE(ExprEqual{}(view.value().z, eview.value().z)) << "z";
}

// ── §5.1 esponenziale ──────────────────────────────────────────────────────

TEST(MeijerGConvertTest, ExpToG) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr x = parse_expr("x", ctx);
    auto res = to_meijerg(ctx, parse_expr("exp(x)", ctx));
    ASSERT_TRUE(res.is_ok());
    // e^x = G^{1,0}_{0,1}(-x | ; 0)
    auto neg_x = ctx.simplify(arena.make<Unary>(UnaryOp::Neg, x));
    ASSERT_TRUE(neg_x.is_ok());
    expect_g_structure(ctx, res.value(), 1, 0, {}, {integer(arena, 0)},
                       neg_x.value());
}

TEST(MeijerGConvertTest, EPowSameAsExp) {
    CASContext ctx;
    ExprPtr x = parse_expr("x", ctx);
    ExprPtr epow = ctx.arena().make<Binary>(BinaryOp::Pow,
        ctx.arena().make<Constant>(MathConstant::E), x);
    auto via_pow = to_meijerg(ctx, epow);
    auto via_call = to_meijerg(ctx, parse_expr("exp(x)", ctx));
    ASSERT_TRUE(via_pow.is_ok());
    ASSERT_TRUE(via_call.is_ok());
    auto s1 = ctx.simplify(via_pow.value());
    auto s2 = ctx.simplify(via_call.value());
    ASSERT_TRUE(s1.is_ok() && s2.is_ok());
    EXPECT_TRUE(ExprEqual{}(s1.value(), s2.value()));
}

// ── §5.2 / §5.3 trigonometriche e iperboliche ──────────────────────────────

TEST(MeijerGConvertTest, SinToG) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    auto res = to_meijerg(ctx, parse_expr("sin(x)", ctx));
    ASSERT_TRUE(res.is_ok());
    // sin x = sqrt(pi) G^{1,0}_{0,2}(x^2/4 | ; 1/2 ; 0)
    auto z = ctx.simplify(parse_expr("x^2/4", ctx));
    ASSERT_TRUE(z.is_ok());
    expect_g_structure(ctx, res.value(), 1, 0, {},
                       {rat(arena, 1, 2), integer(arena, 0)}, z.value());
}

TEST(MeijerGConvertTest, CosToG) {
    CASContext ctx;
    auto res = to_meijerg(ctx, parse_expr("cos(x)", ctx));
    ASSERT_TRUE(res.is_ok());
    auto z = ctx.simplify(parse_expr("x^2/4", ctx));
    ASSERT_TRUE(z.is_ok());
    // cos: b = (0 | 1/2) — group order semantically rigid vs sin.
    const auto* g = find_meijerg(res.value());
    ASSERT_NE(g, nullptr);
    auto view = view_meijerg(*g);
    ASSERT_TRUE(view.is_ok());
    EXPECT_EQ(view.value().m, 1U);
    ASSERT_EQ(view.value().q, 2U);
    const auto* b1 = expr_cast<IntegerLit>(view.value().b[0]);
    ASSERT_NE(b1, nullptr);  // first (m-group) parameter must be 0, not 1/2
    EXPECT_TRUE(b1->value.is_zero());
    EXPECT_TRUE(ExprEqual{}(view.value().z, z.value()));
}

TEST(MeijerGConvertTest, CoshToG_NegatedArgument) {
    CASContext ctx;
    auto res = to_meijerg(ctx, parse_expr("cosh(x)", ctx));
    ASSERT_TRUE(res.is_ok());
    // cosh x = sqrt(pi) G^{1,0}_{0,2}(-x^2/4 | ; 0 ; 1/2), b1=0 branch-safe.
    auto z = ctx.simplify(parse_expr("-(x^2)/4", ctx));
    ASSERT_TRUE(z.is_ok());
    const auto* g = find_meijerg(res.value());
    ASSERT_NE(g, nullptr);
    auto view = view_meijerg(*g);
    ASSERT_TRUE(view.is_ok());
    const auto* b1 = expr_cast<IntegerLit>(view.value().b[0]);
    ASSERT_NE(b1, nullptr);
    EXPECT_TRUE(b1->value.is_zero());
    EXPECT_TRUE(ExprEqual{}(view.value().z, z.value()));
}

TEST(MeijerGConvertTest, SinhToG_RealPrimaryForm) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    auto res = to_meijerg(ctx, parse_expr("sinh(x)", ctx));
    ASSERT_TRUE(res.is_ok());
    // sinh x = (sqrt(pi)/2) x G^{1,0}_{0,2}(-x^2/4 | ; 0 ; -1/2) — §5.3
    // primary REAL form (the -i*sqrt(pi) variant is forbidden, §11).
    auto z = ctx.simplify(parse_expr("-(x^2)/4", ctx));
    ASSERT_TRUE(z.is_ok());
    expect_g_structure(ctx, res.value(), 1, 0, {},
                       {integer(arena, 0), rat(arena, -1, 2)}, z.value());
    // The x prefactor must be present (result != bare G).
    EXPECT_EQ(expr_cast<FuncCall>(res.value()), nullptr);
}

// ── §5.4 logaritmo ─────────────────────────────────────────────────────────

TEST(MeijerGConvertTest, LnOnePlusXToG) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    auto res = to_meijerg(ctx, parse_expr("ln(1 + x)", ctx));
    ASSERT_TRUE(res.is_ok());
    // ln(1+z) = G^{1,2}_{2,2}(z | 1,1 ; 1 ; 0): z = (1+x)-1 = x.
    ExprPtr x = parse_expr("x", ctx);
    expect_g_structure(ctx, res.value(), 1, 2,
                       {integer(arena, 1), integer(arena, 1)},
                       {integer(arena, 1), integer(arena, 0)}, x);
}

// ── §5.5 potenze ───────────────────────────────────────────────────────────

TEST(MeijerGConvertTest, BinomialPowerToG) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    auto res = to_meijerg(ctx, parse_expr("(1 + x)^(-1/3)", ctx));
    ASSERT_TRUE(res.is_ok());
    // (1+z)^{-a} = (1/Gamma(a)) G^{1,1}_{1,1}(z | 1-a ; 0), a = 1/3:
    // 1-a = 2/3, z = x.
    ExprPtr x = parse_expr("x", ctx);
    expect_g_structure(ctx, res.value(), 1, 1, {rat(arena, 2, 3)},
                       {integer(arena, 0)}, x);
}

TEST(MeijerGConvertTest, PolynomialPowerRefused) {
    CASContext ctx;
    auto res = to_meijerg(ctx, parse_expr("(1 + x)^2", ctx));
    ASSERT_TRUE(res.is_error());  // structured Unimplemented, never a wrong G
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
}

TEST(MeijerGConvertTest, PowExpProductToG) {
    CASContext ctx;
    auto input = ctx.simplify(parse_expr("x^(1/3) * exp(-x)", ctx));
    ASSERT_TRUE(input.is_ok());
    auto res = to_meijerg(ctx, input.value());
    ASSERT_TRUE(res.is_ok());
    // z^a e^{-z} = G^{1,0}_{0,1}(z | ; a): bare G, no prefactor.
    const auto* g = expr_cast<FuncCall>(res.value());
    ASSERT_NE(g, nullptr);
    EXPECT_EQ(g->func_id, BuiltinOp::MeijerG);
    auto view = view_meijerg(*g);
    ASSERT_TRUE(view.is_ok());
    EXPECT_EQ(view.value().m, 1U);
    EXPECT_EQ(view.value().n, 0U);
    EXPECT_EQ(view.value().p, 0U);
    ASSERT_EQ(view.value().q, 1U);
    EXPECT_TRUE(ExprEqual{}(view.value().b[0],
                            ctx.simplify(parse_expr("1/3", ctx)).value()));
}

// ── §5.6 / §5.7 / §5.8 ─────────────────────────────────────────────────────

TEST(MeijerGConvertTest, AtanToG_HalfPrefactor) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    auto res = to_meijerg(ctx, parse_expr("arctan(x)", ctx));
    ASSERT_TRUE(res.is_ok());
    // arctan z = (1/2) G^{1,2}_{2,2}(z^2 | 1/2,1 ; 1/2 ; 0) — prefactor 1/2,
    // NOT z/2 (errata Appendice D.2).
    auto z = ctx.simplify(parse_expr("x^2", ctx));
    ASSERT_TRUE(z.is_ok());
    expect_g_structure(ctx, res.value(), 1, 2,
                       {rat(arena, 1, 2), integer(arena, 1)},
                       {rat(arena, 1, 2), integer(arena, 0)}, z.value());
}

TEST(MeijerGConvertTest, AsinToG) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    auto res = to_meijerg(ctx, parse_expr("arcsin(x)", ctx));
    ASSERT_TRUE(res.is_ok());
    auto z = ctx.simplify(parse_expr("-(x^2)", ctx));
    ASSERT_TRUE(z.is_ok());
    expect_g_structure(ctx, res.value(), 1, 2,
                       {rat(arena, 1, 2), rat(arena, 1, 2)},
                       {integer(arena, 0), rat(arena, -1, 2)}, z.value());
}

TEST(MeijerGConvertTest, ErfToG) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    auto res = to_meijerg(ctx, parse_expr("erf(x)", ctx));
    ASSERT_TRUE(res.is_ok());
    auto z = ctx.simplify(parse_expr("x^2", ctx));
    ASSERT_TRUE(z.is_ok());
    expect_g_structure(ctx, res.value(), 1, 1, {integer(arena, 1)},
                       {rat(arena, 1, 2), integer(arena, 0)}, z.value());
}

TEST(MeijerGConvertTest, BesselJToG) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    auto res = to_meijerg(ctx, parse_expr("BesselJ(1/3, x)", ctx));
    ASSERT_TRUE(res.is_ok());
    // J_nu = G^{1,0}_{0,2}(x^2/4 | ; nu/2 ; -nu/2), nu = 1/3.
    auto z = ctx.simplify(parse_expr("x^2/4", ctx));
    ASSERT_TRUE(z.is_ok());
    expect_g_structure(ctx, res.value(), 1, 0, {},
                       {rat(arena, 1, 6), rat(arena, -1, 6)}, z.value());
}

// ── §10.5 ponte pFq (§3.1) ─────────────────────────────────────────────────

TEST(MeijerGConvertTest, Hyp2F1Bridge) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    auto res = to_meijerg(ctx, parse_expr("hyp2F1(1/3, 3/4, 7/5, x)", ctx));
    ASSERT_TRUE(res.is_ok());
    // 2F1 -> G^{1,2}_{2,3}(-x | 1-1/3, 1-3/4 ; 0, 1-7/5):
    // a = (2/3, 1/4), b = (0 | -2/5), m=1, n=p=2.
    auto z = ctx.simplify(parse_expr("-x", ctx));
    ASSERT_TRUE(z.is_ok());
    expect_g_structure(ctx, res.value(), 1, 2,
                       {rat(arena, 2, 3), rat(arena, 1, 4)},
                       {integer(arena, 0), rat(arena, -2, 5)}, z.value());
}

TEST(MeijerGConvertTest, Hyp1F1Bridge) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    auto res = to_meijerg(ctx, parse_expr("hyp1F1(1/3, 7/5, x)", ctx));
    ASSERT_TRUE(res.is_ok());
    auto z = ctx.simplify(parse_expr("-x", ctx));
    ASSERT_TRUE(z.is_ok());
    expect_g_structure(ctx, res.value(), 1, 1, {rat(arena, 2, 3)},
                       {integer(arena, 0), rat(arena, -2, 5)}, z.value());
}

TEST(MeijerGConvertTest, Hyp0F1Bridge) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    auto res = to_meijerg(ctx, parse_expr("hyp0F1(3/2, x)", ctx));
    ASSERT_TRUE(res.is_ok());
    auto z = ctx.simplify(parse_expr("-x", ctx));
    ASSERT_TRUE(z.is_ok());
    expect_g_structure(ctx, res.value(), 1, 0, {},
                       {integer(arena, 0), rat(arena, -1, 2)}, z.value());
}

TEST(MeijerGConvertTest, PfqEmptyListsEqualsExp) {
    // 0F0(;;x) = e^x: the general bridge must agree with the §5.1 table
    // entry (cross-check general path vs fast path, CLAUDE.md cat. 8).
    CASContext ctx;
    ExprPtr x = parse_expr("x", ctx);
    auto via_bridge = pfq_to_meijerg(ctx, {}, {}, x);
    auto via_table = to_meijerg(ctx, parse_expr("exp(x)", ctx));
    ASSERT_TRUE(via_bridge.is_ok());
    ASSERT_TRUE(via_table.is_ok());
    auto s1 = ctx.simplify(via_bridge.value());
    auto s2 = ctx.simplify(via_table.value());
    ASSERT_TRUE(s1.is_ok() && s2.is_ok());
    EXPECT_TRUE(ExprEqual{}(s1.value(), s2.value()));
}

TEST(MeijerGConvertTest, PfqRefusesNonpositiveIntegerUpper) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr x = parse_expr("x", ctx);
    auto res = pfq_to_meijerg(ctx, {integer(arena, -2)},
                              {rat(arena, 7, 5)}, x);
    ASSERT_TRUE(res.is_error());  // polynomial pFq: bridge inapplicable
    EXPECT_EQ(res.error().kind, CASErrorKind::InvalidArgument);
}

TEST(MeijerGConvertTest, PfqRefusesPGreaterThanQPlusOne) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr x = parse_expr("x", ctx);
    auto res = pfq_to_meijerg(ctx,
        {rat(arena, 1, 3), rat(arena, 1, 4), rat(arena, 1, 5)},
        {rat(arena, 7, 5)}, x);
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::InvalidArgument);
}

TEST(MeijerGConvertTest, UnsupportedShapeIsStructuredUnimplemented) {
    CASContext ctx;
    auto res = to_meijerg(ctx, parse_expr("tan(x)", ctx));
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
}

}  // namespace
}  // namespace cas::symbolic
