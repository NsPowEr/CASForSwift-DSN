// A7 Brick 2 — make_meijerg / view_meijerg (Meijer_G_Slater.md §10.2, §7.2).
// Round-trip elementare<->G (§10.1) and the rest of §10 land in later Bricks
// once to_meijerg/from_meijerg exist; this file covers construction only.

#include "cas/meijerg.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

namespace cas::symbolic {
namespace {

[[nodiscard]] ExprPtr sym(AstArena& arena, const char* name) {
    return arena.make<Symbol>(name);
}
[[nodiscard]] ExprPtr rat(AstArena& arena, long long num, long long den) {
    return arena.make<RationalLit>(BigInt(num), BigInt(den));
}
[[nodiscard]] ExprPtr integer(AstArena& arena, long long v) {
    return arena.make<IntegerLit>(BigInt(v));
}

// ── §10.2 Validazione parametri ─────────────────────────────────────────────

TEST(MeijerGConstructTest, RejectsMGreaterThanQ) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr z = sym(arena, "z");
    // m=2, q=1 (b has 1 element) -> m > q.
    auto result = make_meijerg(ctx, 2, 0, {}, {integer(arena, 0)}, z);
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::InvalidArgument);
}

TEST(MeijerGConstructTest, RejectsNGreaterThanP) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr z = sym(arena, "z");
    // n=2, p=1 (a has 1 element) -> n > p.
    auto result = make_meijerg(ctx, 0, 2, {integer(arena, 0)}, {}, z);
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::InvalidArgument);
}

TEST(MeijerGConstructTest, RejectsPoleOverlap_DiffEqualsOne) {
    // erf-shaped node (m=1,n=1,p=1,q=2): a=[1], b=[a1-1=0, 0] deliberately
    // made to overlap: a[0]-b[0] = 1-0 = 1, a positive integer.
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr z = sym(arena, "z");
    auto result = make_meijerg(ctx, 1, 1, {integer(arena, 1)},
                               {integer(arena, 0), integer(arena, 0)}, z);
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::InvalidArgument);
}

TEST(MeijerGConstructTest, AcceptsPoleDiffEqualsZero) {
    // a[0]-b[0] = 0, not a positive integer -> no overlap, construction ok.
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr z = sym(arena, "z");
    auto result = make_meijerg(ctx, 1, 1, {integer(arena, 1)},
                               {integer(arena, 1), integer(arena, 0)}, z);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
}

TEST(MeijerGConstructTest, RejectsSymbolicOverlap_WhenDecidable) {
    // a[0]-b[0] = (x+2) - x = 2, decidable via simplify -> positive integer
    // overlap even though the inputs are not literals.
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr x = sym(arena, "x");
    ExprPtr a0 = arena.make<Binary>(BinaryOp::Add, x, integer(arena, 2));
    auto result = make_meijerg(ctx, 1, 1, {a0}, {x, integer(arena, 0)}, sym(arena, "z"));
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::InvalidArgument);
}

TEST(MeijerGConstructTest, AcceptsUndecidableSymbolicDiff_NoFalseRejection) {
    // a[0]-b[0] = x - y: not decidable as a literal -> accepted (never a
    // false rejection on symbolic parameters, spec §2.2 "solo se decidibile").
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr x = sym(arena, "x");
    ExprPtr y = sym(arena, "y");
    auto result = make_meijerg(ctx, 1, 1, {x}, {y, integer(arena, 0)}, sym(arena, "z"));
    ASSERT_TRUE(result.is_ok()) << result.error().message;
}

TEST(MeijerGConstructTest, RejectsNullArgument) {
    CASContext ctx;
    auto result = make_meijerg(ctx, 0, 0, {}, {}, ExprPtr{});
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::InvalidArgument);
}

TEST(MeijerGConstructTest, RejectsParamCountAboveBudget) {
    CASContext ctx;
    ctx.set_meijerg_max_param_count(3U);
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> a{integer(arena, 1), integer(arena, 2)};
    std::vector<ExprPtr> b{integer(arena, 3), integer(arena, 4)};  // p+q=4 > 3
    auto result = make_meijerg(ctx, 0, 0, a, b, sym(arena, "z"));
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::Unimplemented);
}

// ── §7.2 Canonicalizzazione (structural sharing) ────────────────────────────

TEST(MeijerGConstructTest, CanonicalizesWithinGroup_PointerEquality) {
    // erf-G: m=1,n=1,p=1,q=2, b sub-group b_{m+1..q} = (0) alone (size 1, no
    // permutation possible there) -- exercise the a_{n+1..p} EMPTY case and
    // b_1..b_m size-1 case are trivial; instead build a case with a genuine
    // 2-element sub-group: n=0 (a_1..a_n empty), a_{n+1..p} = (y, x) vs (x, y).
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr x = sym(arena, "x");
    ExprPtr y = sym(arena, "y");
    ExprPtr z = sym(arena, "z");

    auto g_xy = make_meijerg(ctx, 0, 0, {x, y}, {}, z);
    auto g_yx = make_meijerg(ctx, 0, 0, {y, x}, {}, z);
    ASSERT_TRUE(g_xy.is_ok()) << g_xy.error().message;
    ASSERT_TRUE(g_yx.is_ok()) << g_yx.error().message;

    // Same arena -> structural interning means canonicalized-equal nodes are
    // the SAME pointer (Regola 2), not merely structurally equal.
    EXPECT_EQ(g_xy.value(), g_yx.value());
}

TEST(MeijerGConstructTest, DoesNotReorderBetweenGroups) {
    // n=1 (a_1..a_n has 1 elem = x), a_{n+1..p} has 1 elem = y. Even though
    // canonical_compare(x,y) might otherwise interleave them, the two
    // sub-groups must stay in their own args slots (a_n before a_rest).
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr x = sym(arena, "z_high");  // deliberately sorts after y lexically
    ExprPtr y = sym(arena, "a_low");
    ExprPtr z = sym(arena, "z");

    auto built = make_meijerg(ctx, 0, 1, {x, y}, {}, z);
    ASSERT_TRUE(built.is_ok()) << built.error().message;
    const auto* call = expr_cast<FuncCall>(built.value());
    ASSERT_NE(call, nullptr);
    auto view = view_meijerg(*call);
    ASSERT_TRUE(view.is_ok()) << view.error().message;
    // a_1..a_n (n=1) must be exactly {x}; a_{n+1..p} must be exactly {y}.
    EXPECT_EQ(view.value().a[0], x);
    EXPECT_EQ(view.value().a[1], y);
}

// ── view_meijerg round-trip ──────────────────────────────────────────────

TEST(MeijerGConstructTest, ViewRoundTripsConstruction) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr z2 = arena.make<Binary>(BinaryOp::Pow, sym(arena, "z"), integer(arena, 2));
    // erf-G: m=1,n=1,p=1,q=2, a=(1), b=(1/2,0).
    auto built = make_meijerg(ctx, 1, 1, {integer(arena, 1)},
                              {rat(arena, 1, 2), integer(arena, 0)}, z2);
    ASSERT_TRUE(built.is_ok()) << built.error().message;

    const auto* call = expr_cast<FuncCall>(built.value());
    ASSERT_NE(call, nullptr);
    EXPECT_EQ(call->func_id, BuiltinOp::MeijerG);
    EXPECT_EQ(call->args.size(), 4U + 1U + 2U + 1U);

    auto view = view_meijerg(*call);
    ASSERT_TRUE(view.is_ok()) << view.error().message;
    EXPECT_EQ(view.value().m, 1U);
    EXPECT_EQ(view.value().n, 1U);
    EXPECT_EQ(view.value().p, 1U);
    EXPECT_EQ(view.value().q, 2U);
    ASSERT_EQ(view.value().a.size(), 1U);
    ASSERT_EQ(view.value().b.size(), 2U);
    EXPECT_EQ(view.value().z, z2);
}

TEST(MeijerGConstructTest, ViewRejectsNonMeijerGCall) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr call_expr = arena.make<FuncCall>(BuiltinOp::Sin,
        std::vector<ExprPtr>{sym(arena, "x")});
    const auto* call = expr_cast<FuncCall>(call_expr);
    ASSERT_NE(call, nullptr);
    auto view = view_meijerg(*call);
    ASSERT_TRUE(view.is_error());
    EXPECT_EQ(view.error().kind, CASErrorKind::InternalError);
}

}  // namespace
}  // namespace cas::symbolic
