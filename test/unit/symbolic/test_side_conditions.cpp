// A31 fase 1 — side-conditions accumulated by simplify() (test plan §6 of
// Domain_Conditions_Propagation.md). simplify()'s ExprPtr output must stay
// bit-per-bit identical to pre-A31 behaviour; these tests verify that AND
// the new registration channel (last_side_conditions()).

#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/side_conditions.hpp"
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

[[nodiscard]] bool has_condition(const SideConditionSet& set, DomainConditionKind kind,
                                 ExprPtr subject) {
    return set.contains(DomainCondition{kind, subject});
}

[[nodiscard]] bool is_expr_one(ExprPtr expr) {
    const auto* il = expr_cast<IntegerLit>(expr);
    return il != nullptr && il->value == BigInt(1);
}

// ── §6.1 Emission per regola (siti verificati §4.3.1-2) ────────────────────

TEST(SideConditionsTest, XOverXProduct_EmitsNonZero_ViaSimplifyProductFactors) {
    // Product-syntax x*x^-1 bypasses the rewrite_provider Div path (Product
    // nodes never call try_rewrite) and hits simplify_product_factors's
    // exponent-zero drop directly (spec §4.3.1, simplify_arithmetic_chain.cpp).
    CASContext ctx;
    Symbol x{"x"};
    ExprPtr xexpr = ctx.arena().make<Symbol>(x.name);
    ExprPtr inv = ctx.arena().make<Binary>(
        BinaryOp::Pow, xexpr, ctx.arena().make<IntegerLit>(BigInt(-1)));
    ExprPtr prod = ctx.arena().make<Product>(std::vector<ExprPtr>{xexpr, inv});

    auto result = ctx.simplify(prod);
    ASSERT_TRUE(result.is_ok());
    EXPECT_TRUE(is_expr_one(result.value()));  // output unchanged vs pre-A31

    const auto& conds = ctx.last_side_conditions();
    EXPECT_TRUE(has_condition(conds, DomainConditionKind::NonZero, xexpr));
}

TEST(SideConditionsTest, XOverXDiv_EmitsNonZero_ViaGcdRewriteProvider) {
    // Div-syntax x/x reaches the rewrite_provider FIRST (Binary::Div path in
    // simplify_arithmetic.cpp calls try_rewrite before its own case-label),
    // and try_rewrite_algebraic's polynomial-GCD cancellation collapses it
    // (spec §4.3.2, Meccanismo 2 of §1.1 -- a DIFFERENT code path than the
    // Product test above).
    CASContext ctx;
    ExprPtr expr = parse_expr("x / x", ctx);

    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    EXPECT_TRUE(is_expr_one(result.value()));

    Symbol x{"x"};
    ExprPtr xexpr = ctx.arena().make<Symbol>(x.name);
    EXPECT_TRUE(has_condition(ctx.last_side_conditions(),
                              DomainConditionKind::NonZero, xexpr));
}

TEST(SideConditionsTest, AlreadyProvenAssumption_EmitsNothing) {
    // §3.3: if ctx.assumptions() already proves the predicate, it is a FACT
    // of the input, not an assumption the rewrite introduced -- no
    // registration. Controprova of the two tests above.
    CASContext ctx;
    Symbol x{"x"};
    ctx.assumptions().assume_nonzero(x);
    ExprPtr expr = parse_expr("x / x", ctx);

    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    EXPECT_TRUE(is_expr_one(result.value()));
    EXPECT_TRUE(ctx.last_side_conditions().empty());
}

TEST(SideConditionsTest, NoDropNoEmission) {
    // A product with no cancellation emits nothing (sanity: the accumulator
    // isn't polluted by unrelated simplify() traffic).
    CASContext ctx;
    ExprPtr expr = parse_expr("x * y + 1", ctx);
    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    EXPECT_TRUE(ctx.last_side_conditions().empty());
}

// ── §6.2 Dedup / subsumption ────────────────────────────────────────────────

TEST(SideConditionsTest, ExactDuplicate_Dedups) {
    SideConditionSet set;
    Symbol x{"x"};
    AstArena arena;
    ExprPtr xexpr = arena.make<Symbol>(x.name);
    DomainCondition c{DomainConditionKind::NonZero, xexpr};
    EXPECT_TRUE(set.add(c));
    EXPECT_TRUE(set.add(c));
    EXPECT_EQ(set.size(), 1U);
}

TEST(SideConditionsTest, PositiveSubsumesNonZeroAndNonNegative) {
    SideConditionSet set;
    AstArena arena;
    ExprPtr xexpr = arena.make<Symbol>("x");
    EXPECT_TRUE(set.add({DomainConditionKind::NonZero, xexpr}));
    EXPECT_TRUE(set.add({DomainConditionKind::NonNegative, xexpr}));
    EXPECT_EQ(set.size(), 2U);
    EXPECT_TRUE(set.add({DomainConditionKind::Positive, xexpr}));
    // Positive subsumes both weaker conditions on the same subject.
    EXPECT_EQ(set.size(), 1U);
    EXPECT_EQ(set.items().front().kind, DomainConditionKind::Positive);

    // Adding NonZero/NonNegative AFTER Positive is a no-op (already subsumed).
    EXPECT_TRUE(set.add({DomainConditionKind::NonZero, xexpr}));
    EXPECT_EQ(set.size(), 1U);
}

TEST(SideConditionsTest, DeterministicOrder) {
    AstArena arena1, arena2;
    ExprPtr x1 = arena1.make<Symbol>("x");
    ExprPtr y1 = arena1.make<Symbol>("y");
    ExprPtr x2 = arena2.make<Symbol>("x");
    ExprPtr y2 = arena2.make<Symbol>("y");

    SideConditionSet a, b;
    a.add({DomainConditionKind::NonZero, x1});
    a.add({DomainConditionKind::Positive, y1});
    b.add({DomainConditionKind::NonZero, x2});
    b.add({DomainConditionKind::Positive, y2});

    ASSERT_EQ(a.items().size(), b.items().size());
    for (std::size_t i = 0; i < a.items().size(); ++i) {
        EXPECT_EQ(a.items()[i].kind, b.items()[i].kind);
    }
}

// ── §6.3 Cache ───────────────────────────────────────────────────────────

TEST(SideConditionsTest, CacheHit_ReemitsConditions) {
    CASContext ctx;
    Symbol x{"x"};
    ExprPtr xexpr = ctx.arena().make<Symbol>(x.name);
    ExprPtr inv = ctx.arena().make<Binary>(
        BinaryOp::Pow, xexpr, ctx.arena().make<IntegerLit>(BigInt(-1)));
    ExprPtr prod = ctx.arena().make<Product>(std::vector<ExprPtr>{xexpr, inv});

    auto first = ctx.simplify(prod);
    ASSERT_TRUE(first.is_ok());
    EXPECT_TRUE(has_condition(ctx.last_side_conditions(),
                              DomainConditionKind::NonZero, xexpr));

    // A second, DIFFERENT top-level simplify() call that hits the cache for
    // `prod` must still report the condition -- not silently drop it.
    ExprPtr unrelated = ctx.arena().make<IntegerLit>(BigInt(42));
    auto spacer = ctx.simplify(unrelated);
    ASSERT_TRUE(spacer.is_ok());
    EXPECT_TRUE(ctx.last_side_conditions().empty());  // unrelated call: nothing.

    auto second = ctx.simplify(prod);
    ASSERT_TRUE(second.is_ok());
    EXPECT_TRUE(has_condition(ctx.last_side_conditions(),
                              DomainConditionKind::NonZero, xexpr));
}

TEST(SideConditionsTest, AssumptionChange_InvalidatesCacheAndReconditions) {
    CASContext ctx;
    Symbol x{"x"};
    ExprPtr xexpr = ctx.arena().make<Symbol>(x.name);
    ExprPtr inv = ctx.arena().make<Binary>(
        BinaryOp::Pow, xexpr, ctx.arena().make<IntegerLit>(BigInt(-1)));
    ExprPtr prod = ctx.arena().make<Product>(std::vector<ExprPtr>{xexpr, inv});

    auto first = ctx.simplify(prod);
    ASSERT_TRUE(first.is_ok());
    EXPECT_FALSE(ctx.last_side_conditions().empty());

    // F7.0-A4.1: an assumption change bumps Assumptions::revision(), which
    // CASContext::simplify() uses to invalidate simplify_cache_ -- so the
    // NEXT identical call recomputes (and now proves NonZero, so nothing is
    // registered) rather than replaying a stale cached condition.
    ctx.assumptions().assume_nonzero(x);
    auto second = ctx.simplify(prod);
    ASSERT_TRUE(second.is_ok());
    EXPECT_TRUE(ctx.last_side_conditions().empty());
}

TEST(SideConditionsTest, ConsecutiveTopLevelCalls_DoNotCrossPollute) {
    CASContext ctx;
    Symbol x{"x"}, y{"y"};
    ExprPtr xexpr = ctx.arena().make<Symbol>(x.name);
    ExprPtr yexpr = ctx.arena().make<Symbol>(y.name);
    ExprPtr xinv = ctx.arena().make<Binary>(
        BinaryOp::Pow, xexpr, ctx.arena().make<IntegerLit>(BigInt(-1)));
    ExprPtr xprod = ctx.arena().make<Product>(std::vector<ExprPtr>{xexpr, xinv});

    auto r1 = ctx.simplify(xprod);
    ASSERT_TRUE(r1.is_ok());
    EXPECT_TRUE(has_condition(ctx.last_side_conditions(),
                              DomainConditionKind::NonZero, xexpr));

    auto r2 = ctx.simplify(yexpr);  // unrelated: y has no cancellation.
    ASSERT_TRUE(r2.is_ok());
    EXPECT_TRUE(ctx.last_side_conditions().empty());
    EXPECT_FALSE(has_condition(ctx.last_side_conditions(),
                               DomainConditionKind::NonZero, xexpr));
}

// ── §6.5 Meccanismo 1 vs Meccanismo 2 (§1.1) ────────────────────────────────

TEST(SideConditionsTest, DivAndProductSyntax_AgreeOnOutputAndCondition) {
    CASContext ctx1, ctx2;
    Symbol x{"x"};

    ExprPtr div_expr = parse_expr("x / x", ctx1);
    auto div_result = ctx1.simplify(div_expr);
    ASSERT_TRUE(div_result.is_ok());

    ExprPtr xexpr2 = ctx2.arena().make<Symbol>(x.name);
    ExprPtr inv2 = ctx2.arena().make<Binary>(
        BinaryOp::Pow, xexpr2, ctx2.arena().make<IntegerLit>(BigInt(-1)));
    ExprPtr prod_expr = ctx2.arena().make<Product>(std::vector<ExprPtr>{xexpr2, inv2});
    auto prod_result = ctx2.simplify(prod_expr);
    ASSERT_TRUE(prod_result.is_ok());

    EXPECT_TRUE(is_expr_one(div_result.value()));
    EXPECT_TRUE(is_expr_one(prod_result.value()));
    EXPECT_FALSE(ctx1.last_side_conditions().empty());
    EXPECT_FALSE(ctx2.last_side_conditions().empty());
    EXPECT_EQ(ctx1.last_side_conditions().items().front().kind,
              ctx2.last_side_conditions().items().front().kind);
}

TEST(SideConditionsTest, SumOfFractions_DegenerateToNMinusD_EmitsNonZero) {
    // Meccanismo 2 (together()/GCD path) on an aggregated rational function:
    // (x^2 - 1) / (x - 1) collapses to x + 1 via polynomial GCD, no
    // extra condition (denominator degree 1, doesn't hit the =1 special
    // case) -- but (x - 1)/(x - 1) alone (embedded via multiplication) does.
    CASContext ctx;
    ExprPtr expr = parse_expr("(x - 1) / (x - 1)", ctx);
    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    EXPECT_TRUE(is_expr_one(result.value()));
    EXPECT_FALSE(ctx.last_side_conditions().empty());
    EXPECT_EQ(ctx.last_side_conditions().items().front().kind,
              DomainConditionKind::NonZero);
}

// ── GC re-interning (collect_garbage must not leave dangling subjects) ─────

TEST(SideConditionsTest, CollectGarbage_ReinternLiveAndCachedConditions) {
    CASContext ctx;
    Symbol x{"x"};
    ExprPtr xexpr = ctx.arena().make<Symbol>(x.name);
    ExprPtr inv = ctx.arena().make<Binary>(
        BinaryOp::Pow, xexpr, ctx.arena().make<IntegerLit>(BigInt(-1)));
    ExprPtr prod = ctx.arena().make<Product>(std::vector<ExprPtr>{xexpr, inv});

    auto result = ctx.simplify(prod);
    ASSERT_TRUE(result.is_ok());
    ASSERT_FALSE(ctx.last_side_conditions().empty());
    ExprPtr subject_before = ctx.last_side_conditions().items().front().subject;

    ctx.collect_garbage();

    // The live accumulator survives GC with a re-interned (different
    // pointer, same structure) subject -- not a dangling reference into the
    // freed arena.
    ASSERT_FALSE(ctx.last_side_conditions().empty());
    ExprPtr subject_after = ctx.last_side_conditions().items().front().subject;
    EXPECT_NE(subject_after.get(), subject_before.get());
    EXPECT_TRUE(expr_is<Symbol>(subject_after));
    EXPECT_EQ(expr_cast<Symbol>(subject_after)->name, "x");

    // The cached entry (behind simplify_cache_) survives GC too: re-running
    // the identical expression must still be able to serve a cache hit and
    // re-report the condition without touching freed memory (ASan-verified).
    ExprPtr prod2 = ctx.arena().make<Product>(std::vector<ExprPtr>{
        ctx.arena().make<Symbol>(x.name),
        ctx.arena().make<Binary>(BinaryOp::Pow, ctx.arena().make<Symbol>(x.name),
                                 ctx.arena().make<IntegerLit>(BigInt(-1)))});
    auto result2 = ctx.simplify(prod2);
    ASSERT_TRUE(result2.is_ok());
    EXPECT_TRUE(is_expr_one(result2.value()));
    EXPECT_FALSE(ctx.last_side_conditions().empty());
}

}  // namespace
}  // namespace cas::symbolic
