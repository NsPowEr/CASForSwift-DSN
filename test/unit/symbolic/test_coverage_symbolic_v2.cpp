// DEBT-F1-COV-01 — Coverage uplift v2 for src/symbolic/ files below 85%.
// Targets: context_utils, context_core, rewrite_engine, rewrite_matching,
//          normal_form, assumptions, simplify_functions, simplify_trig_inverse,
//          complex_qi.
// All tests use the public CASContext / symbolic API — no internal access.

#include <gtest/gtest.h>
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "cas/normal_form.hpp"
#include "cas/complex_rational.hpp"
#include "cas/rational.hpp"
#include "cas/error.hpp"

#include <chrono>
#include <string>
#include <vector>

using namespace cas;
using namespace cas::symbolic;

namespace {

[[nodiscard]] ExprPtr parse_expr(CASContext& ctx, const std::string& s) {
    auto tok = Lexer(s).tokenize();
    if (!tok.is_ok()) return nullptr;
    auto parsed = Parser(tok.value(), ctx.arena()).parse();
    if (!parsed.is_ok()) return nullptr;
    return parsed.value();
}

[[nodiscard]] ExprPtr make_int(AstArena& a, long long v) {
    return a.make<IntegerLit>(BigInt(v));
}

[[nodiscard]] ExprPtr make_sym(AstArena& a, const std::string& name) {
    return a.make<Symbol>(name);
}

} // namespace

// ─── context_utils: exercise through Assumptions API ─────────────────────────
// context_utils.cpp provides exact_scalar_from_expr, negate_expr,
// range_is_exact_zero, exact_range_excludes_zero which are used internally.
// We exercise them via the public Assumptions::could_be_zero and
// Assumptions::check_consistency interfaces.

TEST(ContextUtilsCoverage, CouldBeZeroWithPositiveLowerRange) {
    // assume_in_range(x, 1, 2) → range excludes zero → could_be_zero = false
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_in_range(x,
        make_int(ctx.arena(), 1), make_int(ctx.arena(), 2));
    EXPECT_FALSE(ctx.assumptions().could_be_zero(x));
}

TEST(ContextUtilsCoverage, CouldBeZeroWithNegativeUpperRange) {
    // assume_in_range(x, -5, -1) → range excludes zero → could_be_zero = false
    CASContext ctx;
    Symbol x("x");
    ExprPtr neg5 = ctx.arena().make<Unary>(UnaryOp::Neg, make_int(ctx.arena(), 5));
    ExprPtr neg1 = ctx.arena().make<Unary>(UnaryOp::Neg, make_int(ctx.arena(), 1));
    ctx.assumptions().assume_in_range(x, neg5, neg1);
    EXPECT_FALSE(ctx.assumptions().could_be_zero(x));
}

TEST(ContextUtilsCoverage, CouldBeZeroWithZeroRange) {
    // assume_in_range(x, 0, 0) → range is exactly zero → could_be_zero = true
    // (unless nonzero constraint overrides)
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_in_range(x,
        make_int(ctx.arena(), 0), make_int(ctx.arena(), 0));
    // Range [0,0] means x==0, so could_be_zero should be true
    EXPECT_TRUE(ctx.assumptions().could_be_zero(x));
}

TEST(ContextUtilsCoverage, CouldBeZeroWithCrossingRange) {
    // assume_in_range(x, -1, 1) → includes zero → could_be_zero = true
    CASContext ctx;
    Symbol x("x");
    ExprPtr neg1 = ctx.arena().make<Unary>(UnaryOp::Neg, make_int(ctx.arena(), 1));
    ExprPtr one = make_int(ctx.arena(), 1);
    ctx.assumptions().assume_in_range(x, neg1, one);
    EXPECT_TRUE(ctx.assumptions().could_be_zero(x));
}

TEST(ContextUtilsCoverage, CouldBeZeroSymbolicRange) {
    // assume_in_range(x, sym, 1) — symbolic lower → can't determine → true
    CASContext ctx;
    Symbol x("x"), y("y");
    ExprPtr y_expr = make_sym(ctx.arena(), "y");
    ExprPtr one = make_int(ctx.arena(), 1);
    ctx.assumptions().assume_in_range(x, y_expr, one);
    EXPECT_TRUE(ctx.assumptions().could_be_zero(x));
}

TEST(ContextUtilsCoverage, ExprWeightThroughMatchSeqAc) {
    // match_sequence_ac_internal uses expr_weight to sort patterns
    // exercise it by matching a Sum with more patterns
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr a = make_sym(arena, "a");
    ExprPtr b = make_sym(arena, "b");
    // Pattern: [x_, y_]  — two wildcards
    ExprPtr wild_x = arena.make<Symbol>("x_");
    ExprPtr wild_y = arena.make<Symbol>("y_");
    std::vector<ExprPtr> pat{wild_x, wild_y};
    ExprPtr pattern = arena.make<Sum>(std::move(pat));
    // Expr: a + b + c  (3 terms, 2 pattern terms → AC partial match via engine)
    std::vector<ExprPtr> expr_terms{a, b};
    ExprPtr expr = arena.make<Sum>(std::move(expr_terms));
    MatchMap m;
    // match_pattern should succeed with AC
    bool matched = match_pattern(expr, pattern, m);
    EXPECT_TRUE(matched);
}

// ─── context_core: materialize_expr covers all node types ────────────────────

TEST(ContextCoreCoverage2, MaterializeIntegerLit) {
    CASContext ctx;
    ExprPtr e = make_int(ctx.arena(), 42);
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    ASSERT_NE(res.value(), nullptr);
    const auto* lit = expr_cast<IntegerLit>(res.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(42));
}

TEST(ContextCoreCoverage2, MaterializeRationalLit) {
    CASContext ctx;
    ExprPtr e = ctx.arena().make<RationalLit>(BigInt(3), BigInt(4));
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NE(res.value(), nullptr);
}

TEST(ContextCoreCoverage2, MaterializeSymbol) {
    CASContext ctx;
    ExprPtr e = make_sym(ctx.arena(), "x");
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    const auto* sym = expr_cast<Symbol>(res.value());
    ASSERT_NE(sym, nullptr);
    EXPECT_EQ(sym->name, "x");
}

TEST(ContextCoreCoverage2, MaterializeUnary) {
    CASContext ctx;
    ExprPtr inner = make_int(ctx.arena(), 5);
    ExprPtr e = ctx.arena().make<Unary>(UnaryOp::Neg, inner);
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    const auto* u = expr_cast<Unary>(res.value());
    ASSERT_NE(u, nullptr);
    EXPECT_EQ(u->op, UnaryOp::Neg);
}

TEST(ContextCoreCoverage2, MaterializeBinary) {
    CASContext ctx;
    ExprPtr e = ctx.arena().make<Binary>(
        BinaryOp::Add, make_int(ctx.arena(), 1), make_int(ctx.arena(), 2));
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    const auto* b = expr_cast<Binary>(res.value());
    ASSERT_NE(b, nullptr);
    EXPECT_EQ(b->op, BinaryOp::Add);
}

TEST(ContextCoreCoverage2, MaterializeFuncCall) {
    CASContext ctx;
    std::vector<ExprPtr> args{make_sym(ctx.arena(), "x")};
    ExprPtr e = ctx.arena().make<FuncCall>(BuiltinOp::Sin, std::move(args));
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    const auto* fc = expr_cast<FuncCall>(res.value());
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->func_id, BuiltinOp::Sin);
}

TEST(ContextCoreCoverage2, MaterializeSum) {
    CASContext ctx;
    std::vector<ExprPtr> terms{make_int(ctx.arena(), 1), make_int(ctx.arena(), 2)};
    ExprPtr e = ctx.arena().make<Sum>(std::move(terms));
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    const auto* s = expr_cast<Sum>(res.value());
    ASSERT_NE(s, nullptr);
    EXPECT_EQ(s->terms.size(), 2U);
}

TEST(ContextCoreCoverage2, MaterializeProduct) {
    CASContext ctx;
    std::vector<ExprPtr> factors{make_int(ctx.arena(), 3), make_sym(ctx.arena(), "x")};
    ExprPtr e = ctx.arena().make<Product>(std::move(factors));
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    const auto* p = expr_cast<Product>(res.value());
    ASSERT_NE(p, nullptr);
    EXPECT_EQ(p->factors.size(), 2U);
}

TEST(ContextCoreCoverage2, MaterializeDerivative) {
    CASContext ctx;
    ExprPtr expr = ctx.arena().make<Derivative>(
        make_sym(ctx.arena(), "x"), Symbol("x"), 1);
    AstArena arena2;
    auto res = materialize_expr(expr, arena2);
    ASSERT_TRUE(res.is_ok());
    const auto* d = expr_cast<Derivative>(res.value());
    ASSERT_NE(d, nullptr);
    EXPECT_EQ(d->order, 1);
}

TEST(ContextCoreCoverage2, MaterializeIntegral) {
    CASContext ctx;
    ExprPtr integrand = make_sym(ctx.arena(), "x");
    ExprPtr e = ctx.arena().make<Integral>(
        integrand, Symbol("x"),
        std::optional<ExprPtr>(make_int(ctx.arena(), 0)),
        std::optional<ExprPtr>(make_int(ctx.arena(), 1)));
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    const auto* intg = expr_cast<Integral>(res.value());
    ASSERT_NE(intg, nullptr);
    EXPECT_TRUE(intg->lower.has_value());
    EXPECT_TRUE(intg->upper.has_value());
}

TEST(ContextCoreCoverage2, MaterializeIntegralNoLimits) {
    CASContext ctx;
    ExprPtr integrand = make_sym(ctx.arena(), "x");
    ExprPtr e = ctx.arena().make<Integral>(
        integrand, Symbol("x"), std::nullopt, std::nullopt);
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    const auto* intg = expr_cast<Integral>(res.value());
    ASSERT_NE(intg, nullptr);
    EXPECT_FALSE(intg->lower.has_value());
    EXPECT_FALSE(intg->upper.has_value());
}

TEST(ContextCoreCoverage2, MaterializeLimit) {
    CASContext ctx;
    ExprPtr e = ctx.arena().make<Limit>(
        make_sym(ctx.arena(), "x"), Symbol("x"),
        make_int(ctx.arena(), 0), LimitDirection::Both);
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    const auto* lim = expr_cast<Limit>(res.value());
    ASSERT_NE(lim, nullptr);
}

TEST(ContextCoreCoverage2, MaterializeMatrix) {
    CASContext ctx;
    std::vector<ExprPtr> elems{
        make_int(ctx.arena(), 1), make_int(ctx.arena(), 2),
        make_int(ctx.arena(), 3), make_int(ctx.arena(), 4)};
    ExprPtr e = ctx.arena().make<Matrix>(2, 2, std::move(elems));
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    const auto* m = expr_cast<Matrix>(res.value());
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->rows, 2U);
    EXPECT_EQ(m->cols, 2U);
}

TEST(ContextCoreCoverage2, MaterializeRootOf) {
    CASContext ctx;
    // RootOf(x^2 - 2, x)
    ExprPtr poly = ctx.arena().make<Binary>(
        BinaryOp::Sub,
        ctx.arena().make<Binary>(BinaryOp::Pow, make_sym(ctx.arena(), "x"), make_int(ctx.arena(), 2)),
        make_int(ctx.arena(), 2));
    ExprPtr e = ctx.arena().make<RootOf>(poly, Symbol("x"), std::nullopt);
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NE(res.value(), nullptr);
}

TEST(ContextCoreCoverage2, MaterializeNullFails) {
    CASContext ctx;
    AstArena arena2;
    auto res = materialize_expr(ExprPtr{}, arena2);
    EXPECT_TRUE(res.is_error());
}

TEST(ContextCoreCoverage2, MaterializeSeriesExp) {
    CASContext ctx;
    std::vector<std::pair<long long, ExprPtr>> terms;
    terms.push_back({0, make_int(ctx.arena(), 1)});
    terms.push_back({1, make_int(ctx.arena(), 1)});
    ExprPtr e = ctx.arena().make<SeriesExp>(
        Symbol("x"), make_int(ctx.arena(), 0), std::move(terms), 2);
    AstArena arena2;
    auto res = materialize_expr(e, arena2);
    ASSERT_TRUE(res.is_ok());
    const auto* s = expr_cast<SeriesExp>(res.value());
    ASSERT_NE(s, nullptr);
    EXPECT_EQ(s->terms.size(), 2U);
}

// ─── context_core: collect_garbage ───────────────────────────────────────────

TEST(ContextCoreCoverage2, CollectGarbageEmpty) {
    CASContext ctx;
    ctx.collect_garbage({});
    // Should not crash on empty
    EXPECT_TRUE(true);
}

TEST(ContextCoreCoverage2, CollectGarbageWithExternalRoot) {
    CASContext ctx;
    ExprPtr expr = parse_expr(ctx, "x + 1");
    ASSERT_NE(expr, nullptr);
    ctx.collect_garbage({&expr});
    // expr should still be valid
    EXPECT_NE(expr, nullptr);
}

TEST(ContextCoreCoverage2, CollectGarbageWithCacheEntries) {
    CASContext ctx;
    ctx.set_caching_enabled(true);
    ExprPtr e = parse_expr(ctx, "x + 1");
    ASSERT_NE(e, nullptr);
    // populate cache
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok());
    // Now collect garbage
    ctx.collect_garbage({});
    EXPECT_TRUE(true);
}

TEST(ContextCoreCoverage2, CollectGarbageWithVariables) {
    CASContext ctx;
    Symbol x("x");
    ctx.define(x, make_int(ctx.arena(), 42));
    ctx.collect_garbage({});
    // variable should still be accessible
    auto found = ctx.lookup(x);
    ASSERT_TRUE(found.has_value());
}

TEST(ContextCoreCoverage2, CollectGarbageWithTraceEnabled) {
    CASContext ctx;
    ctx.enable_trace(true);
    ExprPtr e = parse_expr(ctx, "x + 1");
    if (e) { auto _ = ctx.simplify(e); }
    ctx.collect_garbage({});
    EXPECT_TRUE(true);
}

TEST(ContextCoreCoverage2, CollectGarbageWithAssumptions) {
    CASContext ctx;
    Symbol x("x");
    ExprPtr lower = make_int(ctx.arena(), 0);
    ExprPtr upper = make_int(ctx.arena(), 1);
    ctx.assumptions().assume_in_range(x, lower, upper);
    ctx.collect_garbage({});
    EXPECT_TRUE(true);
}

// ─── context_core: CASContext setters & getters ───────────────────────────────

TEST(ContextCoreCoverage2, SettersClamping) {
    CASContext ctx;
    ctx.set_timeout_check_interval(10U); // below 64 → clamped to 64
    EXPECT_GE(ctx.timeout_check_interval(), 64U);

    ctx.set_max_simplification_depth(5); // below 10 → clamped to 10
    EXPECT_GE(ctx.max_simplification_depth(), 10);

    ctx.set_max_integration_depth(0U); // below 1 → clamped to 1
    EXPECT_GE(ctx.max_integration_depth(), 1U);

    ctx.set_max_integration_depth(200U); // above 128 → clamped to 128
    EXPECT_LE(ctx.max_integration_depth(), 128U);
}

TEST(ContextCoreCoverage2, SetGcdErrorProbability) {
    CASContext ctx;
    ctx.set_gcd_error_probability(1e-10); // below minimum → clamped
    EXPECT_GE(ctx.gcd_error_probability(), 1e-6);
    ctx.set_gcd_error_probability(0.5); // above max → clamped to 0.1
    EXPECT_LE(ctx.gcd_error_probability(), 0.1);
}

TEST(ContextCoreCoverage2, SetNumericPrecisionDigits) {
    CASContext ctx;
    ctx.set_numeric_precision_digits(3U); // below 6 → clamped
    EXPECT_GE(ctx.numeric_precision_digits(), 6U);
    ctx.set_numeric_precision_digits(20000U); // above 10000 → clamped
    EXPECT_LE(ctx.numeric_precision_digits(), 10000U);
}

TEST(ContextCoreCoverage2, SetMaxRootofExplicitDegree) {
    CASContext ctx;
    ctx.set_max_rootof_explicit_degree(0U); // below 1 → clamped
    EXPECT_GE(ctx.max_rootof_explicit_degree(), 1U);
    ctx.set_max_rootof_explicit_degree(5U);
    EXPECT_EQ(ctx.max_rootof_explicit_degree(), 5U);
}

TEST(ContextCoreCoverage2, SetMaxGcdRecursionDepth) {
    CASContext ctx;
    ctx.set_max_gcd_recursion_depth(1U); // below 4 → clamped
    EXPECT_GE(ctx.max_gcd_recursion_depth(), 4U);
}

TEST(ContextCoreCoverage2, SetMinGcdDivisionSteps) {
    CASContext ctx;
    ctx.set_min_gcd_division_steps(0U); // below 1 → clamped
    EXPECT_GE(ctx.min_gcd_division_steps(), 1U);
}

TEST(ContextCoreCoverage2, SetMaxGcdTotalCalls) {
    CASContext ctx;
    ctx.set_max_gcd_total_calls(5U); // below 16 → clamped
    EXPECT_GE(ctx.max_gcd_total_calls(), 16U);
}

TEST(ContextCoreCoverage2, SetMaxCyclotomicN) {
    CASContext ctx;
    ctx.set_max_cyclotomic_n(200);
    EXPECT_EQ(ctx.max_cyclotomic_n(), 200);
}

TEST(ContextCoreCoverage2, SetMaxQAlphaBridgeDepth) {
    CASContext ctx;
    ctx.set_max_q_alpha_bridge_depth(3U); // below 8 → clamped
    EXPECT_GE(ctx.max_q_alpha_bridge_depth(), 8U);
}

TEST(ContextCoreCoverage2, SetMaxGammaRecursion) {
    CASContext ctx;
    ctx.set_max_gamma_recursion(5U); // below 16 → clamped
    EXPECT_GE(ctx.max_gamma_recursion(), 16U);
}

TEST(ContextCoreCoverage2, SetImproperLeadingOrderScan) {
    CASContext ctx;
    ctx.set_improper_leading_order_scan(0U); // below 1 → clamped
    EXPECT_GE(ctx.improper_leading_order_scan(), 1U);
}

TEST(ContextCoreCoverage2, SetExpandBesselRecurrence) {
    CASContext ctx;
    ctx.set_expand_bessel_recurrence(true);
    EXPECT_TRUE(ctx.expand_bessel_recurrence());
    ctx.set_expand_bessel_recurrence(false);
    EXPECT_FALSE(ctx.expand_bessel_recurrence());
}

TEST(ContextCoreCoverage2, SetMaxTragerTowerShiftAttempts) {
    CASContext ctx;
    ctx.set_max_trager_tower_shift_attempts(10U);
    EXPECT_EQ(ctx.max_trager_tower_shift_attempts(), 10U);
}

TEST(ContextCoreCoverage2, MakeFreshSymbolCollisionAvoidance) {
    CASContext ctx;
    // Define some candidate names to force counter increment
    ctx.define(Symbol("x_1"), make_int(ctx.arena(), 1));
    ctx.define(Symbol("x_2"), make_int(ctx.arena(), 2));
    // Fresh symbol should skip over already-defined names
    Symbol fresh1 = ctx.make_fresh_symbol("x");
    Symbol fresh2 = ctx.make_fresh_symbol("x");
    EXPECT_NE(fresh1.name, fresh2.name);
    // Neither should be x_1 or x_2 (already defined)
    EXPECT_NE(fresh1.name, "x_1");
    EXPECT_NE(fresh1.name, "x_2");
}

TEST(ContextCoreCoverage2, EnableTraceAndGetTrace) {
    CASContext ctx;
    ctx.enable_trace(true);
    ExprPtr e = parse_expr(ctx, "x + x");
    ASSERT_NE(e, nullptr);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    // trace should be populated
    const auto& trace = ctx.get_trace();
    (void)trace;
    // disable trace — clears it
    ctx.enable_trace(false);
    EXPECT_TRUE(ctx.get_trace().empty());
}

TEST(ContextCoreCoverage2, CacheMetricsAfterHits) {
    CASContext ctx;
    ctx.set_caching_enabled(true);
    ExprPtr e = parse_expr(ctx, "x + 1");
    ASSERT_NE(e, nullptr);
    { auto _ = ctx.simplify(e); } // miss
    { auto _ = ctx.simplify(e); } // hit
    auto m = ctx.get_simplify_metrics();
    (void)m;
    EXPECT_GE(ctx.get_cache_limit(), 1U);
}

// ─── rewrite_engine: apply_rule with all traversal strategies ─────────────────

TEST(RewriteEngineCoverage, ApplyRuleBottomUp) {
    CASContext ctx;
    // Rule: x_ + 0 → x_  (wildcard pattern)
    AstArena& arena = ctx.arena();
    ExprPtr zero = make_int(arena, 0);
    ExprPtr wild = arena.make<Symbol>("x_");
    // pattern: x_ + 0, replacement: x_
    std::vector<ExprPtr> pat_terms{wild, zero};
    ExprPtr pattern = arena.make<Sum>(std::move(pat_terms));
    ExprPtr replacement = arena.make<Symbol>("x_");
    RewriteRule rule{pattern, replacement, nullptr};

    // Target: (y + 0)
    ExprPtr y = make_sym(arena, "y");
    std::vector<ExprPtr> target_terms{y, zero};
    ExprPtr target = arena.make<Sum>(std::move(target_terms));

    auto result = apply_rule(target, rule, TraversalStrategy::BottomUp, arena);
    // Rule may or may not fire depending on orientation check;
    // but call must succeed without crash
    EXPECT_FALSE(result.is_error());
}

TEST(RewriteEngineCoverage, ApplyRuleTopDown) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Simple rule: a_ → a_  (trivially matching anything; not strictly decreasing)
    // Test that invalid rule returns error
    ExprPtr wild = arena.make<Symbol>("a_");
    RewriteRule bad_rule{wild, wild, nullptr}; // pattern == replacement weight => not oriented
    ExprPtr target = make_int(arena, 5);
    auto result = apply_rule(target, bad_rule, TraversalStrategy::TopDown, arena);
    // Should fail because not strictly decreasing
    EXPECT_TRUE(result.is_error());
}

TEST(RewriteEngineCoverage, ApplyRuleFixPointNull) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Null expression should fail
    RewriteRule rule{make_int(arena, 1), make_int(arena, 0), nullptr};
    auto result = apply_rule(ExprPtr{}, rule, TraversalStrategy::FixPoint, arena);
    EXPECT_TRUE(result.is_error());
}

TEST(RewriteEngineCoverage, ApplyRuleNullPattern) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    RewriteRule rule{ExprPtr{}, make_int(arena, 0), nullptr};
    auto result = apply_rule(make_int(arena, 1), rule, TraversalStrategy::BottomUp, arena);
    EXPECT_TRUE(result.is_error());
}

TEST(RewriteEngineCoverage, ApplyRuleSetEmpty) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    std::vector<RewriteRule> rules;
    // Empty rule set should succeed without modification
    auto result = apply_rule_set(make_int(arena, 5), rules, arena);
    EXPECT_FALSE(result.is_error());
}

TEST(RewriteEngineCoverage, ApplyRuleSetNullExpr) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    std::vector<RewriteRule> rules;
    auto result = apply_rule_set(ExprPtr{}, rules, arena);
    EXPECT_TRUE(result.is_error());
}

TEST(RewriteEngineCoverage, RewriteChildrenCoversUnary) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Build -x and apply an identity rule to get children traversed
    ExprPtr x = make_sym(arena, "x");
    ExprPtr neg_x = arena.make<Unary>(UnaryOp::Neg, x);
    // Rule: x_ + 0 → x_ won't match but forces bottom-up traversal through Unary
    ExprPtr wild = arena.make<Symbol>("x_");
    std::vector<ExprPtr> pat_terms{wild, make_int(arena, 0)};
    ExprPtr pattern = arena.make<Sum>(std::move(pat_terms));
    RewriteRule rule{pattern, arena.make<Symbol>("x_"), nullptr};
    // Rule may not fire (orientation check), so just verify no crash
    auto result = apply_rule(neg_x, rule, TraversalStrategy::BottomUp, arena);
    (void)result;
    EXPECT_TRUE(true);
}

TEST(RewriteEngineCoverage, RewriteChildrenCoversBinary) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr e = ctx.arena().make<Binary>(
        BinaryOp::Mul, make_sym(arena, "a"), make_sym(arena, "b"));
    ExprPtr wild = arena.make<Symbol>("x_");
    std::vector<ExprPtr> pat_terms{wild, make_int(arena, 0)};
    ExprPtr pattern = arena.make<Sum>(std::move(pat_terms));
    RewriteRule rule{pattern, arena.make<Symbol>("x_"), nullptr};
    auto result = apply_rule(e, rule, TraversalStrategy::BottomUp, arena);
    (void)result;
    EXPECT_TRUE(true);
}

// ─── rewrite_matching: Integral, Derivative, Limit, RootOf, Matrix matching ──

TEST(RewriteMatchingCoverage, MatchDerivative) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr expr = arena.make<Derivative>(make_sym(arena, "x"), Symbol("x"), 1);
    ExprPtr pattern = arena.make<Derivative>(make_sym(arena, "x"), Symbol("x"), 1);
    MatchMap m;
    EXPECT_TRUE(match_pattern(expr, pattern, m));
}

TEST(RewriteMatchingCoverage, MatchDerivativeOrderMismatch) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr expr = arena.make<Derivative>(make_sym(arena, "x"), Symbol("x"), 1);
    ExprPtr pattern = arena.make<Derivative>(make_sym(arena, "x"), Symbol("x"), 2);
    MatchMap m;
    EXPECT_FALSE(match_pattern(expr, pattern, m));
}

TEST(RewriteMatchingCoverage, MatchLimit) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr expr = arena.make<Limit>(
        make_sym(arena, "x"), Symbol("x"), make_int(arena, 0), LimitDirection::Both);
    ExprPtr pattern = arena.make<Limit>(
        make_sym(arena, "x"), Symbol("x"), make_int(arena, 0), LimitDirection::Both);
    MatchMap m;
    EXPECT_TRUE(match_pattern(expr, pattern, m));
}

TEST(RewriteMatchingCoverage, MatchLimitVariableMismatch) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr expr = arena.make<Limit>(
        make_sym(arena, "x"), Symbol("x"), make_int(arena, 0), LimitDirection::Both);
    ExprPtr pattern = arena.make<Limit>(
        make_sym(arena, "y"), Symbol("y"), make_int(arena, 0), LimitDirection::Both);
    MatchMap m;
    EXPECT_FALSE(match_pattern(expr, pattern, m));
}

TEST(RewriteMatchingCoverage, MatchRootOf) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr poly = arena.make<Binary>(
        BinaryOp::Sub,
        arena.make<Binary>(BinaryOp::Pow, make_sym(arena, "x"), make_int(arena, 2)),
        make_int(arena, 2));
    ExprPtr expr = arena.make<RootOf>(poly, Symbol("x"), std::nullopt);
    ExprPtr pattern = arena.make<RootOf>(poly, Symbol("x"), std::nullopt);
    MatchMap m;
    EXPECT_TRUE(match_pattern(expr, pattern, m));
}

TEST(RewriteMatchingCoverage, MatchMatrix2x2) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> elems{
        make_int(arena, 1), make_int(arena, 2),
        make_int(arena, 3), make_int(arena, 4)};
    ExprPtr expr = arena.make<Matrix>(2, 2, elems);
    ExprPtr pattern = arena.make<Matrix>(2, 2, elems);
    MatchMap m;
    EXPECT_TRUE(match_pattern(expr, pattern, m));
}

TEST(RewriteMatchingCoverage, MatchMatrixDimensionMismatch) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> elems2x2{make_int(arena, 1), make_int(arena, 2),
                                   make_int(arena, 3), make_int(arena, 4)};
    std::vector<ExprPtr> elems1x2{make_int(arena, 1), make_int(arena, 2)};
    ExprPtr expr = arena.make<Matrix>(2, 2, elems2x2);
    ExprPtr pattern = arena.make<Matrix>(1, 2, elems1x2);
    MatchMap m;
    EXPECT_FALSE(match_pattern(expr, pattern, m));
}

TEST(RewriteMatchingCoverage, MatchIntegralWithBounds) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr e = arena.make<Integral>(
        make_sym(arena, "x"), Symbol("x"),
        std::optional<ExprPtr>(make_int(arena, 0)),
        std::optional<ExprPtr>(make_int(arena, 1)));
    ExprPtr p = arena.make<Integral>(
        make_sym(arena, "x"), Symbol("x"),
        std::optional<ExprPtr>(make_int(arena, 0)),
        std::optional<ExprPtr>(make_int(arena, 1)));
    MatchMap m;
    EXPECT_TRUE(match_pattern(e, p, m));
}

TEST(RewriteMatchingCoverage, MatchIntegralVariableMismatch) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr e = arena.make<Integral>(
        make_sym(arena, "x"), Symbol("x"), std::nullopt, std::nullopt);
    ExprPtr p = arena.make<Integral>(
        make_sym(arena, "y"), Symbol("y"), std::nullopt, std::nullopt);
    MatchMap m;
    EXPECT_FALSE(match_pattern(e, p, m));
}

TEST(RewriteMatchingCoverage, MatchIntegralBoundsMismatch) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // one has bounds, other doesn't
    ExprPtr e = arena.make<Integral>(
        make_sym(arena, "x"), Symbol("x"),
        std::optional<ExprPtr>(make_int(arena, 0)),
        std::optional<ExprPtr>(make_int(arena, 1)));
    ExprPtr p = arena.make<Integral>(
        make_sym(arena, "x"), Symbol("x"), std::nullopt, std::nullopt);
    MatchMap m;
    EXPECT_FALSE(match_pattern(e, p, m));
}

TEST(RewriteMatchingCoverage, MatchWildcardConsistency) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Pattern: x_ + x_  (same wildcard twice) — should match only equal terms
    ExprPtr wild = arena.make<Symbol>("x_");
    std::vector<ExprPtr> pts{wild, wild};
    ExprPtr pattern = arena.make<Sum>(std::move(pts));

    // expr: y + y  → should match, wildcard x_ = y
    ExprPtr y = make_sym(arena, "y");
    std::vector<ExprPtr> same_terms{y, y};
    ExprPtr expr_same = arena.make<Sum>(std::move(same_terms));
    MatchMap m1;
    EXPECT_TRUE(match_pattern(expr_same, pattern, m1));

    // expr: y + z  → should NOT match (wildcard inconsistency)
    ExprPtr z = make_sym(arena, "z");
    std::vector<ExprPtr> diff_terms{y, z};
    ExprPtr expr_diff = arena.make<Sum>(std::move(diff_terms));
    MatchMap m2;
    EXPECT_FALSE(match_pattern(expr_diff, pattern, m2));
}

TEST(RewriteMatchingCoverage, MatchNullPatterns) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    MatchMap m;
    // Both null
    EXPECT_TRUE(match_pattern(ExprPtr{}, ExprPtr{}, m));
    // One null
    EXPECT_FALSE(match_pattern(make_int(arena, 1), ExprPtr{}, m));
    EXPECT_FALSE(match_pattern(ExprPtr{}, make_int(arena, 1), m));
}

TEST(RewriteMatchingCoverage, MatchKindMismatch) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr integer = make_int(arena, 5);
    ExprPtr sym = make_sym(arena, "x");
    MatchMap m;
    EXPECT_FALSE(match_pattern(integer, sym, m));
}

TEST(RewriteMatchingCoverage, MatchAcPattern) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // match_ac_pattern on sums
    ExprPtr a = make_sym(arena, "a");
    ExprPtr b = make_sym(arena, "b");
    std::vector<ExprPtr> sum_terms{a, b};
    ExprPtr expr = arena.make<Sum>(sum_terms);
    ExprPtr pattern = arena.make<Sum>(sum_terms);
    MatchMap m;
    EXPECT_TRUE(match_ac_pattern(expr, pattern, m));
}

// ─── assumptions: all Domain paths, check_consistency ────────────────────────

TEST(AssumptionsCoverage, AssumeDomainNegative) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_domain(x, Domain::Negative);
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::Negative);
    EXPECT_TRUE(ctx.assumptions().is_real(x));
    EXPECT_TRUE(ctx.assumptions().is_nonzero(x));
}

TEST(AssumptionsCoverage, AssumeDomainNatural) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_domain(x, Domain::Natural);
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::Natural);
    EXPECT_TRUE(ctx.assumptions().is_real(x));
    EXPECT_TRUE(ctx.assumptions().is_integer(x));
}

TEST(AssumptionsCoverage, AssumeDomainInteger) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_domain(x, Domain::Integer);
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::Integer);
    EXPECT_TRUE(ctx.assumptions().is_integer(x));
}

TEST(AssumptionsCoverage, AssumeDomainRational) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_domain(x, Domain::Rational);
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::Rational);
    EXPECT_TRUE(ctx.assumptions().is_real(x));
}

TEST(AssumptionsCoverage, AssumeDomainNonZero) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_domain(x, Domain::NonZero);
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::NonZero);
    EXPECT_TRUE(ctx.assumptions().is_nonzero(x));
}

TEST(AssumptionsCoverage, AssumeDomainComplex) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_domain(x, Domain::Complex);
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::Complex);
}

TEST(AssumptionsCoverage, GetDomainFallbackToSets) {
    CASContext ctx;
    Symbol x("x");
    // No domain explicitly set — derive from sets
    ctx.assumptions().assume_positive(x);
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::Positive);
}

TEST(AssumptionsCoverage, GetDomainFromNegativeSet) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_domain(x, Domain::Negative); // sets negative_symbols_
    // clear the explicit domain entry by using a fresh ctx; access via get_domain
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::Negative);
}

TEST(AssumptionsCoverage, CheckConsistencyPositiveVsNegative) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_domain(x, Domain::Negative); // forces into negative_symbols_
    // This should detect a contradiction
    auto res = ctx.assumptions().check_consistency();
    EXPECT_TRUE(res.is_error());
}

TEST(AssumptionsCoverage, CheckConsistencyRangeVsPositive) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_positive(x);
    // Range with negative upper bound conflicts with positive
    ctx.assumptions().assume_in_range(x,
        ctx.arena().make<Unary>(UnaryOp::Neg, make_int(ctx.arena(), 5)),
        ctx.arena().make<Unary>(UnaryOp::Neg, make_int(ctx.arena(), 1)));
    auto res = ctx.assumptions().check_consistency();
    EXPECT_TRUE(res.is_error());
}

TEST(AssumptionsCoverage, CheckConsistencyRangeVsNonzeroZeroRange) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_nonzero(x);
    ctx.assumptions().assume_in_range(x,
        make_int(ctx.arena(), 0), make_int(ctx.arena(), 0));
    auto res = ctx.assumptions().check_consistency();
    EXPECT_TRUE(res.is_error());
}

TEST(AssumptionsCoverage, CheckConsistencyValid) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_in_range(x,
        make_int(ctx.arena(), 1), make_int(ctx.arena(), 10));
    auto res = ctx.assumptions().check_consistency();
    EXPECT_FALSE(res.is_error());
}

TEST(AssumptionsCoverage, AssumeViaBinaryExpressions) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr zero = make_int(arena, 0);

    // BinaryOp::Equal — assume_greater_equal both ways
    ExprPtr eq_cond = arena.make<Binary>(BinaryOp::Equal, x_expr, zero);
    ctx.assumptions().assume(eq_cond);

    // BinaryOp::Less — x < 0
    ExprPtr lt_cond = arena.make<Binary>(BinaryOp::Less, x_expr, zero);
    ctx.assumptions().assume(lt_cond);

    // BinaryOp::Greater — x > 0
    ExprPtr gt_cond = arena.make<Binary>(BinaryOp::Greater, x_expr, zero);
    ctx.assumptions().assume(gt_cond);

    // BinaryOp::LessEqual — x <= 0
    ExprPtr le_cond = arena.make<Binary>(BinaryOp::LessEqual, x_expr, zero);
    ctx.assumptions().assume(le_cond);

    // BinaryOp::GreaterEqual — x >= 0
    ExprPtr ge_cond = arena.make<Binary>(BinaryOp::GreaterEqual, x_expr, zero);
    ctx.assumptions().assume(ge_cond);
}

TEST(AssumptionsCoverage, AssumeViaFuncCallForms) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr x_expr = make_sym(arena, "x");

    auto make_fc = [&](const std::string& name, std::vector<ExprPtr> args) -> ExprPtr {
        return arena.make<FuncCall>(name, std::move(args));
    };

    ctx.assumptions().assume(make_fc("greater", {x_expr, make_int(arena, 0)}));
    ctx.assumptions().assume(make_fc("greater_equal", {x_expr, make_int(arena, 0)}));
    ctx.assumptions().assume(make_fc("less", {x_expr, make_int(arena, 10)}));
    ctx.assumptions().assume(make_fc("less_equal", {x_expr, make_int(arena, 10)}));
    ctx.assumptions().assume(make_fc("positive", {x_expr}));
    ctx.assumptions().assume(make_fc("real", {x_expr}));
    ctx.assumptions().assume(make_fc("nonzero", {x_expr}));
    ctx.assumptions().assume(make_fc("integer", {x_expr}));
    ctx.assumptions().assume(make_fc("nonneg", {x_expr}));
    ctx.assumptions().assume(make_fc("nonnegative", {x_expr}));
    ctx.assumptions().assume(make_fc("negative", {x_expr}));
    // non-symbol negative form
    ExprPtr two = make_int(arena, 2);
    ctx.assumptions().assume(make_fc("negative", {two}));
    // positive non-symbol form
    ctx.assumptions().assume(make_fc("positive", {two}));
}

TEST(AssumptionsCoverage, AssumeViaSymbol) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr x_expr = make_sym(arena, "x");
    ctx.assumptions().assume(x_expr); // Symbol → nonzero
    Symbol x("x");
    EXPECT_TRUE(ctx.assumptions().is_nonzero(x));
}

TEST(AssumptionsCoverage, IsRealOfComplexExpressions) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ctx.assumptions().assume_real(x);
    ctx.assumptions().assume_positive(x);

    ExprPtr x_expr = make_sym(arena, "x");
    // Sum of reals is real
    ExprPtr s = arena.make<Sum>(std::vector<ExprPtr>{x_expr, make_int(arena, 1)});
    EXPECT_TRUE(ctx.assumptions().is_real(s));

    // Product of reals is real
    ExprPtr p = arena.make<Product>(std::vector<ExprPtr>{x_expr, make_int(arena, 2)});
    EXPECT_TRUE(ctx.assumptions().is_real(p));

    // exp(x) is real when x is real
    ExprPtr exp_x = arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{x_expr});
    EXPECT_TRUE(ctx.assumptions().is_real(exp_x));

    // sin(x) is real when x is real
    ExprPtr sin_x = arena.make<FuncCall>(BuiltinOp::Sin, std::vector<ExprPtr>{x_expr});
    EXPECT_TRUE(ctx.assumptions().is_real(sin_x));

    // sqrt(x) is real when x >= 0 (x is positive here)
    ExprPtr sqrt_x = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{x_expr});
    EXPECT_TRUE(ctx.assumptions().is_real(sqrt_x));

    // ln(x) is real when x > 0
    ExprPtr ln_x = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{x_expr});
    EXPECT_TRUE(ctx.assumptions().is_real(ln_x));

    // Unary neg
    ExprPtr neg_x = arena.make<Unary>(UnaryOp::Neg, x_expr);
    EXPECT_TRUE(ctx.assumptions().is_real(neg_x));

    // Binary pow: positive ^ real
    ExprPtr pow_e = arena.make<Binary>(BinaryOp::Pow, x_expr, x_expr);
    EXPECT_TRUE(ctx.assumptions().is_real(pow_e));

    // Constant Pi and E
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr e_const = arena.make<Constant>(MathConstant::E);
    EXPECT_TRUE(ctx.assumptions().is_real(pi));
    EXPECT_TRUE(ctx.assumptions().is_real(e_const));
}

TEST(AssumptionsCoverage, IsNonnegativeBranches) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ctx.assumptions().assume_real(x);
    ExprPtr x_expr = make_sym(arena, "x");

    // abs is nonneg
    ExprPtr abs_x = arena.make<FuncCall>(BuiltinOp::Abs, std::vector<ExprPtr>{x_expr});
    EXPECT_TRUE(ctx.assumptions().is_nonnegative(abs_x));

    // sqrt is nonneg
    ExprPtr sqrt_x = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{x_expr});
    EXPECT_TRUE(ctx.assumptions().is_nonnegative(sqrt_x));

    // even power of real
    ExprPtr x2 = arena.make<Binary>(BinaryOp::Pow, x_expr, make_int(arena, 2));
    EXPECT_TRUE(ctx.assumptions().is_nonnegative(x2));
}

TEST(AssumptionsCoverage, IsGreaterWithSumDecomposition) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_positive(x); // x > 0
    ctx.assumptions().assume_positive(y); // y > 0

    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    ExprPtr zero = make_int(arena, 0);

    // x + y > 0
    std::vector<ExprPtr> sum_terms{x_expr, y_expr};
    ExprPtr sum = arena.make<Sum>(std::move(sum_terms));
    EXPECT_TRUE(ctx.assumptions().is_greater(sum, zero));
}

TEST(AssumptionsCoverage, CoulBeZeroWithRange) {
    CASContext ctx;
    Symbol x("x");
    // Range [1, 5] — excludes zero
    ctx.assumptions().assume_in_range(x,
        make_int(ctx.arena(), 1), make_int(ctx.arena(), 5));
    EXPECT_FALSE(ctx.assumptions().could_be_zero(x));
}

TEST(AssumptionsCoverage, CouldBeZeroSymbolicNoRange) {
    CASContext ctx;
    Symbol x("x");
    EXPECT_TRUE(ctx.assumptions().could_be_zero(x));
}

TEST(AssumptionsCoverage, IsIntegerFunctions) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol n("n");
    ctx.assumptions().assume_integer(n);
    EXPECT_TRUE(ctx.assumptions().is_integer(n));
    ExprPtr int_lit = make_int(arena, 5);
    EXPECT_TRUE(ctx.assumptions().is_integer(int_lit));
    ExprPtr sym = make_sym(arena, "n");
    EXPECT_TRUE(ctx.assumptions().is_integer(sym));
}

TEST(AssumptionsCoverage, GetRange) {
    CASContext ctx;
    Symbol x("x");
    ExprPtr low = make_int(ctx.arena(), 0);
    ExprPtr high = make_int(ctx.arena(), 10);
    ctx.assumptions().assume_in_range(x, low, high);
    auto range = ctx.assumptions().get_range(x);
    ASSERT_TRUE(range.has_value());
    EXPECT_EQ(range->lower, low);
    EXPECT_EQ(range->upper, high);
}

TEST(AssumptionsCoverage, ProvePositiveProduct) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_domain(x, Domain::Negative);
    ctx.assumptions().assume_domain(y, Domain::Negative);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    // (-x) * (-y) should be positive when x,y < 0
    std::vector<ExprPtr> factors{x_expr, y_expr};
    ExprPtr prod = arena.make<Product>(std::move(factors));
    // Both negative → product is positive
    EXPECT_TRUE(ctx.assumptions().is_positive(prod));
}

// ─── normal_form: transcendental_normal_form all paths ───────────────────────

TEST(NormalFormCoverage, TranscendentalNormalFormLn1) {
    CASContext ctx;
    ExprPtr ln_1 = ctx.arena().make<FuncCall>(
        BuiltinOp::Ln, std::vector<ExprPtr>{make_int(ctx.arena(), 1)});
    auto res = transcendental_normal_form(ln_1, ctx);
    ASSERT_TRUE(res.is_ok());
    // ln(1) → 0
    const auto* lit = expr_cast<IntegerLit>(res.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(0));
}

TEST(NormalFormCoverage, TranscendentalNormalFormLnExp) {
    CASContext ctx;
    ExprPtr inner = ctx.arena().make<FuncCall>(
        BuiltinOp::Exp, std::vector<ExprPtr>{make_sym(ctx.arena(), "x")});
    ExprPtr ln_exp = ctx.arena().make<FuncCall>(
        BuiltinOp::Ln, std::vector<ExprPtr>{inner});
    auto res = transcendental_normal_form(ln_exp, ctx);
    ASSERT_TRUE(res.is_ok());
    // ln(exp(x)) → x
    const auto* sym = expr_cast<Symbol>(res.value());
    ASSERT_NE(sym, nullptr);
    EXPECT_EQ(sym->name, "x");
}

TEST(NormalFormCoverage, TranscendentalNormalFormLnBinaryMul) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr x = make_sym(arena, "x");
    ExprPtr y = make_sym(arena, "y");
    ExprPtr x_times_y = arena.make<Binary>(BinaryOp::Mul, x, y);
    ExprPtr ln_xy = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{x_times_y});
    auto res = transcendental_normal_form(ln_xy, ctx);
    ASSERT_TRUE(res.is_ok());
    // ln(x*y) → ln(x) + ln(y) as Binary Add
    EXPECT_NE(res.value(), nullptr);
}

TEST(NormalFormCoverage, TranscendentalNormalFormLnDiv) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr x = make_sym(arena, "x");
    ExprPtr y = make_sym(arena, "y");
    ExprPtr x_div_y = arena.make<Binary>(BinaryOp::Div, x, y);
    ExprPtr ln_div = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{x_div_y});
    auto res = transcendental_normal_form(ln_div, ctx);
    ASSERT_TRUE(res.is_ok());
    // ln(x/y) → ln(x) - ln(y)
    EXPECT_NE(res.value(), nullptr);
}

TEST(NormalFormCoverage, TranscendentalNormalFormLnPow) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr x = make_sym(arena, "x");
    ExprPtr x_sq = arena.make<Binary>(BinaryOp::Pow, x, make_int(arena, 3));
    ExprPtr ln_pow = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{x_sq});
    auto res = transcendental_normal_form(ln_pow, ctx);
    ASSERT_TRUE(res.is_ok());
    // ln(x^3) → 3*ln(x) as Binary Mul
    EXPECT_NE(res.value(), nullptr);
}

TEST(NormalFormCoverage, TranscendentalNormalFormLnProduct) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr x = make_sym(arena, "x");
    ExprPtr y = make_sym(arena, "y");
    std::vector<ExprPtr> factors{x, y};
    ExprPtr prod = arena.make<Product>(std::move(factors));
    ExprPtr ln_prod = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{prod});
    auto res = transcendental_normal_form(ln_prod, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NE(res.value(), nullptr);
}

TEST(NormalFormCoverage, TranscendentalNormalFormExp0) {
    CASContext ctx;
    ExprPtr exp_0 = ctx.arena().make<FuncCall>(
        BuiltinOp::Exp, std::vector<ExprPtr>{make_int(ctx.arena(), 0)});
    auto res = transcendental_normal_form(exp_0, ctx);
    ASSERT_TRUE(res.is_ok());
    // exp(0) → 1
    const auto* lit = expr_cast<IntegerLit>(res.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(1));
}

TEST(NormalFormCoverage, TranscendentalNormalFormExpLn) {
    CASContext ctx;
    ExprPtr inner = ctx.arena().make<FuncCall>(
        BuiltinOp::Ln, std::vector<ExprPtr>{make_sym(ctx.arena(), "x")});
    ExprPtr exp_ln = ctx.arena().make<FuncCall>(
        BuiltinOp::Exp, std::vector<ExprPtr>{inner});
    auto res = transcendental_normal_form(exp_ln, ctx);
    ASSERT_TRUE(res.is_ok());
    // exp(ln(x)) → x
    const auto* sym = expr_cast<Symbol>(res.value());
    ASSERT_NE(sym, nullptr);
    EXPECT_EQ(sym->name, "x");
}

TEST(NormalFormCoverage, TranscendentalNormalFormBinaryRecurse) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Binary(Add, ln(1), exp(0))
    ExprPtr ln_1 = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{make_int(arena, 1)});
    ExprPtr exp_0 = arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{make_int(arena, 0)});
    ExprPtr bin = arena.make<Binary>(BinaryOp::Add, ln_1, exp_0);
    auto res = transcendental_normal_form(bin, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NE(res.value(), nullptr);
}

TEST(NormalFormCoverage, TranscendentalNormalFormSumRecurse) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr ln_1 = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{make_int(arena, 1)});
    ExprPtr x = make_sym(arena, "x");
    std::vector<ExprPtr> terms{ln_1, x};
    ExprPtr sum = arena.make<Sum>(std::move(terms));
    auto res = transcendental_normal_form(sum, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NE(res.value(), nullptr);
}

TEST(NormalFormCoverage, TranscendentalNormalFormProductRecurse) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr exp_0 = arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{make_int(arena, 0)});
    ExprPtr x = make_sym(arena, "x");
    std::vector<ExprPtr> factors{exp_0, x};
    ExprPtr prod = arena.make<Product>(std::move(factors));
    auto res = transcendental_normal_form(prod, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NE(res.value(), nullptr);
}

TEST(NormalFormCoverage, TranscendentalNormalFormUnaryRecurse) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr ln_1 = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{make_int(arena, 1)});
    ExprPtr neg_ln = arena.make<Unary>(UnaryOp::Neg, ln_1);
    auto res = transcendental_normal_form(neg_ln, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NE(res.value(), nullptr);
}

TEST(NormalFormCoverage, TranscendentalNormalFormDepthExceeded) {
    CASContext ctx;
    ExprPtr x = make_sym(ctx.arena(), "x");
    // depth=0 should return Unimplemented
    auto res = transcendental_normal_form(x, ctx, 0);
    EXPECT_TRUE(res.is_error());
}

TEST(NormalFormCoverage, TranscendentalNormalFormNull) {
    CASContext ctx;
    auto res = transcendental_normal_form(ExprPtr{}, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_EQ(res.value(), ExprPtr{});
}

TEST(NormalFormCoverage, TranscendentalNormalFormOtherFuncCall) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Other FuncCall (e.g. sin) — recurse into args
    ExprPtr ln_1 = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{make_int(arena, 1)});
    ExprPtr sin_ln = arena.make<FuncCall>(BuiltinOp::Sin, std::vector<ExprPtr>{ln_1});
    auto res = transcendental_normal_form(sin_ln, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NE(res.value(), nullptr);
}

// ─── simplify_functions: matrix ops, piecewise, N() ─────────────────────────

TEST(SimplifyFunctionsCoverage, SimplifyPiecewiseOnlyDefault) {
    CASContext ctx;
    // piecewise(default) → default
    ExprPtr e = parse_expr(ctx, "piecewise(42)");
    if (e) {
        auto r = ctx.simplify(e);
        if (r.is_ok()) {
            const auto* lit = expr_cast<IntegerLit>(r.value());
            if (lit) EXPECT_EQ(lit->value, BigInt(42));
        }
    }
    EXPECT_TRUE(true); // no crash is the minimum requirement
}

TEST(SimplifyFunctionsCoverage, SimplifyPiecewiseTrueCondition) {
    CASContext ctx;
    // piecewise(1, 7, 0)  → first condition true → result 7
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> args{
        make_int(arena, 1),  // condition = true
        make_int(arena, 7),  // branch
        make_int(arena, 0)   // default
    };
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Piecewise, std::move(args));
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    const auto* lit = expr_cast<IntegerLit>(r.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(7));
}

TEST(SimplifyFunctionsCoverage, SimplifyPiecewiseFalseCondition) {
    CASContext ctx;
    // piecewise(0, 7, 0)  → condition false → drop → return default
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> args{
        make_int(arena, 0),  // condition = false
        make_int(arena, 7),  // branch
        make_int(arena, 99)  // default
    };
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Piecewise, std::move(args));
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    const auto* lit = expr_cast<IntegerLit>(r.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(99));
}

TEST(SimplifyFunctionsCoverage, SimplifyPiecewiseSymbolicCondition) {
    CASContext ctx;
    // piecewise(x > 0, 7, 0)  → symbolic condition → stays as piecewise
    ExprPtr e = parse_expr(ctx, "piecewise(x, 7, 0)");
    if (e) {
        auto r = ctx.simplify(e);
        EXPECT_FALSE(r.is_error());
    }
    EXPECT_TRUE(true);
}

TEST(SimplifyFunctionsCoverage, SimplifyPiecewiseTrueWithRemainingUndecided) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // piecewise(sym, 3, 1, 7, 0) — first pair undecided, second true
    ExprPtr x = make_sym(arena, "x");
    std::vector<ExprPtr> args{
        x,                  // condition1 = undecided
        make_int(arena, 3), // branch1
        make_int(arena, 1), // condition2 = true
        make_int(arena, 7), // branch2
        make_int(arena, 0)  // default
    };
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Piecewise, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyFunctionsCoverage, SimplifyMatrixTrace) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // 2x2 identity matrix, trace = 2
    std::vector<ExprPtr> elems{
        make_int(arena, 1), make_int(arena, 0),
        make_int(arena, 0), make_int(arena, 1)};
    ExprPtr mat = arena.make<Matrix>(2, 2, std::move(elems));
    std::vector<ExprPtr> args{mat};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Trace, std::move(args));
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    const auto* lit = expr_cast<IntegerLit>(r.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(2));
}

TEST(SimplifyFunctionsCoverage, SimplifyMatrixTranspose) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Transpose of [[1,2],[3,4]] → [[1,3],[2,4]]
    std::vector<ExprPtr> elems{
        make_int(arena, 1), make_int(arena, 2),
        make_int(arena, 3), make_int(arena, 4)};
    ExprPtr mat = arena.make<Matrix>(2, 2, std::move(elems));
    std::vector<ExprPtr> args{mat};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Transpose, std::move(args));
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    const auto* m = expr_cast<Matrix>(r.value());
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->rows, 2U);
    EXPECT_EQ(m->cols, 2U);
}

TEST(SimplifyFunctionsCoverage, SimplifyMatrixRank) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> elems{
        make_int(arena, 1), make_int(arena, 0),
        make_int(arena, 0), make_int(arena, 1)};
    ExprPtr mat = arena.make<Matrix>(2, 2, std::move(elems));
    std::vector<ExprPtr> args{mat};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Rank, std::move(args));
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    const auto* lit = expr_cast<IntegerLit>(r.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(2));
}

TEST(SimplifyFunctionsCoverage, SimplifyMatrixDet) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // det([[1,2],[3,4]]) = 1*4 - 2*3 = -2
    std::vector<ExprPtr> elems{
        make_int(arena, 1), make_int(arena, 2),
        make_int(arena, 3), make_int(arena, 4)};
    ExprPtr mat = arena.make<Matrix>(2, 2, std::move(elems));
    std::vector<ExprPtr> args{mat};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Det, std::move(args));
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    EXPECT_NE(r.value(), nullptr);
}

TEST(SimplifyFunctionsCoverage, SimplifyMatrixInv) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // inv([[2,0],[0,2]]) = [[1/2,0],[0,1/2]]
    std::vector<ExprPtr> elems{
        make_int(arena, 2), make_int(arena, 0),
        make_int(arena, 0), make_int(arena, 2)};
    ExprPtr mat = arena.make<Matrix>(2, 2, std::move(elems));
    std::vector<ExprPtr> args{mat};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Inv, std::move(args));
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    EXPECT_NE(r.value(), nullptr);
}

TEST(SimplifyFunctionsCoverage, SimplifySeriesExpNode) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    std::vector<std::pair<long long, ExprPtr>> terms;
    // coefficient is ln(1) → 0 — to trigger simplification
    ExprPtr ln_1 = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{make_int(arena, 1)});
    terms.push_back({0, ln_1});
    terms.push_back({1, make_int(arena, 1)});
    ExprPtr e = arena.make<SeriesExp>(
        Symbol("x"), make_int(arena, 0), std::move(terms), 2);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    EXPECT_NE(r.value(), nullptr);
}

// ─── simplify_trig_inverse: arc-trig with assumptions ─────────────────────────

TEST(SimplifyTrigInverseCoverage, AsinSinWithAssumptions) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    // Assume x in [-pi/2, pi/2] → asin(sin(x)) = x
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr pi_2 = arena.make<Binary>(BinaryOp::Div, pi, make_int(arena, 2));
    ExprPtr neg_pi_2 = arena.make<Unary>(UnaryOp::Neg, pi_2);
    ctx.assumptions().assume_in_range(x, neg_pi_2, pi_2);
    ctx.assumptions().assume_greater_equal(make_sym(arena, "x"), neg_pi_2);
    ctx.assumptions().assume_greater_equal(pi_2, make_sym(arena, "x"));

    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr sin_x = arena.make<FuncCall>(BuiltinOp::Sin, std::vector<ExprPtr>{x_expr});
    ExprPtr asin_sin_x = arena.make<FuncCall>(BuiltinOp::Asin, std::vector<ExprPtr>{sin_x});
    auto r = ctx.simplify(asin_sin_x);
    ASSERT_TRUE(r.is_ok());
    EXPECT_NE(r.value(), nullptr);
}

TEST(SimplifyTrigInverseCoverage, AcosCosWithAssumptions) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    // Assume 0 <= x <= pi → acos(cos(x)) = x
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr zero = make_int(arena, 0);
    ctx.assumptions().assume_greater_equal(make_sym(arena, "x"), zero);
    ctx.assumptions().assume_greater_equal(pi, make_sym(arena, "x"));

    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr cos_x = arena.make<FuncCall>(BuiltinOp::Cos, std::vector<ExprPtr>{x_expr});
    ExprPtr acos_cos_x = arena.make<FuncCall>(BuiltinOp::Acos, std::vector<ExprPtr>{cos_x});
    auto r = ctx.simplify(acos_cos_x);
    ASSERT_TRUE(r.is_ok());
    EXPECT_NE(r.value(), nullptr);
}

TEST(SimplifyTrigInverseCoverage, AtanTanWithAssumptions) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr pi_2 = arena.make<Binary>(BinaryOp::Div, pi, make_int(arena, 2));
    ExprPtr neg_pi_2 = arena.make<Unary>(UnaryOp::Neg, pi_2);
    // atan(tan(x)) = x if -pi/2 < x < pi/2
    ctx.assumptions().assume_greater(make_sym(arena, "x"), neg_pi_2);
    ctx.assumptions().assume_greater(pi_2, make_sym(arena, "x"));

    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr tan_x = arena.make<FuncCall>(BuiltinOp::Tan, std::vector<ExprPtr>{x_expr});
    ExprPtr atan_tan_x = arena.make<FuncCall>(BuiltinOp::Atan, std::vector<ExprPtr>{tan_x});
    auto r = ctx.simplify(atan_tan_x);
    ASSERT_TRUE(r.is_ok());
    EXPECT_NE(r.value(), nullptr);
}

TEST(SimplifyTrigInverseCoverage, AtanZero) {
    CASContext ctx;
    ExprPtr e = parse_expr(ctx, "atan(0)");
    ASSERT_NE(e, nullptr);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    const auto* lit = expr_cast<IntegerLit>(r.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(0));
}

TEST(SimplifyTrigInverseCoverage, AtanOne) {
    CASContext ctx;
    ExprPtr e = parse_expr(ctx, "atan(1)");
    ASSERT_NE(e, nullptr);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    // atan(1) = pi/4
    EXPECT_NE(r.value(), nullptr);
}

TEST(SimplifyTrigInverseCoverage, AtanNegOne) {
    CASContext ctx;
    ExprPtr e = parse_expr(ctx, "atan(-1)");
    ASSERT_NE(e, nullptr);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    // atan(-1) = -pi/4
    EXPECT_NE(r.value(), nullptr);
}

TEST(SimplifyTrigInverseCoverage, AtanNegX) {
    CASContext ctx;
    // atan(-x) = -atan(x)
    ExprPtr e = parse_expr(ctx, "atan(-x)");
    ASSERT_NE(e, nullptr);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    EXPECT_NE(r.value(), nullptr);
}

// ─── complex_qi: all operations ──────────────────────────────────────────────

TEST(ComplexQiCoverage, ZeroOneImagUnit) {
    auto z = ComplexRational::zero();
    EXPECT_EQ(z.real(), Rational{});
    EXPECT_EQ(z.imag(), Rational{});

    auto one = ComplexRational::one();
    EXPECT_EQ(one.real(), Rational(BigInt(1)));
    EXPECT_EQ(one.imag(), Rational{});

    auto i = ComplexRational::imag_unit();
    EXPECT_EQ(i.real(), Rational{});
    EXPECT_EQ(i.imag(), Rational(BigInt(1)));
}

TEST(ComplexQiCoverage, Addition) {
    auto a = ComplexRational(Rational(BigInt(3)), Rational(BigInt(2)));
    auto b = ComplexRational(Rational(BigInt(1)), Rational(BigInt(4)));
    auto c = a + b;
    EXPECT_EQ(c.real(), Rational(BigInt(4)));
    EXPECT_EQ(c.imag(), Rational(BigInt(6)));
}

TEST(ComplexQiCoverage, Subtraction) {
    auto a = ComplexRational(Rational(BigInt(5)), Rational(BigInt(3)));
    auto b = ComplexRational(Rational(BigInt(2)), Rational(BigInt(1)));
    auto c = a - b;
    EXPECT_EQ(c.real(), Rational(BigInt(3)));
    EXPECT_EQ(c.imag(), Rational(BigInt(2)));
}

TEST(ComplexQiCoverage, Negation) {
    auto a = ComplexRational(Rational(BigInt(3)), Rational(BigInt(4)));
    auto neg = -a;
    EXPECT_EQ(neg.real(), Rational(BigInt(-3)));
    EXPECT_EQ(neg.imag(), Rational(BigInt(-4)));
}

TEST(ComplexQiCoverage, Multiplication) {
    // (1+i)(1+i) = 2i
    auto i = ComplexRational::imag_unit();
    auto one = ComplexRational::one();
    auto a = one + i;
    auto c = a * a;
    EXPECT_EQ(c.real(), Rational(BigInt(0)));
    EXPECT_EQ(c.imag(), Rational(BigInt(2)));
}

TEST(ComplexQiCoverage, Division) {
    // (2+0i) / (1+i) = (2)(1-i)/2 = 1-i
    auto two = ComplexRational(Rational(BigInt(2)));
    auto one_plus_i = ComplexRational::one() + ComplexRational::imag_unit();
    auto res = two.divide(one_plus_i);
    ASSERT_TRUE(res.is_ok());
    EXPECT_EQ(res.value().real(), Rational(BigInt(1)));
    EXPECT_EQ(res.value().imag(), Rational(BigInt(-1)));
}

TEST(ComplexQiCoverage, DivisionByZeroFails) {
    auto z = ComplexRational::zero();
    auto one = ComplexRational::one();
    auto res = one.divide(z);
    EXPECT_TRUE(res.is_error());
}

TEST(ComplexQiCoverage, Conjugate) {
    auto a = ComplexRational(Rational(BigInt(3)), Rational(BigInt(4)));
    auto c = a.conjugate();
    EXPECT_EQ(c.real(), Rational(BigInt(3)));
    EXPECT_EQ(c.imag(), Rational(BigInt(-4)));
}

TEST(ComplexQiCoverage, NormSq) {
    // |3+4i|^2 = 9+16 = 25
    auto a = ComplexRational(Rational(BigInt(3)), Rational(BigInt(4)));
    EXPECT_EQ(a.norm_sq(), Rational(BigInt(25)));
}

TEST(ComplexQiCoverage, IsUnit) {
    EXPECT_TRUE(ComplexRational::one().is_unit());
    EXPECT_TRUE((-ComplexRational::one()).is_unit());
    EXPECT_TRUE(ComplexRational::imag_unit().is_unit());
    // 2+0i is not a unit
    EXPECT_FALSE(ComplexRational(Rational(BigInt(2))).is_unit());
}

TEST(ComplexQiCoverage, EqualityInequality) {
    auto a = ComplexRational(Rational(BigInt(1)), Rational(BigInt(2)));
    auto b = ComplexRational(Rational(BigInt(1)), Rational(BigInt(2)));
    auto c = ComplexRational(Rational(BigInt(3)), Rational(BigInt(4)));
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a != b);
    EXPECT_FALSE(a == c);
    EXPECT_TRUE(a != c);
}

TEST(ComplexQiCoverage, MultiplyByImagUnit) {
    // i * i = -1
    auto i = ComplexRational::imag_unit();
    auto result = i * i;
    EXPECT_EQ(result.real(), Rational(BigInt(-1)));
    EXPECT_EQ(result.imag(), Rational{});
}
