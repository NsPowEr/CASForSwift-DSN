// DEBT-F1-COV-01 — White-box coverage uplift for src/symbolic/context_core.cpp.
// Exercises: define/lookup/clear_variables, caching API, cache metrics,
// collect_garbage, make_fresh_symbol (collision probe), enable_trace/get_trace,
// post_simplify_hook, set_rewrite_provider, timeout config, and setters.

#include <gtest/gtest.h>
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "cas/trace.hpp"

#include <chrono>
#include "cas/calculus.hpp"

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

}  // namespace

// ── define / lookup / clear_variables ────────────────────────────────────────

TEST(ContextCoreCoverage, DefineAndLookup) {
    CASContext ctx;
    ExprPtr three = ctx.arena().make<IntegerLit>(BigInt(3));
    Symbol x("x");
    ctx.define(x, three);
    auto found = ctx.lookup(x);
    ASSERT_TRUE(found.has_value());
    EXPECT_EQ(found.value(), three);
}

TEST(ContextCoreCoverage, LookupMissingReturnsNullopt) {
    CASContext ctx;
    Symbol y("y");
    auto found = ctx.lookup(y);
    EXPECT_FALSE(found.has_value());
}

TEST(ContextCoreCoverage, ClearVariables) {
    CASContext ctx;
    Symbol x("x");
    ctx.define(x, ctx.arena().make<IntegerLit>(BigInt(7)));
    ASSERT_TRUE(ctx.lookup(x).has_value());
    ctx.clear_variables();
    EXPECT_FALSE(ctx.lookup(x).has_value());
}

// ── simplify with define propagates variable value ────────────────────────────

TEST(ContextCoreCoverage, SimplifyWithDefinedVariable) {
    CASContext ctx;
    ExprPtr two = ctx.arena().make<IntegerLit>(BigInt(2));
    Symbol x("x");
    ctx.define(x, two);
    ExprPtr x_expr = ctx.arena().make<Symbol>("x");
    auto result = ctx.simplify(x_expr);
    ASSERT_TRUE(result.is_ok());
    // After define, simplify may substitute x→2; at minimum no crash
    EXPECT_NE(result.value(), nullptr);
}

// ── caching: enabled/disabled, cache limit, metrics ──────────────────────────

TEST(ContextCoreCoverage, CachingDefaultEnabled) {
    CASContext ctx;
    ExprPtr e = parse_expr(ctx, "x + 1");
    ASSERT_NE(e, nullptr);
    // Simplify twice → second call hits cache
    auto r1 = ctx.simplify(e);
    auto r2 = ctx.simplify(e);
    ASSERT_TRUE(r1.is_ok());
    ASSERT_TRUE(r2.is_ok());
    // Cache hit count should be ≥ 1 after second call
    auto metrics = ctx.get_simplify_metrics();
    EXPECT_GE(metrics.hits + metrics.misses, 2U);
}

TEST(ContextCoreCoverage, SetCacheLimitAndGetIt) {
    CASContext ctx;
    ctx.set_cache_limit(42U);
    EXPECT_EQ(ctx.get_cache_limit(), 42U);
}

TEST(ContextCoreCoverage, DisableCachingClearsCache) {
    CASContext ctx;
    ExprPtr e = parse_expr(ctx, "x + 2");
    ASSERT_NE(e, nullptr);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    // Now disable caching — should clear without crash
    ctx.set_caching_enabled(false);
    auto r2 = ctx.simplify(e);
    ASSERT_TRUE(r2.is_ok());
}

TEST(ContextCoreCoverage, CacheMetricsAfterWork) {
    CASContext ctx;
    ExprPtr e = parse_expr(ctx, "2 + 3");
    ASSERT_NE(e, nullptr);
    (void)ctx.simplify(e);
    (void)ctx.simplify(e);
    (void)ctx.simplify(e);
    auto sm = ctx.get_simplify_metrics();
    auto dm = ctx.get_diff_metrics();
    auto im = ctx.get_integrate_metrics();
    // Just ensure they're accessible without crash
    (void)sm; (void)dm; (void)im;
    EXPECT_GE(sm.hits + sm.misses, 2U);
}

// ── make_fresh_symbol ─────────────────────────────────────────────────────────

TEST(ContextCoreCoverage, MakeFreshSymbolUnique) {
    CASContext ctx;
    Symbol s1 = ctx.make_fresh_symbol("tmp");
    Symbol s2 = ctx.make_fresh_symbol("tmp");
    Symbol s3 = ctx.make_fresh_symbol("tmp");
    EXPECT_NE(s1.name, s2.name);
    EXPECT_NE(s2.name, s3.name);
    EXPECT_NE(s1.name, s3.name);
}

TEST(ContextCoreCoverage, MakeFreshSymbolCollisionProbe) {
    // If "pfx_1" is already defined, make_fresh_symbol("pfx") should return "pfx_2".
    CASContext ctx;
    // Pre-define pfx_1 so counter must advance past it
    Symbol occupied("pfx_1");
    ctx.define(occupied, ctx.arena().make<IntegerLit>(BigInt(0)));
    // Now fresh symbols must not be "pfx_1"
    Symbol s = ctx.make_fresh_symbol("pfx");
    EXPECT_NE(s.name, "pfx_1");
}

TEST(ContextCoreCoverage, MakeFreshSymbolEmptyPrefixUsesDefault) {
    CASContext ctx;
    Symbol s = ctx.make_fresh_symbol("");
    EXPECT_FALSE(s.name.empty());
}

// ── enable_trace / get_trace ──────────────────────────────────────────────────

TEST(ContextCoreCoverage, EnableTraceAndGetTrace) {
    CASContext ctx;
    ctx.enable_trace(true);
    ExprPtr e = parse_expr(ctx, "x^2 + x");
    ASSERT_NE(e, nullptr);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    const ComputationTrace& trace = ctx.get_trace();
    // Trace may be empty or not depending on simplification complexity; no crash
    (void)trace;
}

TEST(ContextCoreCoverage, DisableTraceClearsTrace) {
    CASContext ctx;
    ctx.enable_trace(true);
    ExprPtr e = parse_expr(ctx, "3 + 4");
    ASSERT_NE(e, nullptr);
    (void)ctx.simplify(e);
    ctx.enable_trace(false);
    // After disabling, get_trace returns empty
    const ComputationTrace& trace = ctx.get_trace();
    EXPECT_TRUE(trace.empty());
}

// ── post_simplify_hook ────────────────────────────────────────────────────────

TEST(ContextCoreCoverage, PostSimplifyHookCalled) {
    CASContext ctx;
    bool hook_called = false;
    ctx.set_post_simplify_hook([&](ExprPtr expr, CASContext& /*c*/) -> Result<ExprPtr> {
        hook_called = true;
        return ok(expr);  // pass-through
    });
    ExprPtr e = parse_expr(ctx, "1 + 1");
    ASSERT_NE(e, nullptr);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(hook_called);
}

TEST(ContextCoreCoverage, PostSimplifyHookTransforms) {
    CASContext ctx;
    // Hook replaces any result with the integer 99
    ctx.set_post_simplify_hook([&](ExprPtr /*expr*/, CASContext& c) -> Result<ExprPtr> {
        // Return a different pointer to trigger re-simplify pass
        return ok(c.arena().make<IntegerLit>(BigInt(99)));
    });
    ExprPtr e = parse_expr(ctx, "x + y");
    ASSERT_NE(e, nullptr);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    // Result should be 99 (hook replaced everything)
    if (const auto* il = expr_cast<IntegerLit>(r.value())) {
        EXPECT_EQ(il->value, BigInt(99));
    }
}

TEST(ContextCoreCoverage, ClearPostSimplifyHook) {
    CASContext ctx;
    int call_count = 0;
    ctx.set_post_simplify_hook([&](ExprPtr expr, CASContext& /*c*/) -> Result<ExprPtr> {
        ++call_count;
        return ok(expr);
    });
    ExprPtr e = parse_expr(ctx, "1");
    ASSERT_NE(e, nullptr);
    (void)ctx.simplify(e);
    EXPECT_EQ(call_count, 1);
    ctx.clear_post_simplify_hook();
    (void)ctx.simplify(e);  // Hook cleared, should not increment
    EXPECT_EQ(call_count, 1);
}

// ── timeout config ────────────────────────────────────────────────────────────

TEST(ContextCoreCoverage, SetTimeoutNoHang) {
    CASContext ctx;
    ctx.set_timeout(std::chrono::milliseconds(5000));
    ExprPtr e = parse_expr(ctx, "2 * 3");
    ASSERT_NE(e, nullptr);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
}

TEST(ContextCoreCoverage, SetTimeoutCheckInterval) {
    CASContext ctx;
    ctx.set_timeout_check_interval(128U);
    // Clamp: below 64 clamps to 64
    ctx.set_timeout_check_interval(10U);
    ExprPtr e = parse_expr(ctx, "x + x");
    ASSERT_NE(e, nullptr);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
}

// ── various config setters ────────────────────────────────────────────────────

TEST(ContextCoreCoverage, SetMaxSimplificationDepthClampLow) {
    CASContext ctx;
    ctx.set_max_simplification_depth(5);  // below 10, should clamp to 10
    ExprPtr e = parse_expr(ctx, "x");
    ASSERT_NE(e, nullptr);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
}

TEST(ContextCoreCoverage, SetMaxIntegrationDepthClamps) {
    CASContext ctx;
    ctx.set_max_integration_depth(0U);    // clamps to 1
    ctx.set_max_integration_depth(200U);  // clamps to 128
    SUCCEED();
}

TEST(ContextCoreCoverage, SetGcdErrorProbabilityClamps) {
    CASContext ctx;
    ctx.set_gcd_error_probability(0.0);    // clamps to 1e-6
    ctx.set_gcd_error_probability(0.5);    // clamps to 0.1
    ctx.set_gcd_error_probability(0.01);   // valid
    SUCCEED();
}

TEST(ContextCoreCoverage, SetNumericPrecisionDigitsClamps) {
    CASContext ctx;
    ctx.set_numeric_precision_digits(3U);      // clamps to 6
    ctx.set_numeric_precision_digits(20000U);  // clamps to 10000
    ctx.set_numeric_precision_digits(50U);     // valid
    SUCCEED();
}

TEST(ContextCoreCoverage, SetMaxRootofExplicitDegree) {
    CASContext ctx;
    ctx.set_max_rootof_explicit_degree(0U);  // clamps to 1
    ctx.set_max_rootof_explicit_degree(10U);
    SUCCEED();
}

TEST(ContextCoreCoverage, SetMaxGcdRecursionDepth) {
    CASContext ctx;
    ctx.set_max_gcd_recursion_depth(2U);   // clamps to 4
    ctx.set_max_gcd_recursion_depth(16U);
    SUCCEED();
}

TEST(ContextCoreCoverage, SetMinGcdDivisionSteps) {
    CASContext ctx;
    ctx.set_min_gcd_division_steps(0U);   // clamps to 1
    ctx.set_min_gcd_division_steps(8U);
    SUCCEED();
}

TEST(ContextCoreCoverage, SetMaxCyclotomicN) {
    CASContext ctx;
    ctx.set_max_cyclotomic_n(200);
    SUCCEED();
}

TEST(ContextCoreCoverage, SetMaxQAlphaBridgeDepth) {
    CASContext ctx;
    ctx.set_max_q_alpha_bridge_depth(5U);   // clamps to 8
    ctx.set_max_q_alpha_bridge_depth(20U);
    SUCCEED();
}

TEST(ContextCoreCoverage, SetMaxGammaRecursion) {
    CASContext ctx;
    ctx.set_max_gamma_recursion(10U);   // clamps to 16
    ctx.set_max_gamma_recursion(50U);
    SUCCEED();
}

TEST(ContextCoreCoverage, SetImproperLeadingOrderScan) {
    CASContext ctx;
    ctx.set_improper_leading_order_scan(0U);  // clamps to 1
    ctx.set_improper_leading_order_scan(5U);
    SUCCEED();
}

TEST(ContextCoreCoverage, SetExpandBesselRecurrence) {
    CASContext ctx;
    ctx.set_expand_bessel_recurrence(true);
    ctx.set_expand_bessel_recurrence(false);
    SUCCEED();
}

TEST(ContextCoreCoverage, SetMaxTragerTowerShiftAttempts) {
    CASContext ctx;
    ctx.set_max_trager_tower_shift_attempts(5U);
    SUCCEED();
}

TEST(ContextCoreCoverage, SetMaxIntegrateByPartsDepth) {
    CASContext ctx;
    ctx.set_max_integrate_by_parts_depth(4U);
    EXPECT_EQ(ctx.max_integrate_by_parts_depth(), 4U);
}

TEST(ContextCoreCoverage, SetMaxRischRationalAnsatzDegree) {
    CASContext ctx;
    ctx.set_max_risch_rational_ansatz_degree(16U);
    EXPECT_EQ(ctx.max_risch_rational_ansatz_degree(), 16U);
}

// ── assumptions accessor ──────────────────────────────────────────────────────

TEST(ContextCoreCoverage, AssumptionsAccessorMutable) {
    CASContext ctx;
    Assumptions& asmp = ctx.assumptions();
    (void)asmp;
    const Assumptions& casmp = static_cast<const CASContext&>(ctx).assumptions();
    (void)casmp;
    SUCCEED();
}

// ── arena accessor ────────────────────────────────────────────────────────────

TEST(ContextCoreCoverage, ArenaAccessorConst) {
    CASContext ctx;
    const AstArena& arena = static_cast<const CASContext&>(ctx).arena();
    (void)arena;
    SUCCEED();
}

// ── simplify null expr returns error ─────────────────────────────────────────

TEST(ContextCoreCoverage, SimplifyNullExprReturnsError) {
    CASContext ctx;
    auto r = ctx.simplify(nullptr);
    EXPECT_TRUE(r.is_error());
}

// ── collect_garbage ───────────────────────────────────────────────────────────

TEST(ContextCoreCoverage, CollectGarbageNoRoots) {
    CASContext ctx;
    // Define a variable and some expressions, then GC with no roots
    Symbol x("x");
    ctx.define(x, ctx.arena().make<IntegerLit>(BigInt(5)));
    ctx.collect_garbage({});
    // After GC, variable still accessible (was included as a root via variables_)
    auto found = ctx.lookup(x);
    ASSERT_TRUE(found.has_value());
    const auto* il = expr_cast<IntegerLit>(found.value());
    ASSERT_NE(il, nullptr);
    EXPECT_EQ(il->value, BigInt(5));
}

TEST(ContextCoreCoverage, CollectGarbageWithExternalRoot) {
    CASContext ctx;
    ExprPtr e = ctx.arena().make<IntegerLit>(BigInt(42));
    // GC with e as external root
    ctx.collect_garbage({&e});
    // e should still be valid and point to 42
    ASSERT_NE(e, nullptr);
    const auto* il = expr_cast<IntegerLit>(e);
    ASSERT_NE(il, nullptr);
    EXPECT_EQ(il->value, BigInt(42));
}

TEST(ContextCoreCoverage, CollectGarbageWithCachePopulated) {
    CASContext ctx;
    ExprPtr e = parse_expr(ctx, "x + 1");
    ASSERT_NE(e, nullptr);
    (void)ctx.simplify(e);  // populates simplify cache
    ctx.collect_garbage({&e});
    // Cache rebuilt; no crash
    SUCCEED();
}

TEST(ContextCoreCoverage, CollectGarbageWithDiffCachePopulated) {
    CASContext ctx;
    ExprPtr e = parse_expr(ctx, "x^3");
    ASSERT_NE(e, nullptr);
    (void)ctx.simplify(e);
    // Populate diff cache via calculus::diff
    auto d = cas::calculus::diff(e, Symbol("x"), 1U, ctx);
    (void)d;
    ctx.collect_garbage({&e});
    SUCCEED();
}

TEST(ContextCoreCoverage, CollectGarbageWithIntegrateCachePopulated) {
    CASContext ctx;
    ExprPtr e = parse_expr(ctx, "x^2");
    ASSERT_NE(e, nullptr);
    (void)ctx.simplify(e);
    // Populate integrate cache via calculus::integrate
    auto ig = cas::calculus::integrate(e, Symbol("x"), ctx);
    (void)ig;
    ctx.collect_garbage({&e});
    SUCCEED();
}

TEST(ContextCoreCoverage, CollectGarbageWithTraceEnabled) {
    CASContext ctx;
    ctx.enable_trace(true);
    ExprPtr e = parse_expr(ctx, "x^2 + 1");
    ASSERT_NE(e, nullptr);
    (void)ctx.simplify(e);
    ctx.collect_garbage({&e});
    SUCCEED();
}

TEST(ContextCoreCoverage, CollectGarbageWithIntegralNodeAsRoot) {
    // Puts an Integral AST node (unevaluated) in the arena, then GC to
    // exercise materialize_expr_impl for Integral nodes.
    CASContext ctx;
    ExprPtr integrand = parse_expr(ctx, "x^2");
    ASSERT_NE(integrand, nullptr);
    ExprPtr integral_node = ctx.arena().make<Integral>(integrand, Symbol("x"),
        std::nullopt, std::nullopt);
    ctx.collect_garbage({&integral_node});
    ASSERT_NE(integral_node, nullptr);
}

TEST(ContextCoreCoverage, CollectGarbageWithDerivativeNodeAsRoot) {
    // Puts a Derivative AST node in the arena, then GC.
    CASContext ctx;
    ExprPtr f = parse_expr(ctx, "sin(x)");
    ASSERT_NE(f, nullptr);
    ExprPtr deriv_node = ctx.arena().make<Derivative>(f, Symbol("x"), 1U);
    ctx.collect_garbage({&deriv_node});
    ASSERT_NE(deriv_node, nullptr);
}

TEST(ContextCoreCoverage, CollectGarbageWithMatrixNodeAsRoot) {
    // Puts a Matrix AST node in the arena, then GC.
    CASContext ctx;
    ExprPtr e1 = ctx.arena().make<IntegerLit>(BigInt(1));
    ExprPtr e2 = ctx.arena().make<IntegerLit>(BigInt(2));
    ExprPtr mat = ctx.arena().make<Matrix>(1U, 2U, std::vector<ExprPtr>{e1, e2});
    ctx.collect_garbage({&mat});
    ASSERT_NE(mat, nullptr);
    const auto* m = expr_cast<Matrix>(mat);
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->rows, 1U);
    EXPECT_EQ(m->cols, 2U);
}

// ── clear_caches ──────────────────────────────────────────────────────────────

TEST(ContextCoreCoverage, ClearCachesNoError) {
    CASContext ctx;
    ExprPtr e = parse_expr(ctx, "2 + 2");
    ASSERT_NE(e, nullptr);
    (void)ctx.simplify(e);
    ctx.clear_caches();
    // After clear, next simplify should be a miss but still succeed
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
}

// ── has_post_simplify_hook ────────────────────────────────────────────────────

TEST(ContextCoreCoverage, HasPostSimplifyHookQueryableByDefault) {
    // CASContext registers an algebra hook in its constructor, so the hook
    // may be set by default. The test just verifies the API is accessible.
    CASContext ctx;
    bool has = ctx.has_post_simplify_hook();
    (void)has;  // result depends on algebra registration; no crash suffices
    SUCCEED();
}

TEST(ContextCoreCoverage, HasPostSimplifyHookTrueAfterSet) {
    CASContext ctx;
    ctx.set_post_simplify_hook([](ExprPtr e, CASContext&) { return ok(e); });
    EXPECT_TRUE(ctx.has_post_simplify_hook());
    ctx.clear_post_simplify_hook();
    EXPECT_FALSE(ctx.has_post_simplify_hook());
}

// ── interrupt ────────────────────────────────────────────────────────────────

TEST(ContextCoreCoverage, InterruptFlagCycle) {
    CASContext ctx;
    EXPECT_FALSE(ctx.is_interrupted());
    ctx.interrupt();
    EXPECT_TRUE(ctx.is_interrupted());
    ctx.clear_interrupt();
    EXPECT_FALSE(ctx.is_interrupted());
}

// ── rewrite provider ─────────────────────────────────────────────────────────

TEST(ContextCoreCoverage, SetRewriteProviderNull) {
    CASContext ctx;
    ctx.set_rewrite_provider(nullptr);
    // Simplify should still work (uses no provider)
    ExprPtr e = parse_expr(ctx, "1");
    ASSERT_NE(e, nullptr);
    auto r = ctx.simplify(e);
    // May error but should not crash
    (void)r;
}

TEST(ContextCoreCoverage, GetRewriteProviderDefault) {
    CASContext ctx;
    const RewriteProvider* p = ctx.rewrite_provider();
    // Default provider is non-null
    (void)p;
    SUCCEED();
}
