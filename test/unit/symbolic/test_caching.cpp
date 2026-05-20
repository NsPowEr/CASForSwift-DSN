#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include <gtest/gtest.h>

namespace cas::symbolic {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

TEST(CASCachingTest, LRUEvictionWorks) {
    CASContext ctx;
    ctx.set_cache_limit(2); // Small limit to test eviction
    
    auto e1 = parse_expr("x+1", ctx.arena()).value();
    auto e2 = parse_expr("x+2", ctx.arena()).value();
    auto e3 = parse_expr("x+3", ctx.arena()).value();
    
    // Fill cache
    (void)ctx.simplify(e1); // Miss
    (void)ctx.simplify(e2); // Miss
    
    auto metrics = ctx.get_simplify_metrics();
    EXPECT_EQ(metrics.misses, 2);
    EXPECT_EQ(metrics.hits, 0);
    EXPECT_EQ(metrics.evictions, 0);
    
    // Add third item, should evict e1 (LRU)
    (void)ctx.simplify(e3); // Miss
    metrics = ctx.get_simplify_metrics();
    EXPECT_EQ(metrics.misses, 3);
    EXPECT_EQ(metrics.evictions, 1);
    
    // e1 should be a miss now
    (void)ctx.simplify(e1); // Miss
    metrics = ctx.get_simplify_metrics();
    EXPECT_EQ(metrics.misses, 4);
    EXPECT_EQ(metrics.hits, 0);
    
    // e3 should be a hit (it was most recent before e1 miss)
    (void)ctx.simplify(e3); // Hit
    metrics = ctx.get_simplify_metrics();
    EXPECT_EQ(metrics.hits, 1);
    
    // e2 should NOT be there anymore (it was evicted when e1 was re-added)
    (void)ctx.simplify(e2); // Miss
    metrics = ctx.get_simplify_metrics();
    EXPECT_EQ(metrics.misses, 5);
    EXPECT_EQ(metrics.hits, 1);
}

TEST(CASCachingTest, MetricsTracking) {
    CASContext ctx;
    ctx.set_cache_limit(10);
    
    auto e1 = parse_expr("a+b", ctx.arena()).value();
    
    (void)ctx.simplify(e1); // Miss
    (void)ctx.simplify(e1); // Hit
    (void)ctx.simplify(e1); // Hit
    
    auto metrics = ctx.get_simplify_metrics();
    EXPECT_EQ(metrics.misses, 1);
    EXPECT_EQ(metrics.hits, 2);
    EXPECT_EQ(metrics.evictions, 0);
}

TEST(CASCachingTest, ConfigurableLimit) {
    CASContext ctx;
    ctx.set_cache_limit(100);
    EXPECT_EQ(ctx.get_cache_limit(), 100);
    
    ctx.set_cache_limit(500);
    EXPECT_EQ(ctx.get_cache_limit(), 500);
}

TEST(CASCachingTest, DisableCachingSkipsLookup) {
    CASContext ctx;
    auto e1 = parse_expr("x*y", ctx.arena()).value();
    
    (void)ctx.simplify(e1); // Miss (1)
    EXPECT_EQ(ctx.get_simplify_metrics().misses, 1);
    
    ctx.set_caching_enabled(false);
    (void)ctx.simplify(e1); // Skipped
    EXPECT_EQ(ctx.get_simplify_metrics().misses, 1);
    
    ctx.set_caching_enabled(true);
    (void)ctx.simplify(e1); // Miss again, because cache was cleared when disabled
    EXPECT_EQ(ctx.get_simplify_metrics().misses, 2);
}

TEST(CASCachingTest, CollectGarbagePreservesCache) {
    CASContext ctx;
    ctx.set_cache_limit(10);
    
    auto e1 = parse_expr("x+y", ctx.arena()).value();
    (void)ctx.simplify(e1); // Miss
    
    EXPECT_EQ(ctx.get_simplify_metrics().misses, 1);
    
    std::vector<ExprPtr*> roots = {&e1};
    ctx.collect_garbage(roots);
    
    (void)ctx.simplify(e1); // Should be a HIT if GC correctly migrated the cache
    EXPECT_EQ(ctx.get_simplify_metrics().hits, 1);
}

} // namespace
} // namespace cas::symbolic
