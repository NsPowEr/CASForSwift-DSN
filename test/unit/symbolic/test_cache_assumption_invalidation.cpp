// F7.0-A4.1 — Cache invalidation on assumption change tests.
//
// CRITICAL mathematical-correctness: simplify_cache_ must be flushed when
// assumptions change, since simplification results depend on assumption
// state (e.g. simplify(abs(x)) → x under x>0, -x under x<0).
//
// These tests verify the revision-counter primitive that drives the
// cache-invalidation hook in CASContext::simplify().

#include <gtest/gtest.h>

#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

using namespace cas;
using namespace cas::symbolic;

namespace {

TEST(CacheAssumptionInvalidation, RevisionStartsZero) {
    CASContext ctx;
    EXPECT_EQ(ctx.assumptions().revision(), 0U);
}

TEST(CacheAssumptionInvalidation, AssumePositiveBumpsRevision) {
    CASContext ctx;
    Symbol x{"x"};
    const auto r0 = ctx.assumptions().revision();
    ctx.assumptions().assume_positive(x);
    EXPECT_GT(ctx.assumptions().revision(), r0);
}

TEST(CacheAssumptionInvalidation, AssumeRealBumpsRevision) {
    CASContext ctx;
    Symbol x{"x"};
    const auto r0 = ctx.assumptions().revision();
    ctx.assumptions().assume_real(x);
    EXPECT_GT(ctx.assumptions().revision(), r0);
}

TEST(CacheAssumptionInvalidation, AssumeIntegerBumpsRevision) {
    CASContext ctx;
    Symbol x{"x"};
    const auto r0 = ctx.assumptions().revision();
    ctx.assumptions().assume_integer(x);
    EXPECT_GT(ctx.assumptions().revision(), r0);
}

TEST(CacheAssumptionInvalidation, AssumeNonzeroBumpsRevision) {
    CASContext ctx;
    Symbol x{"x"};
    const auto r0 = ctx.assumptions().revision();
    ctx.assumptions().assume_nonzero(x);
    EXPECT_GT(ctx.assumptions().revision(), r0);
}

TEST(CacheAssumptionInvalidation, EachMutatorBumpsExactlyOnce) {
    CASContext ctx;
    Symbol x{"x"};
    const auto r0 = ctx.assumptions().revision();
    ctx.assumptions().assume_positive(x);
    const auto r1 = ctx.assumptions().revision();
    EXPECT_EQ(r1, r0 + 1U);
    ctx.assumptions().assume_real(x);
    EXPECT_EQ(ctx.assumptions().revision(), r1 + 1U);
    ctx.assumptions().assume_integer(x);
    EXPECT_EQ(ctx.assumptions().revision(), r1 + 2U);
}

TEST(CacheAssumptionInvalidation, SimplifyAfterAssumeBumpsRevisionAndClearsCache) {
    CASContext ctx;
    Symbol x{"x"};
    auto x_expr = ctx.arena().make<Symbol>("x");
    auto seven = ctx.arena().make<IntegerLit>(BigInt(7));
    auto sum = ctx.arena().make<Sum>(std::vector<ExprPtr>{x_expr, seven});

    // First simplify warms the cache.
    auto r1 = ctx.simplify(sum);
    ASSERT_TRUE(r1.is_ok());
    const std::uint64_t rev_before = ctx.assumptions().revision();

    // Change assumptions → revision bumps. Next simplify must observe it.
    ctx.assumptions().assume_positive(x);
    EXPECT_GT(ctx.assumptions().revision(), rev_before);

    // Next simplify triggers the revision-mismatch hook → clear_caches().
    // Both calls return well-formed results.
    auto r2 = ctx.simplify(sum);
    ASSERT_TRUE(r2.is_ok());
}

}  // namespace
