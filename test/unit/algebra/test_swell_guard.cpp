// F7.0-A3.4 — Swell Guard tests for algebra::expand.
//
// Prevents combinatorial explosion (e.g. (x+y+z+w+v)^100) from saturating
// memory before the OOM-killer fires. Reports a structured Overflow error
// with explicit message rather than crashing the host process.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

using namespace cas;
using namespace cas::symbolic;
using namespace cas::algebra;

namespace {

ExprPtr sum_of_symbols(CASContext& ctx, std::initializer_list<const char*> names) {
    std::vector<ExprPtr> terms;
    for (const char* n : names) {
        terms.push_back(ctx.arena().make<Symbol>(n));
    }
    return ctx.arena().make<Sum>(std::move(terms));
}

ExprPtr pow(CASContext& ctx, ExprPtr base, long long k) {
    return ctx.arena().make<Binary>(
        BinaryOp::Pow,
        base,
        ctx.arena().make<IntegerLit>(BigInt(k)));
}

TEST(SwellGuard, DefaultLimitAllowsSmallExpand) {
    CASContext ctx;
    auto base = sum_of_symbols(ctx, {"x", "y"});  // n = 2
    auto p = pow(ctx, base, 10);                  // 2^10 = 1024 ≤ 100000
    auto res = expand(p, ctx);
    EXPECT_TRUE(res.is_ok());
}

TEST(SwellGuard, AggressiveLimitBlocksMediumExpand) {
    CASContext ctx;
    ctx.set_max_expand_monomials(50);             // tight budget
    auto base = sum_of_symbols(ctx, {"x", "y", "z"});  // n = 3
    auto p = pow(ctx, base, 4);                   // 3^4 = 81 > 50
    auto res = expand(p, ctx);
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::Overflow);
    EXPECT_NE(res.error().message.find("max_expand_monomials"), std::string::npos);
}

TEST(SwellGuard, RaisingLimitAllowsPreviouslyBlocked) {
    CASContext ctx;
    ctx.set_max_expand_monomials(50);
    auto base = sum_of_symbols(ctx, {"x", "y", "z"});
    auto p = pow(ctx, base, 4);
    ASSERT_TRUE(expand(p, ctx).is_error());

    ctx.set_max_expand_monomials(100000);
    auto res = expand(p, ctx);
    EXPECT_TRUE(res.is_ok());
}

TEST(SwellGuard, CatastrophicExpandBlockedFast) {
    CASContext ctx;
    // (a+b+c+d+e+f+g+h+i+j)^100 = 10^100 monomials — must reject in
    // microseconds without touching the allocator.
    auto base = sum_of_symbols(ctx, {"a","b","c","d","e","f","g","h","i","j"});
    auto p = pow(ctx, base, 100);
    auto res = expand(p, ctx);
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::Overflow);
}

}  // namespace
