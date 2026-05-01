#include <gtest/gtest.h>
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"
#include <vector>
#include <string>

using namespace cas::symbolic;
using namespace cas;

TEST(AssumptionsStabilityTest, DeepDeductionChain) {
    CASContext ctx;
    auto& arena = ctx.arena();
    
    // Create x0 > x1 > ... > x10000 > 0
    const int N = 10000;
    std::vector<ExprPtr> symbols;
    for (int i = 0; i <= N; ++i) {
        symbols.push_back(arena.make<Symbol>("x" + std::to_string(i)));
    }
    
    for (int i = 0; i < N; ++i) {
        ctx.assumptions().assume_greater(symbols[i], symbols[i+1]);
    }
    ctx.assumptions().assume_greater(symbols[N], nullptr); // x10000 > 0
    
    // Verify x0 > 0 and x0 > x10000
    // This previously would have caused stack overflow if recursive.
    EXPECT_TRUE(ctx.assumptions().is_positive(symbols[0]));
    EXPECT_TRUE(ctx.assumptions().is_greater(symbols[0], symbols[N]));
}

TEST(AssumptionsStabilityTest, CircularAssumptions) {
    CASContext ctx;
    auto& arena = ctx.arena();
    
    auto x = arena.make<Symbol>("x");
    auto y = arena.make<Symbol>("y");
    auto z = arena.make<Symbol>("z");
    
    // x > y > z > x (Cycle)
    ctx.assumptions().assume_greater(x, y);
    ctx.assumptions().assume_greater(y, z);
    ctx.assumptions().assume_greater(z, x);
    
    // Should not hang/crash.
    EXPECT_TRUE(ctx.assumptions().is_greater(x, z));
    EXPECT_FALSE(ctx.assumptions().is_greater(x, x));
}

TEST(AssumptionsStabilityTest, MixedCycle) {
    CASContext ctx;
    auto& arena = ctx.arena();
    
    auto a = arena.make<Symbol>("a");
    auto b = arena.make<Symbol>("b");
    
    // a >= b, b >= a
    ctx.assumptions().assume_greater_equal(a, b);
    ctx.assumptions().assume_greater_equal(b, a);
    
    EXPECT_TRUE(ctx.assumptions().is_greater_equal(a, b));
    EXPECT_TRUE(ctx.assumptions().is_greater_equal(b, a));
    EXPECT_TRUE(ctx.assumptions().is_greater_equal(a, a));
}
