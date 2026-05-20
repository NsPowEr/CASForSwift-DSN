#include <gtest/gtest.h>
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

using namespace cas::symbolic;
using namespace cas;

TEST(AssumptionsDomainTest, DomainInference) {
    CASContext ctx;
    auto& arena = ctx.arena();
    auto x = arena.make<Symbol>("x");
    
    ctx.assumptions().assume_domain(expr_ref<Symbol>(x), Domain::Positive);
    EXPECT_TRUE(ctx.assumptions().is_positive(x));
    EXPECT_TRUE(ctx.assumptions().is_real(x));
    EXPECT_TRUE(ctx.assumptions().is_nonzero(x));
    EXPECT_FALSE(ctx.assumptions().is_negative(x));
    EXPECT_EQ(ctx.assumptions().get_domain(expr_ref<Symbol>(x)), Domain::Positive);
}

TEST(AssumptionsDomainTest, NaturalDomain) {
    CASContext ctx;
    auto& arena = ctx.arena();
    auto n = arena.make<Symbol>("n");
    
    ctx.assumptions().assume_domain(expr_ref<Symbol>(n), Domain::Natural);
    EXPECT_TRUE(ctx.assumptions().is_integer(n));
    EXPECT_TRUE(ctx.assumptions().is_nonnegative(n));
    EXPECT_TRUE(ctx.assumptions().is_greater_equal(n, nullptr));
}

TEST(AssumptionsInferenceTest, ProductDeduction) {
    CASContext ctx;
    auto& arena = ctx.arena();
    auto x = arena.make<Symbol>("x");
    auto y = arena.make<Symbol>("y");
    
    ctx.assumptions().assume_domain(expr_ref<Symbol>(x), Domain::Negative);
    ctx.assumptions().assume_domain(expr_ref<Symbol>(y), Domain::Negative);
    
    auto prod = arena.make<Product>(std::vector<ExprPtr>{x, y});
    EXPECT_TRUE(ctx.assumptions().is_positive(prod));
    
    auto z = arena.make<Symbol>("z");
    ctx.assumptions().assume_domain(expr_ref<Symbol>(z), Domain::Positive);
    auto prod3 = arena.make<Product>(std::vector<ExprPtr>{x, y, z});
    EXPECT_TRUE(ctx.assumptions().is_positive(prod3));
}

TEST(AssumptionsConsistencyTest, ContradictionPositiveNegative) {
    CASContext ctx;
    auto& arena = ctx.arena();
    auto x = arena.make<Symbol>("x");
    
    ctx.assumptions().assume_positive(expr_ref<Symbol>(x));
    ctx.assumptions().assume_domain(expr_ref<Symbol>(x), Domain::Negative);
    
    auto res = ctx.assumptions().check_consistency();
    EXPECT_TRUE(res.is_error());
}

TEST(AssumptionsConsistencyTest, CyclicStrictInequality) {
    CASContext ctx;
    auto& arena = ctx.arena();
    auto x = arena.make<Symbol>("x");
    auto y = arena.make<Symbol>("y");
    
    ctx.assumptions().assume_greater(x, y);
    ctx.assumptions().assume_greater(y, x);
    
    auto res = ctx.assumptions().check_consistency();
    EXPECT_TRUE(res.is_error());
}

TEST(AssumptionsConsistencyTest, SelfLoopStrict) {
    CASContext ctx;
    auto& arena = ctx.arena();
    auto x = arena.make<Symbol>("x");
    
    ctx.assumptions().assume_greater(x, x);
    
    auto res = ctx.assumptions().check_consistency();
    EXPECT_TRUE(res.is_error());
}

TEST(AssumptionsInferenceTest, TransitiveProperty) {
    CASContext ctx;
    auto& arena = ctx.arena();
    auto a = arena.make<Symbol>("a");
    auto b = arena.make<Symbol>("b");
    auto c = arena.make<Symbol>("c");
    auto d = arena.make<Symbol>("d");
    
    ctx.assumptions().assume_greater(a, b);
    ctx.assumptions().assume_greater_equal(b, c);
    ctx.assumptions().assume_greater(c, d);
    
    EXPECT_TRUE(ctx.assumptions().is_greater(a, d));
    EXPECT_TRUE(ctx.assumptions().is_greater(a, c));
    EXPECT_TRUE(ctx.assumptions().is_greater_equal(a, d));
}
