#include "cas/symbolic.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include <gtest/gtest.h>

namespace cas::symbolic {
namespace {

Result<ExprPtr> parse_and_simplify(const std::string& input, CASContext& ctx) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), ctx.arena());
    auto expr = parser.parse();
    if (expr.is_error()) return expr;
    return ctx.simplify(expr.value());
}

} // namespace

TEST(LinalgIntegrationTest, SimplifiesDeterminant) {
    CASContext ctx;
    auto res = parse_and_simplify("det([[1, 2], [3, 4]])", ctx);
    ASSERT_TRUE(res.is_ok());
    
    auto expected = parse_and_simplify("-2", ctx);
    auto equal = mathematically_equal(res.value(), expected.value(), ctx);
    ASSERT_TRUE(equal.is_ok());
    EXPECT_TRUE(equal.value());
}

TEST(LinalgIntegrationTest, SimplifiesRank) {
    CASContext ctx;
    auto res = parse_and_simplify("rank([[1, 0], [0, 1]])", ctx);
    ASSERT_TRUE(res.is_ok());
    
    auto expected = parse_and_simplify("2", ctx);
    auto equal = mathematically_equal(res.value(), expected.value(), ctx);
    ASSERT_TRUE(equal.is_ok());
    EXPECT_TRUE(equal.value());
}

TEST(LinalgIntegrationTest, SimplifiesTrace) {
    CASContext ctx;
    auto res = parse_and_simplify("trace([[a, b], [c, d]])", ctx);
    ASSERT_TRUE(res.is_ok());
    
    // We expect "a + d"
    // Note: order might change due to term ordering, so we compare mathematically
    auto expected = parse_and_simplify("a + d", ctx);
    auto equal = mathematically_equal(res.value(), expected.value(), ctx);
    ASSERT_TRUE(equal.is_ok());
    EXPECT_TRUE(equal.value());
}

TEST(LinalgIntegrationTest, SimplifiesInverse) {
    CASContext ctx;
    // [[1, 2], [3, 4]]^-1 = [[-2, 1], [1.5, -0.5]]
    // det = -2. adj = [[4, -2], [-3, 1]]. inv = [[-2, 1], [3/2, -1/2]]
    auto res = parse_and_simplify("inv([[1, 2], [3, 4]])", ctx);
    ASSERT_TRUE(res.is_ok());
    
    auto expected = parse_and_simplify("[[-2, 1], [3/2, -1/2]]", ctx);
    auto equal = mathematically_equal(res.value(), expected.value(), ctx);
    ASSERT_TRUE(equal.is_ok());
    EXPECT_TRUE(equal.value());
}

} // namespace cas::symbolic
