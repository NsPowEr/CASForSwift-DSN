#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas::algebra {
namespace {

[[nodiscard]] Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

void expect_math_equal(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    auto expanded_lhs = expand(lhs, ctx);
    auto expanded_rhs = expand(rhs, ctx);
    ASSERT_TRUE(expanded_lhs.is_ok());
    ASSERT_TRUE(expanded_rhs.is_ok());

    auto eq = symbolic::mathematically_equal(expanded_lhs.value(), expanded_rhs.value(), ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value()) << "Expressions are not mathematically equivalent\n"
                            << "LHS: " << debug_print(expanded_lhs.value()) << "\n"
                            << "RHS: " << debug_print(expanded_rhs.value());
}

TEST(MultivariateGcdRobust, BivariateHigherDegree) {
    symbolic::CASContext ctx;
    auto p = parse_expr("(x^2 + y + 1) * (x + y)", ctx.arena());
    auto q = parse_expr("(x^2 + y + 1) * (x - y)", ctx.arena());
    auto expected = parse_expr("x^2 + y + 1", ctx.arena());
    
    auto gcd_res = polynomial_gcd_multivariate(p.value(), q.value(), ctx);
    ASSERT_TRUE(gcd_res.is_ok()) << gcd_res.error().message;
    expect_math_equal(gcd_res.value(), expected.value(), ctx);
}

TEST(MultivariateGcdRobust, TrivariateHigherDegree) {
    symbolic::CASContext ctx;
    auto p = parse_expr("(x^2 + y^2 + z^2) * (x + 1)", ctx.arena());
    auto q = parse_expr("(x^2 + y^2 + z^2) * (y + 1)", ctx.arena());
    auto expected = parse_expr("x^2 + y^2 + z^2", ctx.arena());
    
    auto gcd_res = polynomial_gcd_multivariate(p.value(), q.value(), ctx);
    ASSERT_TRUE(gcd_res.is_ok()) << gcd_res.error().message;
    expect_math_equal(gcd_res.value(), expected.value(), ctx);
}

TEST(MultivariateGcdRobust, QuadrivariateLinear) {
    symbolic::CASContext ctx;
    auto p = parse_expr("(x + y + z + w) * (x - y)", ctx.arena());
    auto q = parse_expr("(x + y + z + w) * (z + w)", ctx.arena());
    auto expected = parse_expr("x + y + z + w", ctx.arena());
    
    auto gcd_res = polynomial_gcd_multivariate(p.value(), q.value(), ctx);
    ASSERT_TRUE(gcd_res.is_ok()) << gcd_res.error().message;
    expect_math_equal(gcd_res.value(), expected.value(), ctx);
}

TEST(MultivariateGcdRobust, LargeCoefficients) {
    symbolic::CASContext ctx;
    auto p = parse_expr("x + 100", ctx.arena());
    auto q = parse_expr("x + 100", ctx.arena());
    auto expected = parse_expr("x + 100", ctx.arena());
    
    auto gcd_res = polynomial_gcd_multivariate(p.value(), q.value(), ctx);
    ASSERT_TRUE(gcd_res.is_ok()) << gcd_res.error().message;
    expect_math_equal(gcd_res.value(), expected.value(), ctx);
}

TEST(MultivariateGcdRobust, NonLinearBivariate) {
    symbolic::CASContext ctx;
    // (x^2 + y^2) * (x + y) and (x^2 + y^2) * (x - y)
    auto p = parse_expr("x^3 + x^2*y + x*y^2 + y^3", ctx.arena());
    auto q = parse_expr("x^3 - x^2*y + x*y^2 - y^3", ctx.arena());
    auto expected = parse_expr("x^2 + y^2", ctx.arena());
    
    auto gcd_res = polynomial_gcd_multivariate(p.value(), q.value(), ctx);
    ASSERT_TRUE(gcd_res.is_ok()) << gcd_res.error().message;
    expect_math_equal(gcd_res.value(), expected.value(), ctx);
}

TEST(MultivariateGcdRobust, IdenticalTrivariateLinear) {
    symbolic::CASContext ctx;
    auto p = parse_expr("x + y + z", ctx.arena());
    auto q = parse_expr("x + y + z", ctx.arena());
    auto expected = parse_expr("x + y + z", ctx.arena());
    
    auto gcd_res = polynomial_gcd_multivariate(p.value(), q.value(), ctx);
    ASSERT_TRUE(gcd_res.is_ok()) << gcd_res.error().message;
    expect_math_equal(gcd_res.value(), expected.value(), ctx);
}

TEST(MultivariateGcdRobust, QuadrivariateHigherDegree) {
    symbolic::CASContext ctx;
    auto p = parse_expr("(x^2 + y^2 + z^2 + w^2) * (x + y)", ctx.arena());
    auto q = parse_expr("(x^2 + y^2 + z^2 + w^2) * (z + w)", ctx.arena());
    auto expected = parse_expr("x^2 + y^2 + z^2 + w^2", ctx.arena());
    
    auto gcd_res = polynomial_gcd_multivariate(p.value(), q.value(), ctx);
    ASSERT_TRUE(gcd_res.is_ok()) << gcd_res.error().message;
    expect_math_equal(gcd_res.value(), expected.value(), ctx);
}

} // namespace
} // namespace cas::algebra
