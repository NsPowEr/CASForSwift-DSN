#include "cas/algebra.hpp"
#include "cas/ast.hpp"
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
    auto eq = symbolic::mathematically_equal(lhs, rhs, ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value()) << "Expressions are not mathematically equivalent";
}

TEST(PolynomialGcd, SubresultantX10Minus1AndX8Minus1) {
    symbolic::CASContext ctx;
    auto p = parse_expr("x^10 - 1", ctx.arena());
    auto q = parse_expr("x^8 - 1", ctx.arena());
    ASSERT_TRUE(p.is_ok());
    ASSERT_TRUE(q.is_ok());
    
    Symbol x("x");
    auto gcd_res = polynomial_gcd(p.value(), q.value(), x, ctx);
    ASSERT_TRUE(gcd_res.is_ok());
    
    auto expected = parse_expr("x^2 - 1", ctx.arena());
    ASSERT_TRUE(expected.is_ok());
    
    expect_math_equal(gcd_res.value(), expected.value(), ctx);
}

TEST(PolynomialGcd, SubresultantWithConstants) {
    symbolic::CASContext ctx;
    auto p = parse_expr("2*x^2 - 2", ctx.arena());
    auto q = parse_expr("4*x - 4", ctx.arena());
    ASSERT_TRUE(p.is_ok());
    ASSERT_TRUE(q.is_ok());
    
    Symbol x("x");
    auto gcd_res = polynomial_gcd(p.value(), q.value(), x, ctx);
    ASSERT_TRUE(gcd_res.is_ok());
    
    auto expected = parse_expr("2*x - 2", ctx.arena());
    ASSERT_TRUE(expected.is_ok());

    expect_math_equal(gcd_res.value(), expected.value(), ctx);
}

TEST(PolynomialGcd, CertifiedBivariateCommonLinearFactor) {
    symbolic::CASContext ctx;
    auto p = parse_expr("x^2 - y^2", ctx.arena());
    auto q = parse_expr("x^2 + 2*x*y + y^2", ctx.arena());
    auto expected = parse_expr("x + y", ctx.arena());
    ASSERT_TRUE(p.is_ok());
    ASSERT_TRUE(q.is_ok());
    ASSERT_TRUE(expected.is_ok());

    auto gcd_res = polynomial_gcd_multivariate(p.value(), q.value(), ctx);
    ASSERT_TRUE(gcd_res.is_ok()) << gcd_res.error().message;
    expect_math_equal(gcd_res.value(), expected.value(), ctx);
}

TEST(PolynomialGcd, CertifiedTrivariateCommonLinearFactor) {
    symbolic::CASContext ctx;
    auto p = parse_expr("(x + y + z) * (x - z)", ctx.arena());
    auto q = parse_expr("(x + y + z) * (y + 2)", ctx.arena());
    auto expected = parse_expr("x + y + z", ctx.arena());
    ASSERT_TRUE(p.is_ok());
    ASSERT_TRUE(q.is_ok());
    ASSERT_TRUE(expected.is_ok());

    auto gcd_res = polynomial_gcd_multivariate(p.value(), q.value(), ctx);
    ASSERT_TRUE(gcd_res.is_ok()) << gcd_res.error().message;
    expect_math_equal(gcd_res.value(), expected.value(), ctx);
}

TEST(PolynomialGcd, DispatcherUsesCertifiedMultivariatePath) {
    symbolic::CASContext ctx;
    auto p = parse_expr("x^2 - y^2", ctx.arena());
    auto q = parse_expr("x^2 + 2*x*y + y^2", ctx.arena());
    auto expected = parse_expr("x + y", ctx.arena());
    ASSERT_TRUE(p.is_ok());
    ASSERT_TRUE(q.is_ok());
    ASSERT_TRUE(expected.is_ok());

    auto gcd_res = polynomial_gcd(p.value(), q.value(), Symbol("x"), ctx);
    ASSERT_TRUE(gcd_res.is_ok()) << gcd_res.error().message;
    expect_math_equal(gcd_res.value(), expected.value(), ctx);
}

} // namespace
} // namespace cas::algebra
