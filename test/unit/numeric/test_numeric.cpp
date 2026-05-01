#include "cas/numeric.hpp"
#include "cas/symbolic.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include <gtest/gtest.h>
#include <cmath>

namespace cas::numeric {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

} // namespace

TEST(NumericTest, EvaluatesBaseExpressions) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("1 + 2 * 3", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    
    auto res = eval(expr.value());
    ASSERT_TRUE(res.is_ok());
    EXPECT_DOUBLE_EQ(res.value(), 7.0);
}

TEST(NumericTest, EvaluatesConstants) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("pi", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    
    auto res = eval(expr.value());
    ASSERT_TRUE(res.is_ok());
    EXPECT_NEAR(res.value(), M_PI, 1e-15);
}

TEST(NumericTest, EvaluatesFunctions) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("sin(pi/2) + exp(1)", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    
    auto res = eval(expr.value());
    ASSERT_TRUE(res.is_ok());
    EXPECT_NEAR(res.value(), 1.0 + std::exp(1.0), 1e-15);
}

TEST(NumericTest, BisectionRootFinding) {
    symbolic::CASContext ctx;
    // f(x) = x^2 - 2, root at sqrt(2) approx 1.4142
    auto expr = parse_expr("x^2 - 2", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    
    auto res = solve_numeric_bisection(expr.value(), "x", 0.0, 2.0);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NEAR(res.value(), std::sqrt(2.0), 1e-8);
}

TEST(NumericTest, NewtonRaphsonRootFinding) {
    symbolic::CASContext ctx;
    // f(x) = x^2 - 2, derivative 2x
    auto expr = parse_expr("x^2 - 2", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    
    auto res = solve_numeric_newton(expr.value(), "x", 1.5, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NEAR(res.value(), std::sqrt(2.0), 1e-9);
}

TEST(NumericTest, AdaptiveSimpsonIntegration) {
    symbolic::CASContext ctx;
    // integrate x^2 from 0 to 1 = 1/3
    auto expr = parse_expr("x^2", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    
    auto res = integrate_numeric(expr.value(), "x", 0.0, 1.0);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NEAR(res.value(), 1.0/3.0, 1e-7);
}

TEST(NumericTest, RK4OdeSolver) {
    symbolic::CASContext ctx;
    // dy/dt = y, y(0)=1 -> y(t) = exp(t). y(1) = e
    auto expr = parse_expr("y", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    
    auto res = solve_ode_rk4(expr.value(), "t", "y", 0.0, 1.0, 1.0, 0.1);
    ASSERT_TRUE(res.is_ok());
    ASSERT_FALSE(res.value().empty());
    EXPECT_NEAR(res.value().back().y, std::exp(1.0), 1e-4);
}

TEST(NumericTest, IntegrationWithSimplifierN) {
    symbolic::CASContext ctx;
    auto tokens = Lexer("N(sqrt(2))").tokenize();
    Parser parser(tokens.value(), ctx.arena());
    auto expr = parser.parse();
    
    auto simplified = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified.is_ok());
    
    const auto* dec = expr_cast<DecimalLit>(simplified.value());
    ASSERT_NE(dec, nullptr);
    EXPECT_NEAR(dec->to_double(), std::sqrt(2.0), 1e-15);
}

} // namespace cas::numeric
