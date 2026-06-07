#include <gtest/gtest.h>
#include "cas/algebra/hensel.hpp"
#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/parser.hpp"
#include "cas/lexer.hpp"
#include "algebra/algebra_internal.hpp"

using namespace cas;
using namespace cas::algebra;
using namespace cas::symbolic;

TEST(MultivariateHenselTest, LinearHenselStepExact) {
    CASContext ctx;
    auto parse_poly = [&](const std::string& str) {
        auto tokens = Lexer(str).tokenize();
        EXPECT_TRUE(tokens.is_ok());
        Parser parser(tokens.value(), ctx.arena());
        auto expr = parser.parse();
        EXPECT_TRUE(expr.is_ok());
        auto poly = parse_multivariate_polynomial(expr.value(), ctx);
        EXPECT_TRUE(poly.is_ok());
        return poly.value();
    };

    auto x = expr_cast<Symbol>(parse_poly("x").to_univariate_coefficients(Symbol("x"), ctx).value()[1]) ? Symbol("x") : Symbol("x"); // hack to just get a symbol

    auto f = parse_poly("x^2 + 2*x + 2");
    auto g0 = parse_poly("x");
    auto h0 = parse_poly("x + 1");
    auto s = parse_poly("-1");
    auto t = parse_poly("1");

    BezoutCoeffs bezout;
    bezout.s = s;
    bezout.t = t;

    Ideal ideal;
    ideal.point.vars.push_back(x);
    ideal.point.values.push_back(BigInt(0));
    ideal.degree = 2;
    
    auto res = linear_hensel_step(f, g0, h0, bezout, ideal, ctx);
    ASSERT_TRUE(res.is_ok());
    
    auto [g1, h1] = res.value();
    
    auto g1_expr = multivariate_to_expr(g1, ctx);
    ASSERT_TRUE(g1_expr.is_ok());
    auto h1_expr = multivariate_to_expr(h1, ctx);
    ASSERT_TRUE(h1_expr.is_ok());

    auto g1_simp = simplify_expr(g1_expr.value(), ctx);
    auto h1_simp = simplify_expr(h1_expr.value(), ctx);
    
    auto exp_g1 = parse_poly("x + 2");
    auto exp_h1 = parse_poly("x");
    
    auto diff_g = subtract_exprs(g1_simp.value(), multivariate_to_expr(exp_g1, ctx).value(), ctx).value();
    auto diff_h = subtract_exprs(h1_simp.value(), multivariate_to_expr(exp_h1, ctx).value(), ctx).value();
    
    EXPECT_TRUE(is_zero_expr(simplify_expr(diff_g, ctx).value()));
    EXPECT_TRUE(is_zero_expr(simplify_expr(diff_h, ctx).value()));
}

TEST(MultivariateHenselTest, BivariateHenselStepExact) {
    CASContext ctx;
    auto parse_poly = [&](const std::string& str) {
        auto tokens = Lexer(str).tokenize();
        EXPECT_TRUE(tokens.is_ok());
        Parser parser(tokens.value(), ctx.arena());
        auto expr = parser.parse();
        EXPECT_TRUE(expr.is_ok());
        auto poly = parse_multivariate_polynomial(expr.value(), ctx);
        EXPECT_TRUE(poly.is_ok());
        return poly.value();
    };

    auto y = Symbol("y");

    // f(x, y) = x^2 - x + x*y - y + y^2
    auto f = parse_poly("x^2 - x + x*y - y + y^2");
    // g0 = x
    auto g0 = parse_poly("x");
    // h0 = x - 1
    auto h0 = parse_poly("x - 1");
    // s = 1, t = -1 (s*g0 + t*h0 = x - (x - 1) = 1)
    auto s = parse_poly("1");
    auto t = parse_poly("-1");

    BezoutCoeffs bezout;
    bezout.s = s;
    bezout.t = t;

    Ideal ideal;
    ideal.point.vars.push_back(y); // eval point in y
    ideal.point.values.push_back(BigInt(0)); // y = 0
    ideal.degree = 2; // mod I^2
    
    auto res = linear_hensel_step(f, g0, h0, bezout, ideal, ctx);
    ASSERT_TRUE(res.is_ok());
    
    auto [g1, h1] = res.value();
    
    auto g1_expr = multivariate_to_expr(g1, ctx);
    ASSERT_TRUE(g1_expr.is_ok());
    auto h1_expr = multivariate_to_expr(h1, ctx);
    ASSERT_TRUE(h1_expr.is_ok());

    auto g1_simp = simplify_expr(g1_expr.value(), ctx);
    auto h1_simp = simplify_expr(h1_expr.value(), ctx);
    
    // We calculated:
    // pg = y - y^2  =>  g1 = x + y - y^2
    // ph = y^2      =>  h1 = x - 1 + y^2
    auto exp_g1 = parse_poly("x + y - y^2");
    auto exp_h1 = parse_poly("x - 1 + y^2");
    
    auto diff_g = subtract_exprs(g1_simp.value(), multivariate_to_expr(exp_g1, ctx).value(), ctx).value();
    auto diff_h = subtract_exprs(h1_simp.value(), multivariate_to_expr(exp_h1, ctx).value(), ctx).value();
    
    EXPECT_TRUE(is_zero_expr(simplify_expr(diff_g, ctx).value()));
    EXPECT_TRUE(is_zero_expr(simplify_expr(diff_h, ctx).value()));
}
