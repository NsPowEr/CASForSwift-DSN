#include <gtest/gtest.h>
#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "algebra/algebra_internal.hpp"
#include <string>
#include <vector>

using namespace cas;
using namespace cas::algebra;

namespace {

MultivariatePolynomial parse_mpoly(const std::string& input, symbolic::CASContext& ctx) {
    Lexer lex_obj(input);
    auto tokens = lex_obj.tokenize();
    EXPECT_TRUE(tokens.is_ok());
    Parser parse_obj(tokens.value(), ctx.arena());
    auto expr = parse_obj.parse();
    EXPECT_TRUE(expr.is_ok());
    auto mpoly_res = parse_multivariate_polynomial(expr.value(), ctx);
    EXPECT_TRUE(mpoly_res.is_ok());
    return mpoly_res.value();
}

bool mpoly_equal(const MultivariatePolynomial& poly, const std::string& expected, symbolic::CASContext& ctx) {
    auto expected_poly = parse_mpoly(expected, ctx);
    auto expected_expr = multivariate_to_expr(expected_poly, ctx).value();
    auto poly_expr = multivariate_to_expr(poly, ctx).value();
    auto e1 = ctx.simplify(expected_expr).value();
    auto e2 = ctx.simplify(poly_expr).value();
    return structural_equal(e1, e2);
}

} // namespace

TEST(MultivariateSquareFree, Basic) {
    symbolic::CASContext ctx;
    
    // Test 1: Primitive polynomial with square-free factors
    // P = (x + y)^2 * (x - y) = x^3 + x^2 y - x y^2 - y^3
    auto poly = parse_mpoly("x^3 + x^2*y - x*y^2 - y^3", ctx);
    
    auto factors_res = square_free_factorize_multivariate(poly, ctx);
    ASSERT_TRUE(factors_res.is_ok());
    auto factors = factors_res.value();
    
    // We expect two factors: (x - y) with multiplicity 1, and (x + y) with multiplicity 2
    ASSERT_EQ(factors.size(), 2);
    
    bool found_1 = false;
    bool found_2 = false;
    for (const auto& f : factors) {
        if (f.multiplicity == 1) {
            EXPECT_TRUE(mpoly_equal(f.factor, "x - y", ctx));
            found_1 = true;
        } else if (f.multiplicity == 2) {
            EXPECT_TRUE(mpoly_equal(f.factor, "x + y", ctx));
            found_2 = true;
        }
    }
    EXPECT_TRUE(found_1 && found_2);
}

TEST(MultivariateSquareFree, WithContent) {
    symbolic::CASContext ctx;
    
    // Test 2: Polynomial with content
    // P = y^2 * (x + 1)^2 * (x + 2) = y^2 * (x^3 + 4x^2 + 5x + 2)
    // = x^3 y^2 + 4 x^2 y^2 + 5 x y^2 + 2 y^2
    auto poly = parse_mpoly("x^3*y^2 + 4*x^2*y^2 + 5*x*y^2 + 2*y^2", ctx);
    
    auto factors_res = square_free_factorize_multivariate(poly, ctx);
    ASSERT_TRUE(factors_res.is_ok());
    auto factors = factors_res.value();
    
    // We expect factors: (x + 2)^1, y^2, (x + 1)^2
    // That means: (x + 2) mult 1
    // y mult 2
    // (x + 1) mult 2
    
    ASSERT_EQ(factors.size(), 3);
    
    int count_mult_1 = 0;
    int count_mult_2 = 0;
    
    for (const auto& f : factors) {
        if (f.multiplicity == 1) {
            EXPECT_TRUE(mpoly_equal(f.factor, "x + 2", ctx));
            count_mult_1++;
        } else if (f.multiplicity == 2) {
            EXPECT_TRUE(mpoly_equal(f.factor, "x + 1", ctx) || mpoly_equal(f.factor, "y", ctx));
            count_mult_2++;
        }
    }
    EXPECT_EQ(count_mult_1, 1);
    EXPECT_EQ(count_mult_2, 2);
}

TEST(MultivariateSquareFree, IntegerContent) {
    symbolic::CASContext ctx;
    
    // Test 3: Integer content
    // P = 12 * (x + y)^3
    auto poly = parse_mpoly("12*x^3 + 36*x^2*y + 36*x*y^2 + 12*y^3", ctx);
    
    auto factors_res = square_free_factorize_multivariate(poly, ctx);
    ASSERT_TRUE(factors_res.is_ok());
    auto factors = factors_res.value();
    
    ASSERT_EQ(factors.size(), 2);
    
    for (const auto& f : factors) {
        if (f.multiplicity == 1) {
            EXPECT_TRUE(mpoly_equal(f.factor, "12", ctx));
        } else if (f.multiplicity == 3) {
            EXPECT_TRUE(mpoly_equal(f.factor, "x + y", ctx));
        } else {
            FAIL() << "Unexpected multiplicity " << f.multiplicity;
        }
    }
}
