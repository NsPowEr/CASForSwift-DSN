#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace cas::algebra {
namespace {

[[nodiscard]] Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

TEST(FactorPolynomialTrager, FactorX4Plus1OverSqrt2) {
    symbolic::CASContext ctx;
    auto poly_res = parse_expr("x^4 + 1", ctx.arena());
    auto ext_res = parse_expr("sqrt(2)", ctx.arena());
    ASSERT_TRUE(poly_res.is_ok());
    ASSERT_TRUE(ext_res.is_ok());
    
    Symbol x("x");
    auto res = factor_polynomial(poly_res.value(), x, ctx, ext_res.value());
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    
    // factor(x^4 + 1, sqrt(2)) -> (x^2 - sqrt(2)x + 1)(x^2 + sqrt(2)x + 1)
    EXPECT_EQ(res.value().factors.size(), 2U);
    
    // Verifichiamo che il prodotto dei fattori (espanso) sia uguale al polinomio originale
    std::vector<ExprPtr> factors;
    for (const auto& f : res.value().factors) {
        for (unsigned int i = 0; i < f.multiplicity; ++i) {
            factors.push_back(f.factor);
        }
    }
    
    ExprPtr product = factors[0];
    for (size_t i = 1; i < factors.size(); ++i) {
        product = ctx.arena().make<Product>(std::vector<ExprPtr>{product, factors[i]});
    }
    
    auto expanded = expand(product, ctx);
    ASSERT_TRUE(expanded.is_ok());
    
    auto simplified_orig = ctx.simplify(poly_res.value());
    ASSERT_TRUE(simplified_orig.is_ok());
    
    // In CAS, structural equality might fail if simplify doesn't handle sqrt(2)^2 = 2 perfectly in expand.
    // But expand() should handle it.
    auto mat_eq = mathematically_equal(expanded.value(), simplified_orig.value(), ctx);
    EXPECT_TRUE(mat_eq.is_ok() && mat_eq.value());
}

TEST(FactorPolynomialTrager, FactorX2Minus2OverSqrt2) {
    symbolic::CASContext ctx;
    auto poly_res = parse_expr("x^2 - 2", ctx.arena());
    auto ext_res = parse_expr("sqrt(2)", ctx.arena());
    ASSERT_TRUE(poly_res.is_ok());
    ASSERT_TRUE(ext_res.is_ok());
    
    Symbol x("x");
    auto res = factor_polynomial(poly_res.value(), x, ctx, ext_res.value());
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    
    // x^2 - 2 -> (x - sqrt(2))(x + sqrt(2))
    EXPECT_EQ(res.value().factors.size(), 2U);
}

} // namespace
} // namespace cas::algebra
