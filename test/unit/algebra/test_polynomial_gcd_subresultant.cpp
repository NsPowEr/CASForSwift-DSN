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
    
    // We can't use structural_equal directly if the order of terms or leading coefficient sign differs
    // But since it's integer poly GCD, it should be normalized.
    // Let's use simplify(gcd - expected) == 0.
    auto diff = ctx.arena().make<Binary>(BinaryOp::Sub, gcd_res.value(), expected.value());
    auto expanded = expand(diff, ctx);
    ASSERT_TRUE(expanded.is_ok());
    auto simp = ctx.simplify(expanded.value());
    ASSERT_TRUE(simp.is_ok());
    
    const auto* il = expr_cast<IntegerLit>(simp.value());
    if (!il) {
        const auto* rl = expr_cast<RationalLit>(simp.value());
        if (rl) {
            EXPECT_TRUE(rl->numerator.is_zero());
        } else {
            FAIL() << "Result is not a constant zero after expansion and simplification.";
        }
    } else {
        EXPECT_TRUE(il->value.is_zero());
    }
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

    auto diff = ctx.arena().make<Binary>(BinaryOp::Sub, gcd_res.value(), expected.value());    auto expanded = expand(diff, ctx);
    ASSERT_TRUE(expanded.is_ok());
    auto simp = ctx.simplify(expanded.value());
    ASSERT_TRUE(simp.is_ok());
    
    const auto* il = expr_cast<IntegerLit>(simp.value());
    if (!il) {
        const auto* rl = expr_cast<RationalLit>(simp.value());
        if (rl) {
            EXPECT_TRUE(rl->numerator.is_zero());
        } else {
            FAIL() << "Result is not a constant zero after expansion and simplification.";
        }
    } else {
        EXPECT_TRUE(il->value.is_zero());
    }
}

} // namespace
} // namespace cas::algebra
