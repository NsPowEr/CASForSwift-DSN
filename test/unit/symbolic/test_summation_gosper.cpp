#include <gtest/gtest.h>
#include "cas/symbolic.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "../../../src/symbolic/summation_gosper.hpp"

using namespace cas;
using namespace cas::symbolic;

class GosperSumTest : public ::testing::Test {
protected:
    CASContext ctx;
    Symbol k{"k"};

    ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
};

TEST_F(GosperSumTest, DISABLED_Polynomial1) {
    auto term = parse("1");
    auto res = gosper_sum(term, k, ctx);
    ASSERT_TRUE(res.is_ok());
    ASSERT_TRUE(res.value().has_value());
    
    auto s = res.value().value();
    
    auto k_plus_1 = parse("k + 1");
    auto s_next = substitute(s, k, k_plus_1, ctx).value();
    // Let's do it manually to be safe
    auto sub_val = ctx.arena().make<Binary>(BinaryOp::Sub, s_next, s);
    auto diff_simp = simplify(sub_val, ctx).value();
    
    EXPECT_TRUE(mathematically_equal(diff_simp, term, ctx).value());
}

TEST_F(GosperSumTest, DISABLED_PolynomialK) {
    auto term = parse("k");
    auto res = gosper_sum(term, k, ctx);
    ASSERT_TRUE(res.is_ok());
    ASSERT_TRUE(res.value().has_value());
    
    auto s = res.value().value();
    
    auto k_plus_1 = parse("k + 1");
    auto s_next = substitute(s, k, k_plus_1, ctx).value();
    auto sub_val = ctx.arena().make<Binary>(BinaryOp::Sub, s_next, s);
    auto diff_simp = simplify(sub_val, ctx).value();
    auto diff_expand = algebra::expand(diff_simp, ctx).value();
    
    EXPECT_TRUE(mathematically_equal(diff_expand, term, ctx).value());
}

TEST_F(GosperSumTest, DISABLED_RationalShift) {
    auto term = parse("1 / (k * (k+1))");
    auto res = gosper_sum(term, k, ctx);
    ASSERT_TRUE(res.is_ok());
    ASSERT_TRUE(res.value().has_value());
    
    auto s = res.value().value();
    
    auto k_plus_1 = parse("k + 1");
    auto s_next = substitute(s, k, k_plus_1, ctx).value();
    auto sub_val = ctx.arena().make<Binary>(BinaryOp::Sub, s_next, s);
    auto diff_simp = simplify(sub_val, ctx).value();
    auto together_diff = algebra::together(diff_simp, ctx).value();
    auto term_simp = algebra::together(simplify(term, ctx).value(), ctx).value();
    
    EXPECT_TRUE(mathematically_equal(together_diff, term_simp, ctx).value());
}

TEST_F(GosperSumTest, DISABLED_NotHypergeometricSummable) {
    auto term = parse("1 / (k^2 + 1)");
    auto res = gosper_sum(term, k, ctx);
    ASSERT_TRUE(res.is_ok());
    ASSERT_FALSE(res.value().has_value());
}