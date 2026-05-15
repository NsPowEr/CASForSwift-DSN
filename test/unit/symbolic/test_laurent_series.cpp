// Tests for laurent_series — rational expansion around a finite center.
// Verifies pole detection, negative-order coefficients (principal part), and
// positive-order Taylor coefficients all in one pass.

#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <string>

namespace cas::test {

namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

}  // namespace

class LaurentSeriesTest : public ::testing::Test {
protected:
    void SetUp() override { ctx = std::make_unique<symbolic::CASContext>(); }

    [[nodiscard]] ExprPtr E(const std::string& src) {
        auto r = parse_expr(src, ctx->arena());
        EXPECT_TRUE(r.is_ok()) << "Parse failed for: " << src
                               << "  err=" << (r.is_error() ? r.error().message : std::string{});
        return r.is_ok() ? r.value() : ExprPtr{};
    }

    // Test that the k-th coefficient (k counted from leading_order, ascending)
    // simplifies to `expected_text`.
    void expect_coeff(const calculus::LaurentExpansion& laurent, std::size_t idx, const std::string& expected_text) {
        ASSERT_LT(idx, laurent.coefficients.size())
            << "Index " << idx << " out of range; size = " << laurent.coefficients.size();
        ExprPtr expected = E(expected_text);
        ASSERT_TRUE(expected);
        auto actual = ctx->simplify(laurent.coefficients[idx]);
        ASSERT_TRUE(actual.is_ok());
        auto exp_s = ctx->simplify(expected);
        ASSERT_TRUE(exp_s.is_ok());
        auto eq = symbolic::mathematically_equal(actual.value(), exp_s.value(), *ctx);
        ASSERT_TRUE(eq.is_ok());
        EXPECT_TRUE(eq.value())
            << "Mismatch coeff idx=" << idx
            << "  got=" << debug_print(actual.value())
            << "  expected=" << debug_print(exp_s.value());
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

TEST_F(LaurentSeriesTest, SimplePoleOrder1) {
    // 1 / (x - 1) around x = 1.  Principal part is just 1/(x-1).
    auto laurent = calculus::laurent_series(E("1/(x-1)"), Symbol("x"), E("1"), 2U, *ctx);
    ASSERT_TRUE(laurent.is_ok()) << laurent.error().message;
    EXPECT_EQ(laurent.value().leading_order, -1);
    EXPECT_EQ(laurent.value().positive_order, 2U);
    ASSERT_EQ(laurent.value().coefficients.size(), 4U);  // m + pos + 1 = 1 + 2 + 1
    expect_coeff(laurent.value(), 0U, "1");
    expect_coeff(laurent.value(), 1U, "0");
    expect_coeff(laurent.value(), 2U, "0");
    expect_coeff(laurent.value(), 3U, "0");
}

TEST_F(LaurentSeriesTest, ProductOfLinearFactorsHasOrder1PoleAtZero) {
    // 1 / (x * (x - 1)) around x = 0.
    // u-series:  1/(x(x-1)) = -1/x - 1 - x - x^2 - ...
    auto laurent = calculus::laurent_series(E("1/(x*(x-1))"), Symbol("x"), E("0"), 3U, *ctx);
    ASSERT_TRUE(laurent.is_ok()) << laurent.error().message;
    EXPECT_EQ(laurent.value().leading_order, -1);
    EXPECT_EQ(laurent.value().positive_order, 3U);
    ASSERT_EQ(laurent.value().coefficients.size(), 5U);
    expect_coeff(laurent.value(), 0U, "-1");  // c_{-1}
    expect_coeff(laurent.value(), 1U, "-1");  // c_0
    expect_coeff(laurent.value(), 2U, "-1");  // c_1
    expect_coeff(laurent.value(), 3U, "-1");  // c_2
    expect_coeff(laurent.value(), 4U, "-1");  // c_3
}

TEST_F(LaurentSeriesTest, DoublePoleWithLinearNumerator) {
    // (x + 2) / (x - 1)^2 around x = 1.
    // Let u = x - 1.  N = u + 3, D = u^2.   Expansion: 3/u^2 + 1/u + 0 + 0 + ...
    auto laurent = calculus::laurent_series(E("(x+2)/(x-1)^2"), Symbol("x"), E("1"), 2U, *ctx);
    ASSERT_TRUE(laurent.is_ok()) << laurent.error().message;
    EXPECT_EQ(laurent.value().leading_order, -2);
    ASSERT_EQ(laurent.value().coefficients.size(), 5U);  // 2 + 2 + 1
    expect_coeff(laurent.value(), 0U, "3");  // c_{-2}
    expect_coeff(laurent.value(), 1U, "1");  // c_{-1}
    expect_coeff(laurent.value(), 2U, "0");  // c_0
    expect_coeff(laurent.value(), 3U, "0");  // c_1
    expect_coeff(laurent.value(), 4U, "0");  // c_2
}

TEST_F(LaurentSeriesTest, AnalyticAtCenterProducesTaylorOnly) {
    // 1 / (1 + x^2) around x = 0  =  1 - x^2 + x^4 - x^6 + ...
    auto laurent = calculus::laurent_series(E("1/(1 + x^2)"), Symbol("x"), E("0"), 4U, *ctx);
    ASSERT_TRUE(laurent.is_ok()) << laurent.error().message;
    EXPECT_EQ(laurent.value().leading_order, 0);
    ASSERT_EQ(laurent.value().coefficients.size(), 5U);  // 0 + 4 + 1
    expect_coeff(laurent.value(), 0U, "1");
    expect_coeff(laurent.value(), 1U, "0");
    expect_coeff(laurent.value(), 2U, "-1");
    expect_coeff(laurent.value(), 3U, "0");
    expect_coeff(laurent.value(), 4U, "1");
}

TEST_F(LaurentSeriesTest, ResidueAtSimplePoleMatchesLaurentMinusOneCoeff) {
    // Cross-check: residue at x = 1 of (x^2 + 1) / ((x - 1)*(x + 2)) is 2/3
    // via Laurent's c_{-1}.  This is the same data residue() already returns.
    auto laurent = calculus::laurent_series(E("(x^2 + 1)/((x - 1)*(x + 2))"), Symbol("x"), E("1"), 0U, *ctx);
    ASSERT_TRUE(laurent.is_ok()) << laurent.error().message;
    EXPECT_EQ(laurent.value().leading_order, -1);
    expect_coeff(laurent.value(), 0U, "2/3");
}

}  // namespace cas::test
