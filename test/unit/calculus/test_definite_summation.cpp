// F5.7 sub-block 0 — Definite hypergeometric summation via Gosper.
//
// The Newton-Leibniz finite-calculus identity
//   Σ_{k=a}^{b} t(k) = S(b+1) − S(a),  where S(k+1) − S(k) = t(k),
// turns Gosper antidifferences into definite-sum closed forms.

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

#include <string>

namespace cas::test {
namespace {

ExprPtr parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    EXPECT_TRUE(tokens.is_ok()) << tokens.error().message;
    Parser parser(tokens.value(), arena);
    auto res = parser.parse();
    EXPECT_TRUE(res.is_ok()) << res.error().message;
    return res.value();
}

bool simplifies_to_zero(ExprPtr a, ExprPtr b, symbolic::CASContext& ctx) {
    auto diff = ctx.arena().make<Binary>(BinaryOp::Sub, a, b);
    auto expanded = algebra::expand(diff, ctx);
    if (expanded.is_error()) return false;
    auto s = ctx.simplify(expanded.value());
    if (s.is_error()) return false;
    if (const auto* il = expr_cast<IntegerLit>(s.value())) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(s.value())) return rl->numerator.is_zero();
    return false;
}

}  // namespace

class DefiniteSummationTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol k{"k"};
    Symbol n{"n"};
};

// Σ_{k=1}^{n} 1 = n.
TEST_F(DefiniteSummationTest, ConstantOne) {
    auto term = parse_expr("1", ctx.arena());
    auto lo = parse_expr("1", ctx.arena());
    auto hi = parse_expr("n", ctx.arena());
    auto res = calculus::sum(term, k, lo, hi, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("n", ctx.arena());
    EXPECT_TRUE(simplifies_to_zero(res.value(), ctx.simplify(expected).value(), ctx));
}

// Σ_{k=1}^{n} k = n·(n+1)/2.
TEST_F(DefiniteSummationTest, ArithmeticSeriesFirstN) {
    auto term = parse_expr("k", ctx.arena());
    auto lo = parse_expr("1", ctx.arena());
    auto hi = parse_expr("n", ctx.arena());
    auto res = calculus::sum(term, k, lo, hi, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("n*(n+1)/2", ctx.arena());
    EXPECT_TRUE(simplifies_to_zero(res.value(), ctx.simplify(expected).value(), ctx));
}

// Σ_{k=1}^{n} k² = n·(n+1)·(2n+1)/6  (closed via Gosper antidifference of k²
// after the F5.7-GOSPER-K2-NORMALIZATION rescaling).
TEST_F(DefiniteSummationTest, SumOfSquares) {
    auto term = parse_expr("k^2", ctx.arena());
    auto lo = parse_expr("1", ctx.arena());
    auto hi = parse_expr("n", ctx.arena());
    auto res = calculus::sum(term, k, lo, hi, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("n*(n+1)*(2*n+1)/6", ctx.arena());
    EXPECT_TRUE(simplifies_to_zero(res.value(), ctx.simplify(expected).value(), ctx));
}

// Σ_{k=1}^{n} 1/(k·(k+1)) = n/(n+1)  (telescoping; Gosper-summable rational).
TEST_F(DefiniteSummationTest, TelescopingRational) {
    auto term = parse_expr("1/(k*(k+1))", ctx.arena());
    auto lo = parse_expr("1", ctx.arena());
    auto hi = parse_expr("n", ctx.arena());
    auto res = calculus::sum(term, k, lo, hi, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("n/(n+1)", ctx.arena());
    // Bring both through `together` first to collapse the partial-fraction
    // representation Gosper produces.
    auto res_tog = algebra::together(res.value(), ctx);
    ASSERT_TRUE(res_tog.is_ok());
    auto exp_tog = algebra::together(ctx.simplify(expected).value(), ctx);
    ASSERT_TRUE(exp_tog.is_ok());
    EXPECT_TRUE(simplifies_to_zero(res_tog.value(), exp_tog.value(), ctx));
}

}  // namespace cas::test
