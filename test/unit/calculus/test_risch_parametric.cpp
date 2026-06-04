// F5.1 / B9-Task#21 — Tests per solve_risch_de_parametric_q + limited_integration_q.
//
// Riferimento: Bronstein, "Symbolic Integration I", §7.1/§7.2.

#include "../../../src/calculus/calculus_internal.hpp"

#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas::calculus {
namespace {

ExprPtr parse_expr(const std::string& s, AstArena& arena) {
    auto t = Lexer(s).tokenize();
    EXPECT_TRUE(t.is_ok()) << s << ": " << t.error().message;
    Parser p(t.value(), arena);
    auto r = p.parse();
    EXPECT_TRUE(r.is_ok()) << s << ": " << r.error().message;
    return r.value();
}

class ParametricRischDeTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};
};

// Caso banale: m=1, f=0, g_1 = 1.  Cerca y polinomio + c con y' = c·1.
// Spazio nullo dim ≥ 1 (es. y = x, c = 1; o y = 0, c = 0 banale).
TEST_F(ParametricRischDeTest, F0_G1Constant_BasisContainsLinearY) {
    auto f = parse_expr("0", ctx.arena());
    auto g1 = parse_expr("1", ctx.arena());
    auto res = solve_risch_de_parametric_q(f, {g1}, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    // Almeno una soluzione (y=x, c=1) deve esistere: c=1 ≠ 0.
    bool found_nontrivial = false;
    for (const auto& sol : res.value()) {
        if (!sol.c.empty() && !sol.c[0].numerator().is_zero()) {
            found_nontrivial = true;
            break;
        }
    }
    EXPECT_TRUE(found_nontrivial)
        << "Expected a basis vector with c_1 ≠ 0 representing y = x, c = 1.";
}

// f = 0, g_1 = 2x.  Allora y' = c · 2x → y = c·x², c arbitrario.
// Base 1-dim (c = 1, y = x²).
TEST_F(ParametricRischDeTest, F0_G2x_YEqualsXSquared) {
    auto f = parse_expr("0", ctx.arena());
    auto g1 = parse_expr("2*x", ctx.arena());
    auto res = solve_risch_de_parametric_q(f, {g1}, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    bool found_match = false;
    for (const auto& sol : res.value()) {
        if (sol.c.size() != 1U) continue;
        // Cerca soluzione (y = c·x², c) con c ≠ 0.
        if (sol.c[0].numerator().is_zero()) continue;
        // y deve essere proporzionale a x² (no termini lineari o costanti).
        // Verifichiamo via differentiation: D(y) = c·2x.
        auto dy = cas::calculus::diff(sol.y, x, 1U, ctx);
        ASSERT_TRUE(dy.is_ok());
        auto expected = parse_expr("2*x", ctx.arena());
        auto c_e = ctx.arena().make<RationalLit>(
            sol.c[0].numerator(), sol.c[0].denominator());
        auto c_g = ctx.arena().make<Binary>(BinaryOp::Mul, c_e, expected);
        auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, dy.value(), c_g);
        auto simp = ctx.simplify(delta);
        ASSERT_TRUE(simp.is_ok());
        const auto* lit = expr_cast<IntegerLit>(simp.value());
        if (lit && lit->value.is_zero()) {
            found_match = true;
            break;
        }
    }
    EXPECT_TRUE(found_match) << "Expected basis vector encoding y' = c · 2x.";
}

// f = 1, g_1 = exp-shifted; ma qui forziamo solo polinomi → testiamo
// f = 1, g_1 = 1, equazione y' + y = c.  Soluzione (y=c, c) sempre.
TEST_F(ParametricRischDeTest, F1_G1Constant_ConstantSolution) {
    auto f = parse_expr("1", ctx.arena());
    auto g1 = parse_expr("1", ctx.arena());
    auto res = solve_risch_de_parametric_q(f, {g1}, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    bool found = false;
    for (const auto& sol : res.value()) {
        if (sol.c.size() != 1U) continue;
        if (sol.c[0].numerator().is_zero()) continue;
        // y deve coincidere con c (costante).
        auto c_e = ctx.arena().make<RationalLit>(
            sol.c[0].numerator(), sol.c[0].denominator());
        auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, sol.y, c_e);
        auto simp = ctx.simplify(delta);
        ASSERT_TRUE(simp.is_ok());
        const auto* lit = expr_cast<IntegerLit>(simp.value());
        if (lit && lit->value.is_zero()) { found = true; break; }
    }
    EXPECT_TRUE(found) << "Expected (y=c, c) constant solution.";
}

// m=2 con g_1 = 1, g_2 = x.  y' + 0·y = c1 + c2·x.
// Soluzione generica: y = c1·x + c2·x²/2.  Base 2-dim.
TEST_F(ParametricRischDeTest, F0_TwoForcings_TwoBasisVectors) {
    auto f = parse_expr("0", ctx.arena());
    auto g1 = parse_expr("1", ctx.arena());
    auto g2 = parse_expr("x", ctx.arena());
    auto res = solve_risch_de_parametric_q(f, {g1, g2}, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_GE(res.value().size(), 2U)
        << "Expected null-space dim ≥ 2 for two independent forcings.";
}

// Limited integration Q[x]: f = 3x² → g = x³ esattamente.
TEST_F(ParametricRischDeTest, LimitedIntegration_PurePolynomial) {
    auto f = parse_expr("3*x^2", ctx.arena());
    auto res = limited_integration_q(f, {}, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("x^3", ctx.arena());
    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, res.value().g, expected);
    auto simp = ctx.simplify(delta);
    ASSERT_TRUE(simp.is_ok());
    const auto* lit = expr_cast<IntegerLit>(simp.value());
    ASSERT_NE(lit, nullptr) << "Expected literal 0, got non-literal residual.";
    EXPECT_TRUE(lit->value.is_zero());
    EXPECT_TRUE(res.value().c.empty());
}

TEST_F(ParametricRischDeTest, LimitedIntegration_LinearWithResidueBasis) {
    auto f = parse_expr("4*x + 1", ctx.arena());
    auto h1 = parse_expr("x", ctx.arena());
    auto res = limited_integration_q(f, {h1}, x, ctx);
    ASSERT_TRUE(res.is_ok());
    // In Q[x] puro la decomposizione f = g' + Σc_i h_i con c_i ≥ 0 lascia
    // c = [0] e g = 2x² + x.
    auto expected_g = parse_expr("2*x^2 + x", ctx.arena());
    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, res.value().g, expected_g);
    auto simp = ctx.simplify(delta);
    ASSERT_TRUE(simp.is_ok());
    const auto* lit = expr_cast<IntegerLit>(simp.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_TRUE(lit->value.is_zero());
    ASSERT_EQ(res.value().c.size(), 1U);
    EXPECT_TRUE(res.value().c[0].numerator().is_zero());
}

}  // namespace
}  // namespace cas::calculus
