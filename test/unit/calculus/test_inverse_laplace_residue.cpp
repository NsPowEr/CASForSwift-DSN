// F5.8 / Task #17 — Tests per inverse Laplace via Bronstein residue.

#include "../../../src/calculus/calculus_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas::calculus {
namespace {

ExprPtr parse_expr(const std::string& src, AstArena& arena) {
    auto t = Lexer(src).tokenize();
    EXPECT_TRUE(t.is_ok()) << src;
    Parser p(t.value(), arena);
    auto r = p.parse();
    EXPECT_TRUE(r.is_ok()) << src;
    return r.value();
}

// Verifica via t = 1: confronta numericamente sostituendo t con valore concreto.
[[nodiscard]] bool same_after_subst_t(
    ExprPtr lhs, ExprPtr rhs, const Symbol& t_sym, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr t_val = arena.make<IntegerLit>(BigInt(0));
    auto lhs_s = ctx.substitute(lhs, t_sym, t_val);
    auto rhs_s = ctx.substitute(rhs, t_sym, t_val);
    if (lhs_s.is_error() || rhs_s.is_error()) return false;
    ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, lhs_s.value(), rhs_s.value());
    auto delta_tog = algebra::together(delta, ctx);
    if (delta_tog.is_error()) return false;
    auto delta_simp = ctx.simplify(delta_tog.value());
    if (delta_simp.is_error()) return false;
    if (const auto* il = expr_cast<IntegerLit>(delta_simp.value()))
        return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(delta_simp.value()))
        return rl->numerator.is_zero();
    return false;
}

class InvLaplaceResidueTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol s{"s"};
    Symbol t{"t"};
};

// L⁻¹{1/s}(t) = 1.  Verify at t=0: both sides = 1.
TEST_F(InvLaplaceResidueTest, OneOverS) {
    auto F = parse_expr("1/s", ctx.arena());
    auto res = inverse_laplace_residue_q(F, s, t, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("1", ctx.arena());
    EXPECT_TRUE(same_after_subst_t(res.value(), expected, t, ctx));
}

// L⁻¹{1/(s-2)}(t) = e^(2t).  Verify at t=0: both sides = 1.
TEST_F(InvLaplaceResidueTest, OneOverSMinus2) {
    auto F = parse_expr("1/(s-2)", ctx.arena());
    auto res = inverse_laplace_residue_q(F, s, t, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("exp(2*t)", ctx.arena());
    EXPECT_TRUE(same_after_subst_t(res.value(), expected, t, ctx));
}

// L⁻¹{1/(s·(s-1))}(t) = e^t - 1.  Verify at t=0: -1 + 1 = 0.
TEST_F(InvLaplaceResidueTest, OneOverSTimesSMinus1) {
    auto F = parse_expr("1/(s*(s-1))", ctx.arena());
    auto res = inverse_laplace_residue_q(F, s, t, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("exp(t) - 1", ctx.arena());
    EXPECT_TRUE(same_after_subst_t(res.value(), expected, t, ctx));
}

// L⁻¹{1/(s²-1)}(t) = (1/2)·(e^t - e^(-t)) = sinh(t).  Verify at t=0: 0.
TEST_F(InvLaplaceResidueTest, OneOverSSquaredMinus1) {
    auto F = parse_expr("1/(s^2-1)", ctx.arena());
    auto res = inverse_laplace_residue_q(F, s, t, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("(1/2)*(exp(t) - exp(-t))", ctx.arena());
    EXPECT_TRUE(same_after_subst_t(res.value(), expected, t, ctx));
}

// Wiring check: inverse_laplace_transform (pattern-based) deve cadere
// sul fallback residue per F(s) NON nelle coppie elementari pre-esistenti.
// Esempio: 1/(s³-s) = 1/(s(s-1)(s+1)) ha 3 polos {0, 1, -1}; pattern
// table classico richiede partial fractions, ma il fallback residue
// chiude direttamente.
TEST_F(InvLaplaceResidueTest, WiringFallback_OneOverSCubedMinusS) {
    auto F = parse_expr("1/(s^3 - s)", ctx.arena());
    auto res = inverse_laplace_transform(F, s, t, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    // L⁻¹{1/(s(s²-1))} = cosh(t) - 1.  Verifica at t=0: cosh(0)-1 = 0.
    auto expected = parse_expr("(1/2)*(exp(t) + exp(-t)) - 1", ctx.arena());
    EXPECT_TRUE(same_after_subst_t(res.value(), expected, t, ctx));
}

}  // namespace
}  // namespace cas::calculus
