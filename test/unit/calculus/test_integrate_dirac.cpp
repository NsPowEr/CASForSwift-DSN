// F5.9-pre — Tests sifting property Dirac per integrazione.

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

[[nodiscard]] bool same_after_simplify(
    ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr delta_e = arena.make<Binary>(BinaryOp::Sub, lhs, rhs);
    auto delta_tog = algebra::together(delta_e, ctx);
    if (delta_tog.is_error()) return false;
    auto delta_simp = ctx.simplify(delta_tog.value());
    if (delta_simp.is_error()) return false;
    if (const auto* il = expr_cast<IntegerLit>(delta_simp.value()))
        return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(delta_simp.value()))
        return rl->numerator.is_zero();
    return false;
}

class IntegrateDiracTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};
};

// ∫ δ(x) · x² dx = (x²)|_{x=0} / 1 = 0.
TEST_F(IntegrateDiracTest, SiftingAtZero) {
    auto e = parse_expr("delta(x) * x^2", ctx.arena());
    auto res = try_integrate_dirac_sifting(e, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("0", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// ∫ δ(x − 3) · x² dx = 3² = 9.
TEST_F(IntegrateDiracTest, SiftingAtThree) {
    auto e = parse_expr("delta(x - 3) * x^2", ctx.arena());
    auto res = try_integrate_dirac_sifting(e, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("9", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// ∫ δ(2x − 4) · (3x + 1) dx = (3·2 + 1) / |2| = 7/2.
//   Polo: 2x − 4 = 0 → x = 2.  c = 2 (positivo).  f(2) = 7.  Risultato 7/2.
TEST_F(IntegrateDiracTest, SiftingScaledArgument) {
    auto e = parse_expr("delta(2*x - 4) * (3*x + 1)", ctx.arena());
    auto res = try_integrate_dirac_sifting(e, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("7/2", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// δ solo (nessun f): ∫ δ(x) dx = 1 / |1| = 1.
TEST_F(IntegrateDiracTest, SiftingDeltaAlone) {
    auto e = parse_expr("delta(x)", ctx.arena());
    auto res = try_integrate_dirac_sifting(e, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("1", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// δ con argomento non lineare: arg = x² − 1 (caso fuori scope).
// Restituisce Unimplemented diagnostico.
TEST_F(IntegrateDiracTest, NonLinearArgumentUnimplemented) {
    auto e = parse_expr("delta(x^2 - 1)", ctx.arena());
    auto res = try_integrate_dirac_sifting(e, x, ctx);
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
}

}  // namespace
}  // namespace cas::calculus
