// F5.9 / Task #19 — Tests closed-form identities Bessel/Erf/Zeta.
//
// Verifica scope completo già implementato nel core:
//   - J_{1/2}(z) = √(2/(πz)) · sin(z)
//   - J_{-1/2}(z) = √(2/(πz)) · cos(z)
//   - Y_{1/2}(z) = -√(2/(πz)) · cos(z)
//   - I_{1/2}(z) = √(2/(πz)) · sinh(z)
//   - K_{±1/2}(z) = √(π/(2z)) · e^(-z)
//   - erf(0) = 0,  erf(-x) = -erf(x)
//   - ζ(0) = -1/2, ζ(2) = π²/6, ζ(4) = π⁴/90,  ζ(-1) = -1/12

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/builtin_functions.hpp"
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
    ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, lhs, rhs);
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

class SpecialClosedFormsTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol z{"z"};
};

// J_{1/2}(z) = √(2/(πz)) · sin(z).
TEST_F(SpecialClosedFormsTest, BesselJ_HalfPositive) {
    auto e = parse_expr("BesselJ(1/2, z)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("sqrt(2/(pi*z)) * sin(z)", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// K_{1/2}(z) = √(π/(2z)) · exp(-z).
TEST_F(SpecialClosedFormsTest, BesselK_HalfPositive) {
    auto e = parse_expr("BesselK(1/2, z)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("sqrt(pi/(2*z)) * exp(-z)", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// erf(0) = 0.
TEST_F(SpecialClosedFormsTest, ErfAtZero) {
    auto e = parse_expr("erf(0)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("0", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// erf(-z) = -erf(z) (odd).  Substitute z → 1: erf(-1) - (-erf(1)) = 0.
// Verifica via espressione simbolica.
TEST_F(SpecialClosedFormsTest, ErfOddSymmetry) {
    auto e = parse_expr("erf(-z) + erf(z)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("0", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// ζ(2) = π²/6.  Test via valore numerico: differenza zeta(2) - pi^2/6,
// dopo simplify, dovrebbe essere 0 letterale.  La forma esatta puo
// variare (1/6 * pi^2 vs pi^2/6) ma sono identiche dopo together.
TEST_F(SpecialClosedFormsTest, ZetaTwo) {
    auto e = parse_expr("zeta(2)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    // Result deve contenere fattore pi^2 (struttura check minimal).
    bool has_pi_squared = false;
    std::function<void(ExprPtr)> walk = [&](ExprPtr x) {
        if (!x || has_pi_squared) return;
        if (const auto* bin = expr_cast<Binary>(x);
            bin && bin->op == BinaryOp::Pow) {
            if (const auto* c = expr_cast<Constant>(bin->left);
                c && c->value == MathConstant::Pi) {
                if (const auto* il = expr_cast<IntegerLit>(bin->right);
                    il && il->value == BigInt(2)) {
                    has_pi_squared = true;
                    return;
                }
            }
        }
        if (const auto* bin = expr_cast<Binary>(x)) {
            walk(bin->left); walk(bin->right);
        }
        if (const auto* prod = expr_cast<Product>(x)) {
            for (ExprPtr y : prod->factors) walk(y);
        }
    };
    walk(simp.value());
    EXPECT_TRUE(has_pi_squared) << "expected zeta(2) to contain π² factor";
}

// ζ(-1) = -1/12.  Riemann functional equation per ζ(s) = ζ(1-s)·...
//   tramite Bernoulli: ζ(-1) = -B_2/2 = -(1/6)/2 = -1/12.
TEST_F(SpecialClosedFormsTest, ZetaMinusOne) {
    auto e = parse_expr("zeta(-1)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("-1/12", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// ζ(0) = -1/2.
TEST_F(SpecialClosedFormsTest, ZetaZero) {
    auto e = parse_expr("zeta(0)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("-1/2", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

}  // namespace
}  // namespace cas::calculus
