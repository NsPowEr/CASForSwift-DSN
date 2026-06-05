// F5.9 / Task #20 — Tests Elliptic integrals K, E, Π, F.
//
// Definizioni (Abramowitz-Stegun §17.4–§17.7, mod-k convention):
//   K(k) = ∫_0^{π/2} dθ / √(1 - k²·sin²θ)        (complete first kind)
//   E(k) = ∫_0^{π/2} √(1 - k²·sin²θ) dθ           (complete second kind)
//   Π(n, k) = ∫_0^{π/2} dθ / ((1 - n·sin²θ)√(1 - k²·sin²θ))
//                                                  (complete third kind)
//   F(φ, k) = ∫_0^φ dθ / √(1 - k²·sin²θ)          (incomplete first kind)

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/calculus.hpp"
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

class EllipticIntegralsTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol k{"k"};
};

// K(0) = π/2.
TEST_F(EllipticIntegralsTest, EllipticKAtZero) {
    auto e = parse_expr("EllipticK(0)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("pi/2", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// E(0) = π/2.
TEST_F(EllipticIntegralsTest, EllipticEAtZero) {
    auto e = parse_expr("EllipticE(0)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("pi/2", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// E(1) = 1.
TEST_F(EllipticIntegralsTest, EllipticEAtOne) {
    auto e = parse_expr("EllipticE(1)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("1", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// F(φ, 0) = φ.
TEST_F(EllipticIntegralsTest, EllipticFAtZeroK) {
    auto e = parse_expr("EllipticF(phi, 0)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("phi", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// Π(0, k) = K(k).
TEST_F(EllipticIntegralsTest, EllipticPiAtNZero) {
    auto e = parse_expr("EllipticPi(0, k)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("EllipticK(k)", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// dK/dk struttura check: deve coinvolgere EllipticE.
TEST_F(EllipticIntegralsTest, EllipticKDerivativeInvolvesE) {
    auto e = parse_expr("EllipticK(k)", ctx.arena());
    auto d = diff(e, k, 1U, ctx);
    ASSERT_TRUE(d.is_ok()) << d.error().message;
    bool has_E = false;
    std::function<void(ExprPtr)> walk = [&](ExprPtr x) {
        if (!x || has_E) return;
        if (const auto* fc = expr_cast<FuncCall>(x);
            fc && fc->func_id == BuiltinOp::EllipticE) {
            has_E = true; return;
        }
        if (const auto* bin = expr_cast<Binary>(x)) {
            walk(bin->left); walk(bin->right); return;
        }
        if (const auto* un = expr_cast<Unary>(x)) walk(un->operand);
        if (const auto* prod = expr_cast<Product>(x)) {
            for (ExprPtr y : prod->factors) walk(y); return;
        }
        if (const auto* sum = expr_cast<Sum>(x)) {
            for (ExprPtr y : sum->terms) walk(y); return;
        }
    };
    walk(d.value());
    EXPECT_TRUE(has_E) << "expected dK/dk to involve EllipticE";
}

// dE/dk = (E(k) - K(k))/k.  Verifica strutturale: contiene sia E che K.
TEST_F(EllipticIntegralsTest, EllipticEDerivativeInvolvesBoth) {
    auto e = parse_expr("EllipticE(k)", ctx.arena());
    auto d = diff(e, k, 1U, ctx);
    ASSERT_TRUE(d.is_ok()) << d.error().message;
    bool has_E = false, has_K = false;
    std::function<void(ExprPtr)> walk = [&](ExprPtr x) {
        if (!x) return;
        if (const auto* fc = expr_cast<FuncCall>(x)) {
            if (fc->func_id == BuiltinOp::EllipticE) has_E = true;
            if (fc->func_id == BuiltinOp::EllipticK) has_K = true;
            return;
        }
        if (const auto* bin = expr_cast<Binary>(x)) {
            walk(bin->left); walk(bin->right); return;
        }
        if (const auto* un = expr_cast<Unary>(x)) walk(un->operand);
        if (const auto* prod = expr_cast<Product>(x)) {
            for (ExprPtr y : prod->factors) walk(y); return;
        }
        if (const auto* sum = expr_cast<Sum>(x)) {
            for (ExprPtr y : sum->terms) walk(y); return;
        }
    };
    walk(d.value());
    EXPECT_TRUE(has_E && has_K);
}

}  // namespace
}  // namespace cas::calculus
