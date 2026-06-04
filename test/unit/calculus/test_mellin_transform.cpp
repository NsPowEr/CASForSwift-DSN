// F5.8 / Task #15 — Tests per Mellin transform.

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

// Semantic check via substitution + simplify: r^(-s)·Γ(s) e Γ(s)/r^s sono
// strutturalmente diversi nel simplifier ma uguali matematicamente; lo
// verifichiamo sostituendo s con valore concreto e comparando.
[[nodiscard]] bool same_after_subst_s(
    ExprPtr lhs, ExprPtr rhs, const Symbol& s_sym, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr s_val = arena.make<IntegerLit>(BigInt(2));
    auto lhs_s = ctx.substitute(lhs, s_sym, s_val);
    auto rhs_s = ctx.substitute(rhs, s_sym, s_val);
    if (lhs_s.is_error() || rhs_s.is_error()) return false;
    return same_after_simplify(lhs_s.value(), rhs_s.value(), ctx);
}

class MellinTransformTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol t{"t"};
    Symbol s{"s"};
};

// M{exp(-t)}(s) = Γ(s).
TEST_F(MellinTransformTest, ExpMinusT) {
    auto f = parse_expr("exp(-t)", ctx.arena());
    auto res = mellin_transform(f, t, s, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("gamma(s)", ctx.arena());
    EXPECT_TRUE(same_after_subst_s(res.value(), expected, s, ctx));
}

// M{exp(-3·t)}(s) = Γ(s)/3^s.
TEST_F(MellinTransformTest, ExpMinus3T) {
    auto f = parse_expr("exp(-3*t)", ctx.arena());
    auto res = mellin_transform(f, t, s, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("gamma(s)/3^s", ctx.arena());
    EXPECT_TRUE(same_after_subst_s(res.value(), expected, s, ctx));
}

// M{t^2 · exp(-t)}(s) = Γ(s+2).
TEST_F(MellinTransformTest, TSquaredExpMinusT) {
    auto f = parse_expr("t^2 * exp(-t)", ctx.arena());
    auto res = mellin_transform(f, t, s, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("gamma(s+2)", ctx.arena());
    EXPECT_TRUE(same_after_subst_s(res.value(), expected, s, ctx));
}

// M{t^3 · exp(-2·t)}(s) = Γ(s+3) / 2^(s+3).
TEST_F(MellinTransformTest, TCubedExpMinus2T) {
    auto f = parse_expr("t^3 * exp(-2*t)", ctx.arena());
    auto res = mellin_transform(f, t, s, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("gamma(s+3)/2^(s+3)", ctx.arena());
    EXPECT_TRUE(same_after_subst_s(res.value(), expected, s, ctx));
}

// Linearità: M{exp(-t) + exp(-2t)}(s) = Γ(s) + Γ(s)/2^s.
TEST_F(MellinTransformTest, LinearityTwoExps) {
    auto f = parse_expr("exp(-t) + exp(-2*t)", ctx.arena());
    auto res = mellin_transform(f, t, s, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("gamma(s) + gamma(s)/2^s", ctx.arena());
    EXPECT_TRUE(same_after_subst_s(res.value(), expected, s, ctx));
}

// (1+t)^(-a) → Γ(s)·Γ(a-s)/Γ(a).
TEST_F(MellinTransformTest, OnePlusTToMinusA) {
    AstArena& arena = ctx.arena();
    Symbol a_sym{"a"};
    ExprPtr a_e = arena.make<Symbol>(a_sym);
    // build (1+t)^(-a)
    ExprPtr one_plus_t = arena.make<Binary>(BinaryOp::Add,
        arena.make<IntegerLit>(BigInt(1)), arena.make<Symbol>(t));
    ExprPtr neg_a = arena.make<Unary>(UnaryOp::Neg, a_e);
    ExprPtr f = arena.make<Binary>(BinaryOp::Pow, one_plus_t, neg_a);
    auto res = mellin_transform(f, t, s, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    // expected: gamma(s)·gamma(a-s)/gamma(a).
    ExprPtr expected = arena.make<Binary>(BinaryOp::Div,
        arena.make<Binary>(BinaryOp::Mul,
            arena.make<FuncCall>(BuiltinOp::Gamma, std::vector<ExprPtr>{arena.make<Symbol>(s)}),
            arena.make<FuncCall>(BuiltinOp::Gamma, std::vector<ExprPtr>{
                arena.make<Binary>(BinaryOp::Sub, a_e, arena.make<Symbol>(s))})),
        arena.make<FuncCall>(BuiltinOp::Gamma, std::vector<ExprPtr>{a_e}));
    EXPECT_TRUE(same_after_subst_s(res.value(), expected, s, ctx));
}

// Pattern non riconosciuto: log(t) → Unimplemented.
TEST_F(MellinTransformTest, LogTUnimplemented) {
    auto f = parse_expr("ln(t)", ctx.arena());
    auto res = mellin_transform(f, t, s, ctx);
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
}

}  // namespace
}  // namespace cas::calculus
