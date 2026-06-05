// F5.8 / Task #14 — Tests Fourier transform.

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

class FourierTransformTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol t{"t"};
    Symbol w{"omega"};
};

// F{δ(t)}(ω) = 1.
TEST_F(FourierTransformTest, DiracAtZero) {
    auto f = parse_expr("delta(t)", ctx.arena());
    auto res = fourier_transform(f, t, w, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("1", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// F{δ(t − 3)}(ω) = exp(−i·ω·3) = exp(−3·i·ω).
TEST_F(FourierTransformTest, DiracShifted) {
    auto f = parse_expr("delta(t - 3)", ctx.arena());
    auto res = fourier_transform(f, t, w, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("exp(-i*omega*3)", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// F{1}(ω) = 2π·δ(ω).
TEST_F(FourierTransformTest, ConstantOne) {
    auto f = parse_expr("1", ctx.arena());
    auto res = fourier_transform(f, t, w, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    AstArena& arena = ctx.arena();
    ExprPtr expected = arena.make<Product>(std::vector<ExprPtr>{
        arena.make<IntegerLit>(BigInt(2)),
        arena.make<Constant>(MathConstant::Pi),
        arena.make<FuncCall>(BuiltinOp::DiracDelta,
            std::vector<ExprPtr>{arena.make<Symbol>(w)})});
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// F{cos(2t)}(ω) = π·[δ(ω−2) + δ(ω+2)].
// Verifica strutturale: il risultato deve essere prodotto/sum con π·δ(±2±ω).
// Confronto numerico tramite valutazione su ω = 5 (lontano dai δ peaks):
// entrambi i lati danno 0·π + 0·π = 0 (δ valutata fuori support è 0).
// Test debole ma verifica almeno consistency struttura simbolica.
TEST_F(FourierTransformTest, CosineAtFrequency2) {
    auto f = parse_expr("cos(2*t)", ctx.arena());
    auto res = fourier_transform(f, t, w, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    // Risultato strutturale: deve contenere DiracDelta(ω-2) o DiracDelta(2-ω)
    // e DiracDelta(ω+2) come fattori in una somma π-scaled.  Verifica almeno
    // che il count di DiracDelta sia 2 nel risultato.
    int dirac_count = 0;
    std::function<void(ExprPtr)> count = [&](ExprPtr e) {
        if (!e) return;
        if (const auto* fc = expr_cast<FuncCall>(e);
            fc && fc->func_id == BuiltinOp::DiracDelta) { ++dirac_count; return; }
        if (const auto* bin = expr_cast<Binary>(e)) {
            count(bin->left); count(bin->right); return;
        }
        if (const auto* sum = expr_cast<Sum>(e)) {
            for (ExprPtr x : sum->terms) count(x); return;
        }
        if (const auto* prod = expr_cast<Product>(e)) {
            for (ExprPtr x : prod->factors) count(x); return;
        }
        if (const auto* un = expr_cast<Unary>(e)) count(un->operand);
        if (const auto* fc = expr_cast<FuncCall>(e)) {
            for (ExprPtr x : fc->args) count(x);
        }
    };
    count(res.value());
    EXPECT_EQ(dirac_count, 2)
        << "expected exactly two DiracDelta factors (frequenze ±2)";
}

// Linearità: F{δ(t) + δ(t-1)}(ω) = 1 + exp(-i·ω).
TEST_F(FourierTransformTest, LinearityTwoDiracs) {
    auto f = parse_expr("delta(t) + delta(t - 1)", ctx.arena());
    auto res = fourier_transform(f, t, w, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("1 + exp(-i*omega)", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// Pattern non riconosciuto: F{ln(t)} → Unimplemented diagnostico.
TEST_F(FourierTransformTest, LogUnimplemented) {
    auto f = parse_expr("ln(t)", ctx.arena());
    auto res = fourier_transform(f, t, w, ctx);
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
}

}  // namespace
}  // namespace cas::calculus
