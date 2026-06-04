// F5.1 / B9-Task#22 — Tests per solve_risch_de_logarithmic_q (Bronstein cap.8).

#include "../../../src/calculus/calculus_internal.hpp"

#include "cas/algebra.hpp"
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

// Verifica round-trip: y' + f·y - g deve semplificare a 0.
[[nodiscard]] bool verify_de(
    ExprPtr y, ExprPtr f, ExprPtr g,
    const Symbol& var, symbolic::CASContext& ctx) {
    auto y_prime = diff(y, var, 1U, ctx);
    if (y_prime.is_error()) return false;
    AstArena& arena = ctx.arena();
    ExprPtr fy  = arena.make<Binary>(BinaryOp::Mul, f, y);
    ExprPtr lhs = arena.make<Binary>(BinaryOp::Add, y_prime.value(), fy);
    ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, lhs, g);
    auto delta_tog = algebra::together(delta, ctx);
    if (delta_tog.is_error()) return false;
    auto simp = ctx.simplify(delta_tog.value());
    if (simp.is_error()) return false;
    if (const auto* il = expr_cast<IntegerLit>(simp.value()))
        return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(simp.value()))
        return rl->numerator.is_zero();
    return false;
}

class RischLogarithmicDeTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};
};

// f = 0, u = x → θ = ln(x), θ' = 1/x.
// g = θ = ln(x): equazione  y' = ln(x).
// Soluzione: y = x·ln(x) - x. Verify D[x·ln(x) - x] = ln(x) + 1 - 1 = ln(x). ✓
TEST_F(RischLogarithmicDeTest, F0_GEqualsTheta_XLogXMinusX) {
    auto f = parse_expr("0", ctx.arena());
    auto g = parse_expr("ln(x)", ctx.arena());
    auto u = parse_expr("x", ctx.arena());
    auto res = solve_risch_de_logarithmic_q(f, g, u, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(verify_de(res.value(), f, g, x, ctx))
        << "round-trip failed for y = " << res.error().message;
}

// f = 0, u = x → θ = ln(x), θ' = 1/x.
// g = 1: y' = 1 → y = x. Trivial polynomial RHS, log-power-0 case.
TEST_F(RischLogarithmicDeTest, F0_G1_YEqualsX) {
    auto f = parse_expr("0", ctx.arena());
    auto g = parse_expr("1", ctx.arena());
    auto u = parse_expr("x", ctx.arena());
    auto res = solve_risch_de_logarithmic_q(f, g, u, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(verify_de(res.value(), f, g, x, ctx));
}

// f = 0, u = x → θ' = 1/x.
// g = θ²= ln(x)²: y' = ln(x)².
// Solution: ∫ln(x)² dx = x·ln(x)² - 2·x·ln(x) + 2x.  Verify round-trip.
TEST_F(RischLogarithmicDeTest, F0_GThetaSquared_QuadraticIntegral) {
    auto f = parse_expr("0", ctx.arena());
    auto g = parse_expr("ln(x)^2", ctx.arena());
    auto u = parse_expr("x", ctx.arena());
    auto res = solve_risch_de_logarithmic_q(f, g, u, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(verify_de(res.value(), f, g, x, ctx));
}

// f = 1, u = x.  Equazione: y' + y = ln(x).
// Costante f=1 in θ, polynomial-in-θ ansatz y(x) + y₁(x)·ln(x) con y₁ poly.
// Coefficienti: i=1: y₁' + y₁ = 1 → y₁ = 1. i=0: y₀' + y₀ = -1/x.
// Quest'ultima ha soluzione y₀ irrazionale (exp-integral); fallback Unimplemented atteso.
TEST_F(RischLogarithmicDeTest, F1_GLogX_FallsBackUnimplemented) {
    auto f = parse_expr("1", ctx.arena());
    auto g = parse_expr("ln(x)", ctx.arena());
    auto u = parse_expr("x", ctx.arena());
    auto res = solve_risch_de_logarithmic_q(f, g, u, x, ctx);
    // Risch DE per y₀ ammette esponenziale-integrale (non-elementare).
    // L'esito atteso è Unimplemented diagnostico (non risultato sbagliato).
    if (res.is_ok()) {
        // Se riesce, deve verificare round-trip — engine ha sorpreso.
        EXPECT_TRUE(verify_de(res.value(), f, g, x, ctx));
    } else {
        EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
    }
}

// f costante non nulla, g = 1, u = x.  Equazione y' + 2·y = 1 (y indip da θ).
// Soluzione: y = 1/2 (costante). Round-trip: 0 + 2·(1/2) = 1. ✓
TEST_F(RischLogarithmicDeTest, F2_G1_YEqualsHalf) {
    auto f = parse_expr("2", ctx.arena());
    auto g = parse_expr("1", ctx.arena());
    auto u = parse_expr("x", ctx.arena());
    auto res = solve_risch_de_logarithmic_q(f, g, u, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(verify_de(res.value(), f, g, x, ctx));
}

// f polinomiale in θ: f = ln(x) (caso non gestito).  Deve restituire
// Unimplemented diagnostico, non risultato sbagliato.
TEST_F(RischLogarithmicDeTest, FPolyInTheta_DiagnosticUnimplemented) {
    auto f = parse_expr("ln(x)", ctx.arena());
    auto g = parse_expr("1", ctx.arena());
    auto u = parse_expr("x", ctx.arena());
    auto res = solve_risch_de_logarithmic_q(f, g, u, x, ctx);
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
}

// ============================================================================
// Caso esponenziale θ = exp(u).
// ============================================================================

class RischExponentialDeTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};
};

// θ = exp(x), θ' = exp(x).  f = 0, g = 1: y' = 1 → y = x.
TEST_F(RischExponentialDeTest, F0_G1_YEqualsX) {
    auto f = parse_expr("0", ctx.arena());
    auto g = parse_expr("1", ctx.arena());
    auto u = parse_expr("x", ctx.arena());
    auto res = solve_risch_de_exponential_q(f, g, u, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(verify_de(res.value(), f, g, x, ctx));
}

// θ = exp(x).  f = 0, g = exp(x).  Grado 1: y_1' + 1·y_1 = 1 → y_1 = 1.
// Grado 0: y_0' = 0 → y_0 = 0.  Risultato: y = exp(x).
TEST_F(RischExponentialDeTest, F0_GExpX_YEqualsExpX) {
    auto f = parse_expr("0", ctx.arena());
    auto g = parse_expr("exp(x)", ctx.arena());
    auto u = parse_expr("x", ctx.arena());
    auto res = solve_risch_de_exponential_q(f, g, u, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(verify_de(res.value(), f, g, x, ctx));
}

// θ = exp(2x), θ' = 2·exp(2x), u' = 2.  f = 0, g = exp(2x).
// Grado 1: y_1' + 2·y_1 = 1 → y_1 = 1/2.
// y = (1/2)·exp(2x).  Round-trip: D[(1/2)·exp(2x)] = exp(2x). ✓
TEST_F(RischExponentialDeTest, F0_GExp2x_YEqualsHalfExp2x) {
    auto f = parse_expr("0", ctx.arena());
    auto g = parse_expr("exp(2*x)", ctx.arena());
    auto u = parse_expr("2*x", ctx.arena());
    auto res = solve_risch_de_exponential_q(f, g, u, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(verify_de(res.value(), f, g, x, ctx));
}

// f polinomiale in θ esponenziale: Unimplemented diagnostico atteso.
TEST_F(RischExponentialDeTest, FPolyInTheta_DiagnosticUnimplemented) {
    auto f = parse_expr("exp(x)", ctx.arena());
    auto g = parse_expr("1", ctx.arena());
    auto u = parse_expr("x", ctx.arena());
    auto res = solve_risch_de_exponential_q(f, g, u, x, ctx);
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
}

// ============================================================================
// Wiring di integrate() — esercita il fast-path che chiama
// solve_risch_de_{logarithmic,exponential}_q quando l'integranda è polinomiale
// in un singolo generatore trascendentale.
// ============================================================================

class RischDeWiringTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};

    [[nodiscard]] bool integrate_and_verify(const std::string& src) {
        AstArena& arena = ctx.arena();
        auto t = Lexer(src).tokenize();
        EXPECT_TRUE(t.is_ok()) << src << ": " << t.error().message;
        Parser p(t.value(), arena);
        auto e_res = p.parse();
        EXPECT_TRUE(e_res.is_ok()) << src;
        ExprPtr e = e_res.value();
        auto r = integrate(e, x, ctx);
        if (r.is_error()) return false;
        auto d_back = diff(r.value(), x, 1U, ctx);
        if (d_back.is_error()) return false;
        // Confronto post-expand_log per gestire identità ln(x²)=2·ln(x).
        ExprPtr Dy_exp = expand_log_args_via_factorization(d_back.value(), x, ctx);
        ExprPtr e_exp  = expand_log_args_via_factorization(e, x, ctx);
        ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, Dy_exp, e_exp);
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
};

// ∫ ln(x) dx = x·ln(x) - x
TEST_F(RischDeWiringTest, IntegralOfLnX) {
    EXPECT_TRUE(integrate_and_verify("ln(x)"));
}

// ∫ ln(x)^2 dx = x·ln(x)^2 - 2x·ln(x) + 2x
TEST_F(RischDeWiringTest, IntegralOfLnXSquared) {
    EXPECT_TRUE(integrate_and_verify("ln(x)^2"));
}

// ∫ exp(2x) dx = (1/2)·exp(2x)
TEST_F(RischDeWiringTest, IntegralOfExp2x) {
    EXPECT_TRUE(integrate_and_verify("exp(2*x)"));
}

// ∫ ln(x²) dx — generatore ln(x²), wiring cap.8 diretto.
// Round-trip valido perché D(2x·log(x²) - ?·x) = ... (verifica via diff).
TEST_F(RischDeWiringTest, IntegralOfLnXSquared_AsXSquaredArg) {
    EXPECT_TRUE(integrate_and_verify("ln(x^2)"));
}

// ∫ (ln(x) + ln(x²)) dx — generatori distinti senza expand,
// dipendenti dopo expand: log(x²) → 2·log(x), quindi integrand = 3·log(x).
// Cap.9 pre-processo + cap.8 wiring → 3·(x·log(x) - x).
TEST_F(RischDeWiringTest, IntegralOfLnX_PlusLnXSquared_FusedViaExpansion) {
    EXPECT_TRUE(integrate_and_verify("ln(x) + ln(x^2)"));
}

// ∫ ln(x·(x+1)) dx — log di prodotto.  Cap.9 espande in log(x)+log(x+1).
// Successivo wiring scopre 2 generatori distinti — skip wiring; ricade
// su Hermite/RT.  Test: integrate non deve crashare (può ritornare
// successo o Unimplemented diagnostico).
TEST_F(RischDeWiringTest, IntegralOfLnXTimesXPlusOne_NoCrash) {
    AstArena& arena = ctx.arena();
    auto t = Lexer("ln(x*(x+1))").tokenize();
    ASSERT_TRUE(t.is_ok());
    Parser p(t.value(), arena);
    auto e = p.parse();
    ASSERT_TRUE(e.is_ok());
    auto r = integrate(e.value(), x, ctx);
    // Esito: ok con verifica round-trip OPPURE Unimplemented esplicito.
    if (r.is_ok()) {
        EXPECT_TRUE(integrate_and_verify("ln(x*(x+1))"));
    } else {
        EXPECT_EQ(r.error().kind, CASErrorKind::Unimplemented);
    }
}

}  // namespace
}  // namespace cas::calculus
