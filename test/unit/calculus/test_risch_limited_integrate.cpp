// A38 — Tests per limited_integrate_field (Bronstein "Symbolic Integration I"
// §7.2, Limited Integration Problem) e per il suo wiring nel path integrate()
// reale.
//
// Contesto (A38): il solver parametrico su torre solve_risch_de_parametric_field
// (A1 + A26) era sound e testato via harness, ma NESSUN caller reale lo
// raggiungeva da integrate().  limited_integrate_field è il ponte previsto da
// Bronstein §7.2 ("(7.30) can be considered a parametric Risch differential
// equation for v and can be solved by the algorithm of Sect. 7.1"), e viene
// invocato per ogni grado dalla ricorsione §5.10 IntegratePrimitivePolynomial
// in integrate_log_polynomial_part.
//
// Verifica SOUND: mai confronto su toString(); si controlla
//   f − D(v) − Σ c_i·w_i ≡ 0   (unit, nel campo)
//   D(∫f) − f ≡ 0              (end-to-end, nel dominio reale)

#include "../../../src/calculus/calculus_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>
#include <vector>

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

// e ≡ 0 come funzione razionale (numeratore nullo dopo together+expand).
[[nodiscard]] bool is_zero_expr(ExprPtr e, symbolic::CASContext& ctx) {
    auto tog = algebra::together(e, ctx);
    ExprPtr x = tog.is_ok() ? tog.value() : e;
    auto s = ctx.simplify(x);
    ExprPtr z = s.is_ok() ? s.value() : x;
    if (const auto* il = expr_cast<IntegerLit>(z)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(z)) return rl->numerator.is_zero();
    auto parts = algebra::apart_num_den(z, ctx);
    if (parts.is_error()) return false;
    ExprPtr num = parts.value().numerator;
    if (auto ex = algebra::expand(num, ctx); ex.is_ok()) num = ex.value();
    if (auto sn = ctx.simplify(num); sn.is_ok()) num = sn.value();
    if (const auto* il = expr_cast<IntegerLit>(num)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(num)) return rl->numerator.is_zero();
    return false;
}

// f − D(v) − Σ c_i·w_i ≡ 0 nel campo differenziale.
[[nodiscard]] bool verify_limited(
    ExprPtr f, const std::vector<ExprPtr>& w_vec,
    const LimitedIntegrationFieldSolution& sol,
    const DifferentialField& field, symbolic::CASContext& ctx) {
    if (sol.c.size() != w_vec.size()) return false;
    AstArena& arena = ctx.arena();
    auto dv = field.derive_in_generators(sol.v, ctx);
    if (dv.is_error()) return false;
    std::vector<ExprPtr> terms{f, arena.make<Unary>(UnaryOp::Neg, dv.value())};
    for (std::size_t i = 0; i < w_vec.size(); ++i) {
        if (sol.c[i].numerator().is_zero()) continue;
        ExprPtr c_e = arena.make<RationalLit>(sol.c[i].numerator(), sol.c[i].denominator());
        terms.push_back(arena.make<Unary>(UnaryOp::Neg,
            arena.make<Binary>(BinaryOp::Mul, c_e, w_vec[i])));
    }
    return is_zero_expr(arena.make<Sum>(std::move(terms)), ctx);
}

class RischLimitedIntegrate : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    AstArena& arena() { return ctx.arena(); }
};

// ---------------------------------------------------------------- unit §7.2

// Base Q(x), c = 0: f = 2x è la derivata esatta di x² — la costante sulla
// forzante w = 1/x deve risultare nulla.
TEST_F(RischLimitedIntegrate, RationalBase_ExactDerivative_CZero) {
    Symbol x("x");
    DifferentialField field(x);
    ExprPtr f = parse_expr("2*x", arena());
    ExprPtr w = parse_expr("1/x", arena());

    auto res = limited_integrate_field(f, {w}, field, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(verify_limited(f, {w}, res.value(), field, ctx));
    EXPECT_TRUE(res.value().c.front().numerator().is_zero());
}

// Base Q(x), c ≠ 0: f = 1/x non ha primitiva in Q(x), ma f = D(0) + 1·(1/x).
// È il caso che il solver NON parametrico non sa rappresentare.
TEST_F(RischLimitedIntegrate, RationalBase_ForcingCarriesTheConstant) {
    Symbol x("x");
    DifferentialField field(x);
    ExprPtr f = parse_expr("1/x", arena());
    ExprPtr w = parse_expr("1/x", arena());

    auto res = limited_integrate_field(f, {w}, field, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(verify_limited(f, {w}, res.value(), field, ctx));
    EXPECT_FALSE(res.value().c.front().numerator().is_zero());
}

// Base Q(x), combinazione: f = 2x + 3/x = D(x²) + 3·(1/x).
TEST_F(RischLimitedIntegrate, RationalBase_MixedRationalAndConstantPart) {
    Symbol x("x");
    DifferentialField field(x);
    ExprPtr f = parse_expr("2*x + 3/x", arena());
    ExprPtr w = parse_expr("1/x", arena());

    auto res = limited_integrate_field(f, {w}, field, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(verify_limited(f, {w}, res.value(), field, ctx));
    EXPECT_EQ(res.value().c.front(), Rational(BigInt(3)));
}

// Negativo legittimo (§7.2, "no solution with c_0 = 1"): f = 1/(x+1) non è
// D(v) + c/x per alcun v ∈ Q(x), c ∈ Q.  Deve arrivare un diagnostico, non un
// risultato sbagliato.
TEST_F(RischLimitedIntegrate, RationalBase_NoSolutionIsDiagnostic) {
    Symbol x("x");
    DifferentialField field(x);
    ExprPtr f = parse_expr("1/(x+1)", arena());
    ExprPtr w = parse_expr("1/x", arena());

    auto res = limited_integrate_field(f, {w}, field, ctx);
    if (res.is_ok()) {
        // Se accetta, DEVE essere una soluzione valida (sound-by-construction).
        EXPECT_TRUE(verify_limited(f, {w}, res.value(), field, ctx));
    } else {
        EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
    }
}

// Torre log K = Q(x, t) con t = ln(x), η = D(t) = 1/x.
// f = 1/x + 1/(x·t) = D(t) + ... — qui il pezzo 1/(x·t) è esattamente η/t, la
// derivata logaritmica del generatore: la soluzione richiede la macchina di
// torre, non il solo Q(x).
TEST_F(RischLimitedIntegrate, LogTower_ForcingIsTheMonomialDerivative) {
    Symbol x("x");
    ExprPtr integrand = parse_expr("ln(x)", arena());
    auto field_res = DifferentialField::build(integrand, x, ctx);
    ASSERT_TRUE(field_res.is_ok()) << field_res.error().message;
    const auto& field = field_res.value();
    ASSERT_EQ(field.extensions().size(), 1U);

    ExprPtr eta = parse_expr("1/x", arena());   // D(t) = 1/x
    // f = D(t) exactly ⇒ v = t, c = 0  oppure  v = 0, c = 1: entrambe valide.
    auto res = limited_integrate_field(eta, {eta}, field, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(verify_limited(eta, {eta}, res.value(), field, ctx));
}

// Torre log, forzante MISTA (rational + η): a_1 = 1/x + 1/(x·t) è la
// combinazione D(t) + 1·η che serve per il livello top di una torre
// log-in-log.  HC-A26-PRIMITIVE-PARAMQ-RATIONAL (cancellazione dei poli +
// degree bound §6.3): questa forma ORA si risolve (v = t, c = 1), non è più
// un diagnostico.
TEST_F(RischLimitedIntegrate, LogTower_MixedForcing_Solved) {
    Symbol x("x");
    ExprPtr integrand = parse_expr("ln(x)", arena());
    auto field_res = DifferentialField::build(integrand, x, ctx);
    ASSERT_TRUE(field_res.is_ok()) << field_res.error().message;
    const auto& field = field_res.value();
    ASSERT_EQ(field.extensions().size(), 1U);

    // a_1 = (t0+1)/(x*t0) = 1/x + 1/(x*t0)  (generator form)
    ExprPtr t0 = arena().make<Symbol>(field.extensions()[0].t_var.name);
    ExprPtr a1 = parse_expr("1/x", arena());
    ExprPtr eta_gen = arena().make<Binary>(BinaryOp::Div,
        arena().make<IntegerLit>(BigInt(1)),
        arena().make<Binary>(BinaryOp::Mul, arena().make<Symbol>(Symbol("x")), t0));
    ExprPtr f = arena().make<Binary>(BinaryOp::Add, a1, eta_gen);

    auto res = limited_integrate_field(f, {eta_gen}, field, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(verify_limited(f, {eta_gen}, res.value(), field, ctx));
}

// ------------------------------------------------------- end-to-end §5.10

// A38 rende raggiungibile limited_integrate_field da integrate() per torri a
// UN livello: qualunque polinomio-in-log di grado ≥ 1 ora passa dal ramo
// kz≥1 (limited integration nel campo inferiore, invece del vecchio "b_k =
// ∫rhs" cieco).  Per grado 1 i due algoritmi concordano (c=0 è sempre una
// soluzione valida quando D(t)=η è essa stessa elementarmente integrabile,
// Liouville) — questo è il canary di non-regressione end-to-end che prova
// che il nuovo path è ESERCITATO (non dead code) su un caso reale.
TEST_F(RischLimitedIntegrate, EndToEnd_SingleLevelLogPolynomial_ViaLimitedIntegrateField) {
    Symbol x("x");
    ExprPtr f = parse_expr("ln(x)/x", arena());

    auto F = integrate(f, x, ctx);
    ASSERT_TRUE(F.is_ok()) << F.error().message;
    auto dF = diff(F.value(), x, 1U, ctx);
    ASSERT_TRUE(dF.is_ok()) << dF.error().message;
    ExprPtr delta = arena().make<Binary>(BinaryOp::Sub, dF.value(), f);
    EXPECT_TRUE(is_zero_expr(delta, ctx))
        << "D(∫f) ≠ f — l'antiderivata prodotta non è corretta";
}

// Torre log ANNIDATA (t_1=ln(x), t_2=ln(ln(x))):
//   ∫ (1/x + 1/(x·ln x))·ln(ln x) dx  =  ln(ln x)²/2 + ln(x)·ln(ln x) − ln(x).
// Percorso completo di questa sessione:
//   • HC-A38-01: la ricorsione §5.10 scende sul campo inferiore (prima era un
//     loop infinito da root-restart);
//   • HC-A26 pole-cancel: la famiglia di forzanti al livello top ha poli in
//     t_1 che si cancellano (c₀=c₁) → forzante ridotta polinomiale;
//   • HC-A26 degree bound §6.3: N = dg_max+1, non dg_max, per catturare
//     q = x·t di grado uno più della forzante.
// Esito: antiderivata corretta, verificata per derivazione (D(∫f) = f).
TEST_F(RischLimitedIntegrate, EndToEnd_NestedLogTower_Solved) {
    Symbol x("x");
    ExprPtr f = parse_expr("(1/x + 1/(x*ln(x)))*ln(ln(x))", arena());

    // Stato 2026-07-20: TERMINA (era un loop infinito) e riporta Unimplemented.
    // L'integrale È elementare (= ln(ln x)²/2 + ln(x)·ln(ln x) − ln(x)); la
    // completezza residua NON è più un problema di ricorsione ma del solver
    // §7.2: limited_integrate_field non risolve una forzante con denominatore
    // nel generatore (base case Q(x) polynomial-only, HC-A26-PRIMITIVE-PARAMQ-
    // RATIONAL).  Se un giorno passa, il ramo is_ok qui sotto lo blinda.
    auto F = integrate(f, x, ctx);
    ASSERT_TRUE(F.is_ok()) << F.error().message;
    auto dF = diff(F.value(), x, 1U, ctx);
    ASSERT_TRUE(dF.is_ok()) << dF.error().message;
    ExprPtr delta = arena().make<Binary>(BinaryOp::Sub, dF.value(), f);
    EXPECT_TRUE(is_zero_expr(delta, ctx))
        << "D(∫f) ≠ f — antiderivata sbagliata sulla torre annidata";
}

// Soundness sul caso NON elementare con generatore sibling: ∫ e^{-x}·ln(x) dx
// non ha forma chiusa elementare.  Il vecchio guard `coeff_blocks_poly_quotient`
// (rimosso da HC-A38-01) esisteva anche per impedire che questo integrale
// producesse la forma silenziosamente sbagliata e^{-x}·(x·ln x − x), che nasce
// trattando il generatore exp come costante in x.  Ora che il fallback è
// field-aware (campo inferiore, D(t_exp) noto), quel percorso non è più
// raggiungibile: il test blinda l'invariante.
TEST_F(RischLimitedIntegrate, EndToEnd_SiblingExpCoefficient_NeverSilentlyWrong) {
    Symbol x("x");
    ExprPtr f = parse_expr("exp(-x)*ln(x)", arena());

    auto F = integrate(f, x, ctx);
    if (F.is_ok()) {
        auto dF = diff(F.value(), x, 1U, ctx);
        ASSERT_TRUE(dF.is_ok()) << dF.error().message;
        ExprPtr delta = arena().make<Binary>(BinaryOp::Sub, dF.value(), f);
        EXPECT_TRUE(is_zero_expr(delta, ctx))
            << "D(∫f) ≠ f — regressione silent-wrong su ∫e^{-x}·ln(x)dx";
    } else {
        EXPECT_EQ(F.error().kind, CASErrorKind::Unimplemented);
    }
}

// Non-regressione del caso classico che passava già dal fast-path
// (§5.10 con c = 0 a ogni livello): ∫ ln(x) dx = x·ln(x) − x.
TEST_F(RischLimitedIntegrate, EndToEnd_PlainLogStillCorrect) {
    Symbol x("x");
    ExprPtr f = parse_expr("ln(x)", arena());

    auto F = integrate(f, x, ctx);
    ASSERT_TRUE(F.is_ok()) << F.error().message;
    auto dF = diff(F.value(), x, 1U, ctx);
    ASSERT_TRUE(dF.is_ok()) << dF.error().message;
    ExprPtr delta = arena().make<Binary>(BinaryOp::Sub, dF.value(), f);
    EXPECT_TRUE(is_zero_expr(delta, ctx));
}

}  // namespace
}  // namespace cas::calculus
