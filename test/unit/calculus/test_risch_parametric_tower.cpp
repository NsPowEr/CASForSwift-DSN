// A26 — Tests per solve_risch_de_parametric_field (parametric Risch DE su torre
// differenziale K(x, t_1, ..., t_n)).  Bronstein, "Symbolic Integration I",
// §6.5 / §7.4 / §8.4 (PolyRischDE parametrico).
//
// Scopo (A26): la funzione `solve_risch_de_parametric_field` era *entry-orphan*
// (raggiungibile solo dalla propria ricorsione `ext_idx-1`, mai invocata
// dall'esterno) e completamente non testata.  Questo file la rende raggiungibile
// e testata in isolamento, verificando la correttezza per back-substitution
// SOUND nel campo differenziale:  D(y) + f·y  ≡  Σ_i c_i·g_i.
//
// Inoltre costruisce deterministicamente un input con deg_t(f) > 0 (il ramo
// "non-cancellation" tuttora Unimplemented) = riproduzione concreta del gap A1,
// che 16 integrandi-torre costruiti a mano non riuscivano a raggiungere.

#include "../../../src/calculus/calculus_internal.hpp"

#include "cas/algebra.hpp"
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

// Verifica SOUND nel campo: D(y) + f·y - Σ c_i·g_i  semplifica a 0.
// Tutti gli argomenti sono in *forma generatore* (polinomi/razionali in x e t_i);
// la derivazione usa field.derive che conosce D(t_i).
[[nodiscard]] bool verify_field_de(
    ExprPtr y, ExprPtr f, const std::vector<ExprPtr>& g_vec,
    const std::vector<Rational>& c, const DifferentialField& field,
    symbolic::CASContext& ctx) {
    if (c.size() != g_vec.size()) return false;
    AstArena& arena = ctx.arena();

    auto dy = field.derive(y, ctx);
    if (dy.is_error()) return false;

    ExprPtr fy = arena.make<Binary>(BinaryOp::Mul, f, y);
    ExprPtr lhs = arena.make<Binary>(BinaryOp::Add, dy.value(), fy);

    ExprPtr rhs = arena.make<IntegerLit>(BigInt(0));
    for (std::size_t i = 0; i < g_vec.size(); ++i) {
        if (c[i].numerator().is_zero()) continue;
        ExprPtr c_e = arena.make<RationalLit>(c[i].numerator(), c[i].denominator());
        ExprPtr term = arena.make<Binary>(BinaryOp::Mul, c_e, g_vec[i]);
        rhs = arena.make<Binary>(BinaryOp::Add, rhs, term);
    }

    ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, lhs, rhs);
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

[[nodiscard]] bool is_nontrivial(const ParametricRischDeQSolution& sol) {
    for (const auto& ci : sol.c)
        if (!ci.numerator().is_zero()) return true;
    return false;
}

class ParametricTowerTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};

    DifferentialField log_tower() {
        // t = log(x),  D(t) = 1/x
        DifferentialExtension ext{ExtensionType::Logarithmic,
                                  parse_expr("x", ctx.arena()), Symbol{"t"}};
        return DifferentialField(x, {ext});
    }

    DifferentialField exp_tower() {
        // t = exp(x),  D(t) = t
        DifferentialExtension ext{ExtensionType::Exponential,
                                  parse_expr("x", ctx.arena()), Symbol{"t"}};
        return DifferentialField(x, {ext});
    }
};

// --- df <= 0, caso esponenziale : il solver DEVE risolvere, soluzioni sound. ---

// Exp tower, f = 0, g = {t} (= exp(x)).  RDE: y' = c·t.  Soluzione: y = c·t.
// Esercita il ramo esponenziale (F_eff = i·u' + f_0 con i>0).
TEST_F(ParametricTowerTest, Exp_F0_GExp_AllSound) {
    DifferentialField field = exp_tower();
    auto f = parse_expr("0", ctx.arena());
    std::vector<ExprPtr> g = {parse_expr("t", ctx.arena())};
    auto res = solve_risch_de_parametric_field(f, g, field.extensions().size(), field, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    bool any_nontrivial = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_field_de(sol.y, f, g, sol.c, field, ctx))
            << "unsound solution in exp case";
        any_nontrivial |= is_nontrivial(sol);
    }
    EXPECT_TRUE(any_nontrivial) << "expected y = c*t solution (c != 0)";
}

// --- Caso primitivo (log): gap noto, ramo incompleto (A26 sub-gap). ---

// Log tower, f = 0, g = {1}.  Esiste matematicamente la soluzione y = c·x, c=1.
// Ma il ramo primitivo della ricorsione genera, dalla correzione i·y·D(t) con
// D(t)=1/x, una forzante RAZIONALE che la base parametric_q (solo Q[x] polinomi)
// non gestisce → ritorna Unimplemented.  Questo NON è un bug silenzioso: il
// contratto REGOLA ZERO (diagnostico esplicito, mai hang/crash/silent-wrong) è
// rispettato.  Il completamento (ParamRischDE su Q(x) razionale, weak-normalizer
// + denominator bound, Bronstein §5.12/§6.5) è la voce ledger HC-A26-PRIMITIVE-
// PARAMQ-RATIONAL.  Quando implementato, questo test passerà a verify_field_de.
TEST_F(ParametricTowerTest, Log_PrimitiveDescent_RationalForcing_CleanDiagnostic) {
    DifferentialField field = log_tower();
    auto f = parse_expr("0", ctx.arena());
    std::vector<ExprPtr> g = {parse_expr("1", ctx.arena())};
    auto res = solve_risch_de_parametric_field(f, g, field.extensions().size(), field, ctx);
    ASSERT_TRUE(res.is_error())
        << "primitive descent currently incomplete: must return a diagnostic, "
           "not a (possibly unsound) solution";
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented)
        << "primitive-case gap must be a clean Unimplemented (no hang/crash/"
           "silent-wrong), got kind=" << static_cast<int>(res.error().kind)
        << " msg=" << res.error().message;
}

// --- df > 0 : ramo non-cancellation (A1).  Deterministicamente raggiunto. ---

// Log tower, f = t (= log(x)), g = {1}.  Dopo il denominator-clearing f_new = t,
// deg_t(f) = 1 > 0 → ramo non-cancellation.  ATTUALMENTE Unimplemented (gap A1).
// Questo test e' la riproduzione concreta che sblocca A1: quando Bronstein §6.5
// sara' implementato, l'asserzione passera' a verify_field_de.
TEST_F(ParametricTowerTest, Log_FEqualsTheta_DfPositive_ReachesNonCancellation) {
    DifferentialField field = log_tower();
    auto f = parse_expr("t", ctx.arena());           // f = log(x), deg_t = 1
    std::vector<ExprPtr> g = {parse_expr("1", ctx.arena())};
    auto res = solve_risch_de_parametric_field(f, g, field.extensions().size(), field, ctx);
    // Il ramo df>0 e' raggiunto: contratto attuale = Unimplemented esplicito,
    // MAI hang ne' risultato sbagliato silenzioso (REGOLA ZERO).
    ASSERT_TRUE(res.is_error())
        << "deg_t(f)>0 branch must currently return Unimplemented, not a solution";
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented)
        << "df>0 must be a diagnostic Unimplemented (Bronstein 6.5 gap, A1), got kind="
        << static_cast<int>(res.error().kind) << " msg=" << res.error().message;
}

}  // namespace
}  // namespace cas::calculus
