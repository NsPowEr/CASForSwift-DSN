// A26/A1 — Tests per solve_risch_de_parametric_field (parametric Risch DE su
// torre differenziale K(x, t_1, ..., t_n)).  Bronstein, "Symbolic Integration
// I", §6.4 SPDE / §6.5 non-cancellation / §7.1 parametric PolyRischDE.
//
// Scopo (A26): la funzione `solve_risch_de_parametric_field` era *entry-orphan*
// (raggiungibile solo dalla propria ricorsione `ext_idx-1`, mai invocata
// dall'esterno) e completamente non testata.  Questo file la rende raggiungibile
// e testata in isolamento, verificando la correttezza per back-substitution
// SOUND nel campo differenziale:  D(y) + f·y  ≡  Σ_i c_i·g_i.
//
// Il ramo deg_t(f) > 0 ("non-cancellation") — gap A1 — è ora IMPLEMENTATO
// (solve_param_poly_risch_de_nocancel1, Bronstein §7.1 ParamPolyRischDENoCancel1):
// vedi i test Log_/Exp_FEqualsTheta_DfPositive_*.

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

// --- Caso primitivo (log): risolto via rational limited integration (A26). ---

// Log tower, f = 0, g = {1}.  La ricorsione primitiva genera, dalla correzione
// i·y·D(t) con D(t)=1/x, una forzante RAZIONALE; il base case ora la risolve via
// solve_param_limited_integration_rational_q (HC-A26-PRIMITIVE-PARAMQ-RATIONAL).
// Soluzione attesa: y = x, c = 1 (D(x) = 1 = c·1).  Ogni soluzione DEVE essere
// sound (back-substitution) e ne deve esistere almeno una non banale.
TEST_F(ParametricTowerTest, Log_F0_G1_PrimitiveDescent_SolvedSound) {
    DifferentialField field = log_tower();
    auto f = parse_expr("0", ctx.arena());
    std::vector<ExprPtr> g = {parse_expr("1", ctx.arena())};
    auto res = solve_risch_de_parametric_field(f, g, field.extensions().size(), field, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    bool any_nontrivial = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_field_de(sol.y, f, g, sol.c, field, ctx))
            << "unsound solution in primitive (log) descent";
        any_nontrivial |= is_nontrivial(sol);
    }
    EXPECT_TRUE(any_nontrivial)
        << "expected the nontrivial solution y = x, c = 1 (∫1 = x is rational)";
}

// --- Direct tests of rational limited integration over Q(x) (A26 helper). ---

// y' = Σ c_i g_i back-substitution check over Q(x) (base field, real diff).
[[nodiscard]] bool verify_base_de(
    ExprPtr y, const std::vector<ExprPtr>& g, const std::vector<Rational>& c,
    const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    auto dy = diff(y, var, 1U, ctx);
    if (dy.is_error()) return false;
    ExprPtr rhs = arena.make<IntegerLit>(BigInt(0));
    for (std::size_t i = 0; i < g.size(); ++i) {
        if (c[i].numerator().is_zero()) continue;
        ExprPtr ce = arena.make<RationalLit>(c[i].numerator(), c[i].denominator());
        rhs = arena.make<Binary>(BinaryOp::Add, rhs, arena.make<Binary>(BinaryOp::Mul, ce, g[i]));
    }
    ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, dy.value(), rhs);
    auto tog = algebra::together(delta, ctx);
    auto s = ctx.simplify(tog.is_ok() ? tog.value() : delta);
    if (s.is_error()) return false;
    if (const auto* il = expr_cast<IntegerLit>(s.value())) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(s.value())) return rl->numerator.is_zero();
    return false;
}

// Pure rational forcing (no log/arctan): every g_i directly integrable.
// g = {2x} → y = x², c = 1.
TEST_F(ParametricTowerTest, RationalLimited_PureRational_DirectAntiderivative) {
    std::vector<ExprPtr> g = {parse_expr("2*x", ctx.arena())};
    auto res = solve_param_limited_integration_rational_q(g, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    bool any_nontrivial = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_base_de(sol.y, g, sol.c, x, ctx)) << "unsound (2x)";
        any_nontrivial |= is_nontrivial(sol);
    }
    EXPECT_TRUE(any_nontrivial) << "expected y = x², c = 1";
}

// Log cancellation: g = {1/x, 1/x} → ∫ both = log(x); c = (1, −1) gives y = 0.
TEST_F(ParametricTowerTest, RationalLimited_LogCancellation_NontrivialC) {
    std::vector<ExprPtr> g = {parse_expr("1/x", ctx.arena()),
                              parse_expr("1/x", ctx.arena())};
    auto res = solve_param_limited_integration_rational_q(g, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    bool found_cancel = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_base_de(sol.y, g, sol.c, x, ctx)) << "unsound cancellation";
        if (sol.c.size() == 2U && !sol.c[0].numerator().is_zero() &&
            !sol.c[1].numerator().is_zero())
            found_cancel = true;
    }
    EXPECT_TRUE(found_cancel) << "expected nontrivial c = (1, −1) with y = 0";
}

// Single surviving log: g = {1/x} → ∫ = log(x), no rational antiderivative for
// c ≠ 0, so the only admissible constant is c = 0 (null space is trivial).
TEST_F(ParametricTowerTest, RationalLimited_SingleLog_NoNontrivialSolution) {
    std::vector<ExprPtr> g = {parse_expr("1/x", ctx.arena())};
    auto res = solve_param_limited_integration_rational_q(g, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_base_de(sol.y, g, sol.c, x, ctx)) << "unsound single-log";
        EXPECT_TRUE(sol.c[0].numerator().is_zero())
            << "c ≠ 0 would require ∫ c/x = c·log(x) to be rational — impossible";
    }
}

// --- Direct tests of full parametric Risch DE over Q(x), f≠0 (A26 fase 3). ---

// Base field Q(x) (no extension): field.derive = d/dx — for verify_field_de.
[[nodiscard]] static DifferentialField base_field_qx(const Symbol& x) {
    return DifferentialField(x, {});
}

// f = 1/x, g = {2}.  y' + y/x = c·2 has y = x, c = 1 (1 + x/x = 2).
TEST_F(ParametricTowerTest, RationalRischDE_FirstOrder_PolynomialSolution) {
    auto f = parse_expr("1/x", ctx.arena());
    std::vector<ExprPtr> g = {parse_expr("2", ctx.arena())};
    auto res = solve_param_risch_de_rational_q(f, g, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    DifferentialField bf = base_field_qx(x);
    bool any_nontrivial = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_field_de(sol.y, f, g, sol.c, bf, ctx)) << "unsound (1/x, 2)";
        any_nontrivial |= is_nontrivial(sol);
    }
    EXPECT_TRUE(any_nontrivial) << "expected y = x, c = 1";
}

// f = −1/x, g = {1/x²}.  y' − y/x = c/x² has the particular solution
// y = −1/(2x), c = 1 (den(y)=x divides lcm(x, x²)=x², in scope for the P/D
// ansatz).  Each returned solution is independently confirmed by the symbolic
// field check (D(y)+f·y ≡ Σ c_i g_i), and the expected nontrivial-c solution
// must be present.
TEST_F(ParametricTowerTest, RationalRischDE_PoleSolution_DenDividesLcm) {
    auto f = parse_expr("-1/x", ctx.arena());
    std::vector<ExprPtr> g = {parse_expr("1/x^2", ctx.arena())};
    auto res = solve_param_risch_de_rational_q(f, g, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    DifferentialField bf = base_field_qx(x);
    bool found_particular = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_field_de(sol.y, f, g, sol.c, bf, ctx)) << "unsound (-1/x, 1/x²)";
        if (!sol.c.empty() && !sol.c[0].numerator().is_zero()) found_particular = true;
    }
    EXPECT_TRUE(found_particular) << "expected the particular solution with c ≠ 0";
}

// f = −2/x, g = {1}.  Homogeneous y' − 2y/x = 0 has rational solution y = x²
// (with c = 0); the cancellation-degree bound (deg = −lc(H)/lc(D)) brings it
// into the ansatz.  The solver must return it inside the solution space.
TEST_F(ParametricTowerTest, RationalRischDE_HomogeneousRationalSolution) {
    auto f = parse_expr("-2/x", ctx.arena());
    std::vector<ExprPtr> g = {parse_expr("1", ctx.arena())};
    auto res = solve_param_risch_de_rational_q(f, g, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    DifferentialField bf = base_field_qx(x);
    bool found_homogeneous = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_field_de(sol.y, f, g, sol.c, bf, ctx)) << "unsound (-2/x, 1)";
        bool c_zero = true;
        for (const auto& ci : sol.c) if (!ci.numerator().is_zero()) c_zero = false;
        if (c_zero) {
            if (const auto* il = expr_cast<IntegerLit>(sol.y); !il || !il->value.is_zero())
                found_homogeneous = true;
        }
    }
    EXPECT_TRUE(found_homogeneous) << "expected homogeneous y = x² (c = 0)";
}

// --- WeakNormalizer (Bronstein 6.1.1): den(y) STRICTLY exceeds lcm. ---
//
// These are the residuo A26 closes.  At a SIMPLE pole of f with positive-integer
// residue n (convention y'+f·y=g ⇒ homogeneous y_h satisfies y_h'/y_h = −f, so a
// residue +n gives y_h = (x−α)^{−n}, a pole of order n), the solution denominator
// is x^n while lcm(den f, den g_i) only carries x^1.  Without WeakNormalizer
// denominator inflation the P/D ansatz cannot represent y and the solution is
// silently absent.  inflate_denominator (Rothstein-Trager residue n at the roots
// of gcd(fn−n·fd', s)) lifts D to x^n, bringing it into scope.  Each returned
// candidate is still independently confirmed by the symbolic field check, so the
// inflation can only affect completeness, never soundness.

// f = 2/x (residue +2 at x=0), g = {1}.  Homogeneous y' + 2y/x = 0 has the
// rational solution y = 1/x² (c = 0): den(y)=x² > lcm(den f=x, den g=1)=x.
// Found ONLY because D is inflated x → x².
TEST_F(ParametricTowerTest, WeakNormalizer_ResidueTwo_HomogeneousPoleBeyondLcm) {
    auto f = parse_expr("2/x", ctx.arena());
    std::vector<ExprPtr> g = {parse_expr("1", ctx.arena())};
    auto res = solve_param_risch_de_rational_q(f, g, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    DifferentialField bf = base_field_qx(x);
    bool found_pole = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_field_de(sol.y, f, g, sol.c, bf, ctx)) << "unsound (2/x, 1)";
        bool c_zero = true;
        for (const auto& ci : sol.c) if (!ci.numerator().is_zero()) c_zero = false;
        // A non-constant rational homogeneous solution (den exceeds lcm).
        if (c_zero) {
            if (const auto* il = expr_cast<IntegerLit>(sol.y); !il || !il->value.is_zero())
                found_pole = true;
        }
    }
    EXPECT_TRUE(found_pole)
        << "expected homogeneous y = 1/x² (den x² > lcm x) via WeakNormalizer inflation";
}

// f = 3/x (residue +3), g = {1}.  Homogeneous y = 1/x³ (c = 0): den order 3,
// needs D inflated x → x³ (two extra factors).  Exercises the n=3 inflation loop.
TEST_F(ParametricTowerTest, WeakNormalizer_ResidueThree_HigherOrderPole) {
    auto f = parse_expr("3/x", ctx.arena());
    std::vector<ExprPtr> g = {parse_expr("1", ctx.arena())};
    auto res = solve_param_risch_de_rational_q(f, g, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    DifferentialField bf = base_field_qx(x);
    bool found_pole = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_field_de(sol.y, f, g, sol.c, bf, ctx)) << "unsound (3/x, 1)";
        bool c_zero = true;
        for (const auto& ci : sol.c) if (!ci.numerator().is_zero()) c_zero = false;
        if (c_zero) {
            if (const auto* il = expr_cast<IntegerLit>(sol.y); !il || !il->value.is_zero())
                found_pole = true;
        }
    }
    EXPECT_TRUE(found_pole)
        << "expected homogeneous y = 1/x³ (den x³ > lcm x) via WeakNormalizer inflation";
}

// Residue at a shifted simple pole: f = 2/(x−1), g = {1}.  y_h = 1/(x−1)².
// Confirms inflation is anchored at the actual pole, not hardcoded to x.
TEST_F(ParametricTowerTest, WeakNormalizer_ResidueTwo_ShiftedPole) {
    auto f = parse_expr("2/(x-1)", ctx.arena());
    std::vector<ExprPtr> g = {parse_expr("1", ctx.arena())};
    auto res = solve_param_risch_de_rational_q(f, g, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    DifferentialField bf = base_field_qx(x);
    bool found_pole = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_field_de(sol.y, f, g, sol.c, bf, ctx)) << "unsound (2/(x-1), 1)";
        bool c_zero = true;
        for (const auto& ci : sol.c) if (!ci.numerator().is_zero()) c_zero = false;
        if (c_zero) {
            if (const auto* il = expr_cast<IntegerLit>(sol.y); !il || !il->value.is_zero())
                found_pole = true;
        }
    }
    EXPECT_TRUE(found_pole)
        << "expected homogeneous y = 1/(x−1)² via WeakNormalizer inflation at x=1";
}

// Irreducible QUADRATIC pole (complex-conjugate poles): f = 4x/(x²+1), residue
// +2.  y_h satisfies y_h'/y_h = −4x/(x²+1) ⇒ y_h = (x²+1)^{−2}.  Proves the
// inflation is generic over squarefree factors (gcd(fn−n·fd', s)), not anchored
// to linear poles: D must be lifted (x²+1) → (x²+1)².
TEST_F(ParametricTowerTest, WeakNormalizer_ResidueTwo_IrreducibleQuadraticPole) {
    auto f = parse_expr("4*x/(x^2+1)", ctx.arena());
    std::vector<ExprPtr> g = {parse_expr("1", ctx.arena())};
    auto res = solve_param_risch_de_rational_q(f, g, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    DifferentialField bf = base_field_qx(x);
    bool found_pole = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_field_de(sol.y, f, g, sol.c, bf, ctx)) << "unsound (4x/(x²+1), 1)";
        bool c_zero = true;
        for (const auto& ci : sol.c) if (!ci.numerator().is_zero()) c_zero = false;
        if (c_zero) {
            if (const auto* il = expr_cast<IntegerLit>(sol.y); !il || !il->value.is_zero())
                found_pole = true;
        }
    }
    EXPECT_TRUE(found_pole)
        << "expected homogeneous y = 1/(x²+1)² via WeakNormalizer at the quadratic pole";
}

// Wiring: full parametric Risch DE reachable through the field-solver base case
// (ext_idx == 0, f ≠ 0 rational).  f = 1/x, g = {2} → y = x, c = 1.
TEST_F(ParametricTowerTest, RationalRischDE_ReachedThroughFieldBaseCase) {
    DifferentialField bf = base_field_qx(x);
    auto f = parse_expr("1/x", ctx.arena());
    std::vector<ExprPtr> g = {parse_expr("2", ctx.arena())};
    auto res = solve_risch_de_parametric_field(f, g, 0U, bf, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    bool any_nontrivial = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_field_de(sol.y, f, g, sol.c, bf, ctx)) << "unsound via base case";
        any_nontrivial |= is_nontrivial(sol);
    }
    EXPECT_TRUE(any_nontrivial) << "f≠0 rational must solve through the field base case";
}

// --- df > 0 : ramo non-cancellation (A1).  Deterministicamente raggiunto. ---

// --- df > 0 : ramo non-cancellation (A1) — IMPLEMENTATO (Bronstein §7.1). ---

// Log tower, f = t (= log(x)), g = {1}.  Dopo il denominator-clearing f_new = t,
// deg_t(f) = 1 > 0 → ramo non-cancellation (ParamPolyRischDENoCancel1).  Qui la
// DE  D(y) + log(x)·y = c·1  NON ha soluzione elementare razionale per c ≠ 0
// (la soluzione omogenea exp(−∫log x) non è nella torre), quindi lo spazio
// soluzione è banale (solo c = 0).  Il solver DEVE ritornare OK (non più
// Unimplemented) con nessuna soluzione non-banale, e ogni soluzione sound.
TEST_F(ParametricTowerTest, Log_FEqualsTheta_DfPositive_NonCancellation_OnlyTrivial) {
    DifferentialField field = log_tower();
    auto f = parse_expr("t", ctx.arena());           // f = log(x), deg_t = 1
    std::vector<ExprPtr> g = {parse_expr("1", ctx.arena())};
    auto res = solve_risch_de_parametric_field(f, g, field.extensions().size(), field, ctx);
    ASSERT_TRUE(res.is_ok()) << "df>0 non-cancellation now implemented: " << res.error().message;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_field_de(sol.y, f, g, sol.c, field, ctx)) << "unsound df>0 solution";
        EXPECT_FALSE(is_nontrivial(sol))
            << "D(y)+log(x)·y = c has no elementary rational solution for c ≠ 0";
    }
}

// Log tower, f = t, g = {t}.  D(y) + log(x)·y = c·log(x) has y = 1, c = 1
// (D(1)+log(x)·1 = log(x)).  POSITIVE non-cancellation solve (h_1 peeled at
// n=0, empty residual → c free).
TEST_F(ParametricTowerTest, Log_FEqualsTheta_DfPositive_ConstantSolution) {
    DifferentialField field = log_tower();
    auto f = parse_expr("t", ctx.arena());
    std::vector<ExprPtr> g = {parse_expr("t", ctx.arena())};
    auto res = solve_risch_de_parametric_field(f, g, field.extensions().size(), field, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    bool any_nontrivial = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_field_de(sol.y, f, g, sol.c, field, ctx)) << "unsound (log t, {t})";
        any_nontrivial |= is_nontrivial(sol);
    }
    EXPECT_TRUE(any_nontrivial) << "expected y = 1, c = 1";
}

// Exp tower, f = t (= exp(x)), g = {t² + t}.  D(y) + exp(x)·y = c·(t²+t) has
// y = t, c = 1  (D(t)+t·t = t + t² = t²+t).  POSITIVE non-cancellation solve
// with N = 1 (two peeling passes n=1,0), exercises D(s·t^n) with D(t)=t.
TEST_F(ParametricTowerTest, Exp_FEqualsTheta_DfPositive_ExpSolution) {
    DifferentialField field = exp_tower();
    auto f = parse_expr("t", ctx.arena());
    std::vector<ExprPtr> g = {parse_expr("t^2 + t", ctx.arena())};
    auto res = solve_risch_de_parametric_field(f, g, field.extensions().size(), field, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    bool any_nontrivial = false;
    for (const auto& sol : res.value()) {
        EXPECT_TRUE(verify_field_de(sol.y, f, g, sol.c, field, ctx)) << "unsound (exp t, {t²+t})";
        any_nontrivial |= is_nontrivial(sol);
    }
    EXPECT_TRUE(any_nontrivial) << "expected y = t (= exp x), c = 1";
}

}  // namespace
}  // namespace cas::calculus
