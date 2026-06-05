#pragma once

#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <string>
#include <vector>
#include <set>

namespace cas::calculus {

struct QuotientView {
    ExprPtr numerator;
    ExprPtr denominator;
};

/// @brief Differentiates transcendental functions.
/// Handles asin, acos, atan, sinh, cosh, tanh.
[[nodiscard]] Result<ExprPtr> differentiate_transcendental(
    BuiltinOp func_id,
    ExprPtr argument,
    const Symbol& var,
    symbolic::CASContext& context);

/// @brief Integration by parts: ∫ u dv = uv - ∫ v du
[[nodiscard]] Result<ExprPtr> integrate_by_parts(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& context);

/// @brief Automated Variable Substitution (u-substitution recognition).
/// Returns nullopt if no substitution was found.
[[nodiscard]] Result<std::optional<ExprPtr>> integrate_by_substitution(
    const ExprPtr& integrand,
    const Symbol& var,
    symbolic::CASContext& ctx);

/// @brief Partial Risch Algorithm for Logarithmic and Exponential extensions.
[[nodiscard]] Result<ExprPtr> integrate_risch(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& context);

/// @brief Weierstrass substitution t = tan(x/2) for ∫ R(sin x, cos x) dx.
[[nodiscard]] Result<ExprPtr> integrate_weierstrass(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& ctx);

/// @brief Finite-difference symbolic derivative formulas (L3-12).
enum class FiniteDiffOrder {
    Forward1,   // O(h)
    Central2,   // O(h²)
    Central4,   // O(h⁴)
};

[[nodiscard]] Result<ExprPtr> numeric_diff(
    ExprPtr expr, const Symbol& var, ExprPtr h,
    FiniteDiffOrder order, symbolic::CASContext& ctx);

/// @brief Laplace transform L{f(t)}(s) elementary-pattern table (L3-07).
[[nodiscard]] Result<ExprPtr> laplace_transform(
    ExprPtr expr, const Symbol& t, const Symbol& s,
    symbolic::CASContext& ctx);

/// @brief Inverse Laplace L⁻¹{F(s)}(t) elementary-pattern table (L3-07).
[[nodiscard]] Result<ExprPtr> inverse_laplace_transform(
    ExprPtr expr, const Symbol& s, const Symbol& t,
    symbolic::CASContext& ctx);

/// @brief Sifting property ∫ δ(arg(t)) · f(t) dt = f(a)/|c| dove a è il
///        polo di δ (arg(a) = 0) e c = arg'(a).  Caso arg lineare in `var`.
///        Restituisce Unimplemented se arg non lineare o δ multipli.
[[nodiscard]] Result<ExprPtr> try_integrate_dirac_sifting(
    ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx);

/// @brief Inverse Laplace via Bronstein residue formula — F5.8 / Task #17.
///        L⁻¹{F(s)}(t) = Σ_k Res_{s = p_k}[F(s)·e^(s·t)] su tutti i poli di F.
///        Più generale del pattern-based inverse_laplace_transform: copre F
///        razionale con polos qualsiasi (anche multipli), usando residue()
///        del sub-engine F5.6.  Poli irrazionali / RootOf → Unimplemented.
[[nodiscard]] Result<ExprPtr> inverse_laplace_residue_q(
    ExprPtr F, const Symbol& s, const Symbol& t,
    symbolic::CASContext& ctx);

/// @brief Mellin transform M{f(t)}(s) = ∫_0^∞ t^(s-1) · f(t) dt — F5.8 / Task #15.
/// Pipeline assiomatica + proprietà:
///   - Linearità su Sum, fattoriale Product scalare;
///   - Coppie base:  exp(-a·t)              → Γ(s)/a^s
///                   t^a · exp(-b·t)        → Γ(s+a)/b^(s+a)
///                   1/(1+t)^a              → Γ(s)·Γ(a-s)/Γ(a)
///   - Proprietà:    f(a·t)                 → a^(-s) · F(s)        (scaling)
///                   t^a · f(t)             → F(s+a)               (modulation)
/// Per pattern non riconosciuti restituisce Unimplemented diagnostico.
[[nodiscard]] Result<ExprPtr> mellin_transform(
    ExprPtr expr, const Symbol& t, const Symbol& s,
    symbolic::CASContext& ctx);

/// @brief Z-transform Z{a_n}(z) = Σ_{n=0}^∞ a_n · z^(−n) — F5.8 / Task #16.
/// Convention: causal unilateral Z-transform.  Pipeline:
///   (1) linearità su Sum + estrazione scalari Product;
///   (2) coppie assiomatiche fondamentali (1, n, n^k, r^n, cos(ωn), sin(ωn),
///       δ[n−k] via z^(−k));
///   (3) proprietà di scaling (modulation r^n·a_n → A(z/r)) e
///       differenziazione (n·a_n → −z·A'(z));
///   (4) Unimplemented diagnostico esplicito per pattern non riconosciuti.
/// Nessun set chiuso silenzioso: ogni coppia base è assioma matematico
/// (Eccezione 2 di CLAUDE.md).
[[nodiscard]] Result<ExprPtr> z_transform(
    ExprPtr expr, const Symbol& n, const Symbol& z,
    symbolic::CASContext& ctx);

/// @brief Risch DE elementare in Q(x) (Bronstein §6.1+§6.2): risolve
///        y' + f·y = g per y razionale in `var`, con f, g razionali in `var`
///        e coefficienti razionali (Q).  Dispatcher unico: tenta prima il
///        ramo polynomial (y ∈ Q[var]), poi il ramo rational (y ∈ Q(var)).
///        Usato come building block per estensione cap.8 (Risch DE su torri
///        logaritmiche/esponenziali).
[[nodiscard]] Result<ExprPtr> solve_risch_de_q(
    ExprPtr f_expr, ExprPtr g_expr, const Symbol& var, symbolic::CASContext& ctx);

/// @brief Risch DE per estensione logaritmica θ = log(u), u ∈ Q(var)
///        (Bronstein §6.4.1, "PolyRischDE_logarithmic").  Risolve
///        y' + f·y = g  con f, g polinomiali in θ a coefficienti in Q(var).
///        Cerca y polinomio in θ a coefficienti in Q(var).
///
///        Algoritmo: per ogni grado i = deg(g) down to 0,
///          y_i' + f_0 · y_i = g_i − (i+1) · θ' · y_{i+1}
///        dove θ' = u'/u ∈ Q(var); back-substitution top-down sui gradi.
///        Ogni equazione di grado i è una Risch DE su Q(var), risolta via
///        solve_risch_de_q.  Restituisce y come Σ y_i · log(u)^i.
[[nodiscard]] Result<ExprPtr> solve_risch_de_logarithmic_q(
    ExprPtr f_expr,
    ExprPtr g_expr,
    ExprPtr u_arg,
    const Symbol& var,
    symbolic::CASContext& ctx);

/// @brief Risch DE per estensione esponenziale θ = exp(u), u ∈ Q(var)
///        (Bronstein §6.4.2, "PolyRischDE_exponential").  Risolve
///        y' + f·y = g  con f, g polinomiali in θ a coefficienti in Q(var).
///        Cerca y polinomio in θ a coefficienti in Q(var).
///
///        Algoritmo: per ciascun grado i (decoupled, a differenza del caso log),
///          y_i' + (i · u' + f_0) · y_i = g_i
///        ogni equazione è una Risch DE indipendente su Q(var), risolta via
///        solve_risch_de_q.  Restituisce y = Σ y_i · exp(u)^i.
[[nodiscard]] Result<ExprPtr> solve_risch_de_exponential_q(
    ExprPtr f_expr,
    ExprPtr g_expr,
    ExprPtr u_arg,
    const Symbol& var,
    symbolic::CASContext& ctx);

/// @brief Cap.9 structure theorem (sub-case log-factorization): decompone
///        log(u(var)) → Σ m_i · log(p_i(var)) + log(content) sfruttando la
///        fattorizzazione di u in Q[var].  Walk dell'AST, applica la
///        riduzione ad ogni occorrenza di FuncCall(Ln, polinomio).  Lascia
///        intatte le log con argomento non polinomiale.
///
///        Caso d'uso: pre-processore per il wiring cap.8 (integrate_risch
///        section 2b-bis); decomporre log argomenti riconduce torre log
///        con generatori multipli ma algebricamente dipendenti a torre
///        con generatore unico, abilitando l'antiderivata esatta.
[[nodiscard]] ExprPtr expand_log_args_via_factorization(
    ExprPtr e, const Symbol& var, symbolic::CASContext& ctx);

/// @brief Parametric Risch DE in Q[x] (Bronstein Symbolic Integration I §7.1).
/// Risolve  y' + f·y = Σ_i c_i · g_i  per y ∈ Q[var] e c_i ∈ Q.
/// Restituisce una base dello spazio nullo del sistema lineare omogeneo:
///   Σ_j (y'+f·y)_j · y_k + Σ_i (-g_i)_j · c_i = 0  per ogni coefficiente j.
/// Ogni soluzione della base è una coppia (y, [c_1,...,c_m]).  La dimensione
/// dello spazio nullo è (numero_unknowns - rank).
///
/// Casi:
///   - Base vuota          → solo soluzione banale (y=0, c=0).  L'integrand
///                           rispetto all'estensione è "non-elementare".
///   - Base con (0, [c])   → combinazione di c_i·g_i è già nulla (no info).
///   - Base con (y, [c])   → identità y' + f·y = Σ c_i·g_i utile per
///                           ricostruzione antiderivata in cap.8/9.
struct ParametricRischDeQSolution {
    ExprPtr y;                              ///< polynomial in `var`
    std::vector<Rational> c;                ///< rational scalars c_1...c_m
};
[[nodiscard]] Result<std::vector<ParametricRischDeQSolution>>
solve_risch_de_parametric_q(
    ExprPtr f_expr,
    const std::vector<ExprPtr>& g_vec,
    const Symbol& var,
    symbolic::CASContext& ctx);

/// @brief Limited Integration in Q[x] (Bronstein §7.2 caso polinomiale).
/// Dato f ∈ Q[var], trova rappresentazione  f = g' + Σ_i c_i · h_i
/// dove g ∈ Q[var] è "elementary" (antiderivata esatta della parte ridotta) e
/// {h_i} è il vettore di forzanti residue non-integrabili in Q[var].  Nel caso
/// polinomiale puro Q[var], la decomposizione è triviale (g = ∫f, residuo
/// vuoto) ma la routine esiste come building block per cap.8 (Risch DE
/// generale su torre trascendentale) e cap.9 (structure theorem).
struct LimitedIntegrationQSolution {
    ExprPtr g;                              ///< polynomial antiderivative
    std::vector<Rational> c;                ///< scalars on residue basis
};
[[nodiscard]] Result<LimitedIntegrationQSolution>
limited_integration_q(
    ExprPtr f_expr,
    const std::vector<ExprPtr>& residue_basis,
    const Symbol& var,
    symbolic::CASContext& ctx);

/// @brief ODE Laplace solver (L3-10).
/// Risolve linear ODE coefficienti costanti: Σ coeffs[k]·y^(k) = forcing(t).
/// initial_conditions[j] = y^(j)(0), j=0..n-1.
[[nodiscard]] Result<ExprPtr> solve_ode_laplace(
    const std::vector<ExprPtr>& coeffs,
    ExprPtr forcing,
    const std::vector<ExprPtr>& initial_conditions,
    const Symbol& t,
    symbolic::CASContext& ctx);

[[nodiscard]] ExprPtr limit_make_integer(AstArena& arena, long long value);
[[nodiscard]] ExprPtr limit_make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs);
[[nodiscard]] bool limit_is_zero(ExprPtr expr);
[[nodiscard]] bool limit_is_one(ExprPtr expr);
[[nodiscard]] bool limit_is_infinity(ExprPtr expr);
[[nodiscard]] std::optional<QuotientView> extract_quotient_view(ExprPtr expr, AstArena& arena);
[[nodiscard]] std::optional<ExprPtr> cancel_common_linear_factor(
    ExprPtr numerator,
    ExprPtr denominator,
    const Symbol& var,
    ExprPtr point,
    AstArena& arena);
[[nodiscard]] std::optional<Result<ExprPtr>> try_polynomial_pole_limit(
    ExprPtr numerator,
    ExprPtr denominator,
    const Symbol& var,
    ExprPtr point,
    LimitDirection dir,
    AstArena& arena);
[[nodiscard]] std::optional<Result<ExprPtr>> try_logarithmic_root_limit(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    LimitDirection dir,
    AstArena& arena);
[[nodiscard]] Result<ExprPtr> try_infinite_limit(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    symbolic::CASContext& ctx);

// Gruntz §3.5 asymptotic-growth comparison.  Returns +1 if `a` dominates `b`
// as var → +∞, −1 if `b` dominates, 0 when both grow at the same rate or the
// ordering cannot be decided structurally.  Recursive: handles arbitrarily
// nested exp / ln towers and falls back on polynomial degree for the same
// rank class.
[[nodiscard]] int compare_growth(
    ExprPtr a,
    ExprPtr b,
    const Symbol& var,
    symbolic::CASContext& ctx);

[[nodiscard]] bool depends_on(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool is_bounded(ExprPtr expr, const Symbol& var);

// Tower height (Gruntz §3.5): maximum nesting of exp/log chain over sub-
// expressions that depend on `var`. Structural AST walk, O(N). Used by
// LimitEngine to derive an adaptive recursion bound from the input
// asymptotic complexity instead of a fixed depth cap.
//
//   tower_height(x)                = 1
//   tower_height(exp(f))           = tower_height(f) + 1 if f depends on var
//   tower_height(ln(f) / log(f))   = max(tower_height(f), 1)
//   tower_height(other functions)  = max over children
//
// Conservative (always ≤ true comparability tower height); a stricter bound
// would require running the MRV algorithm itself.
[[nodiscard]] unsigned int transcendental_tower_depth(ExprPtr expr, const Symbol& var);

[[nodiscard]] Result<ExprPtr> symbolic_sum(
    ExprPtr term,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx);

// L2-05 closure: Laurent expansion via Taylor‑series inversion of the
// denominator.  Splits expr on its top‑level Div node (no rational
// canonicalization), Taylor‑expands numerator and denominator at `center`
// to order positive_order + pole_budget, then divides the two series.
// Handles transcendental denominators with finite‑order zeros at the
// expansion point (e.g. 1/sin(x), 1/(x²·sin(x))) which the rational
// fast path cannot.
[[nodiscard]] Result<LaurentExpansion> laurent_series_general(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr center,
    unsigned int positive_order,
    unsigned int pole_budget,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> residue(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr pole,
    symbolic::CASContext& ctx);

// Hermite reduction on the base field (exact RatPoly-based)
[[nodiscard]] Result<HermiteReduction> hermite_reduction_exact(
    ExprPtr P, ExprPtr Q, const Symbol& var, symbolic::CASContext& ctx);

// MRV Algorithm (Gruntz)
struct MRVCompare {
    bool operator()(ExprPtr lhs, ExprPtr rhs) const noexcept;
};
using MRVSet = std::set<ExprPtr, MRVCompare>;

[[nodiscard]] MRVSet mrv_set(ExprPtr e, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> rewrite_mrv(ExprPtr e, const MRVSet& mrv, ExprPtr w, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> compute_limit_mrv(ExprPtr expr, const Symbol& var, ExprPtr point, symbolic::CASContext& ctx);

}  // namespace cas::calculus
