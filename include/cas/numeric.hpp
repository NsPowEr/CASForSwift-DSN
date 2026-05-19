#pragma once

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include <string>
#include <unordered_map>
#include <vector>
#include <functional>

namespace cas::symbolic { class CASContext; }

namespace cas::numeric {

using NumericEnv = std::unordered_map<std::string, double>;

/**
 * NumericEvaluator: trasforma un AST in un valore floating-point (double).
 * Segue il principio di "Muscoli Matematici":
 * 1. Isolamento: la logica numerica non inquina il core simbolico.
 * 2. Precisione: utilizza le funzioni matematiche standard di sistema (std::math).
 * 3. Fallback: restituisce errore se l'espressione contiene simboli non definiti nell'ambiente.
 */
class NumericEvaluator {
public:
    explicit NumericEvaluator(const NumericEnv& env = {}) : env_(env) {}

    [[nodiscard]] Result<double> evaluate(ExprPtr expr);

private:
    NumericEnv env_;
};

/**
 * Funzione di comodo per valutazioni rapide (double).
 */
[[nodiscard]] Result<double> eval(ExprPtr expr, const NumericEnv& env = {});

/**
 * Valutazione MPFR a precisione arbitraria (L3-01).
 * Ritorna decimal_digits cifre significative come stringa.
 * Esempio: N(pi, 50) -> "3.1415926535897932384626433832795028841971693993751"
 * Implementata in bigfloat_eval.cpp; non richiede bigfloat.hpp in questo header.
 */
[[nodiscard]] Result<std::string> eval_mpfr(ExprPtr expr,
    unsigned int decimal_digits,
    const NumericEnv& env = {});

/**
 * Root finding numerico.
 */
struct RootFindingOptions {
    double tolerance{1e-9};
    std::uint32_t max_iterations{100};
};

[[nodiscard]] Result<double> solve_numeric_bisection(
    ExprPtr expr,
    const std::string& variable,
    double low,
    double high,
    const RootFindingOptions& options = {});

[[nodiscard]] Result<double> solve_numeric_newton(
    ExprPtr expr,
    const std::string& variable,
    double initial_guess,
    symbolic::CASContext& ctx,
    const RootFindingOptions& options = {});

/**
 * Trova tutte le radici di expr(variable)=0 nell'intervallo [low, high].
 * Scansiona num_samples punti equidistanti, localizza sign-changes via bisection,
 * poi polish via Newton-Raphson.  Deduplica radici entro dedup_tolerance.
 * L2-06.
 */
struct MultiRootOptions {
    double low{-10.0};
    double high{10.0};
    std::size_t num_samples{200};
    double dedup_tolerance{1e-6};
    RootFindingOptions root_opts{};
};

[[nodiscard]] Result<std::vector<double>> find_roots_on_interval(
    ExprPtr expr,
    const std::string& variable,
    symbolic::CASContext& ctx,
    const MultiRootOptions& options = {});

// Sturm-sequence based real-root isolation for univariate rational
// polynomials (Sturm 1829). Returns the exact list of real roots of
// `f ∈ Q[x]` lying in `[low, high]`, isolated to within `tol`. Each
// reported root corresponds to a distinct real zero of the squarefree
// part of `f`; multiplicities are not exposed (use squarefree
// decomposition upstream to recover them). Newton polish from the
// midpoint of each isolating interval drives the final tolerance.
//
// `expr` must be parseable as a univariate polynomial in `variable`
// over the rationals; non-polynomial expressions return an error and
// the caller should fall back to grid-based numeric root finding.
[[nodiscard]] Result<std::vector<double>> find_polynomial_roots_sturm(
    ExprPtr expr,
    const std::string& variable,
    symbolic::CASContext& ctx,
    double low,
    double high,
    double tol);

/**
 * Integrazione numerica.
 */
struct IntegrationOptions {
    double tolerance{1e-8};
    std::uint32_t max_depth{20};
};

[[nodiscard]] Result<double> integrate_numeric(
    ExprPtr expr,
    const std::string& variable,
    double a,
    double b,
    const IntegrationOptions& options = {});

/**
 * Risolutore ODE.
 */
struct OdePoint {
    double t;
    double y;
};

[[nodiscard]] Result<std::vector<OdePoint>> solve_ode_rk4(
    ExprPtr expr,
    const std::string& t_var,
    const std::string& y_var,
    double t0,
    double y0,
    double t_end,
    double step_size);

} // namespace cas::numeric
