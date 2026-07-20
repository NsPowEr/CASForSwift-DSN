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
    explicit NumericEvaluator(const NumericEnv& env = {}, std::size_t max_recursion_depth = 256U)
        : env_(env), max_recursion_depth_(max_recursion_depth) {}
    explicit NumericEvaluator(symbolic::CASContext& ctx, const NumericEnv& env = {});

    [[nodiscard]] Result<double> evaluate(ExprPtr expr);

private:
    NumericEnv env_;
    std::size_t max_recursion_depth_{256U};
    std::size_t hyp_1f1_max_terms_{4000U};
    double hyp_1f1_rel_tol_{1e-16};
    std::size_t current_depth_{0U};
    std::unordered_set<ExprPtr, ExprHash> active_nodes_;
};

/**
 * Funzione di comodo per valutazioni rapide (double).
 */
[[nodiscard]] Result<double> eval(ExprPtr expr, const NumericEnv& env = {});
[[nodiscard]] Result<double> eval(ExprPtr expr, symbolic::CASContext& ctx, const NumericEnv& env = {});

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
 * Overload context-aware: usa `ctx.numeric_precision_digits()` come
 * default. L3-03 Float contestuale.
 */
[[nodiscard]] Result<std::string> eval_mpfr(ExprPtr expr,
    symbolic::CASContext& ctx,
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

// F8.0-5.4: rigorous isolating intervals for the real roots of a
// rational polynomial. Same input contract as find_polynomial_roots_sturm,
// but the output is the exact rational interval [low_i, high_i] s.t. each
// interval contains exactly one real root of the squarefree part of `expr`.
// Intervals are returned in ascending order of the contained root.
//
// Each (low_i, high_i) is a pair (BigInt num, BigInt den) for both endpoints,
// enabling construction of RootOf nodes whose identity is decidable by
// interval comparison alone (no need to recompute Sturm at evaluation time).
[[nodiscard]] Result<std::vector<IsolatingBound>> find_polynomial_isolating_intervals(
    ExprPtr expr,
    const std::string& variable,
    symbolic::CASContext& ctx,
    double low,
    double high,
    double tol);

// F8.0-5.3 — see include/cas/numeric_bigfloat.hpp for the MPFR-precision
// Sturm refinement API. (Kept in a separate header to avoid pulling mpfr.h
// into every translation unit that uses numeric.hpp.)

// Lipschitz dyadic refinement (Hansen-style interval narrowing).
// Robust root finder for *transcendental* univariate functions on a
// compact interval [low, high] in which `f` is continuously differentiable.
// Caller supplies a pre-computed symbolic derivative `fp`. Each
// sub-interval is excluded when |f(midpoint)| > L * (width / 2), where L
// is a local 3-point sample of |f'|; sign-change descents and Lipschitz
// exclusion together remove the equispaced-grid blind spot for
// densely-oscillating functions (e.g. sin(1000x)). Tolerance and depth
// derive from the problem: tol = ctx-controlled root precision,
// max_depth = ceil(log2((high-low)/tol)) + a small safety constant.
[[nodiscard]] Result<std::vector<double>> lipschitz_refine_roots(
    ExprPtr f, ExprPtr fp, const std::string& variable,
    double low, double high, double tol, unsigned int max_depth);

/**
 * Integrazione numerica.
 */
enum class IntegrationScheme {
    AdaptiveSimpson,   // 4th-order, low-cost legacy default.
    GaussKronrod,      // F7.3-T1: Gauss-Kronrod 15/7 with nested error estimate.
};

struct IntegrationOptions {
    double tolerance{1e-8};
    std::uint32_t max_depth{20};
    IntegrationScheme scheme{IntegrationScheme::AdaptiveSimpson};
};

[[nodiscard]] Result<double> integrate_numeric(
    ExprPtr expr,
    const std::string& variable,
    double a,
    double b,
    const IntegrationOptions& options = {});

/**
 * F6.D — Adaptive G7/K15 Gauss-Kronrod integration with priority-queue
 * panel splitting (no recursion depth limit). Tolerances and resource
 * caps are driven exclusively by CASContext fields:
 *   - ctx.integration_abs_tol()
 *   - ctx.integration_rel_tol()
 *   - ctx.integration_max_intervals()
 *
 * Convergence rule (Burnikel-style monotone descent):
 *   total_error ≤ max(abs_tol, rel_tol · |total_integral|)
 *
 * Reference: QUADPACK §4 (Piessens et al. 1983).
 */
struct GKIntegrationResult {
    double value{0.0};
    double error_estimate{0.0};
    std::size_t interval_count{0U};
    bool success{false};
};

[[nodiscard]] Result<GKIntegrationResult> integrate_gauss_kronrod_adaptive(
    const std::function<Result<double>(double)>& f,
    double a,
    double b,
    const symbolic::CASContext& ctx);

/**
 * F7.3-T2 — Lagrange polynomial interpolation (barycentric form).
 */
struct InterpolationPoint {
    double x;
    double y;
};

[[nodiscard]] Result<std::vector<double>> lagrange_weights(
    const std::vector<InterpolationPoint>& points);

[[nodiscard]] Result<double> lagrange_evaluate(
    const std::vector<InterpolationPoint>& points,
    double x_eval);

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

// F7.3-B3 — Natural cubic spline interpolation.
// Input: sorted-by-x knot list (x_i, y_i), i = 0..n-1.
// Internal representation: per-segment (a, b, c, d) cubic in s = (x - x_i),
// i.e. S_i(s) = a_i + b_i s + c_i s² + d_i s³ for s ∈ [0, x_{i+1} - x_i].
// Natural boundary: S''(x_0) = S''(x_{n-1}) = 0.
struct CubicSpline {
    std::vector<double> x;
    std::vector<double> a, b, c, d;
};
[[nodiscard]] Result<CubicSpline> build_natural_cubic_spline(
    const std::vector<InterpolationPoint>& knots);
[[nodiscard]] Result<double> cubic_spline_evaluate(
    const CubicSpline& spline, double x_eval);

// F7.3-B3 — Cubic Hermite interpolation.
// Input: (x_i, y_i, dy_i) where dy_i is the prescribed derivative at x_i.
// Each segment is the unique cubic matching (y_i, dy_i, y_{i+1}, dy_{i+1}).
struct HermiteSpline {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> dy;
};
[[nodiscard]] Result<HermiteSpline> build_hermite_spline(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& dy);
[[nodiscard]] Result<double> hermite_evaluate(
    const HermiteSpline& spline, double x_eval);

// F7.3-B3 — Runge-Kutta-Fehlberg 4(5) adaptive ODE solver.
// Solves dy/dt = f(t, y) from (t0, y0) to t_end. Step size adapted to keep
// local truncation error ≤ tol per step. Returns the trajectory.
struct RKF45Options {
    double abs_tol;       // absolute tolerance per step
    double rel_tol;       // relative tolerance per step
    double initial_step;  // initial Δt guess
    double max_step;      // hard cap (0 = unbounded)
    std::uint32_t max_steps;  // safety cap on iteration count
    static RKF45Options default_options() {
        return {1e-6, 1e-6, 0.01, 0.0, 100000U};
    }
};
[[nodiscard]] Result<std::vector<OdePoint>> solve_ode_rkf45(
    ExprPtr expr,
    const std::string& t_var,
    const std::string& y_var,
    double t0,
    double y0,
    double t_end,
    RKF45Options options = RKF45Options::default_options());

} // namespace cas::numeric
