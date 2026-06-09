// F7.2 — Statistics package (HP-Prime parity target).
//
// Provides distributions, sampling, and regression utilities built on top
// of the symbolic engine.  Each distribution exposes:
//   - symbolic `pdf(x; θ)` returning an ExprPtr in the variable `x`,
//   - symbolic `cdf(x; θ)`,
//   - numeric `quantile(p; θ)` (inverse CDF) returning a double.

#pragma once

#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <vector>

namespace cas::statistics {

/// Normal distribution N(μ, σ²) — closed-form pdf via exp, cdf via Erf,
/// quantile via Newton on the residual cdf(z) − p.
[[nodiscard]] ExprPtr normal_pdf(
    ExprPtr x, ExprPtr mu, ExprPtr sigma, symbolic::CASContext& ctx);

[[nodiscard]] ExprPtr normal_cdf(
    ExprPtr x, ExprPtr mu, ExprPtr sigma, symbolic::CASContext& ctx);

[[nodiscard]] Result<double> normal_quantile(double p, double mu, double sigma);

/// F7.2-T2 — Ordinary Least Squares result for univariate linear regression
/// y ≈ intercept + slope · x.
struct LinearRegressionResult {
    double intercept;
    double slope;
    double r_squared;
    double residual_sum_of_squares;
};

[[nodiscard]] Result<LinearRegressionResult> linear_regression(
    const std::vector<double>& x,
    const std::vector<double>& y);

// ── F7.2-T3: discrete distributions ──────────────────────────────────────────

/// Binomial distribution Bin(n, p): PMF P(X = k) = C(n, k) · p^k · (1-p)^(n-k).
/// Returns InvalidArgument if k > n, n < 0, or p ∉ [0, 1].
[[nodiscard]] Result<double> binomial_pmf(long long k, long long n, double p);

/// CDF P(X ≤ k) via cumulative summation of binomial_pmf.  Clamps k to
/// [0, n]; returns 0 for k < 0 and 1 for k ≥ n.
[[nodiscard]] Result<double> binomial_cdf(long long k, long long n, double p);

/// Poisson distribution Pois(λ): PMF P(X = k) = λ^k · exp(-λ) / k!.
/// Returns InvalidArgument if k < 0 or λ < 0.
[[nodiscard]] Result<double> poisson_pmf(long long k, double lambda);

/// CDF P(X ≤ k) via cumulative summation of poisson_pmf.
[[nodiscard]] Result<double> poisson_cdf(long long k, double lambda);

// ── F7.2-T4: continuous distributions (Chi-squared, Student-t, F) ────────────

/// Chi-squared χ²_k PDF.  k ∈ (0, ∞) degrees of freedom; x ≥ 0.
[[nodiscard]] Result<double> chi_squared_pdf(double x, double k);

/// χ²_k CDF computed via the regularised lower incomplete gamma function.
[[nodiscard]] Result<double> chi_squared_cdf(double x, double k);

/// Student-t with `nu` degrees of freedom.  Valid on the whole real line.
[[nodiscard]] Result<double> student_t_pdf(double x, double nu);

/// Student-t CDF via the regularised incomplete Beta function.
[[nodiscard]] Result<double> student_t_cdf(double x, double nu);

/// Fisher-Snedecor F-distribution with (d1, d2) degrees of freedom.
[[nodiscard]] Result<double> f_pdf(double x, double d1, double d2);

/// F-distribution CDF via the regularised incomplete Beta.
[[nodiscard]] Result<double> f_cdf(double x, double d1, double d2);

// ── F7.2-B2: Hypothesis testing ──────────────────────────────────────────

/// Result struct shared by all hypothesis-test entry points.
struct HypothesisTestResult {
    double statistic;   // observed test statistic
    double p_value;     // two-sided p-value (for symmetric tests)
    double df;          // degrees of freedom (0 if N/A)
    bool   reject_at_alpha_005;
};

/// One-sample z-test for the mean against the null mu0 with known sigma.
/// Two-sided alternative; statistic z = (sample_mean - mu0) / (sigma/√n).
[[nodiscard]] Result<HypothesisTestResult> z_test_one_sample(
    const std::vector<double>& sample, double mu0, double sigma);

/// One-sample t-test for the mean against mu0 with unknown variance.
/// Two-sided; df = n - 1; uses sample standard deviation.
[[nodiscard]] Result<HypothesisTestResult> t_test_one_sample(
    const std::vector<double>& sample, double mu0);

/// Welch's two-sample t-test (unequal variances).
[[nodiscard]] Result<HypothesisTestResult> t_test_two_sample(
    const std::vector<double>& s1, const std::vector<double>& s2);

/// Chi-squared goodness-of-fit. Compares observed vs expected counts;
/// df = bins - 1.
[[nodiscard]] Result<HypothesisTestResult> chi_squared_goodness_of_fit(
    const std::vector<double>& observed,
    const std::vector<double>& expected);

/// Two-sample F-test for variance ratio. statistic = s1² / s2².
[[nodiscard]] Result<HypothesisTestResult> f_test_variance(
    const std::vector<double>& s1, const std::vector<double>& s2);

// ── F7.2-B2: Multivariate Ordinary Least Squares ─────────────────────────

struct MultivariateOLSResult {
    std::vector<double> coefficients;          // size p+1: [intercept, β_1, …, β_p]
    std::vector<double> residuals;             // size n
    double r_squared;
    double residual_sum_of_squares;
};

/// y ≈ β_0 + Σ β_j · x_{·,j}.  X is n × p (rows = observations, cols = predictors).
/// Uses the normal equations (X^T X) β = X^T y solved via Gauss-Jordan with
/// partial pivoting on the augmented (p+1) × (p+2) system. Returns
/// InvalidArgument on dimension mismatch or singular (X^T X).
[[nodiscard]] Result<MultivariateOLSResult> multivariate_linear_regression(
    const std::vector<std::vector<double>>& X,
    const std::vector<double>& y);

}  // namespace cas::statistics
