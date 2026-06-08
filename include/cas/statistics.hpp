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

}  // namespace cas::statistics
