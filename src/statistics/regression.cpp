// F7.2-T2 — Ordinary Least Squares (univariate linear regression).
//
// Given paired samples (x_i, y_i), fit  ŷ = β₀ + β₁·x  by minimising
// Σ (y_i − ŷ_i)².  Closed form:
//
//   β₁ = Σ(x_i − x̄)·(y_i − ȳ) / Σ(x_i − x̄)²
//   β₀ = ȳ − β₁·x̄
//   R² = 1 − SS_res / SS_tot       where SS_res = Σ(y_i − ŷ_i)²
//                                         SS_tot = Σ(y_i − ȳ)²
//
// Multivariate OLS β̂ = (Xᵀ·X)⁻¹·Xᵀ·y is left to a follow-up that wires
// matrix_solve from the linalg layer; this scope keeps the public API
// stable for the univariate target documented in PLAN_HP_PRIME_PARITY.md
// F7.2.

#include "cas/statistics.hpp"

#include <cstddef>

namespace cas::statistics {

namespace {
[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    CASError err{};
    err.kind = kind;
    err.message = std::move(message);
    return err;
}
}  // namespace

Result<LinearRegressionResult> linear_regression(
    const std::vector<double>& x,
    const std::vector<double>& y)
{
    const std::size_t n = x.size();
    if (n != y.size() || n < 2U) {
        return fail<LinearRegressionResult>(make_error(
            CASErrorKind::InvalidArgument,
            "linear_regression: x and y must have equal size ≥ 2"));
    }
    double sum_x = 0.0, sum_y = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
    }
    const double mean_x = sum_x / static_cast<double>(n);
    const double mean_y = sum_y / static_cast<double>(n);

    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double dx = x[i] - mean_x;
        const double dy = y[i] - mean_y;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    if (sxx == 0.0) {
        return fail<LinearRegressionResult>(make_error(
            CASErrorKind::Undefined,
            "linear_regression: zero variance in x — slope undefined"));
    }
    const double slope = sxy / sxx;
    const double intercept = mean_y - slope * mean_x;

    double ss_res = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double predicted = intercept + slope * x[i];
        const double r = y[i] - predicted;
        ss_res += r * r;
    }
    const double r_squared = (syy == 0.0) ? 1.0 : 1.0 - ss_res / syy;
    return ok(LinearRegressionResult{
        .intercept = intercept,
        .slope     = slope,
        .r_squared = r_squared,
        .residual_sum_of_squares = ss_res});
}

}  // namespace cas::statistics
