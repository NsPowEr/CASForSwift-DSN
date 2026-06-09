// F7.2-B2 — Multivariate Ordinary Least Squares regression.
//
// Solves the normal equations (X^T X) β = X^T y in pure C++ (no external
// linalg dependency). Augmented Gauss-Jordan with partial pivoting on a
// (p+1) × (p+2) matrix. Suitable for p ≤ ~50; for larger p use the symbolic
// linalg backend (cas::linalg::LU) instead.

#include "cas/statistics.hpp"
#include "cas/error_helpers.hpp"

#include <cmath>
#include <numeric>

namespace cas::statistics {

namespace {

CASError make_error(CASErrorKind kind, std::string message) {
    CASError err{};
    err.kind = kind;
    err.message = std::move(message);
    return err;
}

// Solve A · x = b in-place via Gauss-Jordan with partial pivoting.
// A is n×n, b is n. Returns false on singular.
bool gauss_jordan_solve(std::vector<std::vector<double>>& A,
                        std::vector<double>& b) {
    const std::size_t n = A.size();
    for (std::size_t k = 0; k < n; ++k) {
        // Partial pivot.
        std::size_t piv = k;
        double piv_val = std::abs(A[k][k]);
        for (std::size_t i = k + 1; i < n; ++i) {
            if (std::abs(A[i][k]) > piv_val) {
                piv = i;
                piv_val = std::abs(A[i][k]);
            }
        }
        if (piv_val < 1e-14) return false;
        if (piv != k) {
            std::swap(A[k], A[piv]);
            std::swap(b[k], b[piv]);
        }
        const double diag = A[k][k];
        for (std::size_t j = k; j < n; ++j) A[k][j] /= diag;
        b[k] /= diag;
        for (std::size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            const double f = A[i][k];
            if (f == 0.0) continue;
            for (std::size_t j = k; j < n; ++j) A[i][j] -= f * A[k][j];
            b[i] -= f * b[k];
        }
    }
    return true;
}

}  // namespace

Result<MultivariateOLSResult> multivariate_linear_regression(
    const std::vector<std::vector<double>>& X,
    const std::vector<double>& y) {
    if (X.empty() || y.empty() || X.size() != y.size()) {
        return fail<MultivariateOLSResult>(make_error(
            CASErrorKind::InvalidArgument,
            "multivariate_linear_regression: dimension mismatch or empty"));
    }
    const std::size_t n = X.size();
    const std::size_t p = X[0].size();
    for (const auto& row : X) {
        if (row.size() != p) {
            return fail<MultivariateOLSResult>(make_error(
                CASErrorKind::InvalidArgument,
                "multivariate_linear_regression: ragged X rows"));
        }
    }
    if (n <= p + 1) {
        return fail<MultivariateOLSResult>(make_error(
            CASErrorKind::InvalidArgument,
            "multivariate_linear_regression: need n > p + 1 for valid OLS"));
    }

    // Build augmented matrix [1 | X] : n × (p+1)
    const std::size_t k = p + 1;
    std::vector<std::vector<double>> Z(n, std::vector<double>(k, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        Z[i][0] = 1.0;
        for (std::size_t j = 0; j < p; ++j) Z[i][j + 1] = X[i][j];
    }

    // Normal equations: A = Z^T Z (k × k), b = Z^T y (k).
    std::vector<std::vector<double>> A(k, std::vector<double>(k, 0.0));
    std::vector<double> b(k, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t r = 0; r < k; ++r) {
            b[r] += Z[i][r] * y[i];
            for (std::size_t c = 0; c < k; ++c) {
                A[r][c] += Z[i][r] * Z[i][c];
            }
        }
    }

    if (!gauss_jordan_solve(A, b)) {
        return fail<MultivariateOLSResult>(make_error(
            CASErrorKind::InvalidArgument,
            "multivariate_linear_regression: singular (X^T X) — collinear predictors?"));
    }

    // β = b after solve. Residuals + R².
    std::vector<double> fit(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < k; ++j) fit[i] += Z[i][j] * b[j];
    }
    std::vector<double> resid(n, 0.0);
    double rss = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        resid[i] = y[i] - fit[i];
        rss += resid[i] * resid[i];
    }
    const double y_mean = std::accumulate(y.begin(), y.end(), 0.0) / static_cast<double>(n);
    double tss = 0.0;
    for (double yi : y) {
        const double d = yi - y_mean;
        tss += d * d;
    }
    const double r2 = (tss > 0.0) ? (1.0 - rss / tss) : 1.0;

    return ok(MultivariateOLSResult{
        .coefficients = std::move(b),
        .residuals = std::move(resid),
        .r_squared = r2,
        .residual_sum_of_squares = rss,
    });
}

}  // namespace cas::statistics
