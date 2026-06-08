// F7.3-T2 — Lagrange polynomial interpolation.
//
// Given n+1 distinct sample points (x_i, y_i), the Lagrange interpolating
// polynomial L of degree ≤ n is the unique polynomial satisfying
// L(x_i) = y_i for every i.  We construct it via the barycentric
// (second) form:
//
//   w_i      = 1 / ∏_{j ≠ i} (x_i - x_j)
//   L(x)     = (Σ_i (w_i / (x - x_i)) · y_i) / (Σ_i (w_i / (x - x_i)))
//
// The barycentric formulation is the numerically stable evaluation form
// (Berrut & Trefethen, "Barycentric Lagrange Interpolation", SIAM 2004)
// and avoids the catastrophic cancellation of the classical
// product-of-differences expression when n grows.
//
// API exposed in `cas/numeric.hpp`:
//   - `lagrange_evaluate(points, x_eval)` — evaluate L(x_eval).
//   - `lagrange_weights(points)` — return barycentric weights.
//
// Edge cases handled inline:
//   - Empty `points`         → InvalidArgument.
//   - Duplicate x-coordinate → InvalidArgument (interpolation is ill-defined).
//   - Query at a sample      → return y_i directly (avoids 0/0).

#include "cas/numeric.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace cas::numeric {

namespace {
[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    CASError err{};
    err.kind = kind;
    err.message = std::move(message);
    return err;
}
}  // namespace


Result<std::vector<double>> lagrange_weights(
    const std::vector<InterpolationPoint>& points)
{
    const std::size_t n = points.size();
    if (n == 0U) {
        return fail<std::vector<double>>(make_error(
            CASErrorKind::InvalidArgument,
            "lagrange_weights: at least one sample required"));
    }
    std::vector<double> w(n, 1.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            const double dx = points[i].x - points[j].x;
            if (dx == 0.0) {
                return fail<std::vector<double>>(make_error(
                    CASErrorKind::InvalidArgument,
                    "lagrange_weights: duplicate x-coordinate at index " +
                    std::to_string(i) + " / " + std::to_string(j)));
            }
            w[i] /= dx;
        }
    }
    return ok(std::move(w));
}

Result<double> lagrange_evaluate(
    const std::vector<InterpolationPoint>& points,
    double x_eval)
{
    if (points.empty()) {
        return fail<double>(make_error(
            CASErrorKind::InvalidArgument,
            "lagrange_evaluate: at least one sample required"));
    }
    // Direct return when x_eval coincides with a sample (also avoids the
    // 0/0 form in the barycentric ratio).
    for (const auto& p : points) {
        if (p.x == x_eval) return ok(p.y);
    }
    auto w_res = lagrange_weights(points);
    if (w_res.is_error()) return fail<double>(w_res.error());
    const auto& w = w_res.value();

    double numerator   = 0.0;
    double denominator = 0.0;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const double inv = w[i] / (x_eval - points[i].x);
        numerator   += inv * points[i].y;
        denominator += inv;
    }
    if (denominator == 0.0) {
        return fail<double>(make_error(
            CASErrorKind::Undefined,
            "lagrange_evaluate: degenerate denominator (likely duplicate sample)"));
    }
    return ok(numerator / denominator);
}

}  // namespace cas::numeric
