// F7.3-B3 — Natural cubic spline + cubic Hermite interpolation.

#include "cas/numeric.hpp"

#include <algorithm>
#include <cmath>

namespace cas::numeric {

namespace {

CASError make_error(CASErrorKind kind, std::string message) {
    CASError err{};
    err.kind = kind;
    err.message = std::move(message);
    return err;
}

}  // namespace

Result<CubicSpline> build_natural_cubic_spline(
    const std::vector<InterpolationPoint>& knots) {
    const std::size_t n = knots.size();
    if (n < 2U) {
        return fail<CubicSpline>(make_error(
            CASErrorKind::InvalidArgument,
            "build_natural_cubic_spline: need ≥ 2 knots"));
    }
    // Verify strictly increasing x.
    for (std::size_t i = 1; i < n; ++i) {
        if (!(knots[i].x > knots[i - 1].x)) {
            return fail<CubicSpline>(make_error(
                CASErrorKind::InvalidArgument,
                "build_natural_cubic_spline: x must be strictly increasing"));
        }
    }

    std::vector<double> h(n - 1);
    for (std::size_t i = 0; i < n - 1; ++i) h[i] = knots[i + 1].x - knots[i].x;

    // Solve tridiagonal system A · c = r for c_i (i = 0..n-1) with
    // c_0 = c_{n-1} = 0 (natural BC). c here are S''(x_i) / 2 conceptually,
    // but the canonical CRT form uses M = S''(x_i); we follow the (a,b,c,d)
    // segment form directly. Build the standard system on the M values, then
    // derive a,b,c,d.
    std::vector<double> M(n, 0.0);
    if (n >= 3U) {
        std::vector<double> diag(n - 2, 0.0);
        std::vector<double> upper(n - 3, 0.0);
        std::vector<double> lower(n - 3, 0.0);
        std::vector<double> rhs(n - 2, 0.0);
        for (std::size_t i = 0; i + 2 < n; ++i) {
            diag[i] = 2.0 * (h[i] + h[i + 1]);
            rhs[i] = 6.0 * ((knots[i + 2].y - knots[i + 1].y) / h[i + 1]
                          - (knots[i + 1].y - knots[i].y) / h[i]);
        }
        for (std::size_t i = 0; i + 3 < n; ++i) {
            upper[i] = h[i + 1];
            lower[i] = h[i + 1];
        }
        // Thomas algorithm.
        for (std::size_t i = 1; i + 2 < n; ++i) {
            const double m = lower[i - 1] / diag[i - 1];
            diag[i] -= m * upper[i - 1];
            rhs[i]  -= m * rhs[i - 1];
        }
        std::vector<double> inner(n - 2, 0.0);
        if (n - 2 >= 1) {
            inner[n - 3] = rhs[n - 3] / diag[n - 3];
            for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n) - 4; i >= 0; --i) {
                inner[static_cast<std::size_t>(i)] =
                    (rhs[static_cast<std::size_t>(i)]
                     - upper[static_cast<std::size_t>(i)] * inner[static_cast<std::size_t>(i) + 1])
                    / diag[static_cast<std::size_t>(i)];
            }
        }
        for (std::size_t i = 0; i + 2 < n; ++i) M[i + 1] = inner[i];
    }

    CubicSpline s;
    s.x.resize(n);
    for (std::size_t i = 0; i < n; ++i) s.x[i] = knots[i].x;
    s.a.resize(n - 1);
    s.b.resize(n - 1);
    s.c.resize(n - 1);
    s.d.resize(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        s.a[i] = knots[i].y;
        s.c[i] = M[i] / 2.0;
        s.d[i] = (M[i + 1] - M[i]) / (6.0 * h[i]);
        s.b[i] = (knots[i + 1].y - knots[i].y) / h[i]
               - h[i] * (2.0 * M[i] + M[i + 1]) / 6.0;
    }
    return ok(std::move(s));
}

Result<double> cubic_spline_evaluate(const CubicSpline& spline, double x_eval) {
    if (spline.x.size() < 2U) {
        return fail<double>(make_error(
            CASErrorKind::InvalidArgument, "cubic_spline_evaluate: empty spline"));
    }
    auto it = std::upper_bound(spline.x.begin(), spline.x.end(), x_eval);
    std::size_t i = (it == spline.x.begin())
                      ? 0
                      : static_cast<std::size_t>(it - spline.x.begin() - 1);
    if (i + 1 >= spline.x.size()) i = spline.x.size() - 2U;
    const double s = x_eval - spline.x[i];
    const double v = spline.a[i] + s * (spline.b[i] + s * (spline.c[i] + s * spline.d[i]));
    return ok(v);
}

Result<HermiteSpline> build_hermite_spline(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& dy) {
    if (x.size() != y.size() || x.size() != dy.size() || x.size() < 2U) {
        return fail<HermiteSpline>(make_error(
            CASErrorKind::InvalidArgument,
            "build_hermite_spline: dimension mismatch or n < 2"));
    }
    for (std::size_t i = 1; i < x.size(); ++i) {
        if (!(x[i] > x[i - 1])) {
            return fail<HermiteSpline>(make_error(
                CASErrorKind::InvalidArgument,
                "build_hermite_spline: x must be strictly increasing"));
        }
    }
    return ok(HermiteSpline{ .x = x, .y = y, .dy = dy });
}

Result<double> hermite_evaluate(const HermiteSpline& s, double x_eval) {
    const std::size_t n = s.x.size();
    if (n < 2U) {
        return fail<double>(make_error(
            CASErrorKind::InvalidArgument, "hermite_evaluate: empty spline"));
    }
    auto it = std::upper_bound(s.x.begin(), s.x.end(), x_eval);
    std::size_t i = (it == s.x.begin())
                      ? 0
                      : static_cast<std::size_t>(it - s.x.begin() - 1);
    if (i + 1 >= n) i = n - 2;
    const double h = s.x[i + 1] - s.x[i];
    const double t = (x_eval - s.x[i]) / h;
    // Cubic Hermite basis (h00, h10, h01, h11) on [0,1].
    const double t2 = t * t;
    const double t3 = t2 * t;
    const double h00 =  2.0 * t3 - 3.0 * t2 + 1.0;
    const double h10 =        t3 - 2.0 * t2 + t;
    const double h01 = -2.0 * t3 + 3.0 * t2;
    const double h11 =        t3 -       t2;
    const double v = h00 * s.y[i] + h10 * h * s.dy[i]
                   + h01 * s.y[i + 1] + h11 * h * s.dy[i + 1];
    return ok(v);
}

}  // namespace cas::numeric
