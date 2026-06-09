// F7.3-B3 — Runge-Kutta-Fehlberg 4(5) adaptive ODE solver.
//
// Classic Fehlberg coefficients: produces order-4 and order-5 estimates
// from a single 6-stage tableau; the difference yields the local
// truncation error estimate used to adapt Δt.

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

// Fehlberg 4(5) coefficients (Butcher tableau).
struct Tab {
    static constexpr double a21 = 1.0/4.0;
    static constexpr double a31 = 3.0/32.0,        a32 = 9.0/32.0;
    static constexpr double a41 = 1932.0/2197.0,   a42 = -7200.0/2197.0, a43 = 7296.0/2197.0;
    static constexpr double a51 = 439.0/216.0,     a52 = -8.0,            a53 = 3680.0/513.0,    a54 = -845.0/4104.0;
    static constexpr double a61 = -8.0/27.0,       a62 = 2.0,             a63 = -3544.0/2565.0,  a64 = 1859.0/4104.0,    a65 = -11.0/40.0;
    static constexpr double c2 = 1.0/4.0, c3 = 3.0/8.0, c4 = 12.0/13.0, c5 = 1.0, c6 = 1.0/2.0;
    // Order-4 weights (b4) and order-5 weights (b5).
    static constexpr double b4_1 = 25.0/216.0,  b4_3 = 1408.0/2565.0,  b4_4 = 2197.0/4104.0,  b4_5 = -1.0/5.0;
    static constexpr double b5_1 = 16.0/135.0,  b5_3 = 6656.0/12825.0, b5_4 = 28561.0/56430.0, b5_5 = -9.0/50.0, b5_6 = 2.0/55.0;
};

}  // namespace

Result<std::vector<OdePoint>> solve_ode_rkf45(
    ExprPtr expr,
    const std::string& t_var,
    const std::string& y_var,
    double t0,
    double y0,
    double t_end,
    RKF45Options options) {
    if (!(t_end > t0)) {
        return fail<std::vector<OdePoint>>(make_error(
            CASErrorKind::InvalidArgument, "solve_ode_rkf45: t_end must be > t0"));
    }
    if (!(options.initial_step > 0.0)) {
        return fail<std::vector<OdePoint>>(make_error(
            CASErrorKind::InvalidArgument, "solve_ode_rkf45: initial_step must be > 0"));
    }

    auto rhs = [&](double t, double y) -> Result<double> {
        return eval(expr, {{t_var, t}, {y_var, y}});
    };

    std::vector<OdePoint> traj;
    traj.push_back({t0, y0});

    double t = t0;
    double y = y0;
    double h = options.initial_step;
    if (options.max_step > 0.0 && h > options.max_step) h = options.max_step;

    for (std::uint32_t step = 0; step < options.max_steps; ++step) {
        if (t >= t_end) break;
        if (t + h > t_end) h = t_end - t;

        auto k1r = rhs(t, y);
        if (k1r.is_error()) return fail<std::vector<OdePoint>>(k1r.error());
        const double k1 = h * k1r.value();

        auto k2r = rhs(t + Tab::c2 * h, y + Tab::a21 * k1);
        if (k2r.is_error()) return fail<std::vector<OdePoint>>(k2r.error());
        const double k2 = h * k2r.value();

        auto k3r = rhs(t + Tab::c3 * h, y + Tab::a31 * k1 + Tab::a32 * k2);
        if (k3r.is_error()) return fail<std::vector<OdePoint>>(k3r.error());
        const double k3 = h * k3r.value();

        auto k4r = rhs(t + Tab::c4 * h,
                       y + Tab::a41 * k1 + Tab::a42 * k2 + Tab::a43 * k3);
        if (k4r.is_error()) return fail<std::vector<OdePoint>>(k4r.error());
        const double k4 = h * k4r.value();

        auto k5r = rhs(t + Tab::c5 * h,
                       y + Tab::a51 * k1 + Tab::a52 * k2 + Tab::a53 * k3 + Tab::a54 * k4);
        if (k5r.is_error()) return fail<std::vector<OdePoint>>(k5r.error());
        const double k5 = h * k5r.value();

        auto k6r = rhs(t + Tab::c6 * h,
                       y + Tab::a61 * k1 + Tab::a62 * k2 + Tab::a63 * k3
                         + Tab::a64 * k4 + Tab::a65 * k5);
        if (k6r.is_error()) return fail<std::vector<OdePoint>>(k6r.error());
        const double k6 = h * k6r.value();

        const double y4 = y + Tab::b4_1 * k1 + Tab::b4_3 * k3 + Tab::b4_4 * k4 + Tab::b4_5 * k5;
        const double y5 = y + Tab::b5_1 * k1 + Tab::b5_3 * k3 + Tab::b5_4 * k4 + Tab::b5_5 * k5 + Tab::b5_6 * k6;
        const double err = std::abs(y5 - y4);
        const double tol = options.abs_tol + options.rel_tol * std::max(std::abs(y), std::abs(y5));

        if (err <= tol || h <= 1e-15) {
            t += h;
            y = y5;
            traj.push_back({t, y});
            // Step size adjustment for next step.
            const double factor = (err > 0.0) ? 0.9 * std::pow(tol / err, 0.2) : 5.0;
            h *= std::min(5.0, std::max(0.1, factor));
            if (options.max_step > 0.0 && h > options.max_step) h = options.max_step;
        } else {
            // Reject step, shrink and retry.
            const double factor = 0.9 * std::pow(tol / err, 0.25);
            h *= std::max(0.1, factor);
        }
    }
    if (t < t_end) {
        return fail<std::vector<OdePoint>>(make_error(
            CASErrorKind::Timeout,
            "solve_ode_rkf45: max_steps reached before t_end"));
    }
    return ok(std::move(traj));
}

}  // namespace cas::numeric
