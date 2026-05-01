#include "cas/numeric.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

namespace cas::numeric {

Result<std::vector<OdePoint>> solve_ode_rk4(
    ExprPtr expr,
    const std::string& t_var,
    const std::string& y_var,
    double t0,
    double y0,
    double t_end,
    double initial_step) {

    auto f = [&](double t, double y) -> Result<double> {
        return eval(expr, {{t_var, t}, {y_var, y}});
    };

    std::vector<OdePoint> result;
    double t = t0;
    double y = y0;
    double h = initial_step;
    const double tol = 1e-6;
    const double min_h = 1e-10;
    const double max_h = 0.5;

    result.push_back({t, y});

    while (t < t_end) {
        if (t + h > t_end) h = t_end - t;

        auto k1_res = f(t, y); if (k1_res.is_error()) return fail<std::vector<OdePoint>>(k1_res.error());
        double k1 = k1_res.value();

        auto k2_res = f(t + h/4.0, y + h*k1/4.0); if (k2_res.is_error()) return fail<std::vector<OdePoint>>(k2_res.error());
        double k2 = k2_res.value();

        auto k3_res = f(t + 3.0*h/8.0, y + 3.0*h*k1/32.0 + 9.0*h*k2/32.0); if (k3_res.is_error()) return fail<std::vector<OdePoint>>(k3_res.error());
        double k3 = k3_res.value();

        auto k4_res = f(t + 12.0*h/13.0, y + 1932.0*h*k1/2197.0 - 7200.0*h*k2/2197.0 + 7296.0*h*k3/2197.0); if (k4_res.is_error()) return fail<std::vector<OdePoint>>(k4_res.error());
        double k4 = k4_res.value();

        auto k5_res = f(t + h, y + 439.0*h*k1/216.0 - 8.0*h*k2 + 3680.0*h*k3/513.0 - 845.0*h*k4/4104.0); if (k5_res.is_error()) return fail<std::vector<OdePoint>>(k5_res.error());
        double k5 = k5_res.value();

        auto k6_res = f(t + h/2.0, y - 8.0*h*k1/27.0 + 2.0*h*k2 - 3544.0*h*k3/2565.0 + 1859.0*h*k4/4104.0 - 11.0*h*k5/40.0); if (k6_res.is_error()) return fail<std::vector<OdePoint>>(k6_res.error());
        double k6 = k6_res.value();

        // RK4 step
        double y4 = y + h * (25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0);
        // RK5 step
        double y5 = y + h * (16.0*k1/135.0 + 6656.0*k3/12825.0 + 28561.0*k4/56430.0 - 9.0*k5/50.0 + 2.0*k6/55.0);

        double error = std::abs(y5 - y4);
        double s = 0.84 * std::pow(tol * h / error, 0.25);

        if (error <= tol * h || h <= min_h) {
            t += h;
            y = y5;
            result.push_back({t, y});
        }

        h = std::clamp(h * s, min_h, max_h);
        if (std::isinf(y) || std::isnan(y)) break;
    }

    return ok(std::move(result));
}

} // namespace cas::numeric
