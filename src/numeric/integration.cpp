#include "cas/numeric.hpp"
#include <cmath>

namespace cas::numeric {

namespace {

Result<double> simpson_step(
    double a,
    double b,
    double fa,
    double fb,
    double fm) {
    
    return ok((b - a) / 6.0 * (fa + 4.0 * fm + fb));
}

Result<double> adaptive_simpson(
    const std::function<Result<double>(double)>& f,
    double a,
    double b,
    double tolerance,
    double whole,
    double fa,
    double fb,
    double fm,
    std::uint32_t depth) {

    double mid = (a + b) / 2.0;
    double left_mid = (a + mid) / 2.0;
    double right_mid = (mid + b) / 2.0;

    auto f_lm = f(left_mid);
    if (f_lm.is_error()) return f_lm;
    auto f_rm = f(right_mid);
    if (f_rm.is_error()) return f_rm;

    auto left_res = simpson_step(a, mid, fa, fm, f_lm.value());
    if (left_res.is_error()) return left_res;
    auto right_res = simpson_step(mid, b, fm, fb, f_rm.value());
    if (right_res.is_error()) return right_res;

    double left = left_res.value();
    double right = right_res.value();

    if (depth == 0 || std::abs(left + right - whole) <= 15.0 * tolerance) {
        return ok(left + right + (left + right - whole) / 15.0);
    }

    auto left_total = adaptive_simpson(f, a, mid, tolerance / 2.0, left, fa, fm, f_lm.value(), depth - 1);
    if (left_total.is_error()) return left_total;
    auto right_total = adaptive_simpson(f, mid, b, tolerance / 2.0, right, fm, fb, f_rm.value(), depth - 1);
    if (right_total.is_error()) return right_total;

    return ok(left_total.value() + right_total.value());
}

} // namespace

Result<double> integrate_numeric(
    ExprPtr expr,
    const std::string& variable,
    double a,
    double b,
    const IntegrationOptions& options) {

    auto f = [&](double x) -> Result<double> {
        return eval(expr, {{variable, x}});
    };

    auto fa = f(a); if (fa.is_error()) return fa;
    auto fb = f(b); if (fb.is_error()) return fb;
    auto fm = f((a + b) / 2.0); if (fm.is_error()) return fm;

    auto whole = simpson_step(a, b, fa.value(), fb.value(), fm.value());
    if (whole.is_error()) return whole;

    return adaptive_simpson(f, a, b, options.tolerance, whole.value(), fa.value(), fb.value(), fm.value(), options.max_depth);
}

} // namespace cas::numeric
