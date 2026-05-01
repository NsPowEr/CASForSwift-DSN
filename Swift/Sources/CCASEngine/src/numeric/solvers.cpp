#include "cas/numeric.hpp"
#include "cas/calculus.hpp"
#include "cas/symbolic.hpp"
#include <cmath>

namespace cas::numeric {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

} // namespace

Result<double> solve_numeric_bisection(
    ExprPtr expr,
    const std::string& variable,
    double low,
    double high,
    const RootFindingOptions& options) {

    auto f = [&](double x) -> Result<double> {
        return eval(expr, {{variable, x}});
    };

    auto f_low = f(low);
    if (f_low.is_error()) return f_low;
    auto f_high = f(high);
    if (f_high.is_error()) return f_high;

    if (f_low.value() * f_high.value() > 0) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument, "Function must have different signs at interval boundaries"));
    }

    double a = low;
    double b = high;
    for (std::uint32_t i = 0; i < options.max_iterations; ++i) {
        double mid = (a + b) / 2.0;
        if ((b - a) / 2.0 < options.tolerance) {
            return ok(mid);
        }

        auto f_mid = f(mid);
        if (f_mid.is_error()) return f_mid;

        if (std::abs(f_mid.value()) < 1e-15) return ok(mid);

        if ((f_mid.value() > 0) == (f_low.value() > 0)) {
            a = mid;
            f_low = f_mid;
        } else {
            b = mid;
        }
    }

    return ok((a + b) / 2.0);
}

Result<double> solve_numeric_newton(
    ExprPtr expr,
    const std::string& variable,
    double initial_guess,
    symbolic::CASContext& ctx,
    const RootFindingOptions& options) {

    // Calcolo derivata simbolica: Muscolo Matematico (Integrazione tra moduli)
    auto var_node = Symbol(variable);
    auto derivative_res = calculus::diff(expr, var_node, 1, ctx);
    if (derivative_res.is_error()) return fail<double>(derivative_res.error());
    ExprPtr derivative = derivative_res.value();

    double x = initial_guess;
    for (std::uint32_t i = 0; i < options.max_iterations; ++i) {
        auto fx = eval(expr, {{variable, x}});
        if (fx.is_error()) return fx;
        
        if (std::abs(fx.value()) < options.tolerance) {
            return ok(x);
        }

        auto dfx = eval(derivative, {{variable, x}});
        if (dfx.is_error()) return dfx;

        if (std::abs(dfx.value()) < 1e-15) {
            return fail<double>(make_error(CASErrorKind::Undefined, "Newton-Raphson: Derivative too close to zero"));
        }

        double next_x = x - fx.value() / dfx.value();
        if (std::abs(next_x - x) < options.tolerance) {
            return ok(next_x);
        }
        x = next_x;
    }

    return ok(x);
}

} // namespace cas::numeric
