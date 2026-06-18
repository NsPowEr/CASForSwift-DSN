// algebraic_tower_primitive_absolute.cpp — compute_absolute_resultant_xy implementation.
// Split from algebraic_tower_primitive_internal.hpp (ANTI-MONOLITH rule).
// Included translation units: this file only.
#include "algebraic_tower_primitive_internal.hpp"
#include "algebra_internal.hpp"

namespace cas {
namespace algebra {
namespace primitive_internal {

Result<RatPoly> compute_absolute_resultant_xy(
    const RatPoly& g_y,
    const std::vector<RatPoly>& f_xy,
    const Deadline& deadline) {

    if (g_y.is_zero()) {
        return fail<RatPoly>(CASError{
            CASErrorKind::InvalidArgument,
            "compute_absolute_resultant_xy: zero inner minimal polynomial",
            std::nullopt});
    }
    if (f_xy.empty()) {
        return fail<RatPoly>(CASError{
            CASErrorKind::InvalidArgument,
            "compute_absolute_resultant_xy: empty bivariate polynomial",
            std::nullopt});
    }
    const std::size_t deg_x = f_xy.size() - 1U;
    const std::size_t deg_y_g = g_y.degree();
    const std::size_t deg_R = deg_x * deg_y_g;
    const std::size_t num_pts = deg_R + 1U;

    std::vector<Rational> eval_pts;
    eval_pts.reserve(num_pts);
    for (std::size_t j = 1U; j <= num_pts; ++j) {
        eval_pts.push_back(Rational(BigInt(static_cast<std::int64_t>(j))));
    }

    std::vector<Rational> values;
    values.reserve(num_pts);

    for (const Rational& x_j : eval_pts) {
        if (deadline_exceeded(deadline)) {
            return fail<RatPoly>(CASError{
                CASErrorKind::Unimplemented,
                "compute_absolute_resultant_xy: ctx.timeout() exceeded during "
                "evaluation-interpolation",
                std::nullopt});
        }
        // Build f_y(y) = Σ_i (x_j)^i · f_xy[i](y) ∈ Q[y].
        // Find max degree across f_xy[i].
        std::size_t max_deg_y = 0U;
        for (const auto& c_y : f_xy) {
            if (!c_y.is_zero() && c_y.degree() > max_deg_y) max_deg_y = c_y.degree();
        }
        std::vector<Rational> f_y_coeffs(max_deg_y + 1U, Rational(BigInt(0)));

        Rational x_pow = Rational(BigInt(1));
        for (std::size_t i = 0U; i <= deg_x; ++i) {
            const RatPoly& c_y = f_xy[i];
            for (std::size_t d2 = 0U; d2 < c_y.size(); ++d2) {
                f_y_coeffs[d2] = f_y_coeffs[d2] + x_pow * c_y[d2];
            }
            x_pow = x_pow * x_j;
        }
        RatPoly f_y(f_y_coeffs);
        f_y.normalize([](const Rational& r) { return r.numerator().is_zero(); });
        if (f_y.is_zero()) {
            // f(x_j, y) ≡ 0 ⇒ Res = 0 (degenerate evaluation point).
            values.push_back(Rational(BigInt(0)));
            continue;
        }

        auto res_scalar = resultant_generic<Rational>(
            g_y.coefficients(),
            f_y.coefficients(),
            nullptr,
            deadline);
        if (res_scalar.is_error()) return fail<RatPoly>(res_scalar.error());
        values.push_back(res_scalar.value());
    }

    // Newton interpolation (same divided-differences scheme as compute_shift_resultant).
    std::vector<std::vector<Rational>> table(num_pts, std::vector<Rational>(num_pts, Rational(BigInt(0))));
    for (std::size_t i = 0U; i < num_pts; ++i) table[i][0] = values[i];
    for (std::size_t j = 1U; j < num_pts; ++j) {
        if (deadline_exceeded(deadline)) {
            return fail<RatPoly>(CASError{
                CASErrorKind::Unimplemented,
                "compute_absolute_resultant_xy: ctx.timeout() exceeded during "
                "Newton interpolation",
                std::nullopt});
        }
        for (std::size_t i = j; i < num_pts; ++i) {
            const Rational denom = eval_pts[i] - eval_pts[i - j];
            if (denom.numerator().is_zero()) {
                return fail<RatPoly>(CASError{
                    CASErrorKind::InternalError,
                    "compute_absolute_resultant_xy: duplicate evaluation points",
                    std::nullopt});
            }
            const Rational denom_inv{denom.denominator(), denom.numerator()};
            table[i][j] = (table[i][j - 1U] - table[i - 1U][j - 1U]) * denom_inv;
        }
    }

    std::vector<Rational> result_coeffs(num_pts, Rational(BigInt(0)));
    std::vector<Rational> newton_poly = {Rational(BigInt(1))};
    for (std::size_t j = 0U; j < num_pts; ++j) {
        const Rational coeff_j = table[j][j];
        if (!coeff_j.numerator().is_zero()) {
            for (std::size_t d = 0U; d < newton_poly.size(); ++d) {
                result_coeffs[d] = result_coeffs[d] + coeff_j * newton_poly[d];
            }
        }
        if (j + 1U < num_pts) {
            std::vector<Rational> new_newton(newton_poly.size() + 1U, Rational(BigInt(0)));
            for (std::size_t d = 0U; d < newton_poly.size(); ++d) {
                new_newton[d + 1U] = new_newton[d + 1U] + newton_poly[d];
                new_newton[d] = new_newton[d] - eval_pts[j] * newton_poly[d];
            }
            while (!new_newton.empty() && new_newton.back().numerator().is_zero())
                new_newton.pop_back();
            newton_poly = std::move(new_newton);
        }
    }

    RatPoly result(result_coeffs);
    result.normalize([](const Rational& r) { return r.numerator().is_zero(); });
    return ok(std::move(result));
}

}  // namespace primitive_internal
}  // namespace algebra
}  // namespace cas
