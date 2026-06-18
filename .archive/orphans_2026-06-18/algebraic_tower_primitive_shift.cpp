// algebraic_tower_primitive_shift.cpp — compute_shift_resultant implementation.
// Split from algebraic_tower_primitive_internal.hpp (ANTI-MONOLITH rule).
// Included translation units: this file only.
#include "algebraic_tower_primitive_internal.hpp"
#include "algebra_internal.hpp"

namespace cas {
namespace algebra {
namespace primitive_internal {

Result<RatPoly> compute_shift_resultant(
    const RatPoly& q_prev,
    const RatPoly& m_k,
    const BigInt& s,
    const Deadline& deadline) {

    if (q_prev.is_zero() || m_k.is_zero()) {
        return fail<RatPoly>(CASError{
            CASErrorKind::InvalidArgument,
            "compute_shift_resultant: zero input polynomial",
            std::nullopt});
    }

    const std::size_t deg_R = q_prev.degree() * m_k.degree();
    const std::size_t num_pts = deg_R + 1U;

    std::vector<Rational> eval_pts;
    eval_pts.reserve(num_pts);
    for (std::size_t j = 1U; j <= num_pts; ++j) {
        eval_pts.push_back(Rational(BigInt(static_cast<std::int64_t>(j))));
    }

    std::vector<Rational> values;
    values.reserve(num_pts);

    const Rational s_rat(s);
    const Rational neg_s_rat = Rational(BigInt(0)) - s_rat;
    const std::size_t dq = q_prev.degree();

    // Precompute binomial coefficients C(dq, k).
    std::vector<std::vector<BigInt>> binom_table(dq + 1U, std::vector<BigInt>(dq + 1U, BigInt(0)));
    for (std::size_t row = 0U; row <= dq; ++row) {
        binom_table[row][0] = BigInt(1);
        for (std::size_t col = 1U; col <= row; ++col) {
            binom_table[row][col] = binom_table[row - 1U][col - 1U] + binom_table[row - 1U][col];
        }
    }

    for (const Rational& y_j : eval_pts) {
        if (deadline_exceeded(deadline)) {
            return fail<RatPoly>(CASError{
                CASErrorKind::Unimplemented,
                "compute_shift_resultant: ctx.timeout() exceeded during "
                "evaluation-interpolation",
                std::nullopt});
        }
        // Build Q_j(x) = q_prev(y_j − s·x).
        std::vector<Rational> Q_j_coeffs(q_prev.size(), Rational(BigInt(0)));

        std::vector<Rational> y_pows(q_prev.size());
        std::vector<Rational> neg_s_pows(q_prev.size());
        y_pows[0] = Rational(BigInt(1));
        neg_s_pows[0] = Rational(BigInt(1));
        for (std::size_t p = 1U; p < q_prev.size(); ++p) {
            y_pows[p] = y_pows[p - 1U] * y_j;
            neg_s_pows[p] = neg_s_pows[p - 1U] * neg_s_rat;
        }

        for (std::size_t i = 0U; i < q_prev.size(); ++i) {
            const Rational& ci = q_prev[i];
            if (ci.numerator().is_zero()) continue;
            for (std::size_t k = 0U; k <= i; ++k) {
                const Rational term = ci
                    * Rational(binom_table[i][k])
                    * y_pows[i - k]
                    * neg_s_pows[k];
                Q_j_coeffs[k] = Q_j_coeffs[k] + term;
            }
        }
        RatPoly Q_j(Q_j_coeffs);
        Q_j.normalize([](const Rational& r) { return r.numerator().is_zero(); });

        auto res_scalar = resultant_generic<Rational>(
            m_k.coefficients(),
            Q_j.coefficients(),
            nullptr,
            deadline);
        if (res_scalar.is_error()) return fail<RatPoly>(res_scalar.error());
        values.push_back(res_scalar.value());
    }

    // Newton interpolation: divided differences.
    std::vector<std::vector<Rational>> table(num_pts, std::vector<Rational>(num_pts, Rational(BigInt(0))));
    for (std::size_t i = 0U; i < num_pts; ++i) table[i][0] = values[i];
    for (std::size_t j = 1U; j < num_pts; ++j) {
        if (deadline_exceeded(deadline)) {
            return fail<RatPoly>(CASError{
                CASErrorKind::Unimplemented,
                "compute_shift_resultant: ctx.timeout() exceeded during "
                "Newton interpolation",
                std::nullopt});
        }
        for (std::size_t i = j; i < num_pts; ++i) {
            const Rational denom = eval_pts[i] - eval_pts[i - j];
            if (denom.numerator().is_zero()) {
                return fail<RatPoly>(CASError{
                    CASErrorKind::InternalError,
                    "compute_shift_resultant: duplicate evaluation points",
                    std::nullopt});
            }
            const Rational denom_inv{denom.denominator(), denom.numerator()};
            table[i][j] = (table[i][j - 1U] - table[i - 1U][j - 1U]) * denom_inv;
        }
    }

    // Expand Newton form into monomial basis.
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
