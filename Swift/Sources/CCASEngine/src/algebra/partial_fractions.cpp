#include "cas/algebra.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <vector>

namespace cas::algebra {

[[nodiscard]] ExprPtr make_rational_expr(AstArena& arena, const Rational& value) {
    if (value.denominator() == BigInt(1)) {
        return arena.make<IntegerLit>(value.numerator());
    }
    return arena.make<RationalLit>(value.numerator(), value.denominator());
}

[[nodiscard]] static Rational rational_pow_nonnegative(Rational base, std::size_t exponent) {
    Rational result(BigInt(1));
    while (exponent > 0U) {
        if ((exponent & 1U) != 0U) {
            result *= base;
        }
        exponent >>= 1U;
        if (exponent > 0U) {
            base *= base;
        }
    }
    return result;
}

[[nodiscard]] static Result<std::vector<Rational>> solve_dense_rational_system(
    std::vector<std::vector<Rational>> matrix,
    std::vector<Rational> rhs) {
    const std::size_t size = matrix.size();
    if (rhs.size() != size) {
        return fail<std::vector<Rational>>(make_error(
            CASErrorKind::InvalidArgument,
            "Il sistema lineare razionale ha dimensioni incoerenti"));
    }

    for (std::size_t row = 0; row < size; ++row) {
        if (matrix[row].size() != size) {
            return fail<std::vector<Rational>>(make_error(
                CASErrorKind::InvalidArgument,
                "Il sistema lineare razionale deve essere quadrato"));
        }
    }

    for (std::size_t pivot = 0; pivot < size; ++pivot) {
        std::size_t pivot_row = pivot;
        while (pivot_row < size && matrix[pivot_row][pivot].numerator().is_zero()) {
            ++pivot_row;
        }
        if (pivot_row == size) {
            return fail<std::vector<Rational>>(make_error(
                CASErrorKind::Unimplemented,
                "Il sistema lineare della decomposizione parziale e' singolare"));
        }
        if (pivot_row != pivot) {
            std::swap(matrix[pivot], matrix[pivot_row]);
            std::swap(rhs[pivot], rhs[pivot_row]);
        }

        auto inverse_pivot = checked_divide(Rational(BigInt(1)), matrix[pivot][pivot]);
        if (inverse_pivot.is_error()) {
            return fail<std::vector<Rational>>(inverse_pivot.error());
        }

        for (std::size_t column = pivot; column < size; ++column) {
            matrix[pivot][column] *= inverse_pivot.value();
        }
        rhs[pivot] *= inverse_pivot.value();

        for (std::size_t row = 0; row < size; ++row) {
            if (row == pivot || matrix[row][pivot].numerator().is_zero()) {
                continue;
            }

            const Rational factor = matrix[row][pivot];
            for (std::size_t column = pivot; column < size; ++column) {
                matrix[row][column] -= factor * matrix[pivot][column];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }

    return ok(std::move(rhs));
}

[[nodiscard]] static std::vector<Rational> choose_regular_rational_samples(
    const IntPoly& denominator,
    std::size_t count) {
    std::vector<Rational> samples;
    samples.reserve(count);

    for (std::size_t step = 0; samples.size() < count; ++step) {
        const long long magnitude = static_cast<long long>(step / 2U);
        const long long candidate = (step % 2U == 0U) ? magnitude : -magnitude - 1LL;
        const Rational sample{BigInt(candidate)};
        if (!evaluate_integer_polynomial_at(denominator, sample).numerator().is_zero()) {
            samples.push_back(sample);
        }
    }

    return samples;
}

[[nodiscard]] static Result<ExprPtr> build_factor_power_expr(
    ExprPtr factor,
    unsigned int power,
    symbolic::CASContext& ctx) {
    if (power == 1U) {
        return ok(factor);
    }
    return pow_expr(factor, static_cast<std::size_t>(power), ctx);
}

[[nodiscard]] static Result<ExprPtr> build_quadratic_numerator(
    const Rational& coeff_const,
    const Rational& coeff_x,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    std::vector<ExprPtr> sum_parts;
    if (!coeff_x.numerator().is_zero()) {
        auto x_term = multiply_exprs(
            make_rational_expr(ctx.arena(), coeff_x),
            ctx.arena().make<Symbol>(var.name),
            ctx);
        if (x_term.is_error()) return fail<ExprPtr>(x_term.error());
        sum_parts.push_back(x_term.value());
    }
    if (!coeff_const.numerator().is_zero()) {
        sum_parts.push_back(make_rational_expr(ctx.arena(), coeff_const));
    }
    if (sum_parts.empty()) return ok(make_integer(ctx.arena(), 0));
    if (sum_parts.size() == 1U) return ok(sum_parts[0]);
    return simplify_expr(ctx.arena().make<Sum>(std::move(sum_parts)), ctx);
}

struct PFFactorInfo {
    ExprPtr factor_expr;
    IntPoly coefficients;
    unsigned int multiplicity{1U};
    bool is_quadratic{false};
};

struct PFBasisElement {
    std::size_t factor_index{0U};
    unsigned int power{1U};
    bool is_x_coeff{false};
};

Result<std::vector<ExprPtr>> partial_fractions(
    ExprPtr rational_expr,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    if (!rational_expr) {
        return fail<std::vector<ExprPtr>>(make_error(
            CASErrorKind::InvalidArgument,
            "partial_fractions richiede un'espressione razionale non nulla"));
    }

    auto parts = apart_num_den(rational_expr, ctx);
    if (parts.is_error()) {
        return fail<std::vector<ExprPtr>>(parts.error());
    }

    auto numerator_poly = parse_polynomial(parts.value().numerator, var, ctx);
    if (numerator_poly.is_error()) {
        return fail<std::vector<ExprPtr>>(numerator_poly.error());
    }
    auto denominator_poly = parse_polynomial(parts.value().denominator, var, ctx);
    if (denominator_poly.is_error()) {
        return fail<std::vector<ExprPtr>>(denominator_poly.error());
    }

    if (is_zero_poly(denominator_poly.value())) {
        return fail<std::vector<ExprPtr>>(make_error(
            CASErrorKind::Undefined,
            "partial_fractions ha ricevuto un denominatore nullo"));
    }
    if (is_zero_poly(numerator_poly.value())) {
        return ok(std::vector<ExprPtr>{make_integer(ctx.arena(), 0)});
    }
    if (poly_degree(numerator_poly.value()) >= poly_degree(denominator_poly.value())) {
        return fail<std::vector<ExprPtr>>(make_error(
            CASErrorKind::InvalidArgument,
            "partial_fractions richiede una funzione razionale propria"));
    }

    auto denominator_factors = factor_over_integers(parts.value().denominator, var, ctx);
    if (denominator_factors.is_error()) {
        return fail<std::vector<ExprPtr>>(denominator_factors.error());
    }

    auto numerator_coefficients = poly_to_rational_poly(numerator_poly.value());
    if (numerator_coefficients.is_error()) {
        return fail<std::vector<ExprPtr>>(numerator_coefficients.error());
    }
    auto denominator_coefficients = poly_to_integer_poly(denominator_poly.value());
    if (denominator_coefficients.is_error()) {
        return fail<std::vector<ExprPtr>>(denominator_coefficients.error());
    }

    std::vector<PFFactorInfo> all_factors;
    all_factors.reserve(denominator_factors.value().factors.size());

    for (const PolynomialFactor& pf : denominator_factors.value().factors) {
        auto factor_poly = parse_polynomial(pf.factor, var, ctx);
        if (factor_poly.is_error()) {
            return fail<std::vector<ExprPtr>>(factor_poly.error());
        }
        auto factor_int = poly_to_integer_poly(factor_poly.value());
        if (factor_int.is_error()) {
            return fail<std::vector<ExprPtr>>(factor_int.error());
        }
        const std::size_t factor_deg = factor_int.value().size();
        if (factor_deg != 2U && factor_deg != 3U) {
            return fail<std::vector<ExprPtr>>(make_error(
                CASErrorKind::Unimplemented,
                "partial_fractions supporta solo fattori lineari e quadratici irriducibili"));
        }
        all_factors.push_back(PFFactorInfo{
            .factor_expr = pf.factor,
            .coefficients = std::move(factor_int.value()),
            .multiplicity = pf.multiplicity,
            .is_quadratic = (factor_deg == 3U),
        });
    }

    if (all_factors.empty()) {
        return fail<std::vector<ExprPtr>>(make_error(
            CASErrorKind::Unimplemented,
            "partial_fractions richiede almeno un fattore nel denominatore"));
    }

    std::vector<PFBasisElement> basis;
    for (std::size_t fi = 0; fi < all_factors.size(); ++fi) {
        for (unsigned int pw = 1U; pw <= all_factors[fi].multiplicity; ++pw) {
            if (!all_factors[fi].is_quadratic) {
                basis.push_back(PFBasisElement{fi, pw, false});
            } else {
                basis.push_back(PFBasisElement{fi, pw, false});
                basis.push_back(PFBasisElement{fi, pw, true});
            }
        }
    }

    std::vector<Rational> samples = choose_regular_rational_samples(
        denominator_coefficients.value(),
        basis.size());

    std::vector<std::vector<Rational>> matrix(
        basis.size(),
        std::vector<Rational>(basis.size(), Rational(BigInt(0))));
    std::vector<Rational> rhs;
    rhs.reserve(basis.size());

    for (std::size_t row = 0; row < samples.size(); ++row) {
        const Rational& sample = samples[row];
        const Rational num_val = evaluate_rational_polynomial_at(
            numerator_coefficients.value(), sample);
        const Rational den_val = evaluate_integer_polynomial_at(
            denominator_coefficients.value(), sample);
        auto rhs_entry = checked_divide(num_val, den_val);
        if (rhs_entry.is_error()) {
            return fail<std::vector<ExprPtr>>(rhs_entry.error());
        }
        rhs.push_back(rhs_entry.value());

        for (std::size_t col = 0; col < basis.size(); ++col) {
            const PFBasisElement& be = basis[col];
            const PFFactorInfo& fi = all_factors[be.factor_index];
            const Rational factor_pow = rational_pow_nonnegative(
                evaluate_integer_polynomial_at(fi.coefficients, sample),
                static_cast<std::size_t>(be.power));
            const Rational numerator_factor = be.is_x_coeff ? sample : Rational(BigInt(1));
            auto col_entry = checked_divide(numerator_factor, factor_pow);
            if (col_entry.is_error()) {
                return fail<std::vector<ExprPtr>>(col_entry.error());
            }
            matrix[row][col] = col_entry.value();
        }
    }

    auto solved = solve_dense_rational_system(std::move(matrix), std::move(rhs));
    if (solved.is_error()) {
        return fail<std::vector<ExprPtr>>(solved.error());
    }

    std::vector<ExprPtr> terms;

    std::size_t basis_idx = 0;
    for (std::size_t fi = 0; fi < all_factors.size(); ++fi) {
        for (unsigned int pw = 1U; pw <= all_factors[fi].multiplicity; ++pw) {
            if (!all_factors[fi].is_quadratic) {
                const Rational coeff = solved.value()[basis_idx++];
                if (coeff.numerator().is_zero()) continue;

                auto denom_expr = build_factor_power_expr(all_factors[fi].factor_expr, pw, ctx);
                if (denom_expr.is_error()) return fail<std::vector<ExprPtr>>(denom_expr.error());

                auto term = divide_exprs(
                    make_rational_expr(ctx.arena(), coeff),
                    denom_expr.value(), ctx);
                if (term.is_error()) return fail<std::vector<ExprPtr>>(term.error());
                terms.push_back(term.value());
            } else {
                const Rational coeff_const = solved.value()[basis_idx++];
                const Rational coeff_x    = solved.value()[basis_idx++];
                if (coeff_const.numerator().is_zero() && coeff_x.numerator().is_zero()) continue;

                auto num_expr = build_quadratic_numerator(coeff_const, coeff_x, var, ctx);
                if (num_expr.is_error()) return fail<std::vector<ExprPtr>>(num_expr.error());

                auto denom_expr = build_factor_power_expr(all_factors[fi].factor_expr, pw, ctx);
                if (denom_expr.is_error()) return fail<std::vector<ExprPtr>>(denom_expr.error());

                auto term = divide_exprs(num_expr.value(), denom_expr.value(), ctx);
                if (term.is_error()) return fail<std::vector<ExprPtr>>(term.error());
                terms.push_back(term.value());
            }
        }
    }

    if (terms.empty()) {
        terms.push_back(make_integer(ctx.arena(), 0));
    }
    return ok(std::move(terms));
}

} // namespace cas::algebra
