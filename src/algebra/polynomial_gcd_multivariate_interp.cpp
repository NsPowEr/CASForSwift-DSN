// polynomial_gcd_multivariate_interp.cpp
// Lagrange interpolation and polynomial-value interpolation helpers used by the
// multivariate GCD evaluation/interpolation algorithm.
// Declarations live in polynomial_gcd_multivariate_helpers.hpp.

#include "polynomial_gcd_multivariate_helpers.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <set>
#include <vector>

namespace cas::algebra {

[[nodiscard]] Result<std::vector<Rational>> lagrange_interpolate(
    const std::vector<BigInt>& values,
    const std::vector<BigInt>& points) {
    if (values.size() != points.size()) {
        return fail<std::vector<Rational>>(make_error(
            CASErrorKind::InvalidArgument,
            "polynomial_gcd_multivariate: interpolation points/value size mismatch"));
    }
    if (values.empty()) {
        return ok(std::vector<Rational>{});
    }

    const std::size_t n = values.size();
    std::vector<Rational> result(n, Rational(BigInt(0)));

    for (std::size_t i = 0U; i < n; ++i) {
        if (values[i].is_zero()) {
            continue;
        }

        std::vector<Rational> basis(1U, Rational(BigInt(1)));
        for (std::size_t j = 0U; j < n; ++j) {
            if (j == i) {
                continue;
            }

            const BigInt denominator = points[i] - points[j];
            if (denominator.is_zero()) {
                return fail<std::vector<Rational>>(make_error(
                    CASErrorKind::InvalidArgument,
                    "polynomial_gcd_multivariate: duplicate interpolation points"));
            }

            const Rational scale(BigInt(1), denominator);
            std::vector<Rational> next(basis.size() + 1U, Rational(BigInt(0)));
            for (std::size_t k = 0U; k < basis.size(); ++k) {
                next[k + 1U] = next[k + 1U] + basis[k] * scale;
                next[k] = next[k] - basis[k] * scale * Rational(points[j]);
            }
            basis = std::move(next);
        }

        const Rational value(values[i]);
        if (basis.size() > result.size()) {
            result.resize(basis.size(), Rational(BigInt(0)));
        }
        for (std::size_t k = 0U; k < basis.size(); ++k) {
            result[k] = result[k] + value * basis[k];
        }
    }

    while (!result.empty() && result.back().numerator().is_zero()) {
        result.pop_back();
    }
    return ok(std::move(result));
}

[[nodiscard]] MultivariatePolynomial coefficient_poly_in_var(
    const MultivariatePolynomial& poly,
    const Symbol& var,
    std::size_t degree) {
    std::vector<MultivariateTerm> terms;
    for (const auto& term : poly.terms()) {
        unsigned int exponent = 0U;
        for (const auto& [symbol, factor_exponent] : term.factors) {
            if (symbol.name == var.name) {
                exponent = factor_exponent;
                break;
            }
        }

        if (static_cast<std::size_t>(exponent) != degree) {
            continue;
        }

        std::vector<std::pair<Symbol, unsigned int>> reduced_factors;
        reduced_factors.reserve(term.factors.size());
        for (const auto& [symbol, factor_exponent] : term.factors) {
            if (symbol.name != var.name) {
                reduced_factors.emplace_back(symbol, factor_exponent);
            }
        }

        terms.push_back(MultivariateTerm{
            .coefficient = term.coefficient,
            .factors = std::move(reduced_factors),
        });
    }
    return MultivariatePolynomial(std::move(terms));
}

[[nodiscard]] MultivariatePolynomial multiply_by_variable_power(
    const MultivariatePolynomial& poly,
    const Symbol& var,
    std::size_t power) {
    if (power == 0U || poly.is_zero()) {
        return poly;
    }

    std::vector<MultivariateTerm> terms;
    terms.reserve(poly.terms().size());
    for (const auto& term : poly.terms()) {
        auto factors = term.factors;
        factors.emplace_back(var, static_cast<unsigned int>(power));
        terms.push_back(MultivariateTerm{
            .coefficient = term.coefficient,
            .factors = std::move(factors),
        });
    }
    return MultivariatePolynomial(std::move(terms));
}

[[nodiscard]] Result<MultivariatePolynomial> interpolate_polynomial_values(
    const std::vector<MultivariatePolynomial>& values,
    const std::vector<BigInt>& points,
    const Symbol& interpolation_var) {
    if (values.size() != points.size()) {
        return fail<MultivariatePolynomial>(make_error(
            CASErrorKind::InvalidArgument,
            "polynomial_gcd_multivariate: interpolation polynomial points/value size mismatch"));
    }

    std::set<FactorKey, FactorKeyLess> all_keys;
    std::vector<CoeffMap> value_maps;
    value_maps.reserve(values.size());

    for (const auto& value_poly : values) {
        auto map = to_coeff_map(value_poly);
        for (const auto& [key, coefficient] : map) {
            static_cast<void>(coefficient);
            all_keys.insert(key);
        }
        value_maps.push_back(std::move(map));
    }

    std::vector<MultivariateTerm> terms;

    for (const auto& key : all_keys) {
        std::vector<BigInt> sequence(values.size(), BigInt(0));
        for (std::size_t i = 0U; i < value_maps.size(); ++i) {
            const auto it = value_maps[i].find(key);
            if (it != value_maps[i].end()) {
                sequence[i] = it->second;
            }
        }

        auto interpolation = lagrange_interpolate(sequence, points);
        if (interpolation.is_error()) {
            return fail<MultivariatePolynomial>(interpolation.error());
        }

        const auto& coeffs = interpolation.value();
        for (std::size_t degree = 0U; degree < coeffs.size(); ++degree) {
            const Rational& coefficient = coeffs[degree];
            if (coefficient.numerator().is_zero()) {
                continue;
            }
            if (coefficient.denominator() != BigInt(1)) {
                return fail<MultivariatePolynomial>(make_error(
                    CASErrorKind::InternalError,
                    "polynomial_gcd_multivariate: non-integer coefficient after interpolation"));
            }

            auto factors = key_to_factors(key);
            if (degree > 0U) {
                factors.emplace_back(interpolation_var, static_cast<unsigned int>(degree));
            }

            terms.push_back(MultivariateTerm{
                .coefficient = coefficient.numerator(),
                .factors = std::move(factors),
            });
        }
    }

    return ok(MultivariatePolynomial(std::move(terms)));
}

} // namespace cas::algebra
