#pragma once

#include "cas/algebra.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "polynomial_internal.hpp"

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace cas::algebra {

using FactorKey = std::vector<std::pair<std::string, unsigned int>>;

struct FactorKeyLess {
    [[nodiscard]] bool operator()(const FactorKey& lhs, const FactorKey& rhs) const noexcept;
};

using CoeffMap = std::map<FactorKey, BigInt, FactorKeyLess>;

[[nodiscard]] std::vector<Symbol> collect_all_variables(const MultivariatePolynomial& p,
                                                        const MultivariatePolynomial& q);
[[nodiscard]] std::size_t degree_in_var(const MultivariatePolynomial& poly, const Symbol& var);
[[nodiscard]] CoeffMap to_coeff_map(const MultivariatePolynomial& poly,
                                    const std::optional<std::string>& omit_var = std::nullopt);
[[nodiscard]] bool same_polynomial(const MultivariatePolynomial& lhs, const MultivariatePolynomial& rhs);

[[nodiscard]] MultivariatePolynomial normalize_multivariate_gcd(const MultivariatePolynomial& poly);
[[nodiscard]] bool is_unit_polynomial(const MultivariatePolynomial& poly);

[[nodiscard]] Result<IntPoly> multivariate_single_var_to_intpoly(const MultivariatePolynomial& poly, const Symbol& var);
[[nodiscard]] MultivariatePolynomial intpoly_to_multivariate(const IntPoly& poly, const Symbol& var);

[[nodiscard]] Result<std::vector<Rational>> lagrange_interpolate(
    const std::vector<BigInt>& values,
    const std::vector<BigInt>& points);

[[nodiscard]] MultivariatePolynomial coefficient_poly_in_var(
    const MultivariatePolynomial& poly,
    const Symbol& var,
    std::size_t degree);

[[nodiscard]] MultivariatePolynomial multiply_by_variable_power(
    const MultivariatePolynomial& poly,
    const Symbol& var,
    std::size_t power);

[[nodiscard]] Result<MultivariatePolynomial> interpolate_polynomial_values(
    const std::vector<MultivariatePolynomial>& values,
    const std::vector<BigInt>& points,
    const Symbol& interpolation_var);

[[nodiscard]] Result<std::optional<MultivariatePolynomial>> exact_quotient(
    const MultivariatePolynomial& dividend,
    const MultivariatePolynomial& divisor,
    const std::vector<Symbol>& vars);

} // namespace cas::algebra
