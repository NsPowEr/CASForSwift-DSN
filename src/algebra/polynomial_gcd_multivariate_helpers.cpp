// polynomial_gcd_multivariate_helpers.cpp
// Implementations of FactorKey utilities, polynomial normalization helpers,
// conversion between MultivariatePolynomial and IntPoly, and variable collection.
// Declarations live in polynomial_gcd_multivariate_helpers.hpp.

#include "polynomial_gcd_multivariate_helpers.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <algorithm>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace cas::algebra {

[[nodiscard]] bool FactorKeyLess::operator()(const FactorKey& lhs, const FactorKey& rhs) const noexcept {
    return lhs < rhs;
}

[[nodiscard]] FactorKey make_factor_key(const std::vector<std::pair<Symbol, unsigned int>>& factors,
                                        const std::optional<std::string>& omit_var) {
    FactorKey key;
    key.reserve(factors.size());
    for (const auto& [symbol, exponent] : factors) {
        if (exponent == 0U) {
            continue;
        }
        if (omit_var.has_value() && symbol.name == omit_var.value()) {
            continue;
        }
        key.emplace_back(symbol.name, exponent);
    }

    std::sort(key.begin(), key.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });

    FactorKey merged;
    merged.reserve(key.size());
    for (const auto& [name, exponent] : key) {
        if (!merged.empty() && merged.back().first == name) {
            merged.back().second += exponent;
        } else {
            merged.emplace_back(name, exponent);
        }
    }
    return merged;
}

[[nodiscard]] std::vector<std::pair<Symbol, unsigned int>> key_to_factors(const FactorKey& key) {
    std::vector<std::pair<Symbol, unsigned int>> factors;
    factors.reserve(key.size());
    for (const auto& [name, exponent] : key) {
        if (exponent > 0U) {
            factors.emplace_back(Symbol(name), exponent);
        }
    }
    return factors;
}

[[nodiscard]] std::vector<Symbol> collect_all_variables(const MultivariatePolynomial& p,
                                                        const MultivariatePolynomial& q) {
    auto vars_p = p.variables();
    auto vars_q = q.variables();

    std::vector<Symbol> all_vars = vars_p;
    for (const auto& symbol : vars_q) {
        if (std::find_if(all_vars.begin(), all_vars.end(),
                         [&](const Symbol& candidate) { return candidate.name == symbol.name; }) == all_vars.end()) {
            all_vars.push_back(symbol);
        }
    }

    std::sort(all_vars.begin(), all_vars.end(), [](const Symbol& lhs, const Symbol& rhs) {
        return lhs.name < rhs.name;
    });
    return all_vars;
}

[[nodiscard]] std::size_t degree_in_var(const MultivariatePolynomial& poly, const Symbol& var) {
    std::size_t degree = 0U;
    for (const auto& term : poly.terms()) {
        for (const auto& [symbol, exponent] : term.factors) {
            if (symbol.name == var.name) {
                degree = std::max(degree, static_cast<std::size_t>(exponent));
            }
        }
    }
    return degree;
}

[[nodiscard]] CoeffMap to_coeff_map(const MultivariatePolynomial& poly,
                                    const std::optional<std::string>& omit_var) {
    CoeffMap coefficients;
    for (const auto& term : poly.terms()) {
        coefficients[make_factor_key(term.factors, omit_var)] += term.coefficient;
    }

    for (auto it = coefficients.begin(); it != coefficients.end();) {
        if (it->second.is_zero()) {
            it = coefficients.erase(it);
        } else {
            ++it;
        }
    }
    return coefficients;
}

[[nodiscard]] bool same_polynomial(const MultivariatePolynomial& lhs, const MultivariatePolynomial& rhs) {
    return to_coeff_map(lhs) == to_coeff_map(rhs);
}

[[nodiscard]] MultivariatePolynomial multiply_by_scalar(const MultivariatePolynomial& poly, const BigInt& scalar) {
    if (poly.is_zero() || scalar.is_zero()) {
        return MultivariatePolynomial{};
    }

    std::vector<MultivariateTerm> terms;
    terms.reserve(poly.terms().size());
    for (const auto& term : poly.terms()) {
        terms.push_back(MultivariateTerm{
            .coefficient = term.coefficient * scalar,
            .factors = term.factors,
        });
    }
    return MultivariatePolynomial(std::move(terms));
}

[[nodiscard]] MultivariatePolynomial normalize_multivariate_gcd(const MultivariatePolynomial& poly) {
    if (poly.is_zero()) {
        return poly;
    }

    BigInt content(0);
    for (const auto& term : poly.terms()) {
        content = gcd(content, term.coefficient.abs());
    }

    std::vector<MultivariateTerm> normalized;
    normalized.reserve(poly.terms().size());
    for (const auto& term : poly.terms()) {
        MultivariateTerm reduced = term;
        if (content > BigInt(1)) {
            reduced.coefficient /= content;
        }
        normalized.push_back(std::move(reduced));
    }

    MultivariatePolynomial result(std::move(normalized));
    CoeffMap coeffs = to_coeff_map(result);
    if (!coeffs.empty() && std::prev(coeffs.end())->second.is_negative()) {
        result = multiply_by_scalar(result, BigInt(-1));
    }

    return result;
}

[[nodiscard]] bool is_unit_polynomial(const MultivariatePolynomial& poly) {
    if (poly.terms().size() != 1U) {
        return false;
    }
    const auto& term = poly.terms().front();
    return term.factors.empty() && term.coefficient.abs() == BigInt(1);
}

[[nodiscard]] Result<IntPoly> multivariate_single_var_to_intpoly(const MultivariatePolynomial& poly,
                                                                   const Symbol& var) {
    if (poly.is_zero()) {
        return ok(IntPoly{});
    }

    std::map<std::size_t, BigInt> coefficient_map;
    std::size_t max_degree = 0U;

    for (const auto& term : poly.terms()) {
        std::size_t degree = 0U;
        for (const auto& [symbol, exponent] : term.factors) {
            if (symbol.name == var.name) {
                degree = exponent;
            } else {
                return fail<IntPoly>(make_error(
                    CASErrorKind::Unimplemented,
                    "polynomial_gcd_multivariate: expected univariate polynomial after specialization"));
            }
        }
        coefficient_map[degree] += term.coefficient;
        max_degree = std::max(max_degree, degree);
    }

    IntPoly result;
    result.resize(max_degree + 1U, BigInt(0));
    for (const auto& [degree, coefficient] : coefficient_map) {
        result[degree] = coefficient;
    }
    result.normalize([](const BigInt& value) { return value.is_zero(); });
    return ok(std::move(result));
}

[[nodiscard]] MultivariatePolynomial intpoly_to_multivariate(const IntPoly& poly, const Symbol& var) {
    std::vector<MultivariateTerm> terms;
    terms.reserve(poly.size());
    for (std::size_t degree = 0; degree < poly.size(); ++degree) {
        if (poly[degree].is_zero()) {
            continue;
        }
        std::vector<std::pair<Symbol, unsigned int>> factors;
        if (degree > 0U) {
            factors.emplace_back(var, static_cast<unsigned int>(degree));
        }
        terms.push_back(MultivariateTerm{
            .coefficient = poly[degree],
            .factors = std::move(factors),
        });
    }
    return MultivariatePolynomial(std::move(terms));
}

} // namespace cas::algebra
