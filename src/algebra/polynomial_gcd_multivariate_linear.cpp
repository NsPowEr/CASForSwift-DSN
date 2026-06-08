// polynomial_gcd_multivariate_linear.cpp
// Primitive linear GCD candidate generation: generates all primitive linear
// polynomials in the given variable set that might divide both inputs, used
// as a fast-path before the full eval/interpolation algorithm.
// Declaration lives in polynomial_gcd_multivariate_internal.hpp.

#include "polynomial_gcd_multivariate_internal.hpp"
#include "polynomial_gcd_multivariate_helpers.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <algorithm>
#include <set>
#include <vector>

namespace cas::algebra {

namespace {

[[nodiscard]] BigInt gcd_abs_values(const std::vector<int>& values) {
    BigInt result(0);
    for (int value : values) {
        result = gcd(result, BigInt(value).abs());
    }
    return result;
}

[[nodiscard]] std::optional<MultivariatePolynomial> make_primitive_linear_candidate(
    const std::vector<Symbol>& vars,
    const std::vector<int>& coefficients,
    int constant) {
    bool has_variable_term = false;
    std::vector<int> all_values = coefficients;
    all_values.push_back(constant);
    const BigInt content = gcd_abs_values(all_values);
    if (content != BigInt(1)) {
        return std::nullopt;
    }

    int sign = 1;
    for (int coefficient : coefficients) {
        if (coefficient != 0) {
            sign = coefficient < 0 ? -1 : 1;
            break;
        }
    }

    std::vector<MultivariateTerm> terms;
    for (std::size_t i = 0U; i < vars.size(); ++i) {
        int coefficient = coefficients[i] * sign;
        if (coefficient == 0) {
            continue;
        }
        has_variable_term = true;
        terms.push_back(MultivariateTerm{
            .coefficient = BigInt(coefficient),
            .factors = {{vars[i], 1U}},
        });
    }

    const int normalized_constant = constant * sign;
    if (normalized_constant != 0) {
        terms.push_back(MultivariateTerm{
            .coefficient = BigInt(normalized_constant),
            .factors = {},
        });
    }

    if (!has_variable_term) {
        return std::nullopt;
    }

    return MultivariatePolynomial(std::move(terms));
}

} // namespace

[[nodiscard]] std::vector<MultivariatePolynomial> primitive_linear_candidates(
    const std::vector<Symbol>& vars) {
    if (vars.empty() || vars.size() > 4U) {
        return {};
    }

    std::vector<MultivariatePolynomial> candidates;
    std::set<CoeffMap> seen;
    std::vector<int> coefficients(vars.size(), 0);

    auto emit = [&](int constant) {
        auto candidate = make_primitive_linear_candidate(vars, coefficients, constant);
        if (!candidate.has_value()) {
            return;
        }
        CoeffMap key = to_coeff_map(candidate.value());
        if (seen.insert(key).second) {
            candidates.push_back(std::move(candidate.value()));
        }
    };

    auto walk = [&](auto&& self, std::size_t index) -> void {
        if (index == vars.size()) {
            for (int constant = -3; constant <= 3; ++constant) {
                emit(constant);
            }
            return;
        }
        for (int coefficient = -2; coefficient <= 2; ++coefficient) {
            coefficients[index] = coefficient;
            self(self, index + 1U);
        }
    };

    walk(walk, 0U);
    std::sort(candidates.begin(), candidates.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.terms().size() != rhs.terms().size()) {
            return lhs.terms().size() < rhs.terms().size();
        }
        return to_coeff_map(lhs) < to_coeff_map(rhs);
    });
    return candidates;
}

} // namespace cas::algebra
