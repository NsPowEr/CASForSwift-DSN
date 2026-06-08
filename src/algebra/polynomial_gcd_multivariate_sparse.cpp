// polynomial_gcd_multivariate_sparse.cpp
// SparsePoly arithmetic (to_sparse, monomial ops, exact_quotient) used by the
// multivariate GCD certification step.
// Declarations live in polynomial_gcd_multivariate_internal.hpp.

#include "polynomial_gcd_multivariate_internal.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <map>
#include <optional>
#include <vector>

namespace cas::algebra {

[[nodiscard]] SparsePoly to_sparse(const MultivariatePolynomial& poly, const std::vector<Symbol>& vars) {
    SparsePoly sparse;
    std::map<std::string, std::size_t> indices;
    for (std::size_t i = 0U; i < vars.size(); ++i) {
        indices[vars[i].name] = i;
    }

    for (const auto& term : poly.terms()) {
        Monomial monomial(vars.size(), 0U);
        for (const auto& [symbol, exponent] : term.factors) {
            const auto it = indices.find(symbol.name);
            if (it != indices.end()) {
                monomial[it->second] += exponent;
            }
        }
        sparse[monomial] += term.coefficient;
    }

    for (auto it = sparse.begin(); it != sparse.end();) {
        if (it->second.is_zero()) {
            it = sparse.erase(it);
        } else {
            ++it;
        }
    }

    return sparse;
}

[[nodiscard]] bool monomial_divides(const Monomial& divisor, const Monomial& dividend) {
    for (std::size_t i = 0U; i < divisor.size(); ++i) {
        if (divisor[i] > dividend[i]) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] Monomial monomial_quotient(const Monomial& dividend, const Monomial& divisor) {
    Monomial quotient(dividend.size(), 0U);
    for (std::size_t i = 0U; i < dividend.size(); ++i) {
        quotient[i] = dividend[i] - divisor[i];
    }
    return quotient;
}

void sparse_add_term(SparsePoly& poly, const Monomial& monomial, const BigInt& coefficient) {
    if (coefficient.is_zero()) {
        return;
    }

    poly[monomial] += coefficient;
    if (poly[monomial].is_zero()) {
        poly.erase(monomial);
    }
}

[[nodiscard]] SparsePoly multiply_sparse_by_term(
    const SparsePoly& poly,
    const Monomial& monomial,
    const BigInt& coefficient) {
    SparsePoly result;
    if (coefficient.is_zero()) {
        return result;
    }

    for (const auto& [base_monomial, base_coefficient] : poly) {
        Monomial product(base_monomial.size(), 0U);
        for (std::size_t i = 0U; i < base_monomial.size(); ++i) {
            product[i] = base_monomial[i] + monomial[i];
        }
        sparse_add_term(result, product, base_coefficient * coefficient);
    }
    return result;
}

void sparse_subtract(SparsePoly& lhs, const SparsePoly& rhs) {
    for (const auto& [monomial, coefficient] : rhs) {
        sparse_add_term(lhs, monomial, -coefficient);
    }
}

[[nodiscard]] MultivariatePolynomial sparse_to_multivariate(const SparsePoly& sparse,
                                                             const std::vector<Symbol>& vars) {
    std::vector<MultivariateTerm> terms;
    terms.reserve(sparse.size());

    for (const auto& [monomial, coefficient] : sparse) {
        if (coefficient.is_zero()) {
            continue;
        }

        std::vector<std::pair<Symbol, unsigned int>> factors;
        factors.reserve(vars.size());
        for (std::size_t i = 0U; i < vars.size(); ++i) {
            if (monomial[i] > 0U) {
                factors.emplace_back(vars[i], monomial[i]);
            }
        }

        terms.push_back(MultivariateTerm{
            .coefficient = coefficient,
            .factors = std::move(factors),
        });
    }

    return MultivariatePolynomial(std::move(terms));
}

[[nodiscard]] Result<std::optional<MultivariatePolynomial>> exact_quotient(
    const MultivariatePolynomial& dividend,
    const MultivariatePolynomial& divisor,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx) {
    if (divisor.is_zero()) {
        return fail<std::optional<MultivariatePolynomial>>(make_error(
            CASErrorKind::InvalidArgument,
            "polynomial_gcd_multivariate: exact division by zero divisor"));
    }

    SparsePoly remainder = to_sparse(dividend, vars);
    SparsePoly divisor_sparse = to_sparse(divisor, vars);
    SparsePoly quotient;

    if (divisor_sparse.empty()) {
        return fail<std::optional<MultivariatePolynomial>>(make_error(
            CASErrorKind::InvalidArgument,
            "polynomial_gcd_multivariate: zero sparse divisor"));
    }

    const auto [divisor_lm, divisor_lc] = *std::prev(divisor_sparse.end());
    // HPP-003 FIX 2026-05-28 (CLAUDE.md Cat 2 — no magic constant).
    // Bound derivation: multivariate exact division by repeated leading-monomial
    // cancellation.  Each step eliminates the leading monomial of remainder
    // and subtracts at most (|divisor_sparse|-1) new monomials.  The worst-case
    // total monomial-pair interactions equals |remainder| * |divisor_sparse|.
    // Formula: (remainder.size() + 1) * (divisor_sparse.size() + 1).
    // This replaces the previous (remainder.size()+1)*(vars.size()+1)*16U which
    // used the variable count (irrelevant to step count) and the magic 16U.
    const std::size_t max_steps =
        std::max(ctx.min_gcd_division_steps(),
                 (remainder.size() + 1U) * (divisor_sparse.size() + 1U));
    std::size_t steps = 0U;

    while (!remainder.empty()) {
        if (++steps > max_steps) {
            return fail<std::optional<MultivariatePolynomial>>(make_error(
                CASErrorKind::Unimplemented,
                "polynomial_gcd_multivariate: exact division budget exceeded during certification"));
        }

        const auto [remainder_lm, remainder_lc] = *std::prev(remainder.end());

        if (!monomial_divides(divisor_lm, remainder_lm)) {
            return ok(std::optional<MultivariatePolynomial>{std::nullopt});
        }
        if ((remainder_lc % divisor_lc) != BigInt(0)) {
            return ok(std::optional<MultivariatePolynomial>{std::nullopt});
        }

        const BigInt term_coefficient = remainder_lc / divisor_lc;
        const Monomial term_monomial = monomial_quotient(remainder_lm, divisor_lm);

        sparse_add_term(quotient, term_monomial, term_coefficient);

        SparsePoly subtraction = multiply_sparse_by_term(divisor_sparse, term_monomial, term_coefficient);
        sparse_subtract(remainder, subtraction);

        if (!remainder.empty() && !(std::prev(remainder.end())->first < remainder_lm)) {
            return fail<std::optional<MultivariatePolynomial>>(make_error(
                CASErrorKind::Unimplemented,
                "polynomial_gcd_multivariate: exact division did not reduce leading monomial"));
        }
    }

    return ok(std::optional<MultivariatePolynomial>{sparse_to_multivariate(quotient, vars)});
}

} // namespace cas::algebra
