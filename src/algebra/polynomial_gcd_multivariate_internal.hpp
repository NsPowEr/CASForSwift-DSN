#pragma once

// Internal types shared across polynomial_gcd_multivariate_*.cpp split files.
// NOT part of the public API — do NOT include from include/cas/.

#include "polynomial_gcd_multivariate_helpers.hpp"

#include <map>
#include <optional>
#include <vector>

namespace cas::algebra {

// Local Monomial/SparsePoly types used by sparse arithmetic and exact_quotient.
using Monomial  = std::vector<unsigned int>;
using SparsePoly = std::map<Monomial, BigInt>;

// ---------------------------------------------------------------------------
// Sparse arithmetic helpers (polynomial_gcd_multivariate_sparse.cpp)
// ---------------------------------------------------------------------------
[[nodiscard]] SparsePoly to_sparse(const MultivariatePolynomial& poly,
                                    const std::vector<Symbol>& vars);

[[nodiscard]] bool monomial_divides(const Monomial& divisor, const Monomial& dividend);

[[nodiscard]] Monomial monomial_quotient(const Monomial& dividend, const Monomial& divisor);

void sparse_add_term(SparsePoly& poly, const Monomial& monomial, const BigInt& coefficient);

[[nodiscard]] SparsePoly multiply_sparse_by_term(const SparsePoly& poly,
                                                   const Monomial& monomial,
                                                   const BigInt& coefficient);

void sparse_subtract(SparsePoly& lhs, const SparsePoly& rhs);

[[nodiscard]] MultivariatePolynomial sparse_to_multivariate(const SparsePoly& sparse,
                                                              const std::vector<Symbol>& vars);

[[nodiscard]] Result<std::optional<MultivariatePolynomial>> exact_quotient(
    const MultivariatePolynomial& dividend,
    const MultivariatePolynomial& divisor,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx);

// ---------------------------------------------------------------------------
// Linear-candidate helpers (polynomial_gcd_multivariate_linear.cpp)
// ---------------------------------------------------------------------------
[[nodiscard]] std::vector<MultivariatePolynomial> primitive_linear_candidates(
    const std::vector<Symbol>& vars);

// ---------------------------------------------------------------------------
// Core recursive GCD (polynomial_gcd_multivariate.cpp)
// ---------------------------------------------------------------------------
[[nodiscard]] Result<MultivariatePolynomial> gcd_multivariate_recursive(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx,
    std::size_t depth,
    std::size_t& call_count);

} // namespace cas::algebra
