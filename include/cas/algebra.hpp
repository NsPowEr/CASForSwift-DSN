#pragma once

#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <vector>
#include <utility>

namespace cas {
namespace algebra {

struct PolynomialFactor {
    ExprPtr factor;
    unsigned int multiplicity{1U};
};

struct Factorization {
    ExprPtr content;
    std::vector<PolynomialFactor> factors;
};

struct RationalParts {
    ExprPtr numerator;
    ExprPtr denominator;
};

struct MultivariateTerm {
    BigInt coefficient;
    std::vector<std::pair<Symbol, unsigned int>> factors;
};

class MultivariatePolynomial {
public:
    MultivariatePolynomial() = default;
    explicit MultivariatePolynomial(std::vector<MultivariateTerm> terms);

    [[nodiscard]] bool is_zero() const noexcept;
    [[nodiscard]] std::size_t total_degree() const noexcept;
    [[nodiscard]] std::vector<Symbol> variables() const;
    [[nodiscard]] const std::vector<MultivariateTerm>& terms() const noexcept;

    [[nodiscard]] MultivariatePolynomial operator+(const MultivariatePolynomial& other) const;
    [[nodiscard]] MultivariatePolynomial operator*(const MultivariatePolynomial& other) const;

    [[nodiscard]] Result<MultivariatePolynomial> evaluate_at(
        const Symbol& var,
        ExprPtr value) const;

    [[nodiscard]] Result<ExprPtr> evaluate_at_rational(
        const Symbol& var,
        const Rational& value,
        AstArena& arena) const;

    [[nodiscard]] Result<std::vector<ExprPtr>> to_univariate_coefficients(
        const Symbol& main_var,
        symbolic::CASContext& ctx) const;

private:
    std::vector<MultivariateTerm> terms_;
};

struct SquareFreeFactor {
    ExprPtr factor;
    unsigned int multiplicity{1U};
};

struct SquareFreeFactorization {
    ExprPtr content;
    std::vector<SquareFreeFactor> factors;
};

// Boundary pubblico F3: API stabile, implementazione incrementale dietro Result<T>.
[[nodiscard]] Result<ExprPtr> expand(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> collect(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> together(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<RationalParts> apart_num_den(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> polynomial_gcd(ExprPtr p, ExprPtr q, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> polynomial_gcd_multivariate(ExprPtr p, ExprPtr q, symbolic::CASContext& ctx);
[[nodiscard]] Result<Factorization> factor_over_integers(ExprPtr poly, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] Result<Factorization> factor_polynomial(
    ExprPtr poly,
    const Symbol& var,
    symbolic::CASContext& ctx,
    std::optional<ExprPtr> extension = std::nullopt);
[[nodiscard]] Result<std::vector<ExprPtr>> partial_fractions(
    ExprPtr rational_expr,
    const Symbol& var,
    symbolic::CASContext& ctx);
[[nodiscard]] Result<std::vector<ExprPtr>> solve_polynomial(
    ExprPtr poly,
    const Symbol& var,
    symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> csolve(
    ExprPtr eqs,
    ExprPtr vars,
    symbolic::CASContext& ctx);
[[nodiscard]] Result<SquareFreeFactorization> square_free_factorization(
    ExprPtr poly,
    const Symbol& var,
    symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> polynomial_resultant(
    ExprPtr p,
    ExprPtr q,
    const Symbol& var,
    symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> polynomial_discriminant(
    ExprPtr p,
    const Symbol& var,
    symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> integrate_rational_lrt(
    ExprPtr P,
    ExprPtr Q,
    const Symbol& var,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<std::vector<ExprPtr>> polynomial_groebner(
    const std::vector<ExprPtr>& equations,
    const std::vector<Symbol>& variables,
    symbolic::CASContext& ctx);

}  // namespace algebra
}  // namespace cas
