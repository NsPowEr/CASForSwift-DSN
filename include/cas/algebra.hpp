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
    
    [[nodiscard]] MultivariatePolynomial derivative(const Symbol& var) const;
    [[nodiscard]] BigInt integer_content() const;

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

// F4.2d — Estrai coefficienti di un polinomio univariato in `var`.
// Restituisce vettore [c_0, c_1, ..., c_n] tale che expr = Σ c_k · var^k.
// Coefficienti rimanenti possono dipendere da altri simboli. Fallisce se
// `expr` non riducibile a forma polinomiale in `var` (es. termini con
// var al denominatore o esponenti non interi non-negativi).
[[nodiscard]] Result<std::vector<ExprPtr>> univariate_coefficients(
    ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx);

// F4.2b — Bezout identity per polinomi su Q[x]: dati a, b restituisce
// (g, s, t) tali che s·a + t·b = g = gcd(a, b).  Usato da Smith Q[x].
struct PolynomialBezout {
    ExprPtr gcd;
    ExprPtr s;
    ExprPtr t;
};
[[nodiscard]] Result<PolynomialBezout> polynomial_bezout(
    ExprPtr a, ExprPtr b, const Symbol& var, symbolic::CASContext& ctx);

// F4.2b — divisione esatta a / b in Q[x]. Fallisce se b non divide a esattamente.
[[nodiscard]] Result<ExprPtr> polynomial_exact_divide(
    ExprPtr a, ExprPtr b, const Symbol& var, symbolic::CASContext& ctx);

/**
 * @brief Performs polynomial division with remainder: a = q*b + r.
 * @return struct with quotient q and remainder r.
 */
struct PolynomialDivMod {
    ExprPtr quotient;
    ExprPtr remainder;
};
[[nodiscard]] Result<PolynomialDivMod> polynomial_divmod(
    ExprPtr a, ExprPtr b, const Symbol& var, symbolic::CASContext& ctx);

// F4.2b — grado del polinomio in var (su Q[x], coefficienti razionali ammessi).
[[nodiscard]] Result<std::size_t> polynomial_degree(
    ExprPtr a, const Symbol& var, symbolic::CASContext& ctx);
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
// F3.2 — Wang multivariate factorization over Z (EEZ / Extended Zassenhaus).
// Factors a multivariate polynomial into irreducible factors with multiplicities.
// content holds the integer content; factors hold the irreducible parts.
[[nodiscard]] Result<Factorization> factor_multivariate(
    ExprPtr poly,
    symbolic::CASContext& ctx);
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

// L2-06: fsolve — transcendental equation solver.
// Accepts f(x)=0 or f(x)=g(x) form with one variable.
// Pipeline: (1) try symbolic solve_polynomial; (2) fallback to numeric
// multi-root scan on [low, high] via bisection + Newton polish.
// Returns a column Matrix of root values (exact rational or DecimalLit).
[[nodiscard]] Result<ExprPtr> fsolve(
    ExprPtr equation,
    const Symbol& var,
    symbolic::CASContext& ctx,
    double search_low  = -10.0,
    double search_high =  10.0);
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
