#pragma once

#include "cas/algebra.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include "polynomial_internal.hpp"
#include <optional>
#include <string>
#include <vector>

namespace cas::algebra {

// Utility functions
[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message, std::optional<std::string> hint = std::nullopt);
[[nodiscard]] ExprPtr make_rational_expr(AstArena& arena, const Rational& value);

template <typename T>
[[nodiscard]] Result<T> fail_unimplemented(std::string operation, std::string detail) {
    return fail<T>(make_error(
        CASErrorKind::Unimplemented,
        std::move(operation) + " non e' ancora implementata nel modulo algebra",
        std::move(detail)));
}

[[nodiscard]] ExprPtr make_integer(AstArena& arena, long long value);
[[nodiscard]] bool is_zero_expr(ExprPtr expr);
[[nodiscard]] bool is_one_expr(ExprPtr expr);
[[nodiscard]] bool contains_decimal_literal(ExprPtr expr);

[[nodiscard]] Result<ExprPtr> simplify_expr(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> clone_into_context(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> add_exprs(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> negate_expr(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> multiply_exprs(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> divide_exprs(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> pow_expr(ExprPtr base, std::size_t exponent, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> subtract_exprs(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx);

[[nodiscard]] bool same_generator_expr(
    ExprPtr expr,
    ExprPtr alpha_expr,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<std::vector<ExprPtr>> solve_degree_two_expr(
    const PolyExpr& poly,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<std::vector<ExprPtr>> solve_degree_three_expr(
    const PolyExpr& poly,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<std::vector<ExprPtr>> solve_degree_four_expr(
    const PolyExpr& poly,
    symbolic::CASContext& ctx);

[[nodiscard]] BigInt pow_bigint_nonnegative(BigInt base, unsigned int exponent);

// Polynomial helpers
[[nodiscard]] Result<BigInt> expr_to_integer_coefficient(ExprPtr expr);
[[nodiscard]] BigInt integer_content(const IntPoly& coefficients);
void divide_integer_coefficients_by_scalar(IntPoly& coefficients, const BigInt& scalar);
[[nodiscard]] PolyExpr integer_coefficients_to_poly(const IntPoly& coefficients, AstArena& arena);
[[nodiscard]] Result<ExprPtr> integer_coefficients_to_expr(const IntPoly& coefficients, const Symbol& var, symbolic::CASContext& ctx);

[[nodiscard]] Result<IntegerExponent> parse_integer_exponent(ExprPtr expr);

// Internal polynomial operations
[[nodiscard]] Result<RationalParts> split_num_den(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> expand_expr_impl(ExprPtr expr, symbolic::CASContext& ctx);

// Factorization helpers
void append_factor_with_multiplicity(std::vector<PolynomialFactor>& factors, ExprPtr factor, unsigned int multiplicity = 1U);
[[nodiscard]] IntPoly primitive_integer_poly_local(IntPoly coefficients);

// Internal multivariate helpers
[[nodiscard]] Result<BigInt> expr_to_integer_value_for_multivariate(ExprPtr value);
[[nodiscard]] Result<ExprPtr> build_multivariate_monomial_expr(const MultivariateTerm& term, symbolic::CASContext& ctx);
[[nodiscard]] Result<MultivariatePolynomial> parse_multivariate_polynomial(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> multivariate_to_expr(const MultivariatePolynomial& poly, symbolic::CASContext& ctx);

[[nodiscard]] Result<MultivariatePolynomial> gcd_heuristic(const MultivariatePolynomial& P, const MultivariatePolynomial& Q);
[[nodiscard]] Result<MultivariatePolynomial> gcd_modular(const MultivariatePolynomial& P, const MultivariatePolynomial& Q);
[[nodiscard]] Result<MultivariatePolynomial> gcd_probabilistic(const MultivariatePolynomial& P, const MultivariatePolynomial& Q);
[[nodiscard]] Result<MultivariatePolynomial> gcd_multivariate_eval_interp(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx);

// F3.1 — Brown's modular GCD (GCL §7.4) + Zippel sparse + EZ-GCD with cofactors.
struct GcdWithCofactors {
    MultivariatePolynomial gcd;
    MultivariatePolynomial cofactor_p;  // P / gcd, certified: gcd * cofactor_p == P
    MultivariatePolynomial cofactor_q;  // Q / gcd, certified: gcd * cofactor_q == Q
};
[[nodiscard]] Result<MultivariatePolynomial> gcd_brown(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx);
[[nodiscard]] Result<MultivariatePolynomial> gcd_zippel_sparse(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx);
// Honest name for the legacy Lagrange-over-Z evaluation/interpolation path.
// NOT Brown's modular GCD — it interpolates in Z directly (coefficient growth).
[[nodiscard]] Result<MultivariatePolynomial> gcd_eval_interp_z(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx);
// REAL Brown's modular multivariate GCD (GCL §7.4–7.5).
// Uses modular reduction + Fp-recursive eval/interp + multi-prime CRT lift.
// out_primes_used (optional) receives the prime list actually used (probe).
[[nodiscard]] Result<MultivariatePolynomial> gcd_brown_modular(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx,
    std::vector<BigInt>* out_primes_used);
[[nodiscard]] Result<MultivariatePolynomial> gcd_brown_modular(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx);
// REAL Zippel sparse interpolation GCD (Zippel 1979) — Prony skeleton + Vandermonde.
// out_samples_used (optional) receives the number of evaluation calls (probe).
[[nodiscard]] Result<MultivariatePolynomial> gcd_zippel_prony(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx,
    std::size_t* out_samples_used);
[[nodiscard]] Result<MultivariatePolynomial> gcd_zippel_prony(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx);
[[nodiscard]] Result<GcdWithCofactors> gcd_ez(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx);

// Factorization structures
struct IntegerSquareFreeFactor {
    IntPoly factor;
    unsigned int multiplicity{1U};
};

struct MultivariateSquareFreeFactor {
    MultivariatePolynomial factor;
    unsigned int multiplicity{1U};
};

[[nodiscard]] Result<std::vector<MultivariateSquareFreeFactor>> square_free_factorize_multivariate(
    const MultivariatePolynomial& poly, symbolic::CASContext& ctx);

[[nodiscard]] Result<std::vector<IntegerSquareFreeFactor>> square_free_factorize_integer_poly(const IntPoly& primitive, symbolic::CASContext& ctx);
[[nodiscard]] Result<void> append_integer_factor_component(Factorization& factorization, const IntPoly& component, unsigned int multiplicity, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] BigInt select_factorization_prime(const IntPoly& f);

// F7.5.A4 — Hyperbolic reciprocal/quotient normalisation.
// Rewrites sech/csch/coth/tanh FuncCall nodes to canonical cosh/sinh
// quotient form. Preserves structural sharing when no rewrite fires.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Sech_Csch_Identity.md.
[[nodiscard]] ExprPtr hyperbolic_normalize(ExprPtr expr, AstArena& arena);

// A44: rewrites u! / factorial(u) to gamma(u+1) (exact identity) so the two
// interchangeable spellings become structurally comparable. Applied to both
// sides of mathematically_equal, like hyperbolic_normalize above.
[[nodiscard]] ExprPtr factorial_gamma_normalize(ExprPtr expr, AstArena& arena);

// A43 §4: riduce li/Shi/Chi a Ei ed erfi a erf (identità esatte), così che le
// ortografie interscambiabili della famiglia non elementare diventino
// confrontabili strutturalmente. Applicato ai due lati di mathematically_equal,
// come i due normalizzatori sopra. Si/Ci NON vengono toccati (le loro identità
// verso Ei introdurrebbero un i spurio in un risultato reale — spec §4).
[[nodiscard]] ExprPtr nonelementary_normalize(ExprPtr expr, AstArena& arena);

// F7.5.A1 — Geometric / cyclotomic RootOf expansion (closes
// HC-F75-CYCLOTOMIC-ROOTOF). See src/algebra/algebraic_equal_cyclotomic.cpp.
[[nodiscard]] std::optional<std::vector<ExprPtr>>
enumerate_geometric_rootof(const RootOf& node, symbolic::CASContext& ctx);

// T-025 — genuine cyclotomic Φ_n (composite n) RootOf expansion into its φ(n)
// primitive roots exp(2πi·m/n), gcd(m,n)=1. Complements the geometric enumerator
// (which only covers prime-n Φ_p). See src/algebra/algebraic_equal_cyclotomic.cpp.
[[nodiscard]] std::optional<std::vector<ExprPtr>>
enumerate_cyclotomic_rootof(const RootOf& node, symbolic::CASContext& ctx);

// Returns Some(bool) when the RootOf-specific dispatch (distinct-index
// guard + geometric expansion) can decide the equality outright;
// returns nullopt to fall through to the general mathematically_equal
// pipeline.
[[nodiscard]] std::optional<bool> try_rootof_decision(
    ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx);

// A54 — sospende la raccolta simbolica dei like-term (F1.4) per la durata di
// una costruzione la cui post-condizione è la forma ESPANSA.  F1.4 è l'inversa
// della distribuzione: senza questa sospensione il `simplify` che costruisce i
// Sum di `expand` ri-fattorizza ciò che expand ha appena distribuito, e la
// post-condizione viene violata dalla FORMA dei fattori (basta una base
// condivisa non-Symbol).  Rientrante: ripristina il valore precedente, quindi
// chiamate annidate e chiamanti che l'hanno già sospesa non si disturbano.
class ExpandedFormScope {
public:
    explicit ExpandedFormScope(symbolic::CASContext& ctx)
        : ctx_(ctx), previous_(ctx.symbolic_like_term_factoring()) {
        ctx_.set_symbolic_like_term_factoring(false);
    }
    ~ExpandedFormScope() { ctx_.set_symbolic_like_term_factoring(previous_); }
    ExpandedFormScope(const ExpandedFormScope&) = delete;
    ExpandedFormScope& operator=(const ExpandedFormScope&) = delete;
    ExpandedFormScope(ExpandedFormScope&&) = delete;
    ExpandedFormScope& operator=(ExpandedFormScope&&) = delete;

private:
    symbolic::CASContext& ctx_;
    bool previous_;
};

// A55 — sospende la distribuzione simbolica di un fattore su un `Sum` (Step 8
// di `simplify_product_factors`, `k·(a+b) → k·a+k·b`) per la durata di una
// costruzione la cui post-condizione è la forma COMBINATA `N/D`.  È l'esatto
// complementare di `ExpandedFormScope`: lì si sospendeva la raccolta (inversa
// della distribuzione) perché la post-condizione era la forma espansa; qui si
// sospende la distribuzione (inversa della combinazione) perché la
// post-condizione è una frazione unica.  Senza questa sospensione,
// `divide_exprs(N, D)` in `algebra::together` costruisce `Product(Sum(N),
// Pow(D,-1))`, che Step 8 ridistribuisce immediatamente in una somma di
// frazioni — `together` non "mette insieme" nulla.  Rientrante: ripristina il
// valore precedente.
class CombinedFormScope {
public:
    explicit CombinedFormScope(symbolic::CASContext& ctx)
        : ctx_(ctx), previous_(ctx.symbolic_sum_distribution()) {
        ctx_.set_symbolic_sum_distribution(false);
    }
    ~CombinedFormScope() { ctx_.set_symbolic_sum_distribution(previous_); }
    CombinedFormScope(const CombinedFormScope&) = delete;
    CombinedFormScope& operator=(const CombinedFormScope&) = delete;
    CombinedFormScope(CombinedFormScope&&) = delete;
    CombinedFormScope& operator=(CombinedFormScope&&) = delete;

private:
    symbolic::CASContext& ctx_;
    bool previous_;
};

} // namespace cas::algebra
