#pragma once

// F3.2 — Wang multivariate factorization (EEZ / Extended Zassenhaus).
// Reference: Wang 1978 "An improved multivariate polynomial factoring algorithm";
// Geddes, Czapor, Labahn "Algorithms for Computer Algebra" (GCL) §6.4–6.8.
//
// Internal sparse multivariate polynomial over Z (MPoly): a map from a dense
// exponent vector (indexed by a fixed variable ordering) to a BigInt coefficient.
// All arithmetic is exact (BigInt / Rational). No floating point anywhere.

#include "cas/algebra.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <map>
#include <optional>
#include <vector>

namespace cas::algebra {

// Exponent vector over a fixed variable ordering (index i ↔ vars_[i]).
using Monomial = std::vector<unsigned int>;

// Sparse multivariate polynomial over Z.  Variable ordering is shared by the
// owning WangContext; index 0 is the MAIN variable x_1.
struct MPoly {
    std::map<Monomial, BigInt> terms;  // canonical: no zero coefficients

    [[nodiscard]] bool is_zero() const { return terms.empty(); }
    void prune();  // drop zero-coefficient terms
};

// Shared state for one factorization: the variable ordering (main var first).
struct WangContext {
    std::vector<Symbol> vars;  // vars[0] = main variable x_1
    [[nodiscard]] std::size_t nvars() const { return vars.size(); }
};

// ---- MPoly construction / conversion ----
[[nodiscard]] Result<MPoly> mpoly_from_expr(
    ExprPtr expr, const WangContext& wc, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> mpoly_to_expr(
    const MPoly& p, const WangContext& wc, symbolic::CASContext& ctx);
[[nodiscard]] MPoly mpoly_zero();
[[nodiscard]] MPoly mpoly_constant(const BigInt& c, std::size_t nvars);

// ---- MPoly arithmetic (exact) ----
[[nodiscard]] MPoly mpoly_add(const MPoly& a, const MPoly& b);
[[nodiscard]] MPoly mpoly_sub(const MPoly& a, const MPoly& b);
[[nodiscard]] MPoly mpoly_mul(const MPoly& a, const MPoly& b);
[[nodiscard]] MPoly mpoly_scale(const MPoly& a, const BigInt& s);
// Exact division a / b over Z.  nullopt if not exact (verified, never silent-wrong).
[[nodiscard]] std::optional<MPoly> mpoly_exact_div(const MPoly& a, const MPoly& b);
[[nodiscard]] bool mpoly_equal(const MPoly& a, const MPoly& b);

// ---- Degrees / coefficients ----
[[nodiscard]] unsigned int mpoly_degree_in(const MPoly& p, std::size_t var_index);
[[nodiscard]] BigInt mpoly_integer_content(const MPoly& p);
// Coefficient of var_index^d as an MPoly in the remaining variables (same var space).
[[nodiscard]] MPoly mpoly_coeff_in(const MPoly& p, std::size_t var_index, unsigned int deg);
// Leading coefficient of p w.r.t. var_index, as MPoly in remaining variables.
[[nodiscard]] MPoly mpoly_leading_coeff_in(const MPoly& p, std::size_t var_index);

// Evaluate variable var_index := value (BigInt) → MPoly in remaining var space
// (the evaluated variable's exponent becomes 0; var ordering/length preserved).
[[nodiscard]] MPoly mpoly_eval_var(const MPoly& p, std::size_t var_index, const BigInt& value);

// Collapse an MPoly that is univariate in var_index (all other exps 0) to IntPoly.
// nullopt if it actually depends on another variable.
[[nodiscard]] std::optional<IntPoly> mpoly_to_intpoly(const MPoly& p, std::size_t var_index);
[[nodiscard]] MPoly mpoly_from_intpoly(const IntPoly& f, std::size_t var_index, std::size_t nvars);

// ---- Multivariate Hensel lifting (factor_multivariate_hensel.cpp) ----
// Lift the univariate factors u[i] (in main var, with correct integer LCs already
// imposed) up to factors of `a` in Z[x_1,...,x_n], using ideal-adic Hensel lifting
// over the ideal I = <x_2 - a_2, ..., x_n - a_n>.  evaluation_point[k] is the chosen
// integer for variable index (k+1) (i.e. x_2..x_n).  Returns the lifted multivariate
// factors, or Unimplemented with a diagnostic at any cap / unsupported case.
[[nodiscard]] Result<std::vector<MPoly>> wang_multivariate_hensel(
    const MPoly& a,
    const std::vector<MPoly>& lifted_lc,      // intended leading coeff (MPoly) per factor
    const std::vector<IntPoly>& univar_factors,
    const std::vector<BigInt>& evaluation_point,  // size nvars-1, for x_2..x_n
    const WangContext& wc,
    symbolic::CASContext& ctx);

// ---- Leading-coefficient determination (factor_multivariate_lc.cpp) ----
// Wang's leading-coefficient distribution.  Given the full multivariate poly `a`
// (primitive, squarefree in main var), its univariate image factors, and the
// evaluation point, decide the multivariate leading coefficient assigned to each
// univariate factor.  Returns one MPoly LC per factor (in main-var-free var space)
// and possibly rescales the univariate factors so their integer LC matches the
// evaluated assigned LC.  Returns Unimplemented if the LC cannot be distributed
// uniquely (e.g. ambiguous factor LC matching).
struct LcDistribution {
    std::vector<MPoly> lc;            // assigned multivariate LC per factor
    std::vector<IntPoly> adjusted;    // univariate factors, LC-adjusted
    BigInt overall_constant{BigInt(1)};  // residual integer multiplier
};
[[nodiscard]] Result<LcDistribution> wang_distribute_leading_coeff(
    const MPoly& a,
    const std::vector<IntPoly>& univar_factors,
    const std::vector<BigInt>& evaluation_point,
    const WangContext& wc,
    symbolic::CASContext& ctx);

}  // namespace cas::algebra
