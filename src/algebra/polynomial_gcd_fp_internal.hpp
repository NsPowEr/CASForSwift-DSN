// polynomial_gcd_fp_internal.hpp — shared types/helpers for Fp-recursive multivariate GCD.
// Used by polynomial_gcd_brown_modular.cpp and polynomial_gcd_fp_recursive.cpp.
//
// Not a public header; internal to src/algebra/.

#pragma once

#include "cas/algebra.hpp"
#include "cas/bigint.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <map>
#include <optional>
#include <vector>

namespace cas::algebra::fp_helpers {

using BMonomial   = std::vector<unsigned int>;
using BSparsePoly = std::map<BMonomial, BigInt>;

[[nodiscard]] BigInt pos_mod(const BigInt& a, const BigInt& m);
[[nodiscard]] BigInt bigint_gcd(BigInt a, BigInt b);
[[nodiscard]] BigInt centered_repr(const BigInt& r, const BigInt& M);

[[nodiscard]] BSparsePoly reduce_sparse_mod_p(const BSparsePoly& sp, const BigInt& p);
[[nodiscard]] BigInt      sparse_inf_norm(const BSparsePoly& sp);
[[nodiscard]] std::size_t deg_in_var(const BSparsePoly& sp, std::size_t var_idx);
[[nodiscard]] BSparsePoly eval_var_mod_p(const BSparsePoly& sp, std::size_t var_idx,
                                          const BigInt& val, const BigInt& p);

[[nodiscard]] std::optional<BSparsePoly> univariate_sparse_gcd_fp(
    const BSparsePoly& A, const BSparsePoly& B,
    std::size_t var_idx, const BigInt& p);

[[nodiscard]] std::optional<BSparsePoly> lagrange_interp_fp(
    const std::vector<BigInt>& bs, const std::vector<BSparsePoly>& vals,
    std::size_t var_idx, std::size_t target_deg, const BigInt& p);

[[nodiscard]] std::map<std::size_t, BSparsePoly> layers_by_var(
    const BSparsePoly& sp, std::size_t var_idx);

[[nodiscard]] BSparsePoly reassemble_layers(
    const std::map<std::size_t, BSparsePoly>& layers,
    std::size_t var_idx, const BigInt& p);

[[nodiscard]] BSparsePoly mul_mod_p(
    const BSparsePoly& A, const BSparsePoly& B, const BigInt& p);

[[nodiscard]] std::optional<BSparsePoly> exact_div_fp(
    const BSparsePoly& A, const BSparsePoly& B, const BigInt& p);

// Recursive Brown-in-Fp gcd of A,B in Fp[active_vars].
[[nodiscard]] std::optional<BSparsePoly> sparse_gcd_fp(
    const BSparsePoly& A, const BSparsePoly& B,
    const std::vector<std::size_t>& active_vars,
    const BigInt& p, symbolic::CASContext& ctx, std::size_t depth);

[[nodiscard]] BSparsePoly to_sparse(const MultivariatePolynomial& p,
                                    const std::vector<Symbol>& vars);

[[nodiscard]] MultivariatePolynomial from_sparse(const BSparsePoly& sp,
                                                 const std::vector<Symbol>& vars);

[[nodiscard]] std::vector<Symbol> collect_vars(const MultivariatePolynomial& p,
                                               const MultivariatePolynomial& q);

[[nodiscard]] BigInt multivar_mignotte_bound(const BSparsePoly& P,
                                             const BSparsePoly& Q);

[[nodiscard]] bool divides_sparse_z(
    const BSparsePoly& dividend,
    const BSparsePoly& divisor,
    std::size_t n_vars);

}  // namespace cas::algebra::fp_helpers

namespace cas::algebra {

// Lc-poly-scaling helpers (polynomial_gcd_brown_lc_scaling.cpp).
[[nodiscard]] fp_helpers::BSparsePoly extract_lc_in_var(
    const fp_helpers::BSparsePoly& sp, std::size_t main_var);

[[nodiscard]] fp_helpers::BSparsePoly multiply_sparse_mod_p(
    const fp_helpers::BSparsePoly& A, const fp_helpers::BSparsePoly& B,
    const BigInt& p, std::size_t n_vars);

void scale_by_lc(fp_helpers::BSparsePoly& gp, const BigInt& u, const BigInt& p);

[[nodiscard]] std::optional<BigInt> compute_lc_scalar_ratio(
    const fp_helpers::BSparsePoly& Lp, const fp_helpers::BSparsePoly& lcg,
    const BigInt& p);

[[nodiscard]] fp_helpers::BSparsePoly reduce_lc_bound_mod_p(
    const fp_helpers::BSparsePoly& L, const BigInt& p);

[[nodiscard]] MultivariatePolynomial sub_sparse_to_mv(
    const fp_helpers::BSparsePoly& sp, const std::vector<Symbol>& vars,
    std::size_t main_var);

[[nodiscard]] fp_helpers::BSparsePoly mv_to_sub_sparse(
    const MultivariatePolynomial& p, const std::vector<Symbol>& vars,
    std::size_t main_var);

// Decompose cand into main-var coefficient layers, compute their multivariate
// gcd h in n-1 sub-vars (recursive gcd_brown_modular), exactly divide cand by
// h, and write the quotient into `out`.  Returns true if the division succeeds
// (h | cand exactly) and h is non-trivial; false otherwise.
[[nodiscard]] bool remove_spurious_main_var_factor(
    const fp_helpers::BSparsePoly& cand, const std::vector<Symbol>& vars,
    std::size_t main_var, symbolic::CASContext& ctx,
    fp_helpers::BSparsePoly& out);

// Exact integer division of n-var sparse poly `dividend` by `divisor`.
// Returns true with `quo` populated iff division is exact in Z.
[[nodiscard]] bool exact_divide_sparse_z(
    const fp_helpers::BSparsePoly& dividend,
    const fp_helpers::BSparsePoly& divisor,
    std::size_t n_vars, fp_helpers::BSparsePoly& quo);

}  // namespace cas::algebra
