#pragma once
// Shared internals for the Zippel sparse-interpolation GCD (F3.1 / T-006).
// The per-prime monic core lives in polynomial_gcd_zippel_prony.cpp; the
// multi-prime CRT + Farey rational-reconstruction fallback lives in
// polynomial_gcd_zippel_crt.cpp. Both share the sparse Z-polynomial bridge and
// the exact-division certificate declared here. A dedicated sub-namespace avoids
// clashing with the like-named static helpers in polynomial_gcd_brown_modular.cpp.

#include "cas/algebra.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <map>
#include <vector>

namespace cas::algebra::zippel_detail {

using ZMonomial   = std::vector<unsigned int>;
using ZSparsePoly = std::map<ZMonomial, BigInt>;

// Sparse Z[x_1..x_n] bridge (monomial exponent-vectors aligned to `vars`).
[[nodiscard]] ZSparsePoly to_sparse_z(
    const MultivariatePolynomial& p, const std::vector<Symbol>& vars);
[[nodiscard]] MultivariatePolynomial from_sparse_z(
    const ZSparsePoly& sp, const std::vector<Symbol>& vars);

// Exact certificate: does b divide a in Z[x_1..x_n]? (sparse leading-term division)
[[nodiscard]] bool certify_divides(
    const MultivariatePolynomial& a, const MultivariatePolynomial& b,
    const std::vector<Symbol>& vars);

// Maximality certificate (A37/A40): g_cand is the FULL gcd (not a proper
// divisor) iff its cofactors are coprime, i.e. gcd(P/g, Q/g) = 1. Needed because
// certify_divides above passes for any common divisor, while the per-sample
// Fp-monic step of the Prony core can drop a non-constant LC_{x1}(gcd).
[[nodiscard]] Result<bool> is_maximal_gcd_candidate(
    const MultivariatePolynomial& P, const MultivariatePolynomial& Q,
    const MultivariatePolynomial& g_cand, const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx);

// Sign-normalizes g_cand and returns it only if maximality is certified;
// otherwise Unimplemented (ZIPPEL_PRONY_NOT_MAXIMAL) so the dispatcher falls
// back to gcd_brown_modular. `stage` names the producing path for diagnostics.
[[nodiscard]] Result<MultivariatePolynomial> finish_if_maximal(
    const MultivariatePolynomial& P, const MultivariatePolynomial& Q,
    MultivariatePolynomial g_cand, const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx, std::size_t* out_samples_used,
    std::size_t samples_used, const char* stage);

// One prime of Zippel-Prony: returns the gcd made monic-in-x_1 over F_p, as a
// sparse map with coefficients in [0, p). Errors on an unlucky prime (unstable
// degree, singular Vandermonde, skeleton rejected) so the caller tries the next.
[[nodiscard]] Result<ZSparsePoly> zippel_gcd_modp(
    const ZSparsePoly& spP, const ZSparsePoly& spQ,
    const std::vector<Symbol>& vars, const BigInt& p,
    symbolic::CASContext& ctx, std::size_t* out_samples);

// Multi-prime fallback: CRT-combines monic g mod pᵢ across a sequence of primes
// (starting at first_prime), Farey-reconstructs each coefficient to a rational,
// clears to the primitive integer gcd, and certifies by exact division. Used only
// when the single-prime center-lift fails its certificate (large/non-unit lc).
[[nodiscard]] Result<MultivariatePolynomial> gcd_zippel_prony_crt(
    const ZSparsePoly& spP, const ZSparsePoly& spQ,
    const std::vector<Symbol>& vars, const BigInt& first_prime,
    symbolic::CASContext& ctx, std::size_t* out_samples);

}  // namespace cas::algebra::zippel_detail
