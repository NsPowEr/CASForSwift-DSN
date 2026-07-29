// A6 / CAS-L3-18 — Exact set-resolvent machinery (internal header).
// See galois_setpoly.cpp for the mathematics and certificates.

#pragma once

#include "algebraic_tower_resultant.hpp"  // primitive_internal::Deadline
#include "cas/bigint.hpp"
#include "cas/result.hpp"
#include "polynomial_internal.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace cas::algebra::galois_setpoly {

// R_k(y) = ∏_{|T|=k}(y − Σ_{i∈T} α_i) for monic f = ∏(x − α_i), degree
// C(n,k) — computed exactly via power sums (Newton) and the EGF
// recurrence for elementary symmetric functions in e^{s·α_i}; no
// resultants. Structural certificate: q_0 = C(n,k).
[[nodiscard]] Result<RatPoly> kset_resolvent(
    const RatPoly& f_monic, std::size_t k,
    const primitive_internal::Deadline& deadline = std::nullopt);

// R₂(y) = ∏_{i<j}(y − α_i − α_j): kset_resolvent with k = 2.
[[nodiscard]] Result<RatPoly> two_set_resolvent(
    const RatPoly& f_monic,
    const primitive_internal::Deadline& deadline = std::nullopt);

// R₃(y) = ∏_{i<j<k}(y − α_i − α_j − α_k): kset_resolvent with k = 3.
[[nodiscard]] Result<RatPoly> three_set_resolvent(
    const RatPoly& f_monic,
    const primitive_internal::Deadline& deadline = std::nullopt);

// gcd(p, p′) constant ⇔ squarefree over Q.
[[nodiscard]] Result<bool> is_squarefree_q(const RatPoly& p);

// Minimal-polynomial candidate of β = P(α), P(x) = Σ_{m=1..|c|} c_m·x^m
// (1 ≤ deg P ≤ n−1): g(y) = Res_x(f(x), y − P(x)), monic-normalized.
// The Galois action on the β_i = P(α_i) is the same permutation
// representation as on the α_i (P ∈ Q[x]), so k-set orbit reading is
// preserved whenever g is squarefree of degree n (P injective on roots).
// Completeness of the power basis: the map S ↦ (Σ_{α∈S} α^m)_{m=1..n−1}
// is injective on root subsets (Vandermonde), so for every pair of
// distinct k-subsets some coefficient vector separates their β-sums.
[[nodiscard]] Result<RatPoly> tschirnhaus_general(
    const RatPoly& f_monic, const std::vector<BigInt>& c,
    const primitive_internal::Deadline& deadline = std::nullopt);

// β = α² + c·α: degree-2 special case of tschirnhaus_general.
[[nodiscard]] Result<RatPoly> tschirnhaus_quadratic(
    const RatPoly& f_monic, const BigInt& c,
    const primitive_internal::Deadline& deadline = std::nullopt);

}  // namespace cas::algebra::galois_setpoly
