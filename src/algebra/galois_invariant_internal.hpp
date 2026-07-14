// A6 Brick 3b — exact relative invariants and coset transversals for the
// Stauduhar descent (float-free Fieker-Klüners route; spec Galois_Groups.md
// read per REGOLA 0.1).
//
// Given groups H ≤ G ≤ S_n (BSGS representations), produce a polynomial
// F ∈ Z[x_0..x_{n−1}] — a sum of DISTINCT monomials, closed under H — with
//
//     Stab_G(F) = H        (EXACTLY, certified, never assumed),
//
// together with the [G:H] left-coset representatives σ of H in G, realised
// as the G-orbit of F (σ ↦ σ·F is a bijection cosets → images precisely
// because the stabiliser is exact).
//
// Construction is fully algorithmic, no transcribed invariants:
//   • tier 1 — H-orbit sums of squarefree monomials (k-subset indicator
//     monomials, k = 1..n−1, one candidate per H-orbit; the full orbit
//     e_k is skipped: its stabiliser is all of G);
//   • tier 2 — H-orbit sums of j-tuple monomials x_{t₀}¹x_{t₁}²⋯x_{t_j}^j
//     over ordered tuples of distinct indices, ascending j = 2..n−1: the
//     smallest j whose point stabilisers differ already separates
//     subgroups invisible to set-orbits (e.g. sign-character kernels in
//     wreath nodes) with an orbit-sum FAR smaller than the full Galois
//     monomial. The last level j = n−1 IS the classical Galois resolvent
//     monomial m* = x_1 x_2² ⋯ x_{n−1}^{n−1}, whose S_n-stabiliser is
//     trivial, so Stab_G(Σ_{h∈H} h·m*) = H unconditionally (guaranteed
//     terminator — the search never fails, it only prefers small
//     invariants because the later p-adic precision bound grows with the
//     degree and the term count).
// Exactness of every candidate is CERTIFIED by counting: the G-orbit of F
// is enumerated and accepted iff its size equals [G:H] (Stab_G(F) ⊇ H
// always holds by construction, so equality of counts pins Stab = H).
//
// H is NOT assumed maximal in G — redundant deeper candidates from
// perm_maximal.cpp stay sound (the Brick-2 contract).

#pragma once

#include "perm_bsgs_internal.hpp"

#include "cas/result.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace cas::symbolic {
class CASContext;
}

namespace cas::algebra::galois_invariant {

// Exponent vector on x_0..x_{n−1} (exponents ≤ n−1 ≤ 19 fit easily).
using Monomial = std::vector<std::uint8_t>;

struct RelativeInvariant {
    // The invariant F = Σ monomials (distinct, sorted lexicographically).
    std::vector<Monomial> monomials;
    std::size_t total_degree{0U};  // degree of each monomial's image bound
    // Left-coset representatives of H in G; reps[0] is the identity and
    // reps[i]·F are pairwise distinct polynomials.
    std::vector<permgrp::Perm> coset_reps;
};

// Computes the certified relative invariant and transversal. `max_ops`
// caps the total number of orbit-node expansions across the search (a
// pure anti-runaway belt — termination is proven; on exhaustion returns a
// structured Unimplemented so the caller can raise the budget via
// CASContext). Errors:
//   InvalidArgument — degree mismatch, H ⊄ G, or H = G (index 1: nothing
//                     to descend to);
//   Unimplemented   — ops budget exhausted;
//   InternalError   — a certified impossibility (orbit exceeding [G:H]).
[[nodiscard]] Result<RelativeInvariant> relative_invariant(
    const permgrp::BsgsGroup& G, const permgrp::BsgsGroup& H,
    std::uint64_t max_ops, symbolic::CASContext* ctx);

// The image σ·F of an invariant under a permutation (variables x_i ↦
// x_{σ(i)}), returned sorted — the canonical form used for orbit keys.
[[nodiscard]] std::vector<Monomial> apply_perm_to_invariant(
    const permgrp::Perm& sigma, const std::vector<Monomial>& monomials);

}  // namespace cas::algebra::galois_invariant
