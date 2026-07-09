// A6 / CAS-L3-18 — Exact permutation-group engine for Galois descent (deg ≥ 6).
//
// Foundation for the Stauduhar/Soicher-McKay machinery: permutations of
// {0..n-1} in image form, Lehmer ranking for dense S_n indexing, and
// subgroups represented explicitly as the sorted rank-set of their elements
// (exhaustive closure). Everything here is exact combinatorics — no floats,
// no heuristics, no group tables: subgroups are *generated*, never copied
// from literature (anti-hallucination, CLAUDE.md REGOLA 0.1).
//
// Scale contract: the dense representation targets small degrees (n ≤ 12,
// rank fits in uint32 since 12! < 2^32). Every closure takes an explicit
// `max_order` budget (caller derives it from CASContext) and fails with a
// structured Unimplemented — never silently truncates (Cat-1/Cat-4).
//
// References:
//   • Holt, Eick, O'Brien — "Handbook of Computational Group Theory", §2-4.
//   • Stauduhar 1973, "The determination of Galois groups".
//   • Soicher & McKay 1985, "Computing Galois groups over the rationals".

#pragma once

#include "cas/result.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas::algebra::permgrp {

// A permutation of {0..n-1} in image form: img[i] = σ(i).
using Perm = std::vector<std::uint8_t>;

[[nodiscard]] Perm identity(std::size_t n);
// (a∘b)(i) = a(b(i)) — apply b first, then a.
[[nodiscard]] Perm compose(const Perm& a, const Perm& b);
[[nodiscard]] Perm inverse(const Perm& a);
[[nodiscard]] bool is_valid_perm(const Perm& a);
// Parity: true iff σ is an odd permutation.
[[nodiscard]] bool is_odd(const Perm& a);
// Cycle lengths of σ (including fixed points), sorted descending.
[[nodiscard]] std::vector<std::size_t> cycle_type(const Perm& a);

// n! (n ≤ 20 — fits uint64; asserted).
[[nodiscard]] std::uint64_t factorial_u64(std::size_t n);
// Lehmer rank: bijection S_n → [0, n!) preserving nothing but identity→0.
[[nodiscard]] std::uint32_t lehmer_rank(const Perm& a);
[[nodiscard]] Perm lehmer_unrank(std::uint32_t r, std::size_t n);

// Explicit (dense) subgroup of S_n: sorted Lehmer ranks of all elements.
// Built exclusively via `closure` so the invariant "element set is a group"
// holds by construction.
class PermGroup {
public:
    // Exhaustive closure of `gens` inside S_n. Fails with structured
    // Unimplemented if the generated order would exceed `max_order`
    // (budget owned by the caller — derived from CASContext, never magic).
    [[nodiscard]] static Result<PermGroup> closure(
        std::size_t n, std::vector<Perm> gens, std::uint64_t max_order);

    [[nodiscard]] std::size_t degree() const noexcept { return n_; }
    [[nodiscard]] std::uint64_t order() const noexcept {
        return static_cast<std::uint64_t>(ranks_.size());
    }
    [[nodiscard]] const std::vector<Perm>& generators() const noexcept {
        return gens_;
    }
    // Sorted Lehmer ranks of every element.
    [[nodiscard]] const std::vector<std::uint32_t>& element_ranks()
        const noexcept {
        return ranks_;
    }

    [[nodiscard]] bool contains(const Perm& p) const;
    [[nodiscard]] bool contains_rank(std::uint32_t r) const;
    // True iff every element of `other` lies in *this.
    [[nodiscard]] bool contains_group(const PermGroup& other) const;
    // Structural equality as subgroups of S_n (same element set).
    [[nodiscard]] bool equals(const PermGroup& other) const;

    // Natural action on points {0..n-1}.
    [[nodiscard]] bool is_transitive() const;
    // True iff some element is an odd permutation (i.e. G ⊄ A_n).
    [[nodiscard]] bool has_odd_element() const;

    // Orbit lengths of the induced action on k-subsets of {0..n-1},
    // sorted descending. This is the exact group-side invariant matched
    // against factorization patterns of set-resolvents (Soicher-McKay):
    // the irreducible factors of the k-set resolvent of f have degrees
    // equal to these orbit lengths when the resolvent is squarefree.
    // Requires n ≤ 31 (subset encoded as uint32 bitmask) and 1 ≤ k ≤ n.
    [[nodiscard]] std::vector<std::size_t> orbit_lengths_on_ksubsets(
        std::size_t k) const;

    // Multiset of cycle types over all elements, as a sorted list of
    // (cycle_type, count) pairs. Conjugation-invariant signature: cheap
    // pre-filter before full conjugacy testing in the lattice builder,
    // and the exact reference set for the Dedekind/Frobenius sieve.
    [[nodiscard]] std::vector<std::pair<std::vector<std::size_t>, std::size_t>>
    cycle_type_distribution() const;

private:
    PermGroup() = default;
    std::size_t n_{0U};
    std::vector<Perm> gens_;
    std::vector<std::uint32_t> ranks_;  // sorted
};

// True iff A and B are conjugate subgroups of S_n (σAσ⁻¹ = B for some
// σ ∈ S_n). Exhaustive scan over S_n with cheap invariant pre-filters —
// exact, feasible for the dense-lattice degrees (n ≤ 8).
[[nodiscard]] bool conjugate_in_sn(const PermGroup& a, const PermGroup& b);

// Computational budget for the exhaustive lattice enumeration, owned by
// the caller (derived from CASContext — never a magic literal in here).
struct LatticeBudget {
    // Hard cap on the degree: the dense enumeration is Θ(n!·#classes) in
    // time and Θ(n!) in memory per subgroup; past this the Stauduhar
    // maximal-subgroup descent (A6 follow-up increment) is required.
    std::size_t max_degree{0U};
    // Abort knob: rough upper bound on rank-composition operations. 0 is
    // NOT "unlimited" — the caller must always supply a positive budget.
    std::uint64_t max_ops{0U};
};

// All transitive subgroups of S_n up to conjugacy, enumerated exhaustively:
// BFS over the *full* subgroup lattice by single-generator extension
// (every finite group arises this way from the trivial group), with
// (H,H)-double-coset pruning (⟨H,s⟩ = ⟨H,h₁sh₂⟩) and conjugacy dedup.
// Completeness is by construction — no group tables are consulted
// (CLAUDE.md REGOLA 0.1 anti-hallucination).
// Budget exceeded / degree beyond cap → structured Unimplemented.
[[nodiscard]] Result<std::vector<PermGroup>> transitive_subgroup_classes(
    std::size_t n, const LatticeBudget& budget);

}  // namespace cas::algebra::permgrp
