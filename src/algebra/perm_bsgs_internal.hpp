// A6 / CAS-L3-18 — Scalable permutation group via a verified BSGS
// (base + strong generating set), Brick 1 of the Stauduhar deg ≥ 8 closure.
//
// The dense `PermGroup` (perm_group_internal.hpp) stores the whole element set
// (Θ(|G|) ranks) and the transitive-lattice enumerator stores Θ(n!) bitsets per
// subgroup — both are infeasible past n = 7. A stabilizer chain represents any
// subgroup of S_n by a base B = [b_1,…,b_k] and a strong generating set S in
// polynomial space, giving the exact order and an exact membership test without
// ever materialising the elements.
//
// The chain is built by DETERMINISTIC Schreier-Sims (Holt-Eick-O'Brien
// "Handbook of Computational Group Theory", §4.4): every Schreier generator is
// sifted through the partial chain and, on a non-identity residue, added as a
// new strong generator (extending the base when the residue fixes it). No
// Monte-Carlo shortcut — the resulting BSGS is verified, so order() and
// contains() are exact (never a probabilistic guess; REGOLA ZERO).
//
// Reuses the Perm primitives of perm_group_internal.hpp (identity/compose/
// inverse/is_odd/is_valid_perm) — no duplication. Pure combinatorics, ctx-free.

#pragma once

#include "perm_group_internal.hpp"

#include "cas/result.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace cas::algebra::permgrp {

// Verified base-and-strong-generating-set representation of a subgroup of S_n.
class BsgsGroup {
public:
    // Deterministic Schreier-Sims closure of ⟨gens⟩ inside S_n. Fails structured
    // if a generator is invalid, or n exceeds the u64-order-safe cap (n ≤ 20,
    // matching factorial_u64: |G| ≤ n! ≤ 20! < 2^62). The empty generator set
    // yields the trivial group.
    [[nodiscard]] static Result<BsgsGroup> build(std::size_t n,
                                                 std::vector<Perm> gens);

    [[nodiscard]] std::size_t degree() const noexcept { return n_; }

    // The original generators (⟨gens⟩ = this group).
    [[nodiscard]] const std::vector<Perm>& generators() const noexcept {
        return input_gens_;
    }
    // The strong generating set relative to base_ (a superset of the input
    // generators; every level's generators fix all earlier base points).
    [[nodiscard]] const std::vector<Perm>& strong_generators() const noexcept {
        return strong_gens_;
    }
    [[nodiscard]] const std::vector<std::size_t>& base() const noexcept {
        return base_;
    }

    // |G| = ∏_i |orbit of b_i in G^(i)| (product of transversal sizes). Exact.
    [[nodiscard]] std::uint64_t order() const noexcept;

    // Membership via sifting: p ∈ G ⇔ p strips to the identity along the chain.
    [[nodiscard]] bool contains(const Perm& p) const;

    // Residue of the strip of p through the chain (identity ⇔ p ∈ G). Exposed
    // for tests; contains(p) == (sift(p) is the identity).
    [[nodiscard]] std::optional<Perm> sift(const Perm& p) const;

    // Natural action on {0..n-1}: orbit of point 0 under the strong generators
    // has size n. (For a group built from its generators this is the standard
    // transitivity test.)
    [[nodiscard]] bool is_transitive() const;

    // True iff G ⊄ A_n, i.e. some strong generator is an odd permutation
    // (parity is a homomorphism, so this is exact from the generators).
    [[nodiscard]] bool has_odd_element() const;

private:
    BsgsGroup() = default;

    std::size_t n_{0U};
    std::vector<Perm> input_gens_;
    std::vector<std::size_t> base_;             // b_1..b_k
    std::vector<Perm> strong_gens_;             // full strong generating set
    // level_gens_[i] = strong generators fixing b_0..b_{i-1} pointwise.
    std::vector<std::vector<Perm>> level_gens_;
    // transversal_[i][x] = a coset rep u with b_i^u = x (empty if x ∉ orbit).
    std::vector<std::vector<std::optional<Perm>>> transversal_;
};

}  // namespace cas::algebra::permgrp
