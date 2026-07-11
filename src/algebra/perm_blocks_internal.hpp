// A6 — Block systems of a transitive permutation group, Brick 2 of the
// Stauduhar deg ≥ 8 closure.
//
// A block system of a transitive G ≤ S_n is a G-invariant partition of the
// point set. The *minimal* (finest non-trivial) systems are computed exactly
// with Atkinson's union-find algorithm (M.D. Atkinson, "An algorithm for
// finding the blocks of a permutation group", Math. Comp. 29 (1975)): for a
// fixed point 0 and each candidate companion β ≠ 0, the finest G-congruence
// identifying {0, β} is the union-find closure of that merge under the
// generators. G is primitive iff every such closure is the full point set.
//
// Everything here is exact combinatorics on generators — no group elements
// are enumerated (works on top of BsgsGroup / raw generator sets), no floats,
// no tables. Pure functions, ctx-free like the rest of permgrp.

#pragma once

#include "perm_group_internal.hpp"

#include "cas/result.hpp"

#include <cstddef>
#include <vector>

namespace cas::algebra::permgrp {

// A G-invariant partition of {0..n-1} into num_blocks blocks of equal size
// n / num_blocks (blocks of a transitive group are equicardinal).
// block_of[p] = index of the block containing point p; block indices are
// normalised to first-appearance order (block_of[0] == 0), so two systems
// are equal as partitions iff their block_of vectors are equal.
struct BlockSystem {
    std::vector<std::size_t> block_of;
    std::size_t num_blocks{0U};

    [[nodiscard]] bool operator==(const BlockSystem& other) const {
        return block_of == other.block_of;
    }
};

// All minimal (finest) non-trivial block systems of the transitive group
// ⟨gens⟩ ≤ S_n, deduplicated. Empty result ⇔ the group is primitive.
// Fails structured if a generator is invalid or the group is intransitive
// (blocks of an intransitive action are not equicardinal — out of scope).
[[nodiscard]] Result<std::vector<BlockSystem>> minimal_block_systems(
    std::size_t n, const std::vector<Perm>& gens);

// True iff ⟨gens⟩ is primitive on {0..n-1} (transitive with no non-trivial
// block system). Same preconditions as minimal_block_systems.
[[nodiscard]] Result<bool> is_primitive(std::size_t n,
                                        const std::vector<Perm>& gens);

// Induced action on the blocks of `sys`: for each generator g of the group,
// the permutation ḡ of {0..num_blocks-1} with ḡ(i) = block containing the
// g-image of (any point of) block i. Fails structured if `sys` is not
// actually invariant under some generator (caller bug).
[[nodiscard]] Result<std::vector<Perm>> block_action_gens(
    std::size_t n, const std::vector<Perm>& gens, const BlockSystem& sys);

}  // namespace cas::algebra::permgrp
