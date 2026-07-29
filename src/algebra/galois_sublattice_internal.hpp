// A6 Brick 3.5 — exhaustive MAXIMAL transitive-subgroup classes INSIDE a
// small group H (the descent nodes below the first Stauduhar layer).
//
// After the first descent step the current group H is small enough that
// the universe of the enumeration is [0, |H|) — NOT [0, n!) — and the
// dense single-generator-extension technique of
// galois_transitive_lattice.cpp applies verbatim: every subgroup of H is
// reachable by adding one generator at a time, so completeness holds BY
// CONSTRUCTION (no maximal-subgroup classifications transcribed from the
// literature, REGOLA 0.1). The returned classes are deduplicated up to
// H-conjugacy (exactly the equivalence the Stauduhar test scans: σ ranges
// over the CURRENT group) and then filtered to the MAXIMAL transitive
// ones — since the full lattice is in hand, maximality is exact, and it
// is what keeps the later resolvents small: [H:K] is the resolvent
// degree, and the p-adic precision bound grows with it. Soundness of the
// filter: any proper transitive G_f ≤ H lies in a maximal subgroup M of
// H, and M ⊇ G_f is itself transitive, so scanning the maximal
// transitive classes never misses a descent (and "no class contains G_f"
// certifies G_f = H exactly).
//
// Scale: two independent budgets, both caller-supplied. `max_ops` bounds
// time (element enumeration, Cayley table, double-coset marking, closure,
// conjugate registration, maximality filter). `max_bytes` bounds memory —
// the dense u16 Cayley table alone is 2·|H|² bytes, and the conjugate
// dedup sets add 2·|K| bytes per distinct conjugate — because the biggest
// interior nodes (the 5|2-block wreath family of S_10/A_10, |H| = 28800 /
// 14400) would need a 1.66 / 0.41 GiB table: the below-first-layer walk
// pre-dispatches those (via dense_sublattice_min_bytes) to the structural
// wreath-preimage route of galois_wreath_maximal (A6 Brick 3.75) instead
// of ever attempting the allocation here.
//
// Brick 3.75 additions: the same exhaustive core also serves
//   • maximal_subgroup_classes_in — ALL maximal proper classes (no
//     transitivity filter), used to lift the maximal subgroups of the
//     small sign quotient Q_H ≤ C₂ ≀ S_k through φ (family FA of the
//     structural route);
//   • all_subgroup_classes_in — every class, used by the machine checks
//     that verify the Scott-lemma steps of the structural route on dense
//     ground truth (e.g. the subdirect subgroups of A₅ × A₅).

#pragma once

#include "perm_bsgs_internal.hpp"

#include "cas/result.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace cas::symbolic {
class CASContext;
}

namespace cas::algebra::permgrp {

// Lower bound (exact for the fixed structural allocations) of the memory
// the dense enumeration spends for a node of this order and degree: the
// element store (2·m·n), the Cayley table (2·m²) and the inverse index
// (2·m). The per-class conjugate sets come on top, so a node whose lower
// bound already exceeds the byte budget can NEVER succeed here — the
// below-first-layer walk uses this predicate to route such nodes to the
// structural wreath-preimage generator BEFORE any allocation.
[[nodiscard]] std::uint64_t dense_sublattice_min_bytes(std::uint64_t order,
                                                       std::size_t degree);

// All classes of MAXIMAL transitive PROPER subgroups of H up to
// H-conjugacy, sorted by decreasing order. H itself, intransitive
// subgroups (they can never contain the transitive G_f) and transitive
// classes strictly contained in an H-conjugate of another proper
// transitive class are excluded. Errors:
//   Unimplemented   — ops or byte budget exhausted (raise
//                     CASContext::galois_lattice_max_ops /
//                     galois_sublattice_max_bytes), or |H| beyond the
//                     u16-index universe (2^16 — structural);
//   InvalidArgument — zero budget.
[[nodiscard]] Result<std::vector<BsgsGroup>> transitive_subgroup_classes_in(
    const BsgsGroup& H, std::uint64_t max_ops, std::uint64_t max_bytes,
    symbolic::CASContext* ctx);

// All classes of MAXIMAL PROPER subgroups of H up to H-conjugacy (no
// transitivity filter), sorted by decreasing order. Same budgets/errors.
[[nodiscard]] Result<std::vector<BsgsGroup>> maximal_subgroup_classes_in(
    const BsgsGroup& H, std::uint64_t max_ops, std::uint64_t max_bytes,
    symbolic::CASContext* ctx);

// EVERY subgroup class of H up to H-conjugacy (including the trivial
// group and H itself), sorted by decreasing order. Exhaustive ground
// truth for machine verification of structural constructions; same
// budgets/errors.
[[nodiscard]] Result<std::vector<BsgsGroup>> all_subgroup_classes_in(
    const BsgsGroup& H, std::uint64_t max_ops, std::uint64_t max_bytes,
    symbolic::CASContext* ctx);

}  // namespace cas::algebra::permgrp
