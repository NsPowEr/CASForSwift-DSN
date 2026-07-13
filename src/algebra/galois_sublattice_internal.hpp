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
// 14400) would need a 1.66 / 0.41 GiB table: those await the structural
// wreath-maximal route and must fail STRUCTURED here, before allocating.

#pragma once

#include "perm_bsgs_internal.hpp"

#include "cas/result.hpp"

#include <cstdint>
#include <vector>

namespace cas::symbolic {
class CASContext;
}

namespace cas::algebra::permgrp {

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

}  // namespace cas::algebra::permgrp
