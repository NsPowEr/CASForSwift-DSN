// A6 — On-demand maximal-transitive-subgroup candidates for S_n and A_n,
// Brick 2 of the Stauduhar deg ≥ 8 closure.
//
// For an irreducible f ∈ Q[x] the Galois group G_f is transitive; Stauduhar
// descent from the ambient group (S_n, or A_n when disc(f) is a square)
// tests G_f ≤ H^σ against the maximal subgroups H of the ambient. Since a
// transitive subgroup only ever lies inside a *transitive* maximal, the
// intransitive maximals (S_k × S_{n−k}) are omitted by design.
//
// The candidate lists are ASSEMBLED ALGORITHMICALLY from the structure of n
// (divisors → imprimitive wreaths; n = p^d → affine; n = q+1 → Möbius and
// semilinear; n = (p^d−1)/(p−1) → projective; n = C(m,k) → subset actions)
// plus the even parts (∩ A_n, via Schreier) and, where a class may split in
// A_n, an S_n-conjugate twin. No permutation or group table is transcribed
// (CLAUDE.md REGOLA 0.1).
//
// Soundness contract (what the descent may rely on):
//   • every listed candidate is a proper transitive subgroup of the ambient;
//   • COVERAGE — every maximal transitive subgroup of the ambient is
//     S_n-conjugate to a listed candidate, and every A_n-class of it is
//     A_n-conjugate to one (twin mechanism). Coverage is machine-verified
//     against the exhaustively generated subgroup lattice for n ≤ 7 in the
//     unit tests; for 8 ≤ n ≤ 10 it rests on the published enumerations of
//     transitive groups of small degree (Sims 1970; Butler & McKay, "The
//     transitive groups of degree up to eleven", Comm. Algebra 11 (1983)),
//     cross-validated end-to-end by the Brick-4 polynomial oracle corpus.
//   • a listed candidate need NOT be maximal (a redundant deeper subgroup
//     only costs an extra resolvent test, never a wrong answer).
// Degrees outside [5, 10] fail with a structured Unimplemented: n ≤ 4 has
// exact closed-form treatment elsewhere, n ≥ 11 awaits the coverage
// certification of a later increment.

#pragma once

#include "perm_bsgs_internal.hpp"

#include "cas/result.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace cas::algebra::permgrp {

enum class AmbientGroup : std::uint8_t { Symmetric, Alternating };

struct MaximalCandidate {
    BsgsGroup group;
    // Structural provenance ("even part of S_4 wr S_2", "affine AGL(3,2)",
    // …) — derived from the construction, never a transcribed label.
    std::string provenance;
    // Primitive in its action on {0..n-1} (exact, via block systems).
    bool primitive{false};
};

// The descent candidates below `ambient` in degree n (5 ≤ n ≤ 10).
// For AmbientGroup::Symmetric the list always contains A_n first.
[[nodiscard]] Result<std::vector<MaximalCandidate>>
maximal_transitive_candidates(AmbientGroup ambient, std::size_t n);

}  // namespace cas::algebra::permgrp
