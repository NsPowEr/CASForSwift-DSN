// Kovacic Case 2 (dihedral D∞ subgroup) — SCAFFOLD.
//
// Algorithm body intentionally NOT implemented in this translation unit.
// Returns CASErrorKind::Unimplemented with explicit diagnostic so the
// dispatcher in ode_kovacic.cpp can route through Case 2 cleanly and
// surface a precise "Case 2 not yet implemented" message to the user
// instead of conflating with Case 1 failure.
//
// Algorithm spec (Steps 1-7) is documented in
// .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Kovacic_Case2.md
// (gitignored, agent-local).  Ledger tracking entry: HC-KV-03 in
// HARDCODE_LEDGER.md.
//
// Reference: Kovacic J.J. (1986), "An Algorithm for Solving Second
// Order Linear Homogeneous Differential Equations", J. Symbolic
// Computation 2, pp. 3-43, §3.2 (Case 2).
//
// Implementation gated on:
//   1. Verification of Step 1 necessary-condition predicates against
//      the Kovacic 1986 paper (poles order 1 / 2 / even ≥ 4, plus the
//      ord_∞ classification).
//   2. Verification of the Step 6 polynomial DE
//      (P'' + 3θP' + (3θ² + 3θ' - 4r)P = 0) and the Step 7 quadratic
//      for ω against the paper.
// REGOLA ZERO requires verification before any algorithmic code lands
// (silent wrong-answer on Bessel/Mathieu family is worse than explicit
// Unimplemented).

#include "ode_kovacic_internal.hpp"

namespace cas::calculus {
namespace kovacic_impl {

Result<OmegaPair> case2_omega(
    ExprPtr r, const Symbol& x, symbolic::CASContext& ctx) {
    (void)r;
    (void)x;
    (void)ctx;
    return fail<OmegaPair>(kv_unimpl(
        "Kovacic Case 2 (dihedral D∞ subgroup) not yet implemented. "
        "Algorithm scaffold present; full Steps 1-7 deferred to dedicated "
        "session with access to Kovacic 1986 §3.2. See "
        "HARDCODE_LEDGER.md HC-KV-03 and "
        ".APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Kovacic_Case2.md"));
}

} // namespace kovacic_impl
} // namespace cas::calculus
