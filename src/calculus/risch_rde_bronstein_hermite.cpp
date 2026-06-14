// Risch RDE parametric Hermite reduction — SCAFFOLD.
//
// Algorithm body intentionally NOT implemented in this translation unit.
// Returns CASErrorKind::Unimplemented with explicit diagnostic so the
// dispatcher in risch_rde_bronstein.cpp (solve_risch_de_parametric_field)
// can route through Hermite parametric path cleanly and surface a precise
// "RP-2 Hermite parametric not yet implemented" message instead of falling
// back to trial-constant ansatz (HC-F75-B-TRIAL-CONSTANTS).
//
// Algorithm spec is documented in
// .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Risch_Transcendental_Cap8.md
// (gitignored, agent-local).  Ledger tracking entry: HC-F8-PENDING-17 in
// HARDCODE_LEDGER.md (task ID 17, sub-step RP-2 per PLAN_NEXT_SESSIONS.md
// §Sessione 3 Step 3.1).
//
// Reference: Bronstein M. (2005), "Symbolic Integration I —
// Transcendental Functions", 2nd ed., Springer, §§5.1, 5.2, 6.5
// (parametric problems + PolyRischDE).  Davenport, J. (1981), "On the
// Risch Differential Equation Problem", SIGSAM Bulletin 15(3), 9-12.
//
// Implementation gated on:
//   1. Verification of Bronstein 6.5 PolyRischDE degree bound formula.
//   2. Verification of the Hermite reduction parametric formulation
//      (resolution of unknown coefficients via parametric linear systems
//      across an algebraic basis of the constant field).
//   3. Verification that the weak normalizer / structure theorem feeds
//      the correct denominator factorisation into the parametric solver.
// REGOLA ZERO requires verification before any algorithmic code lands
// (silent wrong-answer on integrate(exp(x²), x), elementary log-rational
// kernels, or differential-tower nested cases is worse than explicit
// Unimplemented).

#include "cas/symbolic.hpp"

namespace cas::calculus {

// Stub entry point for Risch RDE parametric Hermite reduction.
//
// Will be called by solve_risch_de_parametric_field as a higher-quality
// alternative to the existing trial-constant ansatz path when the
// parametric DE coefficients reside in a non-rational differential
// field extension (df > 0 in Bronstein notation).
//
// Currently always returns Unimplemented with the explicit HC-F8-
// PENDING-17 diagnostic so the caller can fall back to the existing
// trial-constants path without conflating diagnostics.
[[nodiscard]] Result<ExprPtr> risch_rde_hermite_parametric_stub(
    const symbolic::CASContext& ctx) {
    (void)ctx;
    return fail<ExprPtr>(CASError{
        .kind = CASErrorKind::Unimplemented,
        .message =
            "Risch RDE parametric Hermite reduction (RP-2) not yet "
            "implemented.  Scaffold present; full implementation deferred "
            "to dedicated session with access to Bronstein 2005 §§5.1, "
            "5.2, 6.5.  See HARDCODE_LEDGER.md HC-F8-PENDING-17 and "
            ".APROJECT_REFERENCES/MISSING_FEATURES_SPECS/"
            "Risch_Transcendental_Cap8.md (sub-step RP-2 per "
            "PLAN_NEXT_SESSIONS.md §Sessione 3 Step 3.1).  Existing "
            "trial-constant path in risch_rde_bronstein.cpp continues "
            "to serve callers; closure of HC-F75-B-TRIAL-CONSTANTS "
            "blocked on this scaffold's full impl.",
        .hint = std::nullopt,
    });
}

} // namespace cas::calculus
