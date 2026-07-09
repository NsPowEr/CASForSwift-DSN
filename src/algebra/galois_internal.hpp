// CAS-L3-18 / F3.6 — Galois module shared helpers.
//
// Provides primitives reused across galois.cpp (degrees 2-4) and
// galois_deg5.cpp (degree 5 Soicher-McKay / Frobenius / Dedekind path).
//
// Anti-monolith split per CLAUDE.md "≤500 LOC per source file":
// adding deg-5 inline to galois.cpp would push it well past 500 LOC.

#pragma once

#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <string>

namespace cas::algebra {

// Identify Galois group of an irreducible monic quintic f ∈ Z[x] via
// Frobenius / Dedekind cycle-type evidence plus discriminant parity.
//
// Algorithm (structural, exact, no floating point):
//   1. Compute disc(f); decide if it is a rational square in Q.
//      • disc square  → group ⊆ A_5  (so ∈ {C5, D5, A5})
//      • disc nonsq.  → group ⊄ A_5  (so ∈ {F20, S5})
//   2. For up to ctx.max_galois_frobenius_primes() small primes p with
//      p ∤ lc(f) and p ∤ disc(f), factor f mod p and read the cycle type
//      (Dedekind's theorem: multiset of irreducible-factor degrees = cycle
//      type of Frobenius_p inside Gal(f/Q)).
//   3. Use the cycle-type membership table to identify the group:
//        C5  : only {1^5, 5^1}
//        D5  : {1^5, 5^1, 2^2·1}
//        F20 : {1^5, 5^1, 2^2·1, 4^1·1}
//        A5  : {1^5, 5^1, 2^2·1, 3^1·1^2}
//        S5  : all of the above + odd patterns {2·1^3, 3·2}
//   4. Cross-check vs disc parity. If insufficient evidence within prime
//      budget, return Result<>::Unimplemented (NOT silent guess).
//
// Public API entry-point — invoked from galois_group() in galois.cpp when
// the input has total degree 5.
[[nodiscard]] Result<std::string> galois_group_quintic_irreducible(
    ExprPtr poly, const Symbol& var, symbolic::CASContext& ctx);

// A6 — identify Galois group of an irreducible sextic f ∈ Q[x] via the
// fully exact pipeline (discriminant parity + Dedekind sieve + 2-set
// resolvent) matched against the exhaustively *generated* transitive
// lattice of S₆ (no transcribed group tables). Ambiguity → structured
// Unimplemented listing the surviving candidates (never a guess).
// Implemented in galois_deg6.cpp.
[[nodiscard]] Result<std::string> galois_group_sextic_irreducible(
    ExprPtr poly, const Symbol& var, symbolic::CASContext& ctx);

// Returns true if rational is a perfect square in Q (i.e. p/q with
// p, q both perfect squares as BigInt and same sign). Shared with deg-4 path.
[[nodiscard]] bool is_rational_square_q(const Rational& r);

// Try to interpret expr as Rational. Returns nullopt otherwise. Recognises
// IntegerLit, RationalLit, and Unary(Neg, ...) wrappers. Shared helper.
[[nodiscard]] std::optional<Rational> as_rational_q(ExprPtr e);

}  // namespace cas::algebra
