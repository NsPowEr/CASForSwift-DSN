// CAS-L3-18 — Galois toolkit MVP for low-degree polynomials.
//
// Identifies Galois group of irreducible polynomial f ∈ Q[x] of
// degree ≤ 4 via discriminant + resolvent analysis (when applicable).
//
// Returned identifier strings (subset of standard group notation):
//   "trivial"  : split completely over Q (e.g. (x-1)(x-2) → Galois grp = 1)
//   "C2"       : cyclic of order 2 (quadratic irreducible)
//   "A3"       : alternating on 3 letters (cyclic of 3, irreducible cubic with square disc)
//   "S3"       : symmetric on 3 letters (irreducible cubic, non-square disc)
//   "C4"       : cyclic of 4 (irreducible quartic with specific resolvent)
//   "V4"       : Klein four (biquadratic without sqrt(2)·sqrt(3) etc.)
//   "D4"       : dihedral of 8 (generic irreducible quartic with sqrt resolvent)
//   "A4"       : alternating on 4 (square disc + irreducible resolvent cubic)
//   "S4"       : symmetric on 4 (generic irreducible quartic)
//   "C5"       : cyclic of order 5
//   "D5"       : dihedral of order 10
//   "F20"      : Frobenius metacyclic group AGL(1,F_5), order 20
//   "A5"       : alternating on 5 letters (order 60, simple non-abelian)
//   "S5"       : symmetric on 5 letters (order 120)
//   "unknown"  : degree > 5 or analysis inconclusive

#pragma once

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <string>

namespace cas::algebra {

// Identify Galois group of monic integer polynomial f ∈ Z[x].
// Degree 2: returns "trivial" (rational roots) or "C2".
// Degree 3: returns "trivial" (split), "A3" (disc square), or "S3".
// Degree 4: returns "trivial", "V4", "C4", "D4", "A4", "S4" via
//           resolvent cubic + discriminant — deferred follow-up.
// Higher degrees: "unknown".
[[nodiscard]] Result<std::string> galois_group(
    ExprPtr poly, const Symbol& var, symbolic::CASContext& ctx);

}  // namespace cas::algebra
