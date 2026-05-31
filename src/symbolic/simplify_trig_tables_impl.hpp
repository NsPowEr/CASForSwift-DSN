#pragma once
// Internal header: Fermat-prime base angle constructors and reference tables.
// Included only by simplify_trig.cpp and simplify_trig_tables.cpp.
// NOT part of the public API.
//
// These functions build exact radical forms for the Gauss-constructible
// "primitive" angles at Fermat primes p=3,5 from which all constructible
// cos(p·π/q) are derived via Chebyshev T_p recursion and half-angle reduction.
//
// Reference: Gauss, Disquisitiones Arithmeticae §VII (1801).

#include "simplify_impl.hpp"

namespace cas::symbolic::detail {

// ── Fermat-prime base angle constructors (p = 5) ─────────────────────────────

// cos(π/5) = (1 + √5) / 4
[[nodiscard]] ExprPtr cos_pi_over_5(AstArena& arena);

// cos(2π/5) = (√5 − 1) / 4
[[nodiscard]] ExprPtr cos_2pi_over_5(AstArena& arena);

// sin(π/5) = √(10 − 2√5) / 4
[[nodiscard]] ExprPtr sin_pi_over_5(AstArena& arena);

// sin(2π/5) = √(10 + 2√5) / 4
[[nodiscard]] ExprPtr sin_2pi_over_5(AstArena& arena);

// ── Reference tables for sin / cos at rational multiples r ∈ [0, 1/2] of π ──

// Table lookup for sin(r·π): returns non-null for r ∈ {0, 1/10, 1/6, 1/5,
// 1/4, 3/10, 1/3, 2/5, 1/2}.  Returns nullptr for any other r.
[[nodiscard]] ExprPtr sin_ref_value_table(Rational ref, AstArena& arena);

// Table lookup for cos(r·π): returns non-null for r ∈ {0, 1/10, 1/6, 1/5,
// 1/4, 3/10, 1/3, 2/5, 1/2}.  Returns nullptr for any other r.
[[nodiscard]] ExprPtr cos_ref_value_table(Rational ref, AstArena& arena);

} // namespace cas::symbolic::detail
