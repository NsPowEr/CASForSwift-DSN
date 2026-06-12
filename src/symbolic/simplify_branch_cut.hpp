#pragma once
// F8.0-6.2 / Task 20 (BC-1..BC-3) — Branch-cut propagation helpers.
//
// Centralised construction of unwinding-number corrections required by
// the Branch_Cut_Propagation.md spec when ctx.strict_branch_cuts() is on
// and the operand domain is not provably real-positive.
//
// All helpers return new arena-interned AST nodes; they do NOT mutate.
// Callers are expected to be the simplifier rules in simplify_exp_log.cpp
// and simplify_arithmetic_power.cpp.
//
// Reference:
//   Kahan W. (1987) "Branch Cuts for Complex Elementary Functions".
//   Corless R.M., Davenport J.H., Jeffrey D.J. (2000) "Unwinding the
//   Branches of the Lambert W Function", SIGSAM Bulletin 34(2), 2–8.

#include "cas/ast.hpp"

namespace cas::symbolic::branch_cut {

// e^(2πi · b · K(a · ln(z)))
//
// Power-of-power correction (Branch_Cut_Propagation.md §2 rule 3):
//   (z^a)^b = z^(a·b) · e^(2πi · b · K(a · ln(z)))
//
// `arena` provides the interning context. The returned expression is
// arena-owned; callers must keep the arena alive for its lifetime.
[[nodiscard]] ExprPtr make_pow_of_pow_correction(
    ExprPtr z,
    ExprPtr a,
    ExprPtr b,
    AstArena& arena);

// -2πi · K(ln(z1) + ln(z2))
//
// Log-of-product correction (Branch_Cut_Propagation.md §2 rule 4):
//   ln(z1 · z2) = ln(z1) + ln(z2) - 2πi · K(ln(z1) + ln(z2))
//
// Returns the additive correction term (negative on the unwinding
// contribution) suitable for appending to ln(z1) + ln(z2).
[[nodiscard]] ExprPtr make_log_product_correction(
    ExprPtr z1,
    ExprPtr z2,
    AstArena& arena);

// -2πi · K(ln(z1) - ln(z2))
//
// Log-of-quotient correction (Branch_Cut_Propagation.md §2 rule 5):
//   ln(z1 / z2) = ln(z1) - ln(z2) - 2πi · K(ln(z1) - ln(z2))
[[nodiscard]] ExprPtr make_log_quotient_correction(
    ExprPtr z1,
    ExprPtr z2,
    AstArena& arena);

// (-1)^K(2·ln(z))
//
// Sqrt-of-square correction (Branch_Cut_Propagation.md §2 rule 1):
//   sqrt(z^2) = z · (-1)^K(2·ln(z))
//
// Returned as Pow(-1, K(2·ln(z))). When K(2·ln(z)) reduces to a known
// integer (e.g. via assumptions on z's argument range) the simplifier
// will collapse Pow(-1, integer) downstream.
[[nodiscard]] ExprPtr make_sqrt_of_square_correction(
    ExprPtr z,
    AstArena& arena);

}  // namespace cas::symbolic::branch_cut
