// algebraic_tower_primitive_nested.hpp — F3.4-DEBT-01 closure API.
// Helper: lift outer RootOf with min-poly in Q(β) → absolute Q-min-poly.
// Included only from algebraic_tower_primitive.cpp and
// algebraic_tower_primitive_nested.cpp.
#pragma once

#include "cas/algebraic_number.hpp"
#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include "algebraic_tower_primitive_internal.hpp"

#include <optional>

namespace cas {
namespace algebra {
namespace primitive_internal {

// Attempt to compute α's absolute monic min-poly over Q given:
//   outer        — RootOf whose polynomial f(x) ∈ Q(β)[x]
//   beta_canon   — canonical ExprPtr for β (used for structural substitution)
//   beta_min_poly — rational min-poly of β (ascending coeffs)
//
// Returns:
//   ok(nullopt)   — β not present in outer's polynomial, OR f(x, β) cannot be
//                   parsed as a bivariate polynomial in (x, y) with rational
//                   coefficients (e.g. depends on a third RootOf γ).
//   ok(coeffs)    — squarefree monic R(x) = Res_y(g(y), f(x, y)) ∈ Q[x].
//   fail(...)     — explicit Unimplemented when R is not squarefree
//                   (irreducible-factor selection over Q(β)[x] not implemented).
[[nodiscard]] Result<std::optional<AlgebraicNumber::CoeffVec>> try_nested_lift_min_poly(
    const RootOf& outer,
    ExprPtr beta_canon,
    const AlgebraicNumber::CoeffVec& beta_min_poly,
    symbolic::CASContext& ctx,
    const Deadline& deadline);

}  // namespace primitive_internal
}  // namespace algebra
}  // namespace cas
