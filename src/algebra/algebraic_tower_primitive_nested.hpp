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

struct ResolvedGen {
    ExprPtr canon;
    const RootOf* node{nullptr};
    AlgebraicNumber::CoeffVec mp;
};

// Attempt to compute α's absolute monic min-poly over Q given:
//   outer        — RootOf whose polynomial f(x) ∈ Q(β₁, ..., β_k)[x]
//   outer_canon  — canonical ExprPtr for outer (for vanishing test)
//   resolved_gens — list of already resolved generators in the tower
[[nodiscard]] Result<std::optional<AlgebraicNumber::CoeffVec>> try_nested_lift_min_poly_multi(
    const RootOf& outer,
    ExprPtr outer_canon,
    const std::vector<ResolvedGen>& resolved_gens,
    symbolic::CASContext& ctx,
    const Deadline& deadline);

// Legacy 1-generator overload
[[nodiscard]] Result<std::optional<AlgebraicNumber::CoeffVec>> try_nested_lift_min_poly(
    const RootOf& outer,
    ExprPtr beta_canon,
    const AlgebraicNumber::CoeffVec& beta_min_poly,
    symbolic::CASContext& ctx,
    const Deadline& deadline);

}  // namespace primitive_internal
}  // namespace algebra
}  // namespace cas

