// A9 / F3.4-DEBT-01 — Primitive element of a sequentially nested tower.
// See algebraic_tower_primitive_chain.cpp for the mathematics and the
// certificate. Included only from algebraic_tower_primitive_nested.cpp.
#pragma once

#include "cas/algebraic_number.hpp"
#include "cas/algebraic_tower_bridge.hpp"
#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <vector>

namespace cas {
namespace algebra {
namespace primitive_internal {

// Fast path for a tower whose levels are each defined by a polynomial that is
// linear in the previous level (nested radicals, iterated square roots, …).
// Takes θ to be the generator of largest degree and expresses every other
// generator as a polynomial in θ, with no shift resultant.
//
// Returns nullopt — never a wrong answer — when the structure is not such a
// chain, when a level is not linear in its parent, or when the exact
// certificate m_i(rep_i) ≡ 0 (mod m_θ) fails for any generator. The caller
// then falls back to the generic Trager merge.
//
// `alphas`, `nodes` and `min_polys` are parallel: nodes[i] is the RootOf that
// alphas[i] canonicalises to, and min_polys[i] its absolute minimal polynomial
// over Q.
[[nodiscard]] Result<std::optional<PrimitiveElementResult>>
try_primitive_element_from_chain(
    const std::vector<ExprPtr>& alphas,
    const std::vector<const RootOf*>& nodes,
    const std::vector<AlgebraicNumber::CoeffVec>& min_polys,
    symbolic::CASContext& ctx);

}  // namespace primitive_internal
}  // namespace algebra
}  // namespace cas
