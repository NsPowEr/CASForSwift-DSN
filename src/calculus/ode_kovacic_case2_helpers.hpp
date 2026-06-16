#pragma once
// ode_kovacic_case2_helpers.hpp — shared between ode_kovacic_case2.cpp and
// ode_kovacic_case2_helpers.cpp.  Splits Step 3 (polynomial P search) and
// Step 4 (ω construction) out of the main TU to respect anti-monolith
// (≤ 500 LOC per file).
//
// Reference: Kovacic J.J. (1986), §4.

#include "ode_kovacic_internal.hpp"
#include <optional>

namespace cas::calculus::kovacic_impl {

// Step 3: search a monic polynomial P ∈ C[x] of degree d satisfying
//   P''' + 3θP'' + (3θ² + 3θ' - 4r)P' + (θ'' + 3θθ' + θ³ - 4rθ - 2r')P = 0
// via undetermined coefficients (linear system).
//
// Returns:
//   ok(P)            — polynomial found.
//   ok(nullopt)      — no polynomial of degree d satisfies the ODE for this (θ, r).
//   Unimplemented    — internal failure (derivative/simplify/csolve error,
//                      degree above ctx.kovacic_case2_max_poly_degree()).
[[nodiscard]] Result<std::optional<ExprPtr>> search_polynomial_P_case2(
    ExprPtr theta, ExprPtr r, long long d,
    const Symbol& x, symbolic::CASContext& ctx);

// Step 4: given θ and the polynomial P returned by search_polynomial_P_case2,
// build ω as the two roots of
//   ω² + φω + (½φ' + ½φ² - r) = 0,  with φ = θ + P'/P.
// Returns OmegaPair{ω₊, ω₋}.
[[nodiscard]] Result<OmegaPair> build_omega_from_phi_case2(
    ExprPtr theta, ExprPtr P, ExprPtr r,
    const Symbol& x, symbolic::CASContext& ctx);

} // namespace cas::calculus::kovacic_impl
