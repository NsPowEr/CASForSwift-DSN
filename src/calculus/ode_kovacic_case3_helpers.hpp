#pragma once
// ode_kovacic_case3_helpers.hpp — shared between ode_kovacic_case3.cpp and
// ode_kovacic_case3_helpers.cpp.  Splits Step 3 (polynomial P search via
// Kovacic recurrence P_n, ..., P_{-1}) and Step 4 (minimal polynomial of ω
// build) out of the main TU to respect anti-monolith (≤ 500 LOC per file).
//
// Reference: Kovacic J.J. (1986), §5.

#include "ode_kovacic_internal.hpp"
#include <optional>

namespace cas::calculus::kovacic_impl {

// Step 3: search a monic polynomial P ∈ C[x] of degree d such that, with
//   P_n      = −P
//   P_{i−1}  = −S·P'_i + ((n−i)·S' − S·θ)·P_i − (n−i)·(i+1)·S²·r·P_{i+1}
// (running from i = n down to i = 0), the polynomial P_{−1} is identically 0.
//
// On success the returned P is the desired monic poly (with the chosen
// coefficients substituted).  ok(nullopt) means no such P exists for the
// chosen (θ, S, r, d, n) tuple.
[[nodiscard]] Result<std::optional<ExprPtr>> search_polynomial_P_case3(
    ExprPtr theta, ExprPtr S, ExprPtr r,
    long long d, unsigned n,
    const Symbol& x, symbolic::CASContext& ctx);

// Step 4: given (θ, S, P, r, n), build the n-th-degree minimal polynomial of
// ω = η'/η as the polynomial in `omega_var`:
//   Σ_{i=0}^n  (S^i · P_i(x) / (n − i)!) · ω^i  =  0
// where P_i are the recurrence quantities of Step 3 (recomputed with the
// substituted P).  Returns the polynomial expression in `omega_var`.
[[nodiscard]] Result<ExprPtr> build_omega_minpoly_case3(
    ExprPtr theta, ExprPtr S, ExprPtr P, ExprPtr r, unsigned n,
    const Symbol& x, const Symbol& omega_var, symbolic::CASContext& ctx);

} // namespace cas::calculus::kovacic_impl
