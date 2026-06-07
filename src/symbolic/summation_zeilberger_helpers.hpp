#pragma once
// Internal helpers shared between summation_zeilberger.cpp and
// summation_zeilberger_helpers.cpp.  Not part of the public CAS API.
//
// Provides:
//   - expand_gamma_int_shifts(e, ctx): rewrites Γ(z+n) → (∏(z+i))·Γ(z) for
//     integer n, enabling the simplifier to cancel matching gamma factors.
//   - cancel_common_factors_in_ratio(r, ctx): explicit structural cancellation
//     of common factors between numerator/denominator (needed because the
//     simplifier defers Γ(z+n) reduction inside Product contexts).
//   - compute_shift_ratio(F, sym, ctx): F(sym→sym+1)/F, simplified with both
//     expansion + cancellation steps applied.

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include <optional>

namespace cas::symbolic::zeilberger_detail {

[[nodiscard]] ExprPtr expand_gamma_int_shifts(
    ExprPtr e, symbolic::CASContext& ctx);

[[nodiscard]] ExprPtr cancel_common_factors_in_ratio(
    ExprPtr ratio, symbolic::CASContext& ctx);

[[nodiscard]] // Compute F(sym+delta)/F(sym) as a rational function in sym.
// Returns {Num, Den} pair where Num/Den is the ratio.
std::optional<std::pair<ExprPtr, ExprPtr>> compute_shift_ratio(
    ExprPtr F, const Symbol& sym, symbolic::CASContext& ctx, long long delta = 1);

[[nodiscard]] bool struct_equal(ExprPtr a, ExprPtr b);

}  // namespace cas::symbolic::zeilberger_detail
