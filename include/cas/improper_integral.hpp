// L2-11: Improper integral convergence classification and Cauchy principal value.
//
// This module provides a convergence diagnostic based on the leading
// Laurent (or asymptotic) order of a rational integrand at the endpoints
// and a Cauchy principal value computation for simple poles strictly
// inside the integration interval.
//
// Conventions:
//   * `leading_order_at_lower` and `leading_order_at_upper` follow the
//     `LaurentExpansion::leading_order` convention at a finite endpoint
//     (negative = pole order).
//   * At +/- infinity, leading_order is the asymptotic order of `expr` as
//     |x| -> infinity (so `1/(1+x^2)` reports -2, `1/x` reports -1,
//     `x` reports +1, etc.).
//   * Diagnostic message is human-readable and explains the verdict.
//
// Restrictions:
//   * Implementation handles rational expressions (consumed via
//     `algebra::apart_num_den` and `calculus::laurent_series`).
//   * Non-rational integrands return `Unknown` with a diagnostic message
//     rather than a hard error.

#pragma once

#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <string>

namespace cas::calculus {

enum class ConvergenceStatus {
    Convergent,
    DivergentAtLowerEnd,
    DivergentAtUpperEnd,
    DivergentAtInteriorPole,
    Unknown,
};

struct ConvergenceReport {
    ConvergenceStatus status;
    int leading_order_at_lower;
    int leading_order_at_upper;
    std::string diagnostic;
};

[[nodiscard]] Result<ConvergenceReport> classify_improper_convergence(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> cauchy_principal_value(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    ExprPtr pole,
    symbolic::CASContext& ctx);

}  // namespace cas::calculus
