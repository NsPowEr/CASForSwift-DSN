#pragma once

#include "cas/ode.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

namespace cas::calculus {

/// Kovacic's algorithm for solving 2nd order linear ODEs with rational coefficients.
/// Returns the general solution as an equation y = C1*y1 + C2*y2.
[[nodiscard]] Result<ExprPtr> solve_ode_kovacic(
    const OdeClassification& classification,
    symbolic::CASContext& ctx);

} // namespace cas::calculus
