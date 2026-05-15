#pragma once

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

namespace cas::calculus {

// Compute the real definite integral
//
//     ∫_{-∞}^{+∞}  P(x) / Q(x)  dx
//
// for a rational expression in `var` whose denominator factors over Q[x]
// into pairwise-coprime irreducible factors of degree 1 (real roots) or
// degree 2 (conjugate complex pair).
//
// Returns success only when:
//   * apart_num_den yields a quotient N/D with deg(D) >= deg(N) + 2,
//   * D has no real roots (every linear factor would be a real pole),
//   * Every irreducible factor of D is at most quadratic.
//
// Otherwise returns CASErrorKind::Unimplemented with a diagnostic message.
//
// The computation uses the standard semicircle-contour residue theorem:
//   I = 2π i  Σ_{α in upper half plane}  Res(N/D, α)
// Each conjugate pair contributes once; for an irreducible
// quadratic factor q(x) = x² + b x + c with discriminant Δ = b² - 4c < 0
// and upper-half root α = -b/2 + i √(-Δ)/2, the residue is expressed in
// Q(α) as r = e + f α, and the contribution to the real integral is
//   2π i · r  →  real part = -π · f · √(-Δ).
[[nodiscard]] Result<ExprPtr> integrate_rational_full_real_line(
    ExprPtr rational_expr,
    const Symbol& var,
    symbolic::CASContext& ctx);

}  // namespace cas::calculus
