#pragma once
// F5.7 — Zeilberger creative telescoping for bivariate hypergeometric sums.
//
// Given F(n,k) proper hypergeometric (F vanishes outside the summation range
// in k), finds a linear recurrence Σ_{i=0}^J p_i(n)·S(n+i) = 0 for the
// definite sum S(n) = Σ_{k=lower}^{upper} F(n,k), then solves it for S(n).
//
// Algorithm: parametric Gosper (Zeilberger 1990) on the operator ansatz
//   Σ p_i(n)·F(n+i,k) = G(n,k+1) − G(n,k),  G = R·F.
//
// Scope: J ≤ ctx.max_zeilberger_order(), deg(p_i) ≤ ctx.max_zeilberger_poly_degree().
// Cases outside this scope produce explicit Unimplemented (never silent failure).

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include <optional>

namespace cas::symbolic {

// Returns the closed-form S(n) for Σ_{k=lower}^{upper} F(n,k),
// or nullopt if F is not proper hypergeometric or the recurrence finder fails.
[[nodiscard]] Result<std::optional<ExprPtr>> zeilberger_sum(
    ExprPtr F,
    const Symbol& n_param,
    const Symbol& k,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx);

}  // namespace cas::symbolic
