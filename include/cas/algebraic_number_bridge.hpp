#pragma once

#include "cas/algebraic_number.hpp"
#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <optional>

namespace cas {
namespace algebra {

// Bridge layer between AST-level RootOf nodes and the field-arithmetic
// AlgebraicNumber class.  Pure algebra; no transcendental handling.
//
// Conventions:
//   - AlgebraicNumber::min_poly() is monic over Q, ascending degree order.
//   - AlgebraicNumber::value() is the coefficient vector of p(x), deg(p) < deg(m),
//     where alpha is the (abstract) generator of Q(alpha) defined by m(x).
//   - The "physical" generator expression (e.g. a RootOf or sqrt node) is supplied
//     separately whenever conversion to/from ExprPtr is needed.

// Extract the canonical monic rational minimal polynomial from a RootOf node.
// Fails if the polynomial is non-rational, constant, or has zero leading coeff.
[[nodiscard]] Result<AlgebraicNumber::CoeffVec> rootof_min_poly(
    const RootOf& root,
    symbolic::CASContext& ctx);

// Construct the AlgebraicNumber alpha itself (value = x) from a RootOf node.
[[nodiscard]] Result<AlgebraicNumber> alpha_from_rootof(
    const RootOf& root,
    symbolic::CASContext& ctx);

// Try to express the ExprPtr `e` as an element of Q(alpha), where the
// generator is identified structurally with `alpha_expr` and the minimal
// polynomial is given by `min_poly` (monic, ascending degree).
//
// Returns the corresponding AlgebraicNumber when `e` is built from rationals,
// the literal `alpha_expr`, and the operations +, -, *, /, integer powers
// (negative power = inversion in Q(alpha)).  Returns std::nullopt when the
// expression contains anything else (symbols, other RootOfs, transcendentals,
// matrices, ...).
//
// Returns an error only on internal arithmetic failure (e.g. division by zero
// or inversion of a non-invertible element when the supplied minimal polynomial
// is not irreducible).
[[nodiscard]] Result<std::optional<AlgebraicNumber>> try_express_in_q_alpha(
    ExprPtr e,
    ExprPtr alpha_expr,
    const AlgebraicNumber::CoeffVec& min_poly,
    symbolic::CASContext& ctx);

// Render an AlgebraicNumber back to canonical ExprPtr form:
//   c_0 + c_1 * alpha + c_2 * alpha^2 + ...
// using `alpha_expr` as the rendered generator and the supplied context's
// arena.  The result is passed through context.simplify().
[[nodiscard]] Result<ExprPtr> algebraic_number_to_expr(
    const AlgebraicNumber& value,
    ExprPtr alpha_expr,
    symbolic::CASContext& ctx);

// Same as algebraic_number_to_expr but without invoking the simplifier.
// Returns the raw Sum of c_k * alpha^k.  Intended for callers that drive
// further simplification themselves (or that need to avoid recursive
// simplify calls inside an active simplification pass).
[[nodiscard]] ExprPtr algebraic_number_to_expr_raw(
    const AlgebraicNumber& value,
    ExprPtr alpha_expr,
    AstArena& arena);

// Apply Q(alpha) reduction on `expr`:
//   1. Detect the unique RootOf appearing in `expr` (if any) with rational
//      minimal polynomial.
//   2. Attempt to express the whole expression in Q(alpha).
//   3. If successful, re-render the canonical Q(alpha) form and simplify.
// Returns the original (already-simplified) expression when no productive
// reduction is possible.  Idempotent in practice.
//
// The caller is responsible for first running ctx.simplify(); this helper
// only performs the algebraic-extension layer.
[[nodiscard]] Result<ExprPtr> try_reduce_in_q_alpha(
    ExprPtr expr,
    symbolic::CASContext& ctx);

// Composite entry point: runs ctx.simplify first, then try_reduce_in_q_alpha,
// then ctx.simplify again if the reduction was productive.  Bounded to one
// extra simplify round to avoid pathological feedback.
[[nodiscard]] Result<ExprPtr> simplify_in_q_alpha(
    ExprPtr expr,
    symbolic::CASContext& ctx);

// Reduce an expression that is a polynomial in `poly_var` whose coefficients
// live in Q(alpha) (where alpha is the unique RootOf appearing in the
// expression with rational minimal polynomial).  For each coefficient of
// (poly_var)^k, the function attempts try_reduce_in_q_alpha; the polynomial
// is then rebuilt with the canonicalised coefficients.  Falls back to the
// input expression if the polynomial-in-x parse fails or no extension is
// detected.
[[nodiscard]] Result<ExprPtr> simplify_polynomial_in_x_over_q_alpha(
    ExprPtr expr,
    const Symbol& poly_var,
    symbolic::CASContext& ctx);

}  // namespace algebra
}  // namespace cas
