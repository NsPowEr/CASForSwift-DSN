#pragma once

#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <string>
#include <vector>
#include <set>

namespace cas::calculus {

struct QuotientView {
    ExprPtr numerator;
    ExprPtr denominator;
};

/// @brief Differentiates transcendental functions.
/// Handles asin, acos, atan, sinh, cosh, tanh.
[[nodiscard]] Result<ExprPtr> differentiate_transcendental(
    BuiltinOp func_id,
    ExprPtr argument,
    const Symbol& var,
    symbolic::CASContext& context);

/// @brief Integration by parts: ∫ u dv = uv - ∫ v du
[[nodiscard]] Result<ExprPtr> integrate_by_parts(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& context);

/// @brief Partial Risch Algorithm for Logarithmic and Exponential extensions.
[[nodiscard]] Result<ExprPtr> integrate_risch(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& context);

/// @brief Weierstrass substitution t = tan(x/2) for ∫ R(sin x, cos x) dx.
[[nodiscard]] Result<ExprPtr> integrate_weierstrass(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& ctx);

[[nodiscard]] ExprPtr limit_make_integer(AstArena& arena, long long value);
[[nodiscard]] ExprPtr limit_make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs);
[[nodiscard]] bool limit_is_zero(ExprPtr expr);
[[nodiscard]] bool limit_is_one(ExprPtr expr);
[[nodiscard]] bool limit_is_infinity(ExprPtr expr);
[[nodiscard]] std::optional<QuotientView> extract_quotient_view(ExprPtr expr, AstArena& arena);
[[nodiscard]] std::optional<ExprPtr> cancel_common_linear_factor(
    ExprPtr numerator,
    ExprPtr denominator,
    const Symbol& var,
    ExprPtr point,
    AstArena& arena);
[[nodiscard]] std::optional<Result<ExprPtr>> try_polynomial_pole_limit(
    ExprPtr numerator,
    ExprPtr denominator,
    const Symbol& var,
    ExprPtr point,
    LimitDirection dir,
    AstArena& arena);
[[nodiscard]] std::optional<Result<ExprPtr>> try_logarithmic_root_limit(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    LimitDirection dir,
    AstArena& arena);
[[nodiscard]] Result<ExprPtr> try_infinite_limit(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    AstArena& arena);

[[nodiscard]] bool depends_on(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool is_bounded(ExprPtr expr, const Symbol& var);

// Tower height (Gruntz §3.5): maximum nesting of exp/log chain over sub-
// expressions that depend on `var`. Structural AST walk, O(N). Used by
// LimitEngine to derive an adaptive recursion bound from the input
// asymptotic complexity instead of a fixed depth cap.
//
//   tower_height(x)                = 1
//   tower_height(exp(f))           = tower_height(f) + 1 if f depends on var
//   tower_height(ln(f) / log(f))   = max(tower_height(f), 1)
//   tower_height(other functions)  = max over children
//
// Conservative (always ≤ true comparability tower height); a stricter bound
// would require running the MRV algorithm itself.
[[nodiscard]] unsigned int transcendental_tower_depth(ExprPtr expr, const Symbol& var);

[[nodiscard]] Result<ExprPtr> symbolic_sum(
    ExprPtr term,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx);

// L2-05 closure: Laurent expansion via Taylor‑series inversion of the
// denominator.  Splits expr on its top‑level Div node (no rational
// canonicalization), Taylor‑expands numerator and denominator at `center`
// to order positive_order + pole_budget, then divides the two series.
// Handles transcendental denominators with finite‑order zeros at the
// expansion point (e.g. 1/sin(x), 1/(x²·sin(x))) which the rational
// fast path cannot.
[[nodiscard]] Result<LaurentExpansion> laurent_series_general(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr center,
    unsigned int positive_order,
    unsigned int pole_budget,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> residue(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr pole,
    symbolic::CASContext& ctx);

// Hermite reduction on the base field (exact RatPoly-based)
[[nodiscard]] Result<HermiteReduction> hermite_reduction_exact(
    ExprPtr P, ExprPtr Q, const Symbol& var, symbolic::CASContext& ctx);

// MRV Algorithm (Gruntz)
struct MRVCompare {
    bool operator()(ExprPtr lhs, ExprPtr rhs) const noexcept;
};
using MRVSet = std::set<ExprPtr, MRVCompare>;

[[nodiscard]] MRVSet mrv_set(ExprPtr e, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> rewrite_mrv(ExprPtr e, const MRVSet& mrv, ExprPtr w, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> compute_limit_mrv(ExprPtr expr, const Symbol& var, ExprPtr point, symbolic::CASContext& ctx);

}  // namespace cas::calculus
