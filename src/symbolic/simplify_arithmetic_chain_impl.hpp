#pragma once
// Internal helpers shared between simplify_arithmetic_chain_liketerm.cpp,
// simplify_arithmetic_chain_sum.cpp, and simplify_arithmetic_chain.cpp.
// NOT part of the public API.

#include "simplify_impl.hpp"

namespace cas::symbolic::detail {

// F1.4 like-term helpers declared here so they can be defined in
// simplify_arithmetic_chain_liketerm.cpp and called from the sum/product files.

/// Decompose expr into (numeric_coeff, sorted symbolic factor list).
/// Returns false only if parsing is structurally impossible (never happens
/// for valid AST nodes — the fallback treats the whole expr as one factor).
[[nodiscard]] bool decompose_term(
    ExprPtr expr,
    Rational& coeff_out,
    std::vector<std::pair<ExprPtr, BigInt>>& factors_out);

/// Build coeff_expr * numeric_coeff * ∏mono_factors.
[[nodiscard]] ExprPtr build_coeff_monomial(
    ExprPtr coeff_expr,
    Rational numeric_coeff,
    const std::vector<std::pair<ExprPtr, BigInt>>& mono_factors,
    AstArena& arena);

/// One pass of pairwise symbolic like-term merging on a flat term list.
/// Returns true if any merge occurred (caller should loop to fixed-point).
[[nodiscard]] bool try_merge_symbolic_like_terms(
    std::vector<ExprPtr>& terms,
    AstArena& arena);

} // namespace cas::symbolic::detail
