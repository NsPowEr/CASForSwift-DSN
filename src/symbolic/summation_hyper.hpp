#pragma once

// F5.7 — Petkovšek's Hyper algorithm (rational characteristic root z ∈ ℚ).
//
// Given a homogeneous linear recurrence with polynomial coefficients
//   Σ_{i=0}^J p[i](n)·y(n+i) = 0          (p[0] ≠ 0, p[J] ≠ 0),
// find the term ratio ρ(n) = y(n+1)/y(n) ∈ ℚ(n) of a hypergeometric solution.
//
// Restricted to a rational characteristic root z (Petkovšek 1992, "A=B" ch.8).
// Candidates with an algebraic (non-rational) z are skipped; if every solution
// would need an algebraic z the result is ok(nullopt) with `needs_algebraic`
// set, so the caller can emit a precise Unimplemented diagnostic instead of a
// silent "no solution".
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Zeilberger_Algorithm.md §A.

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <vector>

namespace cas::symbolic {

struct HyperResult {
    // Term ratio ρ(n) = y(n+1)/y(n) of a found hypergeometric solution, if any.
    std::optional<ExprPtr> ratio;
    // True when a candidate with an algebraic (non-rational) characteristic
    // root z was encountered but skipped — a hypergeometric solution may still
    // exist beyond the rational-z restriction.
    bool needs_algebraic{false};
};

// Returns a verified hypergeometric term ratio for the recurrence, or an empty
// ratio when none exists within the rational-z / configured-degree search.
[[nodiscard]] Result<HyperResult> hyper_term_ratio(
    const std::vector<ExprPtr>& p,
    const Symbol& n,
    symbolic::CASContext& ctx);

}  // namespace cas::symbolic
