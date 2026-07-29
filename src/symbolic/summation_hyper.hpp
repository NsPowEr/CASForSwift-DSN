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

// As above but returns ALL distinct verified hypergeometric term ratios (each a
// rational function in n).  Distinct = not equal as rational functions; this is
// what a closed-form solver needs to span the solution space.
[[nodiscard]] Result<HyperResult> hyper_all_ratios(
    const std::vector<ExprPtr>& p,
    const Symbol& n,
    std::vector<ExprPtr>& ratios_out,
    symbolic::CASContext& ctx);

// Closed form for the sequence satisfying  Σ_{i=0}^J p[i](n)·S(n+i)=0  with the
// J initial values init[t] = S(n0+t), t=0..J−1.  Returns the closed form
// S(n) = Σ c_i h_i(n)  when S is a ℚ-linear combination of the hypergeometric
// solutions of the recurrence (each reconstructed as a Pochhammer/power term),
// verified by extending S through the recurrence beyond the fitted points.
// Returns ok(nullopt) when no such closed form exists — never a silent guess.
[[nodiscard]] Result<std::optional<ExprPtr>> solve_recurrence_closed_form(
    const std::vector<ExprPtr>& p,
    const std::vector<ExprPtr>& init,
    const Symbol& n,
    long long n0,
    symbolic::CASContext& ctx);

// Closed form for the definite sum S(n)=Σ_{k=lower}^{n} F(n,k) whose Zeilberger
// recurrence is `p` (order ≥ 1).  Computes the initial values by direct
// summation, solves the recurrence (solve_recurrence_closed_form), then
// CROSS-VERIFIES the candidate against directly-computed sums beyond the
// fitted range — guarding against telescoping boundary terms.  Requires an
// integer `lower`.  Returns ok(nullopt) when no verified closed form is found.
[[nodiscard]] Result<std::optional<ExprPtr>> sum_closed_form_from_recurrence(
    const std::vector<ExprPtr>& p,
    ExprPtr F,
    const Symbol& n,
    const Symbol& k,
    ExprPtr lower,
    symbolic::CASContext& ctx);

}  // namespace cas::symbolic
