#pragma once
// ode_kovacic_pf_helpers.hpp — partial-fraction extraction helpers shared
// between Kovacic Case 2 (§4) and Case 3 (§5) drivers.
//
// Anti-monolith rationale: case2 + case3 both need to (a) walk the output
// of `algebra::partial_fractions` and recover pole location + power + coeff,
// (b) compute  ord(r, ∞)  and the Laurent coefficient at ∞.  Centralised here
// to avoid copy-paste drift.

#include "ode_kovacic_internal.hpp"
#include <optional>
#include <vector>

namespace cas::calculus::kovacic_impl {

// ─── Generic numeric helpers ────────────────────────────────────────────────

[[nodiscard]] std::optional<Rational> try_get_rational(
    ExprPtr e, symbolic::CASContext& ctx);

[[nodiscard]] std::optional<long long> rational_to_int(const Rational& r);

// ─── PF pole extraction ─────────────────────────────────────────────────────

// Try to interpret a base node as (x - c)-style linear factor.  Returns the
// pole location c (or 0 for raw x) when the node matches; std::nullopt otherwise.
[[nodiscard]] std::optional<ExprPtr> extract_pole_loc(
    ExprPtr base, const Symbol& x, AstArena& a);

// Detect  Pow(base, n)  with n ≤ -1.  Returns (base, |n|) on success.
[[nodiscard]] std::optional<std::pair<ExprPtr, unsigned>> as_neg_pow(ExprPtr e);

// Walk the output of algebra::partial_fractions and collect (pole, power, coeff)
// triples.  Accepted term shapes:
//   (a)  Div(N, D)              — D = (x−c)^k or x^k
//   (b)  Pow(base, −k)          — k ≥ 1, base linear in x
//   (c)  Product([coeffs..., Pow(base, −k)])
[[nodiscard]] std::vector<PFPole> collect_pf_poles(
    const std::vector<ExprPtr>& pf_terms,
    const Symbol& x, symbolic::CASContext& ctx);

// ─── Infinity branch data ───────────────────────────────────────────────────

struct InfinityData {
    long long ord;            // ord(r, ∞) = deg(den) − deg(num)
    Rational  leading_b;      // Coefficient of x^{−2} in Laurent (only set when ord = 2)
    bool      leading_set;    // true ⇔ leading coefficient could be extracted as ℚ
};

[[nodiscard]] std::optional<InfinityData> compute_infinity_data(
    ExprPtr r, const Symbol& x, symbolic::CASContext& ctx);

} // namespace cas::calculus::kovacic_impl
