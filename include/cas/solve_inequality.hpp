// F6.4 — Univariate polynomial inequality solver via Sturm sequence.
//
// Solves p(x) ⊳ 0  with p ∈ Q[x] univariate and ⊳ ∈ {<, ≤, >, ≥}.
// Strategy:
//   1. Isolate every real root of p via Sturm bisection on Q (exact bounds).
//   2. Between consecutive roots evaluate sign(p) at the rational midpoint.
//   3. Emit closed intervals whose sign matches the requested operator,
//      tagging the inclusivity of each endpoint (open vs. closed depending
//      on whether the comparison is strict and whether the endpoint is an
//      actual root).
//
// Domain: real univariate polynomial inequalities only.  Multivariate or
// transcendental inequalities return Unimplemented diagnostic; the caller
// is expected to dispatch upstream (e.g. periodic reduction for sin/cos).
//
// References:
//   - Sturm, "Mémoire sur la résolution des équations numériques", 1829.
//   - Cohen, "A course in computational algebraic number theory", §4.1.
//
// Reference impl: F6.4-T1 (HP-Prime parity, Sturm-based, no CAD).

#pragma once

#include "cas/expr.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <vector>

namespace cas::algebra {

enum class InequalityOp {
    Greater,        // p > 0
    GreaterEqual,   // p ≥ 0
    Less,           // p < 0
    LessEqual,      // p ≤ 0
};

// A half-line interval in Q ∪ {-∞, +∞}.  `lower == std::nullopt`  means -∞;
// `upper == std::nullopt` means +∞.  `*_open` flags carry strictness:
//   - true  → boundary excluded (open),
//   - false → boundary included (closed).
struct InequalityInterval {
    std::optional<Rational> lower;
    std::optional<Rational> upper;
    bool lower_open;
    bool upper_open;
};

// Solve `poly OP 0` where `poly ∈ Q[var]`.
//
// Returns the disjoint, sign-correct list of intervals.  The list is sorted
// in ascending order of `lower` (with -∞ first); intervals are guaranteed
// disjoint and non-empty.
//
// Errors:
//   - InvalidArgument if `poly` is not a Q[var] polynomial.
//   - Unimplemented if `poly` contains symbols outside `var`.
[[nodiscard]] Result<std::vector<InequalityInterval>> solve_inequality_1var(
    ExprPtr poly,
    const Symbol& var,
    InequalityOp op,
    symbolic::CASContext& ctx);

}  // namespace cas::algebra
