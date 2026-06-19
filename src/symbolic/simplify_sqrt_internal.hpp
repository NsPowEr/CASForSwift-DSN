#pragma once
// W9.3 split: sqrt rational/denesting helpers shared between simplify_sqrt.cpp
// (the simplifier sqrt branch) and simplify_sqrt_helpers.cpp (their definitions).
// Extracted from the former monolithic simplify_exp_log.cpp (854 LOC) so each
// translation unit stays under the 500-line anti-monolith limit.

#include "simplify_impl.hpp"

#include <cstddef>
#include <optional>
#include <utility>

namespace cas::symbolic::detail {

// Newton integer sqrt: floor(sqrt(n)) for n >= 0.
[[nodiscard]] BigInt integer_sqrt(const BigInt& n);

// Extract perfect-square factor from n: returns {k, m} with n = k²·m.
[[nodiscard]] std::pair<BigInt, BigInt> extract_square_factor(BigInt n, std::size_t trial_bound);

// sqrt(r) for rational r >= 0: returns k·sqrt(m), k rational, m squarefree.
[[nodiscard]] Result<ExprPtr> simplify_rational_sqrt(const Rational& r, AstArena& arena, std::size_t trial_bound);

// Try to extract a rational sqrt: if r = (p/q)² returns p/q, else nullopt.
[[nodiscard]] std::optional<Rational> try_rational_sqrt(const Rational& r);

// Borodin-Fagin-Hopcroft-Tompa (1985) denesting of sqrt(a + b·sqrt(c)).
[[nodiscard]] std::optional<ExprPtr> try_denest_borodin_fagin(ExprPtr radicand, AstArena& arena);

}  // namespace cas::symbolic::detail
