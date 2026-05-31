#pragma once
// Internal header: Chebyshev trig generator.
// Included only by simplify_trig.cpp and simplify_trig_chebyshev.cpp.
// NOT part of the public API.
//
// Algorithm reference:
//   Gauss, Disquisitiones Arithmeticae §VII (1801): constructibility of
//   regular n-gons ↔ n = 2^k · ∏ distinct Fermat primes.
//   The minimal polynomial of 2cos(2π/n) over Q is the n-th "cosine minimal
//   polynomial" Ψ_n(y), obtained from Φ_n(x) by the substitution y = x+1/x
//   (Lang, Algebra §IV.6; also TAOCP App. B §B5 "algebraic numbers").

#include "simplify_impl.hpp"
#include "../algebra/polynomial_internal.hpp"

namespace cas::symbolic::detail {

// Minimal polynomial of 2cos(π/q) over Q.
// = Ψ_{2q}(y) obtained from Φ_{2q}(x) via y = x + x^{-1}.
// Degree = φ(2q)/2 = φ(q) for odd q; φ(2q)/2 for even q.
// Returns empty IntPoly if q ≤ 0 or too large (> kCosPolyMaxQ).
[[nodiscard]] algebra::IntPoly min_poly_2cos_pi_q(int q);

// Build RootOf(Ψ_{2q}(y), _t) ÷ 2 = cos(π/q) as an ExprPtr.
// root_index is the index (0-based) of the real root closest to cos(π/q)
// in the sorted real root list — for q ≥ 2 this is always index 0
// (largest positive real root in [0,1]).
[[nodiscard]] ExprPtr build_rootof_cos_pi_q(int q, AstArena& arena);

// Maximum denominator supported by min_poly_2cos_pi_q.
constexpr int kCosPolyMaxQ = 500;

} // namespace cas::symbolic::detail
