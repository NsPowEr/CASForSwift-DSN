// A6 Brick 3a — shared internals between galois_padic_roots.cpp (mod-p
// layer: Φ selection, residue-field roots) and galois_padic_split.cpp
// (Newton lifting, Frobenius, certification, public entry points).
// Not part of the module API — include galois_padic_internal.hpp for that.

#pragma once

#include "galois_padic_internal.hpp"

#include "cas/bigint.hpp"
#include "cas/result.hpp"

#include <cstddef>
#include <random>
#include <vector>

namespace cas::symbolic {
class CASContext;
}

namespace cas::algebra::galois_padic::detail {

// Interrupt + deadline poll for the long-running loops.
[[nodiscard]] Result<void> poll(symbolic::CASContext* ctx,
                                const primitive_internal::Deadline& dl);

// g reduced mod p and made monic over the residue field K (k == 1).
[[nodiscard]] Result<RingPoly> to_field_poly(const PadicRing& K,
                                             const IntPoly& g);

// A monic integer Φ of degree L irreducible mod p: a degree-L factor of
// f mod p when available, else the first hit of the exhaustive canonical
// sweep certified by Rabin (terminates: irreducibles of every degree
// exist over F_p and the sweep covers the whole finite space).
[[nodiscard]] Result<IntPoly> find_phi(
    const BigInt& p, const std::vector<IntPoly>& factors, std::size_t L,
    symbolic::CASContext* ctx, const primitive_internal::Deadline& dl);

// All deg(g) roots of the monic irreducible-mod-p factor g inside
// K = GF(p^L): linear read-off / Frobenius orbit of t when g is the
// modulus / Cantor-Zassenhaus split (Las Vegas — randomness affects time
// only; every returned root is verified downstream by the certificates).
[[nodiscard]] Result<std::vector<RingElem>> roots_in_field(
    const PadicRing& K, const IntPoly& g, std::mt19937& rng,
    symbolic::CASContext* ctx, const primitive_internal::Deadline& dl);

[[nodiscard]] IntPoly derivative(const IntPoly& f);

// Hensel/Newton: lift x with g(x) ≡ 0 mod p^{m0}, g'(x) a unit, to
// precision k by doubling: x ← x − g(x)·g'(x)⁻¹ mod p^{min(2m,k)}.
[[nodiscard]] Result<RingElem> newton_lift(const IntPoly& g,
                                           const IntPoly& gprime,
                                           const PadicRing& base_ring,
                                           RingElem x, std::size_t m0,
                                           std::size_t k);

}  // namespace cas::algebra::galois_padic::detail
