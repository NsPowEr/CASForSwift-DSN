// A6 — Exact p-adic splitting engine, Brick 3a of the Stauduhar deg ≥ 8
// closure (float-free Fieker-Klüners descent).
//
// Represents the unramified Galois ring GR(p^k, L) = (Z/p^k)[t]/(Φ) with Φ
// monic of degree L and irreducible mod p, and lifts the n roots of a monic
// squarefree-mod-p integer polynomial f into it by Hensel/Newton iteration
// (quadratic convergence; f'(root) is a unit because p ∤ disc f). All
// arithmetic is exact BigInt — no float enters at any point (REGOLA 1 and
// the approved float-free course correction over the MPFR route of
// Galois_Groups.md; spec read per REGOLA 0.1).
//
// Every lifted splitting is CERTIFIED before being returned (InternalError
// on violation, never a silent wrong result):
//   • f(r_i) ≡ 0 mod p^k for every root;
//   • the r_i are pairwise distinct mod p (squarefreeness transported);
//   • Newton's identities: e_1 and e_n of the roots match −f_{n−1} and
//     (−1)^n f_0 mod p^k;
//   • the Frobenius permutation (matching σ(r_i) among the roots, σ the
//     canonical lift of x ↦ x^p) has cycle type equal to the factor-degree
//     multiset of f mod p (Dedekind).
//
// The randomized equal-degree splitting (Cantor-Zassenhaus) draws from
// ctx->rng() when available, else from an input-derived seed — never a
// fixed literal seed (REGOLA hardcode cat. 6); it is a Las Vegas step:
// randomness affects time only, never the (verified) output.

#pragma once

#include "algebraic_tower_resultant.hpp"  // primitive_internal::Deadline
#include "perm_group_internal.hpp"        // Perm
#include "polynomial_internal.hpp"        // IntPoly

#include "cas/bigint.hpp"
#include "cas/result.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace cas::symbolic {
class CASContext;
}

namespace cas::algebra::galois_padic {

// Element of GR(p^k, L): coefficient vector of length L on the power basis
// 1, t, …, t^{L−1}, each coefficient reduced into [0, p^k).
using RingElem = std::vector<BigInt>;

// The ring GR(p^k, L) = (Z/p^k)[t]/(Φ). Φ is fixed once (an integer monic
// polynomial irreducible mod p), so the same Φ presents the ring at every
// precision k — raising precision never changes coordinates already valid.
class PadicRing {
public:
    // Requires: phi monic with integer coefficients, deg phi ≥ 1, k ≥ 1.
    // (Irreducibility mod p is the caller's certified responsibility.)
    [[nodiscard]] static Result<PadicRing> make(BigInt p, std::size_t k,
                                                IntPoly phi);

    [[nodiscard]] const BigInt& prime() const noexcept { return p_; }
    [[nodiscard]] std::size_t precision() const noexcept { return k_; }
    [[nodiscard]] const BigInt& modulus() const noexcept { return pk_; }
    [[nodiscard]] const IntPoly& phi() const noexcept { return phi_; }
    [[nodiscard]] std::size_t ext_degree() const noexcept { return ext_; }

    [[nodiscard]] RingElem zero() const;
    [[nodiscard]] RingElem one() const;
    [[nodiscard]] RingElem from_int(const BigInt& v) const;
    // t^j for j < L (a power-basis vector).
    [[nodiscard]] RingElem basis_power(std::size_t j) const;

    [[nodiscard]] RingElem add(const RingElem& a, const RingElem& b) const;
    [[nodiscard]] RingElem sub(const RingElem& a, const RingElem& b) const;
    [[nodiscard]] RingElem neg(const RingElem& a) const;
    [[nodiscard]] RingElem mul(const RingElem& a, const RingElem& b) const;
    // a^e for e ≥ 0 by binary powering.
    [[nodiscard]] RingElem pow(const RingElem& a, const BigInt& e) const;

    [[nodiscard]] bool is_zero(const RingElem& a) const;
    [[nodiscard]] bool equal(const RingElem& a, const RingElem& b) const;
    // a mod p == b mod p (agreement at precision 1).
    [[nodiscard]] bool equal_mod_p(const RingElem& a, const RingElem& b) const;

    // Multiplicative inverse; error iff a ≡ 0 mod p is not a unit (i.e. the
    // reduction of a mod p is zero — for k = 1, L prime field this is the
    // usual field inverse). Newton-lifted from the mod-p inverse.
    [[nodiscard]] Result<RingElem> inv(const RingElem& a) const;

    // f(x) for integer f by Horner (coefficients reduced into the ring).
    [[nodiscard]] RingElem eval_int_poly(const IntPoly& f,
                                         const RingElem& x) const;

    // The symmetric integer residue of a, iff a lies in the prime subring
    // Z/p^k (all t-coordinates ≥ 1 vanish): value in (−p^k/2, p^k/2].
    // nullopt certifies a ∉ Z (an integer embeds with pure coordinates).
    [[nodiscard]] std::optional<BigInt> integer_residue(const RingElem& a)
        const;

    // Same ring at precision k2 (identical p and Φ).
    [[nodiscard]] Result<PadicRing> with_precision(std::size_t k2) const;

private:
    PadicRing() = default;
    BigInt p_;
    std::size_t k_{0U};
    BigInt pk_;
    IntPoly phi_;
    std::size_t ext_{0U};
};

// Dense polynomial over the ring; index = degree on y.
using RingPoly = std::vector<RingElem>;

// Ring-polynomial helpers shared by root finding and the Stauduhar
// resolvent (Brick 3c). Field-only operations state their k = 1 contract.
void rp_normalize(const PadicRing& R, RingPoly& a);
[[nodiscard]] RingPoly rp_mul(const PadicRing& R, const RingPoly& a,
                              const RingPoly& b);
// Remainder of a modulo monic-leading b (leading coeff must be a unit).
[[nodiscard]] Result<RingPoly> rp_rem(const PadicRing& R, RingPoly a,
                                      const RingPoly& b);
// Monic gcd over the residue FIELD — requires k == 1.
[[nodiscard]] Result<RingPoly> rp_gcd_monic(const PadicRing& R, RingPoly a,
                                            RingPoly b);
// base^e mod m (binary powering over the ring).
[[nodiscard]] Result<RingPoly> rp_powmod(const PadicRing& R, RingPoly base,
                                         const BigInt& e, const RingPoly& m);
[[nodiscard]] RingElem rp_eval(const PadicRing& R, const RingPoly& a,
                               const RingElem& x);

// A prime p with f squarefree mod p (⇔ p ∤ disc f), together with the monic
// irreducible factors of f mod p and L = lcm of their degrees. Chosen among
// the first `prime_budget` unramified primes to minimise L (tie → smaller p):
// L is the extension degree every later ring operation pays for (O(L²)).
struct SplittingPrime {
    BigInt p;
    std::vector<IntPoly> factors;  // monic, irreducible mod p
    std::size_t ext_degree{0U};    // L = lcm of the factor degrees
};
[[nodiscard]] Result<SplittingPrime> choose_splitting_prime(
    const IntPoly& f_monic, std::size_t prime_budget,
    symbolic::CASContext* ctx,
    const primitive_internal::Deadline& deadline = std::nullopt);

// The n roots of f lifted into GR(p^k, L), with the Frobenius action.
struct PadicSplitting {
    PadicRing ring;
    IntPoly f;                    // the monic integer model
    std::vector<RingElem> roots;  // n roots, pairwise distinct mod p
    permgrp::Perm frobenius;      // i ↦ j with σ(r_i) = r_j
};

// Builds the certified splitting at precision k (see file header for the
// certificate list). `sp` must come from choose_splitting_prime for f.
[[nodiscard]] Result<PadicSplitting> build_padic_splitting(
    const IntPoly& f_monic, const SplittingPrime& sp, std::size_t precision_k,
    symbolic::CASContext* ctx,
    const primitive_internal::Deadline& deadline = std::nullopt);

// The same splitting (same root ORDER — descent state depends on it) at
// higher precision k2, re-certified.
[[nodiscard]] Result<PadicSplitting> raise_splitting_precision(
    const PadicSplitting& s, std::size_t k2, symbolic::CASContext* ctx,
    const primitive_internal::Deadline& deadline = std::nullopt);

}  // namespace cas::algebra::galois_padic
