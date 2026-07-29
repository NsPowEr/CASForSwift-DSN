// A6 Brick 3c — the certified Stauduhar descent step and the first-layer
// driver (float-free Fieker-Klüners route; spec Galois_Groups.md read per
// REGOLA 0.1, with the approved exact p-adic course correction over the
// spec's MPFR sketch).
//
// Mathematics of one step (G the current group with G_f ≤ G certified,
// H ≤ G a candidate with exact invariant F and transversal {σ_i} from
// Brick 3b, roots lifted in GR(p^k, L) from Brick 3a):
//
//   • v_i = (σ_i·F)(r_1..r_n). The resolvent R(y) = ∏_i (y − v_i) has
//     G-invariant integer-polynomial coefficients, hence — G_f ≤ G — its
//     coefficients are RATIONAL algebraic integers: exact members of Z.
//   • Archimedean certificate: every complex embedding of v_i is bounded
//     by B = #terms(F)·ρ^{deg F} with ρ = 1 + max|f_j| (Cauchy). So
//     |coeff_j(R)| ≤ C(m,j)·B^j ≤ (1+B)^m, and computing at p-adic
//     precision p^k > 2·(1+B)^m makes the symmetric lift of each
//     coefficient EXACT. Any impure coefficient or violated bound is an
//     InternalError (the theory guarantees purity — never silenced).
//   • v_i pairwise distinct mod p^k ⇒ pairwise distinct ⇒ R squarefree.
//   • If v_i has a pure residue c and R(c) = 0 exactly over Z, then c IS
//     v_i (distinct roots of the squarefree R cannot share a residue), so
//     by Stauduhar's theorem G_f ≤ σ_i·H·σ_i⁻¹: certified descent. If no
//     coset passes, G_f is contained in NO conjugate of H — also exact.
//   • Residue collisions (possible p-adic coincidence or a genuinely
//     non-squarefree resolvent) surface as Collision: the caller first
//     doubles the precision, then sweeps Tschirnhaus transforms on the
//     moment curve (β = P_t(α) evaluated on the SAME lifted roots — the
//     root order, hence the descent state, is preserved) with the sweep
//     length derived from the degree bound of the pairwise separating
//     polynomials, never a magic constant. Exhaustion → structured
//     Unimplemented (REGOLA ZERO: no guess).
//
// First layer (degrees 5..10): ambient = S_n, or A_n when disc(f) is a
// rational square; candidates from perm_maximal.cpp (Brick 2, coverage
// contract). If NO transitive maximal candidate contains G_f, then
// G_f = ambient EXACTLY (coverage + maximality). Below the first layer
// (Brick 3.5) the walk of stauduhar_identify repeats the certified step
// on the maximal transitive classes of each interior node — dense
// exhaustive (galois_sublattice) when the node fits the byte budget,
// structural wreath-preimage (galois_wreath_maximal, Brick 3.75, covers
// the 5|2-block degree-10 nodes of order 28800/14400) otherwise — until
// no class contains G_f, which certifies G_f = current node, exactly
// (HC-F8-PENDING-09 stays open for the Brick-4 naming and wiring).

#pragma once

#include "galois_invariant_internal.hpp"
#include "galois_padic_internal.hpp"
#include "perm_bsgs_internal.hpp"

#include "cas/result.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>

namespace cas::symbolic {
class CASContext;
}

namespace cas::algebra::galois_stauduhar {

enum class StepStatus : std::uint8_t {
    Descended,      // G_f ≤ conjugated_subgroup (certified)
    NotContained,   // G_f in no G-conjugate of H (certified)
    NeedPrecision,  // splitting precision below the derived bound
    Collision,      // residues collided at this precision — caller retries
};

struct StepOutcome {
    StepStatus status{StepStatus::Collision};
    // Set iff status == Descended.
    std::optional<permgrp::Perm> conjugator;
    std::optional<permgrp::BsgsGroup> conjugated_subgroup;
    // Set iff status == NeedPrecision: the derived bound the splitting
    // must reach (the step itself never mutates the splitting — the
    // caller owns the precision source, which under a Tschirnhaus model
    // is the ORIGINAL splitting, the only one Newton can safely refine).
    std::size_t required_precision{0U};
};

// One certified Stauduhar test of H below G. `split.f` must be the monic
// integer model whose roots index G and H.
[[nodiscard]] Result<StepOutcome> stauduhar_step(
    const permgrp::BsgsGroup& G, const permgrp::BsgsGroup& H,
    const galois_invariant::RelativeInvariant& inv,
    const galois_padic::PadicSplitting& split, symbolic::CASContext* ctx,
    const primitive_internal::Deadline& deadline = std::nullopt);

// A certified containment G_f ≤ subgroup (= σ·H·σ⁻¹ for the tested H).
struct DescentHit {
    permgrp::Perm conjugator;
    permgrp::BsgsGroup subgroup;
};

// Certified quadratic-character descent for an INDEX-2 candidate H of G
// (index 2 ⇒ H ◁ G ⇒ no conjugator): when the character of G/H is
// realised by a ±δ alternating form (per-block/top/global Vandermonde
// products over a minimal block system of G — see
// galois_stauduhar_quadratic.cpp), the resolvent is y² − v² with
// v = δ(r) evaluated as a PRODUCT: the |H|-monomial orbit-sum invariant
// (unavoidably ≥ |H|/2 terms for sign-character kernels) is bypassed.
//   ok(true)    — certified G_f ≤ H;
//   ok(false)   — certified G_f ⊄ H;
//   ok(nullopt) — no δ-form realises the character, or the matched form
//                 degenerates on these roots: the caller MUST fall back
//                 to the general invariant route (always sound).
// `base` may come back at a higher (certified) precision.
[[nodiscard]] Result<std::optional<bool>> quadratic_character_descent(
    const permgrp::BsgsGroup& G, const permgrp::BsgsGroup& H,
    const IntPoly& f_monic, galois_padic::PadicSplitting& base,
    symbolic::CASContext& ctx,
    const primitive_internal::Deadline& deadline = std::nullopt);

// Runs the full retry protocol for ONE candidate H below the current
// group: identity model first, then precision raises on demand
// (NeedPrecision — applied to `base`, the ORIGINAL splitting, the only
// one Newton can safely refine) and the derived Tschirnhaus sweep on
// certified multiple integer roots (Collision). Returns the hit, or
// nullopt for a certified NotContained. `base` may come back at a
// higher precision (kept for the following candidates). Structured
// Unimplemented when the derived sweep bound is exhausted.
[[nodiscard]] Result<std::optional<DescentHit>> test_candidate_with_retries(
    const permgrp::BsgsGroup& current, const permgrp::BsgsGroup& H,
    const galois_invariant::RelativeInvariant& inv, const IntPoly& f_monic,
    galois_padic::PadicSplitting& base, symbolic::CASContext& ctx,
    const primitive_internal::Deadline& deadline = std::nullopt);

struct FirstLayerDescent {
    bool disc_square{false};
    // True: G_f equals the ambient group EXACTLY (certified by the
    // coverage contract: no maximal transitive candidate contains it).
    bool is_ambient{false};
    std::string ambient_name;  // structural: "S<n>" / "A<n>"
    // Set iff !is_ambient: a certified proper containment G_f ≤ subgroup
    // (one descent step below the ambient; identification continues in
    // the next brick).
    std::optional<permgrp::BsgsGroup> subgroup;
    std::string provenance;  // Brick-2 provenance of the descended node
    std::optional<galois_padic::PadicSplitting> splitting;
};

// Runs the first descent layer for a monic irreducible squarefree
// f ∈ Z[x] of degree 5..10. Preconditions checked structurally
// (monicity, squarefreeness via disc ≠ 0); irreducibility is the
// caller's certified responsibility (as in galois.cpp).
[[nodiscard]] Result<FirstLayerDescent> stauduhar_first_layer(
    const IntPoly& f_monic, symbolic::CASContext& ctx,
    const primitive_internal::Deadline& deadline = std::nullopt);

// S_n (alternating = false) or A_n (true) as a verified BsgsGroup, with
// the exact order cross-checked. Shared by the first layer and the
// below-first-layer walk.
[[nodiscard]] Result<permgrp::BsgsGroup> ambient_bsgs(std::size_t n,
                                                      bool alternating);

// The exact Galois group, certified end to end (A6 Brick 3.5).
struct GaloisIdentification {
    bool disc_square{false};
    std::string ambient_name;  // structural: "S<n>" / "A<n>"
    // G_f on the root indexing fixed by the p-adic splitting. EXACT: each
    // descent step is a certified Stauduhar containment, and the terminal
    // node is certified by exhausting its maximal transitive classes (any
    // proper transitive subgroup lies inside one of them).
    permgrp::BsgsGroup group;
    std::size_t descent_steps{0U};  // 0 ⟺ G_f = ambient
    std::string first_layer_provenance;  // empty ⟺ G_f = ambient
};

// Full Stauduhar identification for a monic irreducible squarefree
// f ∈ Z[x] of degree 5..10: first layer (Brick-2 maximal candidates),
// then the below-first-layer walk on the maximal transitive classes of
// each interior node (galois_sublattice dense route, or the structural
// wreath-preimage route of galois_wreath_maximal for the 5|2-block
// degree-10 nodes beyond the dense byte budget — A6 Brick 3.75).
[[nodiscard]] Result<GaloisIdentification> stauduhar_identify(
    const IntPoly& f_monic, symbolic::CASContext& ctx,
    const primitive_internal::Deadline& deadline = std::nullopt);

}  // namespace cas::algebra::galois_stauduhar
