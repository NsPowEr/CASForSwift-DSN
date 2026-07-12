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
// G_f = ambient EXACTLY (coverage + maximality). One certified descent
// step below the ambient is returned as a containment; continuing the
// walk below the first layer needs the maximals of the inner nodes —
// that is the next brick, surfaced as a structured Unimplemented by the
// full identification driver (HC-F8-PENDING-09 stays open).

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

}  // namespace cas::algebra::galois_stauduhar
