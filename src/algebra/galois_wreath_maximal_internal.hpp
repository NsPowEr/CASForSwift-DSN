// A6 Brick 3.75 — structural maximal transitive subgroups of the
// wreath-preimage descent nodes that exceed the dense sublattice budget
// (the 5|2-block family of degree 10: S₅ ≀ S₂ of order 28800 and its two
// transitive index-2 subgroups of order 14400).
//
// Setting. Let W = S_d ≀ S_k in its imprimitive action on n = d·k points
// and φ : W → Q = C₂ ≀ S_k the sign-quotient homomorphism (per-block
// parity of the in-block map, induced block permutation); ker φ = A_d^k.
// A walk node H is a WREATH-PREIMAGE node when H = φ⁻¹(Q_H) for
// Q_H = φ(H) — certified exactly by |H| = |Q_H|·(d!/2)^k on H's minimal
// block system (unique for such nodes: A_d, d ≥ 5, is primitive on each
// block, so the d-block system is the finest and any other invariant
// partition coarsens it).
//
// Theorem (coverage; hypotheses guarded structurally below). Let d ≥ 5,
// d ≠ 6, k prime, and let some odd prime p | d satisfy p > k. Every
// maximal transitive subgroup M of H = φ⁻¹(Q_H) is H-conjugate to a
// member of the union of three structurally generated families:
//
//   FA — φ⁻¹(Q_m) for Q_m a maximal subgroup class of Q_H. Covers every
//        M ⊇ ker φ (bijection with the maximal classes of Q_H, computed
//        EXHAUSTIVELY by the dense enumeration on |Q_H| ≤ 2^k·k!), and
//        every M with M·ker φ < H (such an M lies below a proper full
//        preimage, hence below an FA member).
//   FB — (K ≀ S_k)^t ∩ H for K a maximal transitive shell of S_d (the
//        Brick-2 candidates, A_d included) and t a right-transversal rep
//        of H in W. Covers every remaining M whose base intersection
//        D = M ∩ A_d^k has PROPER coordinate projections D_i: aligning
//        the D_i by a base conjugation w (possible because m ∈ M
//        conjugates D_i onto D_{σ(i)}) puts M inside N_{S_d}(D₁) ≀ S_k;
//        N_{S_d}(D₁) is proper (D₁ ◁ S_d would force D₁ ∈ {1, A_d} by
//        d ≥ 5), and if it is intransitive so is M (every in-block part
//        of M lies in it). Writing w = h·t with h ∈ H realises M, up to
//        H-conjugacy, inside (shell ≀ S_k)^t ∩ H.
//   FD — N_W(Δ^t) ∩ H for Δ = {(a, …, a)} the full diagonal of A_d^k.
//        Covers every remaining M with all projections D_i = A_d and
//        D proper subdirect: by Scott's lemma (A_d non-abelian simple)
//        D is a product of full diagonals over a partition of the k
//        coordinates; the partition is invariant under the transitive
//        block action of M, so k prime forces the FULL diagonal Δ_ψ,
//        with ψ realised by conjugation (Aut(A_d) = S_d for d ≠ 6 —
//        the d = 6 guard), i.e. Δ_ψ = Δ^w with w ∈ W, and M ≤ N_W(Δ^w);
//        w = h·t as in FB. N_W(Δ) = {(c, …, c)·σ} because
//        C_{S_d}(A_d) = 1 for d ≥ 4.
//        The last case, D = 1, cannot yield a transitive M: the in-block
//        image of the block stabiliser of M is a section of M/(M ∩ ker φ)
//        ≅ Q_H of order dividing 2^k·k!, while transitivity on the block
//        demands order divisible by the guarded prime p | d with p > k,
//        and p ∤ 2^k·k!.
//
// References: L. L. Scott, "Representations in characteristic p", Proc.
// Symp. Pure Math. 37 (1980), Lemma p. 328 (subdirect subgroups of
// powers of a non-abelian simple group are products of diagonals);
// M. Aschbacher, L. Scott, "Maximal subgroups of finite groups",
// J. Algebra 92 (1985) (the M ⊇ ker φ / M·ker φ = H supplement split).
// The Scott step and the FB normalizer step are additionally
// machine-verified on dense ground truth in the unit tests (all maximal
// classes of A₅ × A₅; normalizers of the subgroup classes of A₅ inside
// S₅), and the (5,2) family lists are cross-checked end-to-end by the
// walk on real polynomials.
//
// The right-transversal twins are PROVEN sufficient: any w ∈ W splits as
// w = h·t (h ∈ H, t in the right transversal of H in W), and
// (shell)^{h·t} ∩ H = h·((shell)^t ∩ H)·h⁻¹ is H-conjugate to the listed
// (shell)^t ∩ H. Listed candidates need not be maximal and two listed
// candidates may be H-conjugate — both only cost redundant resolvent
// tests in the walk, never a wrong answer (the walk's termination
// certificate needs COVERAGE only). Intersections with H are computed
// exactly through φ: (S)^t ∩ H = preimage inside S^t of
// φ(S)^{φ(t)} ∩ Q_H, assembled from the even part of the shell (kernel
// side) plus lifted members of the small sign-image intersection — no
// element enumeration of H ever happens.

#pragma once

#include "perm_bsgs_internal.hpp"

#include "cas/result.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace cas::symbolic {
class CASContext;
}

namespace cas::algebra::permgrp {

// The certified wreath-preimage shape of a node H (see file header).
struct WreathPreimageShape {
    std::size_t d{0U};  // block size
    std::size_t k{0U};  // number of blocks
    // Standard-coordinates alignment: standard point j·d + x (x-th point
    // of standard block j) ↦ node point align[j·d + x].
    Perm align;
    // H's generators conjugated into standard coordinates.
    std::vector<Perm> std_gens;
    // Q_H = φ(H) ≤ C₂ ≀ S_k realised on 2k points (pairs {2i, 2i+1}).
    BsgsGroup q_image;
};

// Detects whether H is a wreath-preimage node on one of its minimal
// block systems (|H| = |φ(H)|·(d!/2)^k). Returns nullopt when it is not
// (including intransitive H). Never enumerates H's elements.
[[nodiscard]] Result<std::optional<WreathPreimageShape>>
detect_wreath_preimage(const BsgsGroup& H);

// The FA ∪ FB ∪ FD candidate classes for a certified wreath-preimage
// node (see the coverage theorem in the file header), in NODE
// coordinates, transitive proper subgroups of H only, deduplicated,
// direct containments pruned, sorted by decreasing order. Structured
// Unimplemented when a theorem hypothesis fails (d < 5, d = 6,
// composite k, or no odd prime p | d with p > k) — each names the
// failing hypothesis. Budgets as in galois_sublattice.
[[nodiscard]] Result<std::vector<BsgsGroup>>
wreath_preimage_maximal_transitive(const BsgsGroup& H,
                                   const WreathPreimageShape& shape,
                                   std::uint64_t max_ops,
                                   std::uint64_t max_bytes,
                                   symbolic::CASContext* ctx);

// The walk's single entry point for the maximal transitive classes of a
// descent node: the dense exhaustive route when the node fits the byte
// budget and the u16 universe (identical behaviour to Brick 3.5), else
// the structural wreath-preimage route. A node beyond the dense budgets
// that is not a certified wreath-preimage node fails with a structured
// Unimplemented naming both routes.
[[nodiscard]] Result<std::vector<BsgsGroup>> node_maximal_transitive_classes(
    const BsgsGroup& H, std::uint64_t max_ops, std::uint64_t max_bytes,
    symbolic::CASContext* ctx);

}  // namespace cas::algebra::permgrp
