// A6 Brick 3.75 — implementation of the structural wreath-preimage
// maximal-transitive generator. Coverage theorem, hypotheses and the
// φ-arithmetic used everywhere here are documented in
// galois_wreath_maximal_internal.hpp.

#include "galois_wreath_maximal_internal.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include "galois_sublattice_internal.hpp"
#include "perm_blocks_internal.hpp"
#include "perm_construct_internal.hpp"
#include "perm_group_internal.hpp"
#include "perm_maximal_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::algebra::permgrp {

namespace {

[[nodiscard]] bool is_prime_size(std::size_t v) {
    if (v < 2U) return false;
    for (std::size_t q = 2U; q * q <= v; ++q) {
        if (v % q == 0U) return false;
    }
    return true;
}

// An odd prime p | d with p > k (excludes the D = 1 supplement case).
[[nodiscard]] bool has_odd_prime_factor_above(std::size_t d, std::size_t k) {
    std::size_t m = d;
    for (std::size_t q = 2U; q * q <= m; ++q) {
        while (m % q == 0U) {
            if (q % 2U == 1U && q > k) return true;
            m /= q;
        }
    }
    return m > 1U && m % 2U == 1U && m > k;
}

[[nodiscard]] Perm conj_by(const Perm& t, const Perm& g) {
    return compose(t, compose(g, inverse(t)));
}

[[nodiscard]] std::vector<Perm> conj_all(const Perm& t,
                                         const std::vector<Perm>& gens) {
    std::vector<Perm> out;
    out.reserve(gens.size());
    for (const Perm& g : gens) out.push_back(conj_by(t, g));
    return out;
}

// φ of a standard-coordinates block-preserving g ∈ S_d ≀ S_k, realised on
// the 2k points of C₂ ≀ S_k (pairs {2i, 2i+1}): q(2j + b) = 2σ(j) + (b ⊕
// ε_j) with σ the block permutation and ε_j the parity of the in-block
// map of block j. This is a homomorphism under compose(a, b) = a ∘ b.
[[nodiscard]] Result<Perm> phi_image(const Perm& g, std::size_t d,
                                     std::size_t k) {
    Perm q(2U * k);
    Perm in_block(d);
    for (std::size_t j = 0U; j < k; ++j) {
        const std::size_t target = static_cast<std::size_t>(g[j * d]) / d;
        for (std::size_t x = 0U; x < d; ++x) {
            const std::size_t img = g[j * d + x];
            if (img / d != target) {
                return fail<Perm>(CASError{
                    .kind = CASErrorKind::InternalError,
                    .message = "phi_image: generator does not preserve the "
                               "standard block system"});
            }
            in_block[x] = static_cast<std::uint8_t>(img - target * d);
        }
        const std::size_t eps = is_odd(in_block) ? 1U : 0U;
        q[2U * j] = static_cast<std::uint8_t>(2U * target + eps);
        q[2U * j + 1U] = static_cast<std::uint8_t>(2U * target + (1U - eps));
    }
    return ok(std::move(q));
}

// A lift of q = (σ, ε⃗) ∈ C₂ ≀ S_k back to S_d ≀ S_k: block j maps onto
// block σ(j) through the identity (ε_j = 0) or through `decorator`
// (ε_j = 1; decorator must be an ODD element of the intended shell so
// φ(lift) = q and the lift stays inside the shell).
[[nodiscard]] Perm lift_q(const Perm& q, std::size_t d, std::size_t k,
                          const Perm& decorator) {
    Perm w(d * k);
    for (std::size_t j = 0U; j < k; ++j) {
        const std::size_t sigma_j = static_cast<std::size_t>(q[2U * j]) / 2U;
        const bool eps = (static_cast<std::size_t>(q[2U * j]) % 2U) == 1U;
        for (std::size_t x = 0U; x < d; ++x) {
            const std::size_t y = eps ? decorator[x] : x;
            w[j * d + x] = static_cast<std::uint8_t>(sigma_j * d + y);
        }
    }
    return w;
}

// g ∈ S_d embedded in block j (identity elsewhere).
[[nodiscard]] Perm embed_in_block(const Perm& g, std::size_t j,
                                  std::size_t d, std::size_t k) {
    Perm w = identity(d * k);
    for (std::size_t x = 0U; x < d; ++x) {
        w[j * d + x] = static_cast<std::uint8_t>(j * d + g[x]);
    }
    return w;
}

// g ∈ S_d applied diagonally in every block.
[[nodiscard]] Perm embed_diagonal(const Perm& g, std::size_t d,
                                  std::size_t k) {
    Perm w(d * k);
    for (std::size_t j = 0U; j < k; ++j) {
        for (std::size_t x = 0U; x < d; ++x) {
            w[j * d + x] = static_cast<std::uint8_t>(j * d + g[x]);
        }
    }
    return w;
}

// τ ∈ S_k as a sign-free element of C₂ ≀ S_k on 2k points.
[[nodiscard]] Perm q_top(const Perm& tau, std::size_t k) {
    Perm q(2U * k);
    for (std::size_t j = 0U; j < k; ++j) {
        q[2U * j] = static_cast<std::uint8_t>(2U * tau[j]);
        q[2U * j + 1U] = static_cast<std::uint8_t>(2U * tau[j] + 1U);
    }
    return q;
}

// BFS element enumeration of a SMALL permutation group (≤ 2^k·k! here).
[[nodiscard]] std::vector<Perm> closure_perms(std::size_t n,
                                              const std::vector<Perm>& gens) {
    std::vector<Perm> elems{identity(n)};
    for (std::size_t i = 0U; i < elems.size(); ++i) {
        for (const Perm& g : gens) {
            Perm nxt = compose(g, elems[i]);
            if (std::find(elems.begin(), elems.end(), nxt) == elems.end()) {
                elems.push_back(std::move(nxt));
            }
        }
    }
    return elems;
}

// The first odd strong generator of K, if any (parity is a homomorphism,
// so K has odd elements iff a strong generator is odd).
[[nodiscard]] const Perm* find_odd(const BsgsGroup& k_group) {
    for (const Perm& g : k_group.strong_generators()) {
        if (is_odd(g)) return &g;
    }
    return nullptr;
}

[[nodiscard]] Result<std::vector<BsgsGroup>> hypothesis_fail(
    const std::string& what) {
    return fail<std::vector<BsgsGroup>>(CASError{
        .kind = CASErrorKind::Unimplemented,
        .message = "wreath_preimage_maximal_transitive: " + what});
}

}  // namespace

Result<std::optional<WreathPreimageShape>> detect_wreath_preimage(
    const BsgsGroup& H) {
    using Out = std::optional<WreathPreimageShape>;
    const std::size_t n = H.degree();
    if (!H.is_transitive()) return ok(Out{});
    auto systems = minimal_block_systems(n, H.generators());
    if (systems.is_error()) return fail<Out>(systems.error());
    for (const BlockSystem& sys : systems.value()) {
        const std::size_t k = sys.num_blocks;
        const std::size_t d = n / k;
        if (d < 2U || k < 2U || d > 20U) continue;
        // Alignment: standard point j·d + x ↦ x-th smallest point of
        // block j (block indices are first-appearance normalised).
        std::vector<std::vector<std::uint8_t>> blocks(k);
        for (std::size_t p = 0U; p < n; ++p) {
            blocks[sys.block_of[p]].push_back(static_cast<std::uint8_t>(p));
        }
        Perm align(n);
        for (std::size_t j = 0U; j < k; ++j) {
            for (std::size_t x = 0U; x < d; ++x) {
                align[j * d + x] = blocks[j][x];
            }
        }
        const Perm align_inv = inverse(align);
        std::vector<Perm> std_gens;
        std_gens.reserve(H.generators().size());
        std::vector<Perm> q_gens;
        q_gens.reserve(H.generators().size());
        for (const Perm& g : H.generators()) {
            Perm gs = compose(align_inv, compose(g, align));
            auto q = phi_image(gs, d, k);
            if (q.is_error()) return fail<Out>(q.error());
            q_gens.push_back(std::move(q.value()));
            std_gens.push_back(std::move(gs));
        }
        auto q_built = BsgsGroup::build(2U * k, std::move(q_gens));
        if (q_built.is_error()) return fail<Out>(q_built.error());
        // Full-preimage certificate: |H| = |Q_H|·(d!/2)^k, u64-checked.
        const std::uint64_t half_fact = factorial_u64(d) / 2ULL;
        std::uint64_t pow = 1ULL;
        bool overflow = false;
        for (std::size_t i = 0U; i < k; ++i) {
            if (half_fact != 0ULL && pow > UINT64_MAX / half_fact) {
                overflow = true;
                break;
            }
            pow *= half_fact;
        }
        if (overflow) continue;
        const std::uint64_t q_order = q_built.value().order();
        if (q_order != 0ULL && pow > UINT64_MAX / q_order) continue;
        if (H.order() != q_order * pow) continue;
        return ok(Out{WreathPreimageShape{d, k, std::move(align),
                                          std::move(std_gens),
                                          std::move(q_built.value())}});
    }
    return ok(Out{});
}

Result<std::vector<BsgsGroup>> wreath_preimage_maximal_transitive(
    const BsgsGroup& H, const WreathPreimageShape& shape,
    std::uint64_t max_ops, std::uint64_t max_bytes,
    symbolic::CASContext* ctx) {
    const std::size_t d = shape.d;
    const std::size_t k = shape.k;
    const std::size_t n = d * k;
    if (d < 5U) {
        return hypothesis_fail(
            "A_d is not non-abelian simple for d < 5 (Scott hypothesis); "
            "such nodes are covered by the dense sublattice route");
    }
    if (d == 6U) {
        return hypothesis_fail(
            "d = 6 needs the exceptional Out(A_6) diagonal classes — a "
            "later A6 increment");
    }
    if (!is_prime_size(k)) {
        return hypothesis_fail(
            "composite k admits partitioned-diagonal families — a later "
            "A6 increment");
    }
    if (!has_odd_prime_factor_above(d, k)) {
        return hypothesis_fail(
            "no odd prime p | d with p > k: the trivial-intersection "
            "supplement case is not excluded for this (d, k)");
    }

    const BsgsGroup& q_h = shape.q_image;
    const Perm t01 = [&] {
        Perm t = identity(d);
        t[0] = 1U;
        t[1] = 0U;
        return t;
    }();

    // ── right transversal of H in W, computed in the sign quotient ────────
    auto w_full = wreath_gens(2U, k);
    if (w_full.is_error()) {
        return fail<std::vector<BsgsGroup>>(w_full.error());
    }
    const std::vector<Perm> q_elems = closure_perms(2U * k, w_full.value());
    std::vector<Perm> transversal_q;
    for (const Perm& q : q_elems) {
        bool fresh = true;
        for (const Perm& r : transversal_q) {
            // Same right coset Q_H·t ⟺ q·r⁻¹ ∈ Q_H.
            if (q_h.contains(compose(q, inverse(r)))) {
                fresh = false;
                break;
            }
        }
        if (fresh) transversal_q.push_back(q);
    }

    // ── assemble the raw candidate generator sets (standard coords) ────────
    std::vector<std::vector<Perm>> raw;

    // FA — lifts of the maximal subgroup classes of Q_H (kernel A_d^k).
    std::vector<Perm> kernel_full;
    for (std::size_t j = 0U; j < k; ++j) {
        for (const Perm& g : alternating_gens(d)) {
            kernel_full.push_back(embed_in_block(g, j, d, k));
        }
    }
    auto q_max = maximal_subgroup_classes_in(q_h, max_ops, max_bytes, ctx);
    if (q_max.is_error()) {
        return fail<std::vector<BsgsGroup>>(q_max.error());
    }
    for (const BsgsGroup& qm : q_max.value()) {
        std::vector<Perm> gens = kernel_full;
        for (const Perm& g : qm.generators()) {
            gens.push_back(lift_q(g, d, k, t01));
        }
        raw.push_back(std::move(gens));
    }

    // Shared shell machinery for FB and FD: candidate = (shell)^t ∩ H =
    // t·(kernel side)·t⁻¹ plus t-conjugated lifts of the members e of the
    // shell's sign image with t·e·t⁻¹ ∈ Q_H (all exact through φ).
    auto push_shell = [&](const std::vector<Perm>& shell_kernel,
                          const std::vector<Perm>& img_gens,
                          const Perm& decorator) {
        const std::vector<Perm> img_elems = closure_perms(2U * k, img_gens);
        const Perm q_id = identity(2U * k);
        for (const Perm& tq : transversal_q) {
            const Perm tq_inv = inverse(tq);
            const Perm w_t = lift_q(tq, d, k, t01);
            std::vector<Perm> gens = conj_all(w_t, shell_kernel);
            for (const Perm& e : img_elems) {
                if (e == q_id) continue;
                if (!q_h.contains(compose(tq, compose(e, tq_inv)))) continue;
                gens.push_back(conj_by(w_t, lift_q(e, d, k, decorator)));
            }
            raw.push_back(std::move(gens));
        }
    };

    // FB — (K ≀ S_k)^t ∩ H for the Brick-2 maximal transitive shells of
    // S_d (A_d included, listed first by the Brick-2 contract).
    auto shells = maximal_transitive_candidates(AmbientGroup::Symmetric, d);
    if (shells.is_error()) {
        return fail<std::vector<BsgsGroup>>(shells.error());
    }
    std::vector<Perm> tops_q;
    for (const Perm& tau : symmetric_gens(k)) tops_q.push_back(q_top(tau, k));
    for (const MaximalCandidate& shell : shells.value()) {
        auto k_even = even_part_gens(d, shell.group.generators());
        if (k_even.is_error()) {
            return fail<std::vector<BsgsGroup>>(k_even.error());
        }
        std::vector<Perm> shell_kernel;
        for (std::size_t j = 0U; j < k; ++j) {
            for (const Perm& g : k_even.value()) {
                shell_kernel.push_back(embed_in_block(g, j, d, k));
            }
        }
        const Perm* odd_k = find_odd(shell.group);
        std::vector<Perm> img_gens = tops_q;
        if (odd_k != nullptr) {
            // A single block-0 sign flip: with the tops it generates the
            // full sign image C₂ ≀ S_k of an odd shell.
            Perm flip0 = identity(2U * k);
            flip0[0] = 1U;
            flip0[1] = 0U;
            img_gens.push_back(std::move(flip0));
        }
        push_shell(shell_kernel, img_gens,
                   odd_k != nullptr ? *odd_k : t01);
    }

    // FD — N_W(Δ^t) ∩ H, Δ the full diagonal of A_d^k. Kernel side
    // diag(A_d); sign image ⟨uniform flip, tops⟩ ≅ C₂ × S_k; lifts are
    // diagonal decorations, so they stay inside N_W(Δ) = {(c, …, c)·σ}.
    {
        std::vector<Perm> diag_kernel;
        for (const Perm& g : alternating_gens(d)) {
            diag_kernel.push_back(embed_diagonal(g, d, k));
        }
        Perm all_flip(2U * k);
        for (std::size_t j = 0U; j < k; ++j) {
            all_flip[2U * j] = static_cast<std::uint8_t>(2U * j + 1U);
            all_flip[2U * j + 1U] = static_cast<std::uint8_t>(2U * j);
        }
        std::vector<Perm> img_gens = tops_q;
        img_gens.push_back(std::move(all_flip));
        push_shell(diag_kernel, img_gens, t01);
    }

    // ── verify, filter, dedup, prune, sort, map back to node coords ────────
    auto h_std = BsgsGroup::build(n, shape.std_gens);
    if (h_std.is_error()) {
        return fail<std::vector<BsgsGroup>>(h_std.error());
    }
    const std::uint64_t h_order = H.order();
    std::vector<BsgsGroup> kept;
    for (std::vector<Perm>& gens : raw) {
        auto built = BsgsGroup::build(n, std::move(gens));
        if (built.is_error()) {
            return fail<std::vector<BsgsGroup>>(built.error());
        }
        BsgsGroup cand = std::move(built.value());
        if (cand.order() >= h_order) {
            if (cand.order() > h_order) {
                return fail<std::vector<BsgsGroup>>(CASError{
                    .kind = CASErrorKind::InternalError,
                    .message = "wreath_preimage_maximal_transitive: "
                               "candidate exceeds the node order"});
            }
            continue;  // H itself: not a descent
        }
        if (!cand.is_transitive()) continue;
        for (const Perm& g : cand.generators()) {
            if (!h_std.value().contains(g)) {
                return fail<std::vector<BsgsGroup>>(CASError{
                    .kind = CASErrorKind::InternalError,
                    .message = "wreath_preimage_maximal_transitive: "
                               "candidate escapes the node (φ arithmetic "
                               "violated)"});
            }
        }
        bool duplicate = false;
        for (const BsgsGroup& o : kept) {
            if (o.order() != cand.order()) continue;
            bool inside = true;
            for (const Perm& g : cand.generators()) {
                if (!o.contains(g)) {
                    inside = false;
                    break;
                }
            }
            if (inside) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) kept.push_back(std::move(cand));
    }
    // Direct containment prune (same coordinates): a candidate lying in a
    // strictly bigger one is reachable through it later in the walk.
    std::vector<bool> drop(kept.size(), false);
    for (std::size_t i = 0U; i < kept.size(); ++i) {
        for (std::size_t j = 0U; j < kept.size() && !drop[i]; ++j) {
            if (i == j || drop[j]) continue;
            if (kept[j].order() <= kept[i].order() ||
                kept[j].order() % kept[i].order() != 0ULL) {
                continue;
            }
            bool inside = true;
            for (const Perm& g : kept[i].generators()) {
                if (!kept[j].contains(g)) {
                    inside = false;
                    break;
                }
            }
            if (inside) drop[i] = true;
        }
    }
    std::vector<BsgsGroup> out;
    for (std::size_t i = 0U; i < kept.size(); ++i) {
        if (drop[i]) continue;
        std::vector<Perm> gens;
        gens.reserve(kept[i].generators().size());
        for (const Perm& g : kept[i].generators()) {
            gens.push_back(conj_by(shape.align, g));
        }
        auto back = BsgsGroup::build(n, std::move(gens));
        if (back.is_error()) return fail<std::vector<BsgsGroup>>(back.error());
        if (back.value().order() != kept[i].order()) {
            return fail<std::vector<BsgsGroup>>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "wreath_preimage_maximal_transitive: alignment "
                           "conjugation changed a candidate's order"});
        }
        out.push_back(std::move(back.value()));
    }
    std::sort(out.begin(), out.end(),
              [](const BsgsGroup& a, const BsgsGroup& b) {
                  return a.order() > b.order();
              });
    return ok(std::move(out));
}

Result<std::vector<BsgsGroup>> node_maximal_transitive_classes(
    const BsgsGroup& H, std::uint64_t max_ops, std::uint64_t max_bytes,
    symbolic::CASContext* ctx) {
    if (H.order() <= 65536ULL &&
        dense_sublattice_min_bytes(H.order(), H.degree()) <= max_bytes) {
        return transitive_subgroup_classes_in(H, max_ops, max_bytes, ctx);
    }
    auto shape = detect_wreath_preimage(H);
    if (shape.is_error()) {
        return fail<std::vector<BsgsGroup>>(shape.error());
    }
    if (!shape.value().has_value()) {
        return fail<std::vector<BsgsGroup>>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message =
                "node_maximal_transitive_classes: node exceeds the dense "
                "sublattice budgets (raise "
                "CASContext::galois_sublattice_max_bytes) and is not a "
                "certified wreath-preimage node — structural maximals for "
                "this shape await a later A6 increment"});
    }
    return wreath_preimage_maximal_transitive(H, *shape.value(), max_ops,
                                              max_bytes, ctx);
}

}  // namespace cas::algebra::permgrp
