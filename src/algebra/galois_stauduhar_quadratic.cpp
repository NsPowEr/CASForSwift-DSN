// A6 Brick 3.75 — certified quadratic-character descent for index-2
// candidates (the δ-form fast path of the Stauduhar walk).
//
// Mathematics. An index-2 subgroup H of the current group G is NORMAL,
// so it is the kernel of the quadratic character χ : G → {±1} with
// χ(g) = −1 ⟺ g ∉ H, and the Stauduhar test needs no conjugator:
// G_f ≤ H ⟺ G_f fixes any χ-relative invariant δ (g·δ = χ(g)·δ), i.e.
// ⟺ the p-adic value v = δ(r₁..r_n) is a RATIONAL integer. The resolvent
// is R(y) = y² − v² with v² a G-invariant integer: the Brick-3 purity
// certificate applies verbatim with m = 2 — v pure with residue c and
// 2c ≢ 0 (automatic for c ≠ 0 at |2c| < p^k) certifies the descent, v
// impure at the derived precision certifies non-containment.
//
// δ-forms. Instead of the |H|-monomial orbit-sum invariant (whose
// evaluation dominated the walk: a sign-character kernel admits NO small
// monomial orbit sum — the seed stabiliser must avoid every same-block
// transposition, forcing ≥ |H|/2 terms), χ is matched against ±1
// combinations of ALTERNATING PRODUCTS evaluated as products:
//   δ_b   = ∏_{i<j ∈ block b} (x_i − x_j)   (per-block Vandermonde),
//   δ_top = ∏_{b<b'} (S_b − S_{b'}),  S_b = Σ_{i ∈ block b} x_i,
//   δ_Ω   = ∏_{i<j} (x_i − x_j)             (the discriminant root),
// over each minimal block system of G (δ_Ω needs none). Candidate forms:
// ∏_b δ_b, (∏_b δ_b)·δ_top, δ_top, δ_Ω and — for two blocks — the
// diagonal-sign forms δ₀ ± δ₁. A form matches iff g·form = χ(g)·form for
// every generator of G (the action permutes blocks and multiplies by the
// exact in-block/top/global signs; characters agree on generators ⟹
// they agree on G). No match → the caller falls back to the general
// monomial-invariant route, which is a guaranteed terminator — this
// module is an exact fast path, never a semantic fork (CLAUDE.md cat. 8:
// fast path with an always-present general route).
//
// Degeneracy. A pure ∏_b δ_b form never vanishes mod p (its square
// divides disc(f)^k-free products of distinct-root differences and p is
// unramified), so the test always concludes. Forms involving δ_top or
// δ₀ ± δ₁ CAN vanish on the specific roots (equal block sums, opposite
// block Vandermondes): v ≡ 0 at the certified precision is then reported
// as "no match" (nullopt) and the general route decides — sound, never
// silent.

#include "galois_stauduhar_internal.hpp"

#include "cas/bigint.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include "galois_internal.hpp"
#include "perm_blocks_internal.hpp"
#include "polynomial_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace cas::algebra::galois_stauduhar {

namespace {

using galois_padic::PadicRing;
using galois_padic::PadicSplitting;
using galois_padic::RingElem;
using permgrp::BlockSystem;
using permgrp::BsgsGroup;
using permgrp::is_odd;
using permgrp::Perm;
using primitive_internal::Deadline;

// One additive term of a δ-form: coef·(∏_{b ∈ mask} δ_b)·δ_top^t·δ_Ω^o.
struct DeltaTerm {
    int coef{1};             // ±1
    std::uint32_t mask{0U};  // per-block Vandermonde factors
    std::uint8_t top{0U};    // block-sum Vandermonde power (0/1)
    std::uint8_t omega{0U};  // full Vandermonde power (0/1)

    [[nodiscard]] bool operator<(const DeltaTerm& o) const {
        if (mask != o.mask) return mask < o.mask;
        if (top != o.top) return top < o.top;
        if (omega != o.omega) return omega < o.omega;
        return coef < o.coef;
    }
    [[nodiscard]] bool operator==(const DeltaTerm& o) const {
        return coef == o.coef && mask == o.mask && top == o.top &&
               omega == o.omega;
    }
};
using DeltaForm = std::vector<DeltaTerm>;  // kept sorted

// The exact signs with which g acts on the δ-building-blocks of `blocks`
// (points of each block sorted ascending — the sign convention every
// evaluation below shares).
struct GeneratorSigns {
    std::vector<std::size_t> block_image;  // b ↦ σ(b)
    std::vector<int> block_sign;           // sgn of the induced in-block map
    int top_sign{1};                       // sgn of σ on the blocks
    int omega_sign{1};                     // sgn of g on all points
};

[[nodiscard]] Result<GeneratorSigns> generator_signs(
    const Perm& g, const std::vector<std::vector<std::uint8_t>>& blocks) {
    const std::size_t k = blocks.size();
    GeneratorSigns out;
    out.block_image.resize(k);
    out.block_sign.resize(k, 1);
    Perm sigma(k);
    for (std::size_t b = 0U; b < k; ++b) {
        const std::size_t d = blocks[b].size();
        // Image block: locate g(first point of b).
        const std::uint8_t img0 = g[blocks[b][0]];
        std::size_t tb = k;
        for (std::size_t c = 0U; c < k; ++c) {
            if (std::binary_search(blocks[c].begin(), blocks[c].end(),
                                   img0)) {
                tb = c;
                break;
            }
        }
        if (tb == k) {
            return fail<GeneratorSigns>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "quadratic_character_descent: point escapes "
                           "the block system"});
        }
        // Induced map in sorted-position coordinates.
        Perm in_block(d);
        for (std::size_t x = 0U; x < d; ++x) {
            const std::uint8_t img = g[blocks[b][x]];
            const auto it = std::lower_bound(blocks[tb].begin(),
                                             blocks[tb].end(), img);
            if (it == blocks[tb].end() || *it != img) {
                return fail<GeneratorSigns>(CASError{
                    .kind = CASErrorKind::InternalError,
                    .message = "quadratic_character_descent: block system "
                               "not invariant"});
            }
            in_block[x] = static_cast<std::uint8_t>(
                std::distance(blocks[tb].begin(), it));
        }
        out.block_image[b] = tb;
        out.block_sign[b] = is_odd(in_block) ? -1 : 1;
        sigma[b] = static_cast<std::uint8_t>(tb);
    }
    out.top_sign = is_odd(sigma) ? -1 : 1;
    out.omega_sign = is_odd(g) ? -1 : 1;
    return ok(std::move(out));
}

// g·form: permute the block masks, multiply the collected signs.
[[nodiscard]] DeltaForm apply_signs(const GeneratorSigns& s,
                                    const DeltaForm& form) {
    DeltaForm out;
    out.reserve(form.size());
    for (const DeltaTerm& t : form) {
        DeltaTerm r = t;
        r.mask = 0U;
        for (std::size_t b = 0U; b < s.block_image.size(); ++b) {
            if ((t.mask >> b) & 1U) {
                r.mask |= 1U << s.block_image[b];
                r.coef *= s.block_sign[b];
            }
        }
        if (t.top != 0U) r.coef *= s.top_sign;
        if (t.omega != 0U) r.coef *= s.omega_sign;
        out.push_back(r);
    }
    std::sort(out.begin(), out.end());
    return out;
}

[[nodiscard]] DeltaForm negate(DeltaForm f) {
    for (DeltaTerm& t : f) t.coef = -t.coef;
    std::sort(f.begin(), f.end());
    return f;
}

// χ_δ(g) ∈ {+1, −1}, or nullopt when g does not scale the form by ±1.
[[nodiscard]] std::optional<int> form_character(const GeneratorSigns& s,
                                                const DeltaForm& form) {
    const DeltaForm img = apply_signs(s, form);
    if (img == form) return 1;
    if (img == negate(form)) return -1;
    return std::optional<int>{};
}

// |δ-form(r)| bound from Cauchy: |r_i| ≤ ρ = 1 + max|f_j|. `n` is the
// full point count (δ_Ω needs it even when the form carries no blocks).
[[nodiscard]] BigInt form_bound(const IntPoly& f, const DeltaForm& form,
                                const std::vector<std::vector<std::uint8_t>>&
                                    blocks,
                                std::size_t n) {
    BigInt maxc(0);
    for (const auto& c : f.coefficients()) {
        BigInt a = c;
        if (a.is_negative()) a = -a;
        if (a > maxc) maxc = a;
    }
    const BigInt rho = maxc + BigInt(1);
    const std::size_t k = blocks.size();
    BigInt total(0);
    for (const DeltaTerm& t : form) {
        BigInt term(1);
        for (std::size_t b = 0U; b < k; ++b) {
            if (((t.mask >> b) & 1U) == 0U) continue;
            const std::size_t d = blocks[b].size();
            term = term * bigint_pow_nonnegative(BigInt(2) * rho,
                                                 d * (d - 1U) / 2U);
        }
        if (t.top != 0U) {
            // |S_b − S_b'| ≤ 2·d·ρ with d = max block size.
            std::size_t dmax = 0U;
            for (const auto& b : blocks) dmax = std::max(dmax, b.size());
            term = term *
                   bigint_pow_nonnegative(
                       BigInt(2) * BigInt(static_cast<std::int64_t>(dmax)) *
                           rho,
                       k * (k - 1U) / 2U);
        }
        if (t.omega != 0U) {
            term = term * bigint_pow_nonnegative(BigInt(2) * rho,
                                                 n * (n - 1U) / 2U);
        }
        total = total + term;
    }
    return total;
}

// Smallest k with p^k > 2·B² (covers the purity certificate of both v
// and the G-invariant integer v²). Bit-length over-approximation.
[[nodiscard]] std::size_t precision_for(const BigInt& p, const BigInt& B) {
    const std::size_t bits = 2U * (B + BigInt(1)).bit_length() + 3U;
    const std::size_t denom = p.bit_length() > 1U ? p.bit_length() - 1U : 1U;
    return bits / denom + 1U;
}

[[nodiscard]] RingElem eval_form(const PadicSplitting& split,
                                 const DeltaForm& form,
                                 const std::vector<std::vector<std::uint8_t>>&
                                     blocks) {
    const PadicRing& R = split.ring;
    const std::size_t k = blocks.size();
    // Per-block Vandermonde and block sums.
    std::vector<RingElem> vand(k, R.one());
    std::vector<RingElem> sums(k, R.zero());
    for (std::size_t b = 0U; b < k; ++b) {
        for (std::size_t i = 0U; i < blocks[b].size(); ++i) {
            sums[b] = R.add(sums[b], split.roots[blocks[b][i]]);
            for (std::size_t j = i + 1U; j < blocks[b].size(); ++j) {
                vand[b] = R.mul(
                    vand[b], R.sub(split.roots[blocks[b][i]],
                                   split.roots[blocks[b][j]]));
            }
        }
    }
    RingElem top = R.one();
    for (std::size_t b = 0U; b < k; ++b) {
        for (std::size_t c = b + 1U; c < k; ++c) {
            top = R.mul(top, R.sub(sums[b], sums[c]));
        }
    }
    RingElem omega = R.one();
    const std::size_t n = split.roots.size();
    for (std::size_t i = 0U; i < n; ++i) {
        for (std::size_t j = i + 1U; j < n; ++j) {
            omega = R.mul(omega, R.sub(split.roots[i], split.roots[j]));
        }
    }
    RingElem acc = R.zero();
    for (const DeltaTerm& t : form) {
        RingElem term = R.one();
        for (std::size_t b = 0U; b < k; ++b) {
            if ((t.mask >> b) & 1U) term = R.mul(term, vand[b]);
        }
        if (t.top != 0U) term = R.mul(term, top);
        if (t.omega != 0U) term = R.mul(term, omega);
        acc = (t.coef > 0) ? R.add(acc, term) : R.sub(acc, term);
    }
    return acc;
}

}  // namespace

Result<std::optional<bool>> quadratic_character_descent(
    const BsgsGroup& G, const BsgsGroup& H, const IntPoly& f_monic,
    PadicSplitting& base, symbolic::CASContext& ctx, const Deadline& dl) {
    using Out = std::optional<bool>;
    if (G.order() != 2ULL * H.order()) {
        return fail<Out>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "quadratic_character_descent: [G:H] must be 2"});
    }
    // χ on the generators of G (index 2 ⇒ H ◁ G ⇒ χ is a character).
    std::vector<int> chi;
    chi.reserve(G.generators().size());
    for (const Perm& g : G.generators()) {
        chi.push_back(H.contains(g) ? 1 : -1);
    }

    const std::size_t n = G.degree();
    // Candidate (blocks, form) pairs: the block-free δ_Ω first, then the
    // per-system forms over each minimal block system of G.
    std::vector<std::vector<std::uint8_t>> no_blocks;
    std::vector<std::pair<std::vector<std::vector<std::uint8_t>>, DeltaForm>>
        candidates;
    candidates.emplace_back(no_blocks,
                            DeltaForm{DeltaTerm{1, 0U, 0U, 1U}});
    auto systems = permgrp::minimal_block_systems(n, G.generators());
    if (systems.is_error()) return fail<Out>(systems.error());
    for (const BlockSystem& sys : systems.value()) {
        const std::size_t k = sys.num_blocks;
        if (k > 32U) continue;  // mask width (n ≤ 20 in this engine)
        std::vector<std::vector<std::uint8_t>> blocks(k);
        for (std::size_t p = 0U; p < n; ++p) {
            blocks[sys.block_of[p]].push_back(static_cast<std::uint8_t>(p));
        }
        const std::uint32_t all =
            (k == 32U) ? 0xFFFFFFFFU : ((1U << k) - 1U);
        candidates.emplace_back(blocks,
                                DeltaForm{DeltaTerm{1, all, 0U, 0U}});
        candidates.emplace_back(blocks,
                                DeltaForm{DeltaTerm{1, all, 1U, 0U}});
        candidates.emplace_back(blocks,
                                DeltaForm{DeltaTerm{1, 0U, 1U, 0U}});
        if (k == 2U) {
            candidates.emplace_back(
                blocks, DeltaForm{DeltaTerm{1, 1U, 0U, 0U},
                                  DeltaTerm{1, 2U, 0U, 0U}});
            candidates.emplace_back(
                blocks, DeltaForm{DeltaTerm{-1, 2U, 0U, 0U},
                                  DeltaTerm{1, 1U, 0U, 0U}});
        }
    }

    for (auto& [blocks, form] : candidates) {
        std::sort(form.begin(), form.end());
        // Signs of every generator on this system (δ_Ω: k = 0 blocks).
        bool match = true;
        for (std::size_t gi = 0U; gi < G.generators().size() && match;
             ++gi) {
            auto s = generator_signs(G.generators()[gi], blocks);
            if (s.is_error()) return fail<Out>(s.error());
            const auto c = form_character(s.value(), form);
            match = c.has_value() && *c == chi[gi];
        }
        if (!match) continue;
        // Certified evaluation at the derived precision.
        const BigInt bound = form_bound(f_monic, form, blocks, n);
        const std::size_t need =
            precision_for(base.ring.prime(), bound);
        if (need > base.ring.precision()) {
            auto raised = galois_padic::raise_splitting_precision(
                base, need, &ctx, dl);
            if (raised.is_error()) return fail<Out>(raised.error());
            base = std::move(raised.value());
        }
        const RingElem v = eval_form(base, form, blocks);
        if (base.ring.is_zero(v)) {
            // Degenerate on these roots (possible only for δ_top / sum
            // forms): let the general invariant route decide.
            continue;
        }
        // The subdiscriminant criterion: D = δ(r)² is a G-invariant
        // ALGEBRAIC integer fixed by G_f ≤ G, hence a rational integer,
        // with |D| ≤ B² < p^k/2 — so its symmetric residue IS D exactly
        // (an impure v² would mean the certificate chain is broken).
        // δ(r) = ±√D, and δ(r) ∈ Z ⟺ D is a perfect square; δ(r) ∈ Z
        // is fixed by all of G_f while every g ∉ H sends δ to −δ ≠ δ
        // (D ≠ 0), so:  D square ⟺ G_f ≤ ker χ = H. Exact, and no
        // residue collision can ever arise (the resolvent has the two
        // distinct roots ±√D).
        const auto vsq = base.ring.integer_residue(base.ring.mul(v, v));
        if (!vsq.has_value()) {
            return fail<Out>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "quadratic_character_descent: G-invariant v² "
                           "is not a pure integer — certificate violated"});
        }
        const BigInt& d_val = *vsq;
        if (d_val.is_zero()) {
            // v ≢ 0 but v² ≡ 0: the value sits at positive p-valuation
            // (only the δ_top / sum forms can do this): general route.
            continue;
        }
        BigInt abs_d = d_val;
        if (abs_d.is_negative()) abs_d = -abs_d;
        if (abs_d > bound * bound) {
            return fail<Out>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "quadratic_character_descent: v² beyond the "
                           "Archimedean bound — certificate violated"});
        }
        return ok(Out{is_rational_square_q(Rational(d_val))});
    }
    return ok(Out{});  // no δ-form realises χ — general route
}

}  // namespace cas::algebra::galois_stauduhar
