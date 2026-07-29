// A6 Brick 3b — implementation of the certified relative invariants.
// See galois_invariant_internal.hpp for the contract and the two-tier
// candidate strategy (k-subset orbit sums, then j-tuple orbit sums of
// growing length, whose last level is the guaranteed Galois resolvent
// monomial).

#include "galois_invariant_internal.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include "perm_group_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>

namespace cas::algebra::galois_invariant {

namespace {

using permgrp::BsgsGroup;
using permgrp::compose;
using permgrp::identity;
using permgrp::Perm;

// x_i ↦ x_{σ(i)}: the image exponent vector has out[σ(i)] = m[i].
[[nodiscard]] Monomial apply_perm(const Perm& sigma, const Monomial& m) {
    Monomial out(m.size(), 0U);
    for (std::size_t i = 0U; i < m.size(); ++i) out[sigma[i]] = m[i];
    return out;
}

// Budget poll: one unit per orbit-node expansion (anti-runaway belt only).
[[nodiscard]] Result<void> spend(std::uint64_t& ops, std::uint64_t max_ops,
                                 symbolic::CASContext* ctx) {
    if (ctx) {
        if (auto chk = ctx->check_interrupt(); chk.is_error()) return chk;
    }
    if (++ops > max_ops) {
        return fail<void>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "relative_invariant: ops budget exhausted (raise "
                       "the configured max_ops)"});
    }
    return ok();
}

// The H-orbit of a monomial as a sorted set (BFS over the generators).
[[nodiscard]] Result<std::set<Monomial>> monomial_orbit(
    const std::vector<Perm>& gens, const Monomial& m0, std::uint64_t& ops,
    std::uint64_t max_ops, symbolic::CASContext* ctx) {
    std::set<Monomial> orbit{m0};
    std::vector<Monomial> queue{m0};
    while (!queue.empty()) {
        const Monomial cur = std::move(queue.back());
        queue.pop_back();
        for (const auto& g : gens) {
            if (auto s = spend(ops, max_ops, ctx); s.is_error()) {
                return fail<std::set<Monomial>>(s.error());
            }
            Monomial img = apply_perm(g, cur);
            if (orbit.insert(img).second) queue.push_back(std::move(img));
        }
    }
    return ok(std::move(orbit));
}

// Enumerates the G-orbit of the invariant F0. Returns the coset
// representatives iff the orbit size equals `expected` ([G:H]); nullopt
// iff the orbit is strictly smaller (stabiliser bigger than H — candidate
// rejected). An orbit exceeding `expected` is mathematically impossible
// (Stab_G(F0) ⊇ H) and trips an InternalError.
[[nodiscard]] Result<std::optional<std::vector<Perm>>> g_orbit_transversal(
    const BsgsGroup& G, const std::vector<Monomial>& F0,
    std::uint64_t expected, std::uint64_t& ops, std::uint64_t max_ops,
    symbolic::CASContext* ctx) {
    std::map<std::vector<Monomial>, Perm> visited;
    const Perm id = identity(G.degree());
    visited.emplace(F0, id);
    std::vector<const std::vector<Monomial>*> queue{&visited.begin()->first};
    std::vector<Perm> reps{id};
    while (!queue.empty()) {
        const std::vector<Monomial>* cur = queue.back();
        queue.pop_back();
        const Perm rep = visited.at(*cur);  // copy: map may rehash below
        for (const auto& g : G.generators()) {
            if (auto s = spend(ops, max_ops, ctx); s.is_error()) {
                return fail<std::optional<std::vector<Perm>>>(s.error());
            }
            std::vector<Monomial> img = apply_perm_to_invariant(g, *cur);
            auto [it, inserted] =
                visited.emplace(std::move(img), compose(g, rep));
            if (inserted) {
                if (visited.size() > expected) {
                    return fail<std::optional<std::vector<Perm>>>(CASError{
                        .kind = CASErrorKind::InternalError,
                        .message = "relative_invariant: G-orbit of an "
                                   "H-invariant exceeded [G:H]"});
                }
                reps.push_back(it->second);
                queue.push_back(&it->first);
            }
        }
    }
    if (visited.size() < expected) {
        return ok(std::optional<std::vector<Perm>>{});  // Stab_G(F) ⊋ H
    }
    return ok(std::optional<std::vector<Perm>>{std::move(reps)});
}

}  // namespace

std::vector<Monomial> apply_perm_to_invariant(
    const Perm& sigma, const std::vector<Monomial>& monomials) {
    std::vector<Monomial> out;
    out.reserve(monomials.size());
    for (const auto& m : monomials) out.push_back(apply_perm(sigma, m));
    std::sort(out.begin(), out.end());
    return out;
}

Result<RelativeInvariant> relative_invariant(const BsgsGroup& G,
                                             const BsgsGroup& H,
                                             std::uint64_t max_ops,
                                             symbolic::CASContext* ctx) {
    const std::size_t n = G.degree();
    if (H.degree() != n) {
        return fail<RelativeInvariant>(
            CASError{.kind = CASErrorKind::InvalidArgument,
                     .message = "relative_invariant: degree mismatch"});
    }
    for (const auto& h : H.generators()) {
        if (!G.contains(h)) {
            return fail<RelativeInvariant>(
                CASError{.kind = CASErrorKind::InvalidArgument,
                         .message = "relative_invariant: H is not a "
                                    "subgroup of G"});
        }
    }
    const std::uint64_t og = G.order();
    const std::uint64_t oh = H.order();
    if (oh == 0U || og % oh != 0U) {
        return fail<RelativeInvariant>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "relative_invariant: |H| does not divide |G| "
                       "(Lagrange violated — corrupt BSGS)"});
    }
    const std::uint64_t index = og / oh;
    if (index == 1U) {
        return fail<RelativeInvariant>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "relative_invariant: H = G (index 1, nothing to "
                       "descend to)"});
    }
    std::uint64_t ops = 0U;

    // ── tier 1: k-subset orbit sums, one candidate per H-orbit ─────────────
    for (std::size_t k = 1U; k < n; ++k) {
        std::set<Monomial> covered;
        // Canonical k-subset sweep (lex on indicator vectors).
        std::vector<std::size_t> idx(k);
        for (std::size_t i = 0U; i < k; ++i) idx[i] = i;
        const std::size_t all = [&] {
            // C(n, k) — exact, small (n ≤ 20 in this engine).
            std::size_t r = 1U;
            for (std::size_t i = 0U; i < k; ++i) r = r * (n - i) / (i + 1U);
            return r;
        }();
        while (true) {
            Monomial seed(n, 0U);
            for (const std::size_t i : idx) seed[i] = 1U;
            if (covered.find(seed) == covered.end()) {
                auto orb =
                    monomial_orbit(H.generators(), seed, ops, max_ops, ctx);
                if (orb.is_error()) {
                    return fail<RelativeInvariant>(orb.error());
                }
                covered.insert(orb.value().begin(), orb.value().end());
                // The full orbit is the symmetric e_k — its stabiliser is
                // all of G; only proper orbits can separate H.
                if (orb.value().size() < all) {
                    std::vector<Monomial> F0(orb.value().begin(),
                                             orb.value().end());
                    auto tr = g_orbit_transversal(G, F0, index, ops,
                                                  max_ops, ctx);
                    if (tr.is_error()) {
                        return fail<RelativeInvariant>(tr.error());
                    }
                    if (tr.value().has_value()) {
                        return ok(RelativeInvariant{std::move(F0), k,
                                                    std::move(*tr.value())});
                    }
                }
            }
            // Next k-subset in lex order (position i may hold at most
            // n − k + i).
            bool advanced = false;
            std::size_t i = k;
            while (i-- > 0U) {
                if (idx[i] < n - k + i) {
                    ++idx[i];
                    for (std::size_t j = i + 1U; j < k; ++j) {
                        idx[j] = idx[j - 1U] + 1U;
                    }
                    advanced = true;
                    break;
                }
            }
            if (!advanced) break;
        }
    }

    // ── tier 2: j-tuple orbit sums, j = 2 .. n−1 ───────────────────────────
    // Seeds x_{t₀}¹·x_{t₁}²⋯x_{t_{j−1}}^j over the ordered j-tuples of
    // distinct indices, ascending j: k-subsets (tier 1) cannot separate a
    // subgroup whose set-orbits coincide with G's (e.g. the sign-character
    // kernels inside a wreath node), while a partial tuple already breaks
    // the symmetry at the smallest j where the point stabilisers differ —
    // with an orbit-sum FAR smaller than the full Galois monomial (the
    // p-adic precision bound of the later resolvent grows with both the
    // degree and the coefficient size of F, so invariant size is resolvent
    // cost). j = n−1 IS the classical Galois resolvent monomial: its
    // S_n-stabiliser is trivial, so Stab_G(Σ_{h∈H} h·m*) = H
    // unconditionally — the tier is a guaranteed terminator.
    for (std::size_t j = 2U; j < n; ++j) {
        std::uint64_t tuple_count = 1U;
        for (std::size_t i = 0U; i < j; ++i) tuple_count *= (n - i);
        std::set<Monomial> covered;
        std::vector<std::size_t> t(j);
        for (std::size_t i = 0U; i < j; ++i) t[i] = i;
        while (true) {
            Monomial seed(n, 0U);
            for (std::size_t i = 0U; i < j; ++i) {
                seed[t[i]] = static_cast<std::uint8_t>(i + 1U);
            }
            if (covered.find(seed) == covered.end()) {
                auto orb =
                    monomial_orbit(H.generators(), seed, ops, max_ops, ctx);
                if (orb.is_error()) return fail<RelativeInvariant>(orb.error());
                covered.insert(orb.value().begin(), orb.value().end());
                // Full tuple orbit ⇒ the sum is symmetric (Stab ⊇ S_n on
                // the tuple pattern): it cannot separate H from G.
                if (static_cast<std::uint64_t>(orb.value().size()) <
                    tuple_count) {
                    std::vector<Monomial> F0(orb.value().begin(),
                                             orb.value().end());
                    auto tr = g_orbit_transversal(G, F0, index, ops,
                                                  max_ops, ctx);
                    if (tr.is_error()) return fail<RelativeInvariant>(tr.error());
                    if (tr.value().has_value()) {
                        return ok(RelativeInvariant{std::move(F0),
                                                    j * (j + 1U) / 2U,
                                                    std::move(*tr.value())});
                    }
                }
            }
            // Next ordered tuple of distinct indices in lex order.
            bool advanced = false;
            std::size_t i = j;
            while (i-- > 0U) {
                std::vector<bool> used(n, false);
                for (std::size_t a = 0U; a < i; ++a) used[t[a]] = true;
                std::size_t v = t[i] + 1U;
                while (v < n && used[v]) ++v;
                if (v < n) {
                    t[i] = v;
                    used[v] = true;
                    std::size_t w = 0U;
                    for (std::size_t a = i + 1U; a < j; ++a) {
                        while (used[w]) ++w;
                        t[a] = w;
                        used[w] = true;
                    }
                    advanced = true;
                    break;
                }
            }
            if (!advanced) break;
        }
    }
    // Unreachable: the j = n−1 pass contains the Galois resolvent
    // monomial, whose exactness certificate is a theorem.
    return fail<RelativeInvariant>(CASError{
        .kind = CASErrorKind::InternalError,
        .message = "relative_invariant: Galois resolvent monomial failed "
                   "the exactness certificate — impossible"});
}

}  // namespace cas::algebra::galois_invariant
