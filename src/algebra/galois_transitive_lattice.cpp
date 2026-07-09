// A6 / CAS-L3-18 — Exhaustive transitive-subgroup lattice of S_n (n small).
//
// Enumerates ALL subgroups of S_n up to conjugacy by breadth-first
// single-generator extension starting from the trivial group, then filters
// the transitive ones. Every finite group is reachable this way (add one
// generator at a time), so completeness holds *by construction*: no
// transitive-group tables from the literature are consulted, eliminating
// transcription/hallucination risk (CLAUDE.md REGOLA 0.1). Classical counts
// (e.g. 5 transitive classes for n=5, 16 for n=6, 7 for n=7) are asserted
// in the unit tests as an independent cross-check, not used by the engine.
//
// Complexity: Θ(n!) memory per subgroup (dense bitset over Lehmer ranks)
// and roughly Θ(#classes · n! ) rank-compositions with the (H,H)-double-
// coset pruning ⟨H,s⟩ = ⟨H, h₁·s·h₂⟩. Practical for n ≤ 7; the caller-owned
// LatticeBudget bounds both degree and operations, failing structured —
// never truncating silently.

#include "perm_group_internal.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas::algebra::permgrp {

namespace {

// Dense bitset over the Lehmer-rank space [0, n!).
class RankSet {
public:
    explicit RankSet(std::size_t universe)
        : words_((universe + 63U) / 64U, 0ULL) {}
    [[nodiscard]] bool test(std::uint32_t r) const {
        return ((words_[r >> 6U] >> (r & 63U)) & 1ULL) != 0ULL;
    }
    // Returns true if the bit was newly set.
    bool set(std::uint32_t r) {
        const std::uint64_t mask = 1ULL << (r & 63U);
        const bool fresh = (words_[r >> 6U] & mask) == 0ULL;
        words_[r >> 6U] |= mask;
        return fresh;
    }

private:
    std::vector<std::uint64_t> words_;
};

// Cache of all n! permutations indexed by Lehmer rank, so rank-level
// composition avoids repeated unranking.
struct SymmetricGroupCache {
    std::size_t n{0U};
    std::vector<Perm> perms;  // perms[r] = lehmer_unrank(r, n)
};

[[nodiscard]] SymmetricGroupCache build_cache(std::size_t n) {
    SymmetricGroupCache c;
    c.n = n;
    const std::uint64_t order = factorial_u64(n);
    c.perms.reserve(static_cast<std::size_t>(order));
    for (std::uint64_t r = 0U; r < order; ++r) {
        c.perms.push_back(lehmer_unrank(static_cast<std::uint32_t>(r), n));
    }
    return c;
}

[[nodiscard]] std::uint32_t compose_ranks(const SymmetricGroupCache& c,
                                          std::uint32_t a, std::uint32_t b) {
    return lehmer_rank(compose(c.perms[a], c.perms[b]));
}

// Conjugation-invariant fingerprint used to avoid the exhaustive conjugacy
// scan for obviously-distinct subgroups: (order, sorted multiset of element
// cycle types).
[[nodiscard]] std::vector<std::size_t> signature_of(const PermGroup& g) {
    std::vector<std::size_t> sig;
    sig.push_back(static_cast<std::size_t>(g.order()));
    for (const auto& [ct, cnt] : g.cycle_type_distribution()) {
        sig.push_back(cnt);
        sig.push_back(ct.size());
        for (const std::size_t l : ct) sig.push_back(l);
    }
    return sig;
}

}  // namespace

bool conjugate_in_sn(const PermGroup& a, const PermGroup& b) {
    if (a.degree() != b.degree()) return false;
    if (a.order() != b.order()) return false;
    if (a.equals(b)) return true;
    if (signature_of(a) != signature_of(b)) return false;
    const std::size_t n = a.degree();
    const std::uint64_t order = factorial_u64(n);
    for (std::uint64_t sr = 0U; sr < order; ++sr) {
        const Perm sigma = lehmer_unrank(static_cast<std::uint32_t>(sr), n);
        const Perm sigma_inv = inverse(sigma);
        bool all_in = true;
        for (const auto& g : a.generators()) {
            const Perm conj = compose(sigma, compose(g, sigma_inv));
            if (!b.contains(conj)) {
                all_in = false;
                break;
            }
        }
        // Generators of A map into B and |A| = |B| ⇒ σAσ⁻¹ = B.
        if (all_in) return true;
    }
    return false;
}

Result<std::vector<PermGroup>> transitive_subgroup_classes(
    std::size_t n, const LatticeBudget& budget) {
    if (n < 1U || n > budget.max_degree) {
        return fail<std::vector<PermGroup>>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "transitive_subgroup_classes: degree exceeds the "
                       "configured dense-lattice cap (A6: degrees past the "
                       "cap require the Stauduhar maximal-subgroup descent "
                       "increment)"});
    }
    if (budget.max_ops == 0U) {
        return fail<std::vector<PermGroup>>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "transitive_subgroup_classes: caller must supply a "
                       "positive operations budget (no implicit unlimited "
                       "mode)"});
    }
    const std::uint64_t order = factorial_u64(n);
    std::uint64_t ops_left = budget.max_ops;
    auto spend = [&ops_left](std::uint64_t cost) -> bool {
        if (ops_left < cost) return false;
        ops_left -= cost;
        return true;
    };
    const SymmetricGroupCache cache = build_cache(n);

    auto budget_fail = []() {
        return fail<std::vector<PermGroup>>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "transitive_subgroup_classes: operations budget "
                       "exhausted — raise the caller-side CASContext budget "
                       "(A6 lattice enumeration)"});
    };

    // BFS over subgroup classes. Level 0: the trivial group.
    auto trivial = PermGroup::closure(n, {}, order);
    if (trivial.is_error())
        return fail<std::vector<PermGroup>>(trivial.error());
    std::vector<PermGroup> classes{std::move(trivial.value())};
    std::vector<std::vector<std::size_t>> class_sigs{signature_of(classes[0])};
    std::size_t next_to_expand = 0U;

    while (next_to_expand < classes.size()) {
        // NOTE: copy, not reference — `classes` reallocates inside the loop.
        const PermGroup h = classes[next_to_expand];
        ++next_to_expand;
        if (h.order() == order) continue;  // S_n itself: nothing above.

        // Elements of H as ranks, for double-coset marking.
        const std::vector<std::uint32_t>& h_ranks = h.element_ranks();
        RankSet visited(static_cast<std::size_t>(order));
        for (const std::uint32_t r : h_ranks) visited.set(r);

        for (std::uint32_t s = 0U; s < order; ++s) {
            if (visited.test(s)) continue;
            // Mark the whole double coset H·s·H: ⟨H, x⟩ is the same
            // subgroup for every x in it.
            for (const std::uint32_t h1 : h_ranks) {
                const std::uint32_t h1s = compose_ranks(cache, h1, s);
                for (const std::uint32_t h2 : h_ranks) {
                    visited.set(compose_ranks(cache, h1s, h2));
                }
            }
            if (!spend(static_cast<std::uint64_t>(h_ranks.size()) *
                       static_cast<std::uint64_t>(h_ranks.size()))) {
                return budget_fail();
            }

            // Extend: K = ⟨H, s⟩.
            std::vector<Perm> gens = h.generators();
            gens.push_back(cache.perms[s]);
            auto k_res = PermGroup::closure(n, std::move(gens), order);
            if (k_res.is_error())
                return fail<std::vector<PermGroup>>(k_res.error());
            PermGroup k = std::move(k_res.value());
            if (!spend(k.order())) return budget_fail();

            // Dedup up to conjugacy (signature pre-filter inside).
            bool known = false;
            const auto k_sig = signature_of(k);
            for (std::size_t i = 0U; i < classes.size(); ++i) {
                if (class_sigs[i] != k_sig) continue;
                if (!spend(order)) return budget_fail();
                if (conjugate_in_sn(classes[i], k)) {
                    known = true;
                    break;
                }
            }
            if (!known) {
                classes.push_back(std::move(k));
                class_sigs.push_back(k_sig);
            }
        }
    }

    std::vector<PermGroup> transitive;
    for (auto& g : classes) {
        if (g.is_transitive()) transitive.push_back(std::move(g));
    }
    std::sort(transitive.begin(), transitive.end(),
              [](const PermGroup& a, const PermGroup& b) {
                  return a.order() < b.order();
              });
    return ok(std::move(transitive));
}

}  // namespace cas::algebra::permgrp
