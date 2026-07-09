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
// Three things keep it practical (n = 7 runs in ~23 s under ASan, n = 6 in
// ~0.5 s — down from 740 s / 11 s before these were added):
//   • a dense Cayley table over Lehmer ranks, built without allocating a Perm
//     per entry via a 3-bits-per-image perfect-hash rank LUT;
//   • a rank-only closure, so a candidate ⟨H,s⟩ is tested for novelty before
//     any PermGroup is materialised (n=7: 16237 candidates, 96 materialised);
//   • an exact `seen` set holding the element set of every subgroup conjugate
//     to an accepted class, turning the per-candidate O(n!) conjugacy scan
//     (74% of the old runtime) into one hash lookup.
//
// (H,H)-double-coset pruning ⟨H,s⟩ = ⟨H, h₁·s·h₂⟩ bounds the candidate count.
// Memory is Θ(n!) per subgroup (dense bitset over Lehmer ranks). The caller-
// owned LatticeBudget bounds both degree and operations, failing structured —
// never truncating silently.

#include "perm_group_internal.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <unordered_set>
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

// Dense Cayley table of S_n over Lehmer ranks: mul[a·order + b] = rank(σ_a ∘ σ_b).
// The double-coset marking below performs Θ(n!·|H|) rank compositions summed
// over the whole lattice; recomputing lehmer_rank (Θ(n²)) inside that loop is
// the dominant cost. Precomputing the table trades one Θ((n!)²·n²) pass and
// (n!)²·2 bytes of memory (n=7 → 50 MB, the enumeration cap) for O(1) lookups.
// Ranks fit in uint16 while n! ≤ 65536 (n ≤ 8), but the table is only built
// when (n!)² stays within kMaxTableEntries — n = 7 (50 MB) is the largest that
// qualifies; past that the enumeration falls back to on-the-fly ranking (and
// is bounded by the caller's ops budget anyway).
class CayleyTable {
public:
    static constexpr std::size_t kMaxTableEntries =
        static_cast<std::size_t>(5040U) * 5040U;

    explicit CayleyTable(const SymmetricGroupCache& c)
        : order_(static_cast<std::size_t>(factorial_u64(c.n))) {
        if (order_ * order_ > kMaxTableEntries) return;  // stay unbuilt
        const std::size_t n = c.n;
        // Perfect-hash rank lookup: a permutation of n ≤ 8 points packs into
        // 3 bits per image, so key(p) = Σ p[i]·8^i indexes a dense table of
        // its Lehmer rank. Using it here removes both the Θ(n²) lehmer_rank
        // and the heap-allocating `compose` (Perm is a std::vector) from the
        // (n!)²-iteration build below — together they dominated it.
        auto pack_key = [n](const Perm& p) {
            std::size_t key = 0U;
            for (std::size_t i = 0U; i < n; ++i) {
                key |= static_cast<std::size_t>(p[i]) << (3U * i);
            }
            return key;
        };
        std::vector<std::uint16_t> rank_lut(std::size_t{1} << (3U * n), 0U);
        for (std::size_t r = 0U; r < order_; ++r) {
            rank_lut[pack_key(c.perms[r])] = static_cast<std::uint16_t>(r);
        }
        mul_.resize(order_ * order_);
        for (std::size_t a = 0U; a < order_; ++a) {
            const Perm& pa = c.perms[a];
            for (std::size_t b = 0U; b < order_; ++b) {
                const Perm& pb = c.perms[b];
                // (pa ∘ pb)[i] = pa[pb[i]] — packed directly, no allocation.
                std::size_t key = 0U;
                for (std::size_t i = 0U; i < n; ++i) {
                    key |= static_cast<std::size_t>(pa[pb[i]]) << (3U * i);
                }
                mul_[a * order_ + b] = rank_lut[key];
            }
        }
    }
    [[nodiscard]] bool built() const { return !mul_.empty(); }
    [[nodiscard]] std::uint32_t compose_rank(const SymmetricGroupCache& c,
                                             std::uint32_t a,
                                             std::uint32_t b) const {
        if (mul_.empty()) return compose_ranks(c, a, b);
        return mul_[static_cast<std::size_t>(a) * order_ + b];
    }

private:
    std::size_t order_;
    std::vector<std::uint16_t> mul_;
};

// Rank-only closure of ⟨gen_ranks⟩ inside S_n: BFS by left-multiplication using
// the Cayley table and a dense bitset, allocating no Perm objects. Returns the
// sorted element ranks. This is the same set PermGroup::closure would produce
// (identical BFS), but it lets the lattice test a candidate for "already known"
// before paying for a full PermGroup — only genuinely new classes are
// materialised.
[[nodiscard]] std::vector<std::uint32_t> closure_ranks(
    const CayleyTable& cayley, const SymmetricGroupCache& cache,
    const std::vector<std::uint32_t>& gen_ranks, std::size_t order,
    std::uint32_t id_rank) {
    RankSet seen(order);
    std::vector<std::uint32_t> elems;
    elems.push_back(id_rank);
    seen.set(id_rank);
    for (std::size_t i = 0U; i < elems.size(); ++i) {
        const std::uint32_t cur = elems[i];
        for (const std::uint32_t g : gen_ranks) {
            const std::uint32_t nxt = cayley.compose_rank(cache, g, cur);
            if (seen.set(nxt)) elems.push_back(nxt);
        }
    }
    std::sort(elems.begin(), elems.end());
    return elems;
}

// Hash for a sorted element-rank vector (a subgroup, identified by its set).
struct RankVecHash {
    [[nodiscard]] std::size_t operator()(
        const std::vector<std::uint32_t>& v) const noexcept {
        std::size_t h = 1469598103934665603ULL;
        for (const std::uint32_t x : v) {
            h ^= static_cast<std::size_t>(x);
            h *= 1099511628211ULL;
        }
        return h;
    }
};

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
    const CayleyTable cayley(cache);

    auto budget_fail = []() {
        return fail<std::vector<PermGroup>>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "transitive_subgroup_classes: operations budget "
                       "exhausted — raise the caller-side CASContext budget "
                       "(A6 lattice enumeration)"});
    };

    const std::uint32_t id_rank = lehmer_rank(identity(n));

    // Rank of σ⁻¹ for every σ, so conjugation of an element rank r by σ is the
    // two Cayley lookups  σ·r·σ⁻¹.
    std::vector<std::uint32_t> inv_rank(static_cast<std::size_t>(order));
    for (std::uint64_t r = 0U; r < order; ++r) {
        inv_rank[static_cast<std::size_t>(r)] =
            lehmer_rank(inverse(cache.perms[static_cast<std::size_t>(r)]));
    }

    // `seen` holds the element-rank set of EVERY subgroup conjugate to an
    // already-accepted class. A candidate K is then known ⟺ its (sorted)
    // element-rank set is in `seen` — an exact O(|K|) hash lookup with full
    // vector comparison, replacing the per-candidate O(n!) conjugacy scan that
    // dominated the enumeration (74% of the n=6 run). Registering a class costs
    // n!·|K| Cayley lookups once; membership is then free for every later
    // candidate. Memory: Σ_classes (n!/|N(K)|)·|K| ≤ #classes·n! ranks.
    std::unordered_set<std::vector<std::uint32_t>, RankVecHash> seen;
    auto register_conjugates = [&](const PermGroup& k) {
        const std::vector<std::uint32_t>& ranks = k.element_ranks();
        std::vector<std::uint32_t> conj(ranks.size());
        for (std::uint64_t sr = 0U; sr < order; ++sr) {
            const std::uint32_t s_r = static_cast<std::uint32_t>(sr);
            const std::uint32_t si = inv_rank[static_cast<std::size_t>(sr)];
            for (std::size_t i = 0U; i < ranks.size(); ++i) {
                conj[i] = cayley.compose_rank(
                    cache, cayley.compose_rank(cache, s_r, ranks[i]), si);
            }
            std::sort(conj.begin(), conj.end());
            seen.insert(conj);
        }
    };

    // BFS over subgroup classes. Level 0: the trivial group.
    auto trivial = PermGroup::closure(n, {}, order);
    if (trivial.is_error())
        return fail<std::vector<PermGroup>>(trivial.error());
    std::vector<PermGroup> classes{std::move(trivial.value())};
    register_conjugates(classes[0]);
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
        std::vector<std::uint32_t> h_gen_ranks;
        h_gen_ranks.reserve(h.generators().size() + 1U);
        for (const auto& g : h.generators()) h_gen_ranks.push_back(lehmer_rank(g));

        for (std::uint32_t s = 0U; s < order; ++s) {
            if (visited.test(s)) continue;
            // Mark the whole double coset H·s·H: ⟨H, x⟩ is the same
            // subgroup for every x in it.
            for (const std::uint32_t h1 : h_ranks) {
                const std::uint32_t h1s = cayley.compose_rank(cache, h1, s);
                for (const std::uint32_t h2 : h_ranks) {
                    visited.set(cayley.compose_rank(cache, h1s, h2));
                }
            }
            if (!spend(static_cast<std::uint64_t>(h_ranks.size()) *
                       static_cast<std::uint64_t>(h_ranks.size()))) {
                return budget_fail();
            }

            // Extend: K = ⟨H, s⟩. Take the rank-only closure first — the vast
            // majority of candidates re-generate an already-known class, and
            // those never need a materialised PermGroup.
            std::vector<std::uint32_t> k_gen_ranks = h_gen_ranks;
            k_gen_ranks.push_back(s);
            const std::vector<std::uint32_t> k_ranks =
                closure_ranks(cayley, cache, k_gen_ranks,
                              static_cast<std::size_t>(order), id_rank);
            if (!spend(k_ranks.size())) return budget_fail();

            // Dedup up to conjugacy: exact set membership in `seen`.
            const bool known = seen.find(k_ranks) != seen.end();
            if (!known) {
                std::vector<Perm> gens = h.generators();
                gens.push_back(cache.perms[s]);
                auto k_res = PermGroup::closure(n, std::move(gens), order);
                if (k_res.is_error())
                    return fail<std::vector<PermGroup>>(k_res.error());
                PermGroup k = std::move(k_res.value());
                if (!spend(order * k.order())) return budget_fail();
                register_conjugates(k);
                classes.push_back(std::move(k));
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
