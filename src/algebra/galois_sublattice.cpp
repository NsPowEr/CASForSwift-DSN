// A6 Brick 3.5 — implementation of the maximal transitive-subgroup classes
// inside a small group H. See galois_sublattice_internal.hpp for the
// contract. The three optimisations of galois_transitive_lattice.cpp are
// reused in index space: a dense Cayley table over element indices,
// index-only closures (a BsgsGroup is materialised only for the surviving
// maximal classes), and an exact `seen` set holding the element-index set
// of every subgroup H-conjugate to an accepted class. Both budgets (ops =
// time, bytes = memory) are spent BEFORE the corresponding allocation or
// scan — the table is never allocated past the byte budget.

#include "galois_sublattice_internal.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include "perm_group_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cas::algebra::permgrp {

namespace {

struct PermHash {
    [[nodiscard]] std::size_t operator()(const Perm& p) const noexcept {
        std::size_t h = 1469598103934665603ULL;
        for (const std::uint8_t x : p) {
            h ^= static_cast<std::size_t>(x);
            h *= 1099511628211ULL;
        }
        return h;
    }
};

struct IdxVecHash {
    [[nodiscard]] std::size_t operator()(
        const std::vector<std::uint16_t>& v) const noexcept {
        std::size_t h = 1469598103934665603ULL;
        for (const std::uint16_t x : v) {
            h ^= static_cast<std::size_t>(x);
            h *= 1099511628211ULL;
        }
        return h;
    }
};

// Dense bitset over element indices [0, |H|).
class IdxSet {
public:
    explicit IdxSet(std::size_t universe)
        : words_((universe + 63U) / 64U, 0ULL) {}
    [[nodiscard]] bool test(std::size_t r) const {
        return ((words_[r >> 6U] >> (r & 63U)) & 1ULL) != 0ULL;
    }
    bool set(std::size_t r) {
        const std::uint64_t mask = 1ULL << (r & 63U);
        const bool fresh = (words_[r >> 6U] & mask) == 0ULL;
        words_[r >> 6U] |= mask;
        return fresh;
    }

private:
    std::vector<std::uint64_t> words_;
};

// Index-only closure of ⟨gen_idxs⟩ inside H via the Cayley table.
[[nodiscard]] std::vector<std::uint16_t> closure_idx(
    const std::vector<std::uint16_t>& mul, std::size_t m,
    const std::vector<std::uint16_t>& gen_idxs, std::uint16_t id_idx) {
    IdxSet seen(m);
    std::vector<std::uint16_t> elems{id_idx};
    seen.set(id_idx);
    for (std::size_t i = 0U; i < elems.size(); ++i) {
        const std::uint16_t cur = elems[i];
        for (const std::uint16_t g : gen_idxs) {
            const std::uint16_t nxt =
                mul[static_cast<std::size_t>(g) * m + cur];
            if (seen.set(nxt)) elems.push_back(nxt);
        }
    }
    std::sort(elems.begin(), elems.end());
    return elems;
}

// A discovered subgroup class: BFS generators and sorted element indices.
struct NodeClass {
    std::vector<std::uint16_t> gen_idxs;
    std::vector<std::uint16_t> elem_idxs;
};

}  // namespace

Result<std::vector<BsgsGroup>> transitive_subgroup_classes_in(
    const BsgsGroup& H, std::uint64_t max_ops, std::uint64_t max_bytes,
    symbolic::CASContext* ctx) {
    if (max_ops == 0U || max_bytes == 0U) {
        return fail<std::vector<BsgsGroup>>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "transitive_subgroup_classes_in: caller must supply "
                       "positive ops and byte budgets"});
    }
    const std::uint64_t order64 = H.order();
    if (order64 > 65536ULL) {
        return fail<std::vector<BsgsGroup>>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "transitive_subgroup_classes_in: |H| exceeds the "
                       "u16-index universe (2^16); such a node needs the "
                       "structural wreath-maximal route"});
    }
    const std::size_t m = static_cast<std::size_t>(order64);
    const std::size_t n = H.degree();
    std::uint64_t ops_left = max_ops;
    std::uint64_t bytes_left = max_bytes;
    auto spend = [](std::uint64_t& left, std::uint64_t cost) -> bool {
        if (left < cost) return false;
        left -= cost;
        return true;
    };
    auto ops_fail = []() {
        return fail<std::vector<BsgsGroup>>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "transitive_subgroup_classes_in: operations budget "
                       "exhausted — raise "
                       "CASContext::galois_lattice_max_ops"});
    };
    auto bytes_fail = []() {
        return fail<std::vector<BsgsGroup>>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "transitive_subgroup_classes_in: memory budget "
                       "exhausted — this node needs the structural "
                       "wreath-maximal route, or raise "
                       "CASContext::galois_sublattice_max_bytes"});
    };
    auto poll = [&]() -> bool {
        if (!ctx) return true;
        return !ctx->check_interrupt().is_error();
    };

    // ── enumerate H's elements (BFS closure) and index them ────────────────
    // Memory: m perms of n bytes each, twice (vector + hash key copies).
    if (!spend(bytes_left, 2ULL * order64 * n)) return bytes_fail();
    std::vector<Perm> elems{identity(n)};
    std::unordered_map<Perm, std::uint16_t, PermHash> index_of;
    {
        std::unordered_set<Perm, PermHash> seen_elems{elems[0]};
        for (std::size_t i = 0U; i < elems.size(); ++i) {
            for (const auto& g : H.generators()) {
                Perm nxt = compose(g, elems[i]);
                if (seen_elems.insert(nxt).second) {
                    elems.push_back(std::move(nxt));
                }
            }
            if (!spend(ops_left, H.generators().size())) return ops_fail();
        }
        if (elems.size() != m) {
            return fail<std::vector<BsgsGroup>>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "transitive_subgroup_classes_in: element "
                           "enumeration disagrees with the BSGS order"});
        }
        // Canonical (lexicographic) element order — indices are stable.
        std::sort(elems.begin(), elems.end());
        for (std::size_t i = 0U; i < m; ++i) {
            index_of.emplace(elems[i], static_cast<std::uint16_t>(i));
        }
    }
    const std::uint16_t id_idx = index_of.at(identity(n));

    // ── Cayley table over indices: mul[a·m + b] = idx(elems[a] ∘ elems[b]) ──
    // The byte spend precedes the allocation: an oversized node fails
    // structured here instead of attempting a multi-GiB table.
    if (!spend(bytes_left, 2ULL * order64 * order64 + 2ULL * order64)) {
        return bytes_fail();
    }
    if (!spend(ops_left, order64 * order64)) return ops_fail();
    std::vector<std::uint16_t> mul(m * m);
    std::vector<std::uint16_t> inv_idx(m);
    {
        Perm scratch(n);
        for (std::size_t a = 0U; a < m; ++a) {
            if (!poll()) return ops_fail();
            const Perm& pa = elems[a];
            for (std::size_t b = 0U; b < m; ++b) {
                const Perm& pb = elems[b];
                for (std::size_t i = 0U; i < n; ++i) scratch[i] = pa[pb[i]];
                mul[a * m + b] = index_of.at(scratch);
            }
            inv_idx[a] = index_of.at(inverse(pa));
        }
    }

    // ── seen: element-index set of every H-conjugate of accepted classes ───
    std::unordered_set<std::vector<std::uint16_t>, IdxVecHash> seen;
    bool byte_exhausted = false;
    auto register_conjugates =
        [&](const std::vector<std::uint16_t>& k_idxs) -> bool {
        if (!spend(ops_left,
                   order64 * static_cast<std::uint64_t>(k_idxs.size()))) {
            return false;
        }
        std::vector<std::uint16_t> conj(k_idxs.size());
        for (std::size_t s = 0U; s < m; ++s) {
            const std::size_t si = inv_idx[s];
            for (std::size_t i = 0U; i < k_idxs.size(); ++i) {
                const std::uint16_t sk = mul[s * m + k_idxs[i]];
                conj[i] = mul[static_cast<std::size_t>(sk) * m + si];
            }
            std::sort(conj.begin(), conj.end());
            if (seen.insert(conj).second &&
                !spend(bytes_left,
                       2ULL * static_cast<std::uint64_t>(conj.size()))) {
                byte_exhausted = true;
                return false;
            }
        }
        return true;
    };
    auto budget_fail = [&]() {
        return byte_exhausted ? bytes_fail() : ops_fail();
    };

    // ── BFS over subgroup classes from the trivial group ────────────────────
    std::vector<NodeClass> classes{
        NodeClass{{}, std::vector<std::uint16_t>{id_idx}}};
    if (!register_conjugates(classes[0].elem_idxs)) return budget_fail();
    std::size_t next_to_expand = 0U;
    while (next_to_expand < classes.size()) {
        // Copy: `classes` reallocates inside the loop.
        const NodeClass h = classes[next_to_expand];
        ++next_to_expand;
        if (h.elem_idxs.size() == m) continue;  // H itself: nothing above
        IdxSet visited(m);
        for (const std::uint16_t r : h.elem_idxs) visited.set(r);
        for (std::size_t s = 0U; s < m; ++s) {
            if (visited.test(s)) continue;
            if (!poll()) return ops_fail();
            // Mark the double coset K·s·K: ⟨K, x⟩ is the same subgroup for
            // every x in it.
            for (const std::uint16_t h1 : h.elem_idxs) {
                const std::uint16_t h1s =
                    mul[static_cast<std::size_t>(h1) * m + s];
                for (const std::uint16_t h2 : h.elem_idxs) {
                    visited.set(mul[static_cast<std::size_t>(h1s) * m + h2]);
                }
            }
            if (!spend(ops_left,
                       static_cast<std::uint64_t>(h.elem_idxs.size()) *
                           h.elem_idxs.size())) {
                return ops_fail();
            }
            std::vector<std::uint16_t> k_gens = h.gen_idxs;
            k_gens.push_back(static_cast<std::uint16_t>(s));
            std::vector<std::uint16_t> k_elems =
                closure_idx(mul, m, k_gens, id_idx);
            if (!spend(ops_left, k_elems.size())) return ops_fail();
            if (seen.find(k_elems) == seen.end()) {
                if (!register_conjugates(k_elems)) return budget_fail();
                classes.push_back(
                    NodeClass{std::move(k_gens), std::move(k_elems)});
            }
        }
    }

    // ── keep the PROPER transitive classes, largest first ──────────────────
    std::vector<NodeClass> trans;
    for (auto& c : classes) {
        if (c.elem_idxs.size() == m) continue;  // proper only
        // Transitivity: orbit of 0 under the class's generators.
        std::vector<bool> orb(n, false);
        orb[0] = true;
        std::vector<std::size_t> queue{0U};
        std::size_t reached = 1U;
        while (!queue.empty()) {
            const std::size_t pt = queue.back();
            queue.pop_back();
            for (const std::uint16_t g : c.gen_idxs) {
                const std::size_t img = elems[g][pt];
                if (!orb[img]) {
                    orb[img] = true;
                    ++reached;
                    queue.push_back(img);
                }
            }
        }
        if (reached == n) trans.push_back(std::move(c));
    }
    std::sort(trans.begin(), trans.end(), [](const auto& a, const auto& b) {
        return a.elem_idxs.size() > b.elem_idxs.size();
    });

    // ── maximality filter: drop K contained in a conjugate of a bigger
    //    proper transitive class (any maximal subgroup above a transitive
    //    K is itself transitive, so bigger TRANSITIVE classes suffice) ────
    std::vector<bool> dominated(trans.size(), false);
    for (std::size_t big = 0U; big < trans.size(); ++big) {
        if (dominated[big]) continue;  // domination is transitive
        IdxSet in_big(m);
        for (const std::uint16_t r : trans[big].elem_idxs) in_big.set(r);
        for (std::size_t small = big + 1U; small < trans.size(); ++small) {
            if (dominated[small]) continue;
            const std::size_t osm = trans[small].elem_idxs.size();
            const std::size_t obg = trans[big].elem_idxs.size();
            if (osm == obg || obg % osm != 0U) continue;  // Lagrange
            if (!poll()) return ops_fail();
            if (!spend(ops_left,
                       order64 * trans[small].gen_idxs.size())) {
                return ops_fail();
            }
            // K ⊆ s·BIG·s⁻¹ ⟺ every generator g of K has s⁻¹gs ∈ BIG.
            for (std::size_t s = 0U; s < m; ++s) {
                bool all_in = true;
                for (const std::uint16_t g : trans[small].gen_idxs) {
                    const std::uint16_t gs =
                        mul[static_cast<std::size_t>(g) * m + s];
                    const std::uint16_t sgs = mul[
                        static_cast<std::size_t>(inv_idx[s]) * m + gs];
                    if (!in_big.test(sgs)) {
                        all_in = false;
                        break;
                    }
                }
                if (all_in) {
                    dominated[small] = true;
                    break;
                }
            }
        }
    }

    std::vector<BsgsGroup> out;
    for (std::size_t i = 0U; i < trans.size(); ++i) {
        if (dominated[i]) continue;
        std::vector<Perm> gens;
        gens.reserve(trans[i].gen_idxs.size());
        for (const std::uint16_t g : trans[i].gen_idxs) {
            gens.push_back(elems[g]);
        }
        auto b = BsgsGroup::build(n, std::move(gens));
        if (b.is_error()) return fail<std::vector<BsgsGroup>>(b.error());
        if (b.value().order() != trans[i].elem_idxs.size()) {
            return fail<std::vector<BsgsGroup>>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "transitive_subgroup_classes_in: BSGS order "
                           "disagrees with the index closure"});
        }
        out.push_back(std::move(b.value()));
    }
    return ok(std::move(out));
}

}  // namespace cas::algebra::permgrp
