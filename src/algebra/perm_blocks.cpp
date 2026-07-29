// A6 — Block-system computation (see perm_blocks_internal.hpp).
// Atkinson 1975 union-find closure; exact, generator-level, no element
// enumeration. Conventions from perm_group_internal.hpp: point action
// pt^g = g[pt].

#include "perm_blocks_internal.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace cas::algebra::permgrp {

namespace {

// Path-halving union-find over {0..n-1}.
class UnionFind {
public:
    explicit UnionFind(std::size_t n) : parent_(n) {
        for (std::size_t i = 0U; i < n; ++i) parent_[i] = i;
    }
    [[nodiscard]] std::size_t find(std::size_t x) {
        while (parent_[x] != x) {
            parent_[x] = parent_[parent_[x]];
            x = parent_[x];
        }
        return x;
    }
    // Returns false if already in the same class.
    bool unite(std::size_t a, std::size_t b) {
        a = find(a);
        b = find(b);
        if (a == b) return false;
        if (a > b) std::swap(a, b);  // keep the smaller root (0 stays a root)
        parent_[b] = a;
        return true;
    }

private:
    std::vector<std::size_t> parent_;
};

// Finest G-congruence on {0..n-1} identifying {0, beta} (Atkinson): merge
// 0~beta, then propagate: whenever classes of u and v merge, also merge the
// classes of g[u] and g[v] for every generator g. The worklist holds merged
// pairs still to propagate; each successful merge decreases the class count,
// so at most n-1 merges occur and termination is guaranteed.
[[nodiscard]] BlockSystem congruence_closure(std::size_t n,
                                             const std::vector<Perm>& gens,
                                             std::size_t beta) {
    UnionFind uf(n);
    std::vector<std::pair<std::size_t, std::size_t>> work;
    uf.unite(0U, beta);
    work.emplace_back(0U, beta);
    while (!work.empty()) {
        const auto [u, v] = work.back();
        work.pop_back();
        for (const Perm& g : gens) {
            const std::size_t gu = g[u];
            const std::size_t gv = g[v];
            if (uf.unite(gu, gv)) work.emplace_back(gu, gv);
        }
    }
    BlockSystem sys;
    sys.block_of.assign(n, n);  // n = "unassigned" sentinel
    std::size_t next_id = 0U;
    std::vector<std::size_t> id_of_root(n, n);
    for (std::size_t p = 0U; p < n; ++p) {
        const std::size_t r = uf.find(p);
        if (id_of_root[r] == n) id_of_root[r] = next_id++;
        sys.block_of[p] = id_of_root[r];
    }
    sys.num_blocks = next_id;
    return sys;
}

[[nodiscard]] Result<void> validate_transitive_gens(
    std::size_t n, const std::vector<Perm>& gens) {
    if (n == 0U) {
        return fail<void>(
            CASError{.kind = CASErrorKind::InvalidArgument,
                     .message = "block systems: degree must be positive"});
    }
    for (const Perm& g : gens) {
        if (g.size() != n || !is_valid_perm(g)) {
            return fail<void>(CASError{
                .kind = CASErrorKind::InvalidArgument,
                .message = "block systems: generator is not a valid "
                           "permutation of the requested degree"});
        }
    }
    // Transitivity: orbit of 0 under gens covers all points.
    std::vector<bool> seen(n, false);
    std::vector<std::size_t> stack{0U};
    seen[0] = true;
    std::size_t count = 1U;
    while (!stack.empty()) {
        const std::size_t pt = stack.back();
        stack.pop_back();
        for (const Perm& g : gens) {
            const std::size_t im = g[pt];
            if (!seen[im]) {
                seen[im] = true;
                ++count;
                stack.push_back(im);
            }
        }
    }
    if (count != n) {
        return fail<void>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "block systems: group is intransitive — block theory "
                       "here is defined for transitive actions only"});
    }
    return ok();
}

}  // namespace

Result<std::vector<BlockSystem>> minimal_block_systems(
    std::size_t n, const std::vector<Perm>& gens) {
    auto valid = validate_transitive_gens(n, gens);
    if (valid.is_error()) return fail<std::vector<BlockSystem>>(valid.error());

    // f(β) = finest system identifying {0, β}. Every minimal system arises as
    // f(β) for each non-zero β of its 0-block, so scanning all β finds them
    // all; f(β) is minimal iff f(γ) == f(β) for every other γ in its 0-block
    // (any such f(γ) refines f(β) by minimality of the closure).
    std::vector<BlockSystem> all(n);  // all[β] used for β ≥ 1
    for (std::size_t beta = 1U; beta < n; ++beta) {
        all[beta] = congruence_closure(n, gens, beta);
    }
    std::vector<BlockSystem> out;
    for (std::size_t beta = 1U; beta < n; ++beta) {
        const BlockSystem& sys = all[beta];
        if (sys.num_blocks <= 1U) continue;  // trivial coarse partition
        bool minimal = true;
        bool first_rep = true;
        for (std::size_t gamma = 1U; gamma < n && minimal; ++gamma) {
            if (sys.block_of[gamma] != 0U) continue;  // not in the 0-block
            if (gamma < beta) { first_rep = false; break; }  // dedup: keep smallest rep
            if (!(all[gamma] == sys)) minimal = false;
        }
        if (minimal && first_rep) out.push_back(sys);
    }
    return ok(std::move(out));
}

Result<bool> is_primitive(std::size_t n, const std::vector<Perm>& gens) {
    auto valid = validate_transitive_gens(n, gens);
    if (valid.is_error()) return fail<bool>(valid.error());
    if (n == 1U) return ok(true);
    for (std::size_t beta = 1U; beta < n; ++beta) {
        if (congruence_closure(n, gens, beta).num_blocks > 1U) return ok(false);
    }
    return ok(true);
}

Result<std::vector<Perm>> block_action_gens(std::size_t n,
                                            const std::vector<Perm>& gens,
                                            const BlockSystem& sys) {
    if (sys.block_of.size() != n || sys.num_blocks == 0U) {
        return fail<std::vector<Perm>>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "block_action_gens: block system does not match the "
                       "degree"});
    }
    // Representative point of each block (first occurrence).
    std::vector<std::size_t> rep(sys.num_blocks, n);
    for (std::size_t p = 0U; p < n; ++p) {
        const std::size_t b = sys.block_of[p];
        if (b >= sys.num_blocks) {
            return fail<std::vector<Perm>>(
                CASError{.kind = CASErrorKind::InvalidArgument,
                         .message = "block_action_gens: block index out of "
                                    "range"});
        }
        if (rep[b] == n) rep[b] = p;
    }
    std::vector<Perm> out;
    out.reserve(gens.size());
    for (const Perm& g : gens) {
        Perm img(sys.num_blocks);
        for (std::size_t b = 0U; b < sys.num_blocks; ++b) {
            img[b] = static_cast<std::uint8_t>(sys.block_of[g[rep[b]]]);
        }
        // Invariance check: every point of block b must land in block img[b].
        for (std::size_t p = 0U; p < n; ++p) {
            if (sys.block_of[g[p]] !=
                static_cast<std::size_t>(img[sys.block_of[p]])) {
                return fail<std::vector<Perm>>(CASError{
                    .kind = CASErrorKind::InvalidArgument,
                    .message = "block_action_gens: partition is not invariant "
                               "under a generator (not a block system)"});
            }
        }
        if (!is_valid_perm(img)) {
            return fail<std::vector<Perm>>(CASError{
                .kind = CASErrorKind::InvalidArgument,
                .message = "block_action_gens: induced block map is not a "
                           "permutation (not a block system)"});
        }
        out.push_back(std::move(img));
    }
    return ok(std::move(out));
}

}  // namespace cas::algebra::permgrp
