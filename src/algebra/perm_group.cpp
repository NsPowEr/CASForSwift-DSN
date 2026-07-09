// A6 / CAS-L3-18 — Exact permutation-group engine (see perm_group_internal.hpp).

#include "perm_group_internal.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"

#include <algorithm>
#include <cassert>
#include <functional>
#include <cstddef>
#include <cstdint>
#include <map>
#include <queue>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cas::algebra::permgrp {

Perm identity(std::size_t n) {
    Perm p(n);
    for (std::size_t i = 0U; i < n; ++i) p[i] = static_cast<std::uint8_t>(i);
    return p;
}

Perm compose(const Perm& a, const Perm& b) {
    assert(a.size() == b.size());
    Perm r(a.size());
    for (std::size_t i = 0U; i < a.size(); ++i) r[i] = a[b[i]];
    return r;
}

Perm inverse(const Perm& a) {
    Perm r(a.size());
    for (std::size_t i = 0U; i < a.size(); ++i) r[a[i]] = static_cast<std::uint8_t>(i);
    return r;
}

bool is_valid_perm(const Perm& a) {
    std::vector<bool> seen(a.size(), false);
    for (std::size_t i = 0U; i < a.size(); ++i) {
        if (a[i] >= a.size() || seen[a[i]]) return false;
        seen[a[i]] = true;
    }
    return true;
}

bool is_odd(const Perm& a) {
    // Parity = Σ_cycles (len − 1) mod 2 = (n − #cycles) mod 2.
    std::vector<bool> seen(a.size(), false);
    std::size_t cycles = 0U;
    for (std::size_t i = 0U; i < a.size(); ++i) {
        if (seen[i]) continue;
        ++cycles;
        for (std::size_t j = i; !seen[j]; j = a[j]) seen[j] = true;
    }
    return ((a.size() - cycles) % 2U) == 1U;
}

std::vector<std::size_t> cycle_type(const Perm& a) {
    std::vector<bool> seen(a.size(), false);
    std::vector<std::size_t> cyc;
    for (std::size_t i = 0U; i < a.size(); ++i) {
        if (seen[i]) continue;
        std::size_t len = 0U;
        for (std::size_t j = i; !seen[j]; j = a[j]) {
            seen[j] = true;
            ++len;
        }
        cyc.push_back(len);
    }
    std::sort(cyc.begin(), cyc.end(), std::greater<std::size_t>());
    return cyc;
}

std::uint64_t factorial_u64(std::size_t n) {
    assert(n <= 20U && "factorial_u64: n! overflows uint64 for n > 20");
    std::uint64_t f = 1U;
    for (std::size_t i = 2U; i <= n; ++i) f *= static_cast<std::uint64_t>(i);
    return f;
}

std::uint32_t lehmer_rank(const Perm& a) {
    // rank = Σ_i c_i · (n−1−i)!  where c_i = #{j > i : a[j] < a[i]}.
    const std::size_t n = a.size();
    assert(factorial_u64(n) <= UINT32_MAX && "lehmer_rank: degree too large");
    std::uint64_t r = 0U;
    for (std::size_t i = 0U; i < n; ++i) {
        std::size_t c = 0U;
        for (std::size_t j = i + 1U; j < n; ++j) {
            if (a[j] < a[i]) ++c;
        }
        r += static_cast<std::uint64_t>(c) * factorial_u64(n - 1U - i);
    }
    return static_cast<std::uint32_t>(r);
}

Perm lehmer_unrank(std::uint32_t r, std::size_t n) {
    std::vector<std::uint8_t> pool;
    pool.reserve(n);
    for (std::size_t i = 0U; i < n; ++i) pool.push_back(static_cast<std::uint8_t>(i));
    Perm out;
    out.reserve(n);
    std::uint64_t rem = r;
    for (std::size_t i = 0U; i < n; ++i) {
        const std::uint64_t f = factorial_u64(n - 1U - i);
        const std::size_t idx = static_cast<std::size_t>(rem / f);
        rem %= f;
        out.push_back(pool[idx]);
        pool.erase(pool.begin() + static_cast<std::ptrdiff_t>(idx));
    }
    return out;
}

Result<PermGroup> PermGroup::closure(
    std::size_t n, std::vector<Perm> gens, std::uint64_t max_order) {
    for (const auto& g : gens) {
        if (g.size() != n || !is_valid_perm(g)) {
            return fail<PermGroup>(CASError{
                .kind = CASErrorKind::InvalidArgument,
                .message = "PermGroup::closure: generator is not a valid "
                           "permutation of the requested degree"});
        }
    }
    // BFS closure over left-multiplication by generators. The generated set
    // is closed under products of generators and contains the identity; in a
    // finite ambient group that set is exactly ⟨gens⟩ (inverses arise as
    // powers).
    std::unordered_set<std::uint32_t> seen;
    std::queue<Perm> frontier;
    const Perm id = identity(n);
    seen.insert(lehmer_rank(id));
    frontier.push(id);
    while (!frontier.empty()) {
        Perm cur = std::move(frontier.front());
        frontier.pop();
        for (const auto& g : gens) {
            Perm nxt = compose(g, cur);
            const std::uint32_t rk = lehmer_rank(nxt);
            if (seen.insert(rk).second) {
                if (static_cast<std::uint64_t>(seen.size()) > max_order) {
                    return fail<PermGroup>(CASError{
                        .kind = CASErrorKind::Unimplemented,
                        .message = "PermGroup::closure: generated group order "
                                   "exceeds the configured budget (A6 lattice "
                                   "enumeration cap) — raise the caller-side "
                                   "CASContext budget to proceed"});
                }
                frontier.push(std::move(nxt));
            }
        }
    }
    PermGroup out;
    out.n_ = n;
    out.gens_ = std::move(gens);
    out.ranks_.assign(seen.begin(), seen.end());
    std::sort(out.ranks_.begin(), out.ranks_.end());
    return ok(std::move(out));
}

bool PermGroup::contains(const Perm& p) const {
    if (p.size() != n_) return false;
    return contains_rank(lehmer_rank(p));
}

bool PermGroup::contains_rank(std::uint32_t r) const {
    return std::binary_search(ranks_.begin(), ranks_.end(), r);
}

bool PermGroup::contains_group(const PermGroup& other) const {
    if (other.n_ != n_) return false;
    return std::includes(ranks_.begin(), ranks_.end(),
                         other.ranks_.begin(), other.ranks_.end());
}

bool PermGroup::equals(const PermGroup& other) const {
    return n_ == other.n_ && ranks_ == other.ranks_;
}

bool PermGroup::is_transitive() const {
    if (n_ == 0U) return false;
    std::vector<bool> in_orbit(n_, false);
    std::vector<std::size_t> stack{0U};
    in_orbit[0] = true;
    std::size_t count = 1U;
    while (!stack.empty()) {
        const std::size_t pt = stack.back();
        stack.pop_back();
        for (const auto& g : gens_) {
            const std::size_t im = g[pt];
            if (!in_orbit[im]) {
                in_orbit[im] = true;
                ++count;
                stack.push_back(im);
            }
        }
    }
    return count == n_;
}

bool PermGroup::has_odd_element() const {
    // G ⊆ A_n iff every generator is even (parity is a homomorphism).
    for (const auto& g : gens_) {
        if (is_odd(g)) return true;
    }
    return false;
}

std::vector<std::size_t> PermGroup::orbit_lengths_on_ksubsets(
    std::size_t k) const {
    assert(n_ <= 31U && k >= 1U && k <= n_ &&
           "orbit_lengths_on_ksubsets: subset bitmask requires n <= 31");
    // Enumerate all k-subsets as bitmasks, BFS orbits under the generators.
    std::vector<std::uint32_t> subsets;
    // Gosper-style enumeration of k-bit masks below 2^n.
    std::uint32_t mask = (k == 0U) ? 0U : ((1U << k) - 1U);
    const std::uint32_t limit = 1U << n_;
    while (mask < limit) {
        subsets.push_back(mask);
        const std::uint32_t c = mask & (~mask + 1U);
        const std::uint32_t r = mask + c;
        if (r >= limit) break;
        mask = (((r ^ mask) >> 2U) / c) | r;
    }
    std::unordered_set<std::uint32_t> visited;
    std::vector<std::size_t> lengths;
    auto apply = [&](const Perm& g, std::uint32_t s) -> std::uint32_t {
        std::uint32_t out = 0U;
        for (std::size_t i = 0U; i < n_; ++i) {
            if ((s >> i) & 1U) out |= (1U << g[i]);
        }
        return out;
    };
    for (const std::uint32_t start : subsets) {
        if (visited.count(start) != 0U) continue;
        std::size_t len = 0U;
        std::vector<std::uint32_t> stack{start};
        visited.insert(start);
        while (!stack.empty()) {
            const std::uint32_t s = stack.back();
            stack.pop_back();
            ++len;
            for (const auto& g : gens_) {
                const std::uint32_t im = apply(g, s);
                if (visited.insert(im).second) stack.push_back(im);
            }
        }
        lengths.push_back(len);
    }
    std::sort(lengths.begin(), lengths.end(), std::greater<std::size_t>());
    return lengths;
}

std::vector<std::pair<std::vector<std::size_t>, std::size_t>>
PermGroup::cycle_type_distribution() const {
    std::map<std::vector<std::size_t>, std::size_t> dist;
    for (const std::uint32_t r : ranks_) {
        dist[cycle_type(lehmer_unrank(r, n_))] += 1U;
    }
    return {dist.begin(), dist.end()};
}

}  // namespace cas::algebra::permgrp
