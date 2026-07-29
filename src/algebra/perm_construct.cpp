// A6 — Algorithmic permutation-group constructors (see
// perm_construct_internal.hpp). Conventions from perm_group_internal.hpp:
// image form, compose(a,b)[i] = a[b[i]], point action pt^g = g[pt].

#include "perm_construct_internal.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas::algebra::permgrp {

namespace {

// σ with σ(pts[i]) = pts[i+1 mod len], identity elsewhere.
[[nodiscard]] Perm cycle_perm(std::size_t n,
                              const std::vector<std::size_t>& pts) {
    Perm p = identity(n);
    for (std::size_t i = 0U; i < pts.size(); ++i) {
        p[pts[i]] = static_cast<std::uint8_t>(pts[(i + 1U) % pts.size()]);
    }
    return p;
}

}  // namespace

std::vector<Perm> symmetric_gens(std::size_t n) {
    std::vector<Perm> out;
    if (n < 2U) return out;
    out.push_back(cycle_perm(n, {0U, 1U}));
    if (n >= 3U) {
        std::vector<std::size_t> full(n);
        for (std::size_t i = 0U; i < n; ++i) full[i] = i;
        out.push_back(cycle_perm(n, full));
    }
    return out;
}

std::vector<Perm> alternating_gens(std::size_t n) {
    std::vector<Perm> out;
    if (n < 3U) return out;
    out.push_back(cycle_perm(n, {0U, 1U, 2U}));
    if (n == 3U) return out;
    if (n % 2U == 1U) {
        std::vector<std::size_t> full(n);  // n-cycle, even for odd n
        for (std::size_t i = 0U; i < n; ++i) full[i] = i;
        out.push_back(cycle_perm(n, full));
    } else {
        std::vector<std::size_t> tail(n - 1U);  // (n−1)-cycle on 1..n−1, even
        for (std::size_t i = 0U; i + 1U < n; ++i) tail[i] = i + 1U;
        out.push_back(cycle_perm(n, tail));
    }
    return out;
}

Result<std::vector<Perm>> wreath_gens(std::size_t a, std::size_t b) {
    if (a < 1U || b < 1U || a * b > 255U) {
        return fail<std::vector<Perm>>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "wreath_gens: need a,b >= 1 with ab within the Perm "
                       "image bound (255)"});
    }
    const std::size_t n = a * b;
    std::vector<Perm> out;
    // S_a acting on block 0 (its conjugates under the block cycle give the
    // copies on every other block, so ⟨out⟩ is the full wreath product).
    if (a >= 2U) out.push_back(cycle_perm(n, {0U, 1U}));
    if (a >= 3U) {
        std::vector<std::size_t> blk(a);
        for (std::size_t i = 0U; i < a; ++i) blk[i] = i;
        out.push_back(cycle_perm(n, blk));
    }
    // S_b permuting the blocks pointwise: swap of blocks 0,1 and b-cycle.
    if (b >= 2U) {
        Perm swap01 = identity(n);
        for (std::size_t t = 0U; t < a; ++t) {
            swap01[t] = static_cast<std::uint8_t>(a + t);
            swap01[a + t] = static_cast<std::uint8_t>(t);
        }
        out.push_back(std::move(swap01));
    }
    if (b >= 3U) {
        Perm bcycle(n);
        for (std::size_t j = 0U; j < b; ++j) {
            for (std::size_t t = 0U; t < a; ++t) {
                bcycle[j * a + t] =
                    static_cast<std::uint8_t>(((j + 1U) % b) * a + t);
            }
        }
        out.push_back(std::move(bcycle));
    }
    return ok(std::move(out));
}

Result<std::vector<Perm>> even_part_gens(std::size_t n,
                                         const std::vector<Perm>& gens) {
    for (const Perm& g : gens) {
        if (g.size() != n || !is_valid_perm(g)) {
            return fail<std::vector<Perm>>(CASError{
                .kind = CASErrorKind::InvalidArgument,
                .message = "even_part_gens: generator is not a valid "
                           "permutation of the requested degree"});
        }
    }
    const Perm* first_odd = nullptr;
    for (const Perm& g : gens) {
        if (is_odd(g)) { first_odd = &g; break; }
    }
    if (first_odd == nullptr) return ok(std::vector<Perm>(gens));

    // Schreier generators for ker(parity) with transversal {e, t}:
    //   g even → g and t·g·t⁻¹ ;  g odd → g·t⁻¹ and t·g.
    const Perm& t = *first_odd;
    const Perm t_inv = inverse(t);
    std::vector<Perm> out;
    auto push_nonid = [&out](Perm p) {
        for (std::size_t i = 0U; i < p.size(); ++i) {
            if (p[i] != static_cast<std::uint8_t>(i)) {
                out.push_back(std::move(p));
                return;
            }
        }
    };
    for (const Perm& g : gens) {
        if (is_odd(g)) {
            push_nonid(compose(g, t_inv));
            push_nonid(compose(t, g));
        } else {
            push_nonid(g);
            push_nonid(compose(t, compose(g, t_inv)));
        }
    }
    return ok(std::move(out));
}

std::uint64_t binomial_u64(std::size_t m, std::size_t k) {
    if (k > m) return 0U;
    if (k > m - k) k = m - k;
    std::uint64_t r = 1U;
    for (std::size_t i = 1U; i <= k; ++i) {
        // Exact at every step: r accumulates C(m-k+i, i)·i!/i! — the division
        // is exact because C(m-k+i, i) = r·(m-k+i)/i is an integer.
        r = r * (m - k + i) / i;
    }
    return r;
}

Result<std::vector<Perm>> ksubset_action_gens(std::size_t m, std::size_t k,
                                              const std::vector<Perm>& gens) {
    if (m < 2U || m > 20U || k < 1U || k >= m) {
        return fail<std::vector<Perm>>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "ksubset_action_gens: need 2 <= m <= 20 and "
                       "1 <= k <= m-1"});
    }
    const std::uint64_t count = binomial_u64(m, k);
    if (count > 255U) {
        return fail<std::vector<Perm>>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "ksubset_action_gens: C(m,k) exceeds the Perm image "
                       "bound (255) — needs the wider point type of a later "
                       "increment"});
    }
    for (const Perm& g : gens) {
        if (g.size() != m || !is_valid_perm(g)) {
            return fail<std::vector<Perm>>(CASError{
                .kind = CASErrorKind::InvalidArgument,
                .message = "ksubset_action_gens: generator is not a valid "
                           "permutation of degree m"});
        }
    }
    // Enumerate the k-subsets as bitmasks in increasing numeric order; the
    // sorted vector doubles as the index (binary search for lookup).
    std::vector<std::uint32_t> masks;
    masks.reserve(static_cast<std::size_t>(count));
    for (std::uint32_t mask = 0U; mask < (1U << m); ++mask) {
        if (static_cast<std::size_t>(__builtin_popcount(mask)) == k) {
            masks.push_back(mask);
        }
    }
    std::vector<Perm> out;
    out.reserve(gens.size());
    for (const Perm& g : gens) {
        Perm img(masks.size());
        for (std::size_t idx = 0U; idx < masks.size(); ++idx) {
            std::uint32_t image_mask = 0U;
            for (std::size_t pt = 0U; pt < m; ++pt) {
                if ((masks[idx] >> pt) & 1U) image_mask |= 1U << g[pt];
            }
            const auto it =
                std::lower_bound(masks.begin(), masks.end(), image_mask);
            img[idx] = static_cast<std::uint8_t>(
                static_cast<std::size_t>(it - masks.begin()));
        }
        out.push_back(std::move(img));
    }
    return ok(std::move(out));
}

}  // namespace cas::algebra::permgrp
