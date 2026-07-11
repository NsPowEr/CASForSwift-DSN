// A6 / CAS-L3-18 — BsgsGroup: deterministic Schreier-Sims (see
// perm_bsgs_internal.hpp). Conventions inherited from perm_group_internal.hpp:
//   • Perm in image form, img[i] = σ(i).
//   • compose(a,b)[i] = a[b[i]]  (apply b, then a).
//   • point action  pt^g := g[pt].
// A transversal rep u_x for base point b satisfies u_x[b] = x (b^{u_x} = x).

#include "perm_bsgs_internal.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace cas::algebra::permgrp {

namespace {

[[nodiscard]] bool is_identity_perm(const Perm& p) {
    for (std::size_t i = 0U; i < p.size(); ++i) {
        if (p[i] != static_cast<std::uint8_t>(i)) return false;
    }
    return true;
}

// Orbit of `base_point` under `gens`, with a Schreier-tree transversal:
// transversal[x] = u with u[base_point] = x (nullopt if x ∉ orbit).
struct OrbitData {
    std::vector<std::optional<Perm>> transversal;  // indexed by point, size n
    std::vector<std::size_t> orbit;                // discovery order
};

[[nodiscard]] OrbitData orbit_transversal(std::size_t n, std::size_t base_point,
                                          const std::vector<Perm>& gens) {
    OrbitData od;
    od.transversal.assign(n, std::nullopt);
    od.transversal[base_point] = identity(n);
    od.orbit.push_back(base_point);
    for (std::size_t idx = 0U; idx < od.orbit.size(); ++idx) {
        const std::size_t x = od.orbit[idx];
        const Perm ux = od.transversal[x].value();  // copy: od grows in-loop
        for (const Perm& g : gens) {
            const std::size_t y = g[x];
            if (!od.transversal[y].has_value()) {
                od.transversal[y] = compose(g, ux);  // (g∘ux)[base] = g[x] = y
                od.orbit.push_back(y);
            }
        }
    }
    return od;
}

}  // namespace

Result<BsgsGroup> BsgsGroup::build(std::size_t n, std::vector<Perm> gens) {
    if (n < 1U || n > 20U) {
        return fail<BsgsGroup>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "BsgsGroup::build: degree outside the u64-order-safe "
                       "range [1,20] (|G| <= n! must fit uint64)"});
    }
    for (const auto& g : gens) {
        if (g.size() != n || !is_valid_perm(g)) {
            return fail<BsgsGroup>(CASError{
                .kind = CASErrorKind::InvalidArgument,
                .message = "BsgsGroup::build: generator is not a valid "
                           "permutation of the requested degree"});
        }
    }

    // Working chain state.
    std::vector<std::size_t> base;
    std::vector<Perm> strong;
    std::vector<std::vector<Perm>> lgens;  // lgens[i] = strong gens fixing b_0..b_{i-1}
    std::vector<OrbitData> ods;            // ods[i] = orbit/transversal of b_i under lgens[i]

    auto fixes_prefix = [&](const Perm& s, std::size_t lvl) {
        for (std::size_t t = 0U; t < lvl; ++t) {
            if (s[base[t]] != static_cast<std::uint8_t>(base[t])) return false;
        }
        return true;
    };

    // Recompute lgens + transversals for the whole chain from (base, strong).
    auto rebuild = [&]() {
        lgens.assign(base.size(), {});
        ods.assign(base.size(), {});
        for (std::size_t i = 0U; i < base.size(); ++i) {
            for (const Perm& s : strong) {
                if (fixes_prefix(s, i)) lgens[i].push_back(s);
            }
            ods[i] = orbit_transversal(n, base[i], lgens[i]);
        }
    };

    // Strip h through the chain. Returns (fall-out level, residue). The residue
    // is the identity iff h ∈ ⟨strong⟩ = the group; a fall-out level < base
    // size means h moves that base point out of its transversal orbit.
    auto strip = [&](Perm h) -> std::pair<std::size_t, Perm> {
        for (std::size_t i = 0U; i < base.size(); ++i) {
            const std::size_t x = h[base[i]];
            if (!ods[i].transversal[x].has_value()) return {i, std::move(h)};
            h = compose(inverse(ods[i].transversal[x].value()), h);
        }
        return {base.size(), std::move(h)};
    };

    // Add a strong generator, extending the base if it fixes every base point,
    // then rebuild the chain.
    auto add_strong = [&](const Perm& r) {
        strong.push_back(r);
        bool moves_a_base = false;
        for (const std::size_t b : base) {
            if (r[b] != static_cast<std::uint8_t>(b)) { moves_a_base = true; break; }
        }
        if (!moves_a_base) {
            for (std::size_t p = 0U; p < n; ++p) {
                if (r[p] != static_cast<std::uint8_t>(p)) { base.push_back(p); break; }
            }
        }
        rebuild();
    };

    for (const Perm& g : gens) {
        if (!is_identity_perm(g)) add_strong(g);
    }

    // Verify the chain: every Schreier generator must strip to the identity.
    // On a non-identity residue, add it and restart from its fall-out level.
    // Deterministic termination: each addition strictly enlarges a stabilizer
    // (orders bounded by n!), so only finitely many occur.
    long li = static_cast<long>(base.size()) - 1;
    while (li >= 0) {
        const std::size_t lvl = static_cast<std::size_t>(li);
        // orbit/transversal are copied: the loop condition re-reads orbit.size()
        // and add_strong() below calls rebuild(), which reassigns ods and would
        // dangle a reference into it. Li stays a reference — it is only iterated
        // up to the first add_strong(), which breaks out immediately.
        const std::vector<std::size_t> orbit = ods[lvl].orbit;
        const std::vector<std::optional<Perm>> transversal = ods[lvl].transversal;
        const std::vector<Perm>& Li = lgens[lvl];
        bool complete = true;
        for (std::size_t oi = 0U; oi < orbit.size() && complete; ++oi) {
            const std::size_t x = orbit[oi];
            const Perm ux = transversal[x].value();
            for (const Perm& s : Li) {
                const std::size_t y = s[x];  // y ∈ orbit (closed under Li)
                const Perm uy = transversal[y].value();
                // Schreier generator u_y^{-1} ∘ s ∘ u_x fixes base[lvl].
                Perm sg = compose(inverse(uy), compose(s, ux));
                auto [j, r] = strip(std::move(sg));
                if (!is_identity_perm(r)) {
                    add_strong(r);  // may extend base + rebuild ods/lgens
                    li = static_cast<long>(j);
                    complete = false;
                    break;
                }
            }
        }
        if (complete) --li;
    }

    BsgsGroup out;
    out.n_ = n;
    out.input_gens_ = std::move(gens);
    out.base_ = std::move(base);
    out.strong_gens_ = std::move(strong);
    out.level_gens_ = std::move(lgens);
    out.transversal_.reserve(ods.size());
    for (auto& od : ods) out.transversal_.push_back(std::move(od.transversal));
    return ok(std::move(out));
}

std::uint64_t BsgsGroup::order() const noexcept {
    std::uint64_t o = 1U;
    for (const auto& tr : transversal_) {
        std::uint64_t sz = 0U;
        for (const auto& u : tr) {
            if (u.has_value()) ++sz;
        }
        o *= sz;
    }
    return o;
}

std::optional<Perm> BsgsGroup::sift(const Perm& p) const {
    if (p.size() != n_) return std::nullopt;
    Perm h = p;
    for (std::size_t i = 0U; i < base_.size(); ++i) {
        const std::size_t x = h[base_[i]];
        if (!transversal_[i][x].has_value()) return h;  // residue (non-identity)
        h = compose(inverse(transversal_[i][x].value()), h);
    }
    return h;
}

bool BsgsGroup::contains(const Perm& p) const {
    auto r = sift(p);
    return r.has_value() && is_identity_perm(*r);
}

bool BsgsGroup::is_transitive() const {
    if (n_ == 0U) return false;
    std::vector<bool> in_orbit(n_, false);
    std::vector<std::size_t> stack{0U};
    in_orbit[0] = true;
    std::size_t count = 1U;
    while (!stack.empty()) {
        const std::size_t pt = stack.back();
        stack.pop_back();
        for (const Perm& g : input_gens_) {
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

bool BsgsGroup::has_odd_element() const {
    // G ⊆ A_n iff every generator is even (parity is a homomorphism).
    for (const Perm& g : input_gens_) {
        if (is_odd(g)) return true;
    }
    return false;
}

}  // namespace cas::algebra::permgrp
