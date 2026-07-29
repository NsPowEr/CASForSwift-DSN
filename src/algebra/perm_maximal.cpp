// A6 — Assembly of the maximal-transitive candidate lists (see
// perm_maximal_internal.hpp). Families are generated from the arithmetic of
// n, verified structurally (BSGS order, transitivity, parity, primitivity)
// and filtered: ambient-sized groups are dropped, even/odd parity routes a
// candidate to the A_n / S_n list, direct containments are pruned, and even
// candidates without a proof that their S_n-class stays one A_n-class get a
// transposition-conjugate twin (covers class splitting in A_n).

#include "perm_maximal_internal.hpp"

#include "perm_blocks_internal.hpp"
#include "perm_construct_fields_internal.hpp"
#include "perm_construct_internal.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace cas::algebra::permgrp {

namespace {

// A generator set with its structural provenance, before BSGS verification.
struct RawFamily {
    std::vector<Perm> gens;
    std::string provenance;
    // True when the construction PROVES an odd element of S_n normalises the
    // group (then its S_n-class is a single A_n-class and no twin is needed):
    // either the group is normal in an odd family member, or it is the even
    // part of an odd group (which normalises it, having it as index-2).
    bool no_split_certificate{false};
};

[[nodiscard]] Perm conj(const Perm& t, const Perm& g) {
    return compose(t, compose(g, inverse(t)));
}

[[nodiscard]] std::vector<Perm> conj_all(const Perm& t,
                                         const std::vector<Perm>& gens) {
    std::vector<Perm> out;
    out.reserve(gens.size());
    for (const Perm& g : gens) out.push_back(conj(t, g));
    return out;
}

// n = p^e with p prime? (exact trial factorisation)
[[nodiscard]] bool as_prime_power(std::size_t n, std::uint32_t& p,
                                  std::size_t& e) {
    if (n < 2U) return false;
    std::size_t m = n;
    std::uint32_t q = 2U;
    while (m % q != 0U) ++q;
    e = 0U;
    while (m % q == 0U) {
        m /= q;
        ++e;
    }
    if (m != 1U) return false;
    p = q;
    return true;
}

// Möbius/semilinear layer on P¹(F_q), q = n−1: PSL(2,q) plus every subgroup
// generated over it by the multiplier x↦γx and the Frobenius x↦x^p (these
// generate PΓL(2,q)/PSL(2,q), so all intermediate subgroups of the small
// abelian quotient are covered; duplicates are pruned later).
void moebius_families(const Gf& gf, std::vector<RawFamily>& out) {
    const std::vector<Perm> psl = psl2_point_gens(gf);
    const Perm gmul = gamma_mult_point_perm(gf);
    const Perm frob = frobenius_point_perm(gf);
    const Perm gfrob = compose(gmul, frob);
    // PΓL(2,q) contains every listed subgroup as a normal subgroup of its
    // abelian-over-PSL layer; if PΓL has an odd element the classes cannot
    // split in A_n.
    bool pgammal_odd = is_odd(gmul) || is_odd(frob);
    for (const Perm& g : psl) pgammal_odd = pgammal_odd || is_odd(g);

    const std::string base = "(2," + std::to_string(gf.q) + ") on P1(F_" +
                             std::to_string(gf.q) + ")";
    auto push = [&](std::vector<Perm> gens, const std::string& name) {
        out.push_back(RawFamily{std::move(gens), name, pgammal_odd});
    };
    push(psl, "PSL" + base);
    std::vector<Perm> pgl = psl;
    pgl.push_back(gmul);
    push(std::move(pgl), "PGL" + base);
    if (gf.e >= 2U) {
        std::vector<Perm> psigmal = psl;
        psigmal.push_back(frob);
        push(std::move(psigmal), "PSL.<frob>" + base);
        std::vector<Perm> twisted = psl;
        twisted.push_back(gfrob);
        push(std::move(twisted), "PSL.<gamma*frob>" + base);
        std::vector<Perm> pgammal = psl;
        pgammal.push_back(gmul);
        pgammal.push_back(frob);
        push(std::move(pgammal), "PGammaL" + base);
    }
}

}  // namespace

Result<std::vector<MaximalCandidate>> maximal_transitive_candidates(
    AmbientGroup ambient, std::size_t n) {
    if (n < 5U || n > 10U) {
        return fail<std::vector<MaximalCandidate>>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message =
                "maximal_transitive_candidates: degree outside [5,10] — "
                "n <= 4 has closed-form Galois treatment, n >= 11 awaits "
                "the coverage certification of a later A6 increment"});
    }

    // ---- assemble the raw families from the arithmetic of n ----
    std::vector<RawFamily> raw;

    // Imprimitive wreaths S_a ≀ S_b for every divisor split n = a·b.
    for (std::size_t a = 2U; a < n; ++a) {
        if (n % a != 0U) continue;
        const std::size_t b = n / a;
        if (b < 2U) continue;
        auto w = wreath_gens(a, b);
        if (w.is_error())
            return fail<std::vector<MaximalCandidate>>(w.error());
        raw.push_back(RawFamily{std::move(w.value()),
                                "S_" + std::to_string(a) + " wr S_" +
                                    std::to_string(b) + " (" +
                                    std::to_string(b) + " blocks of " +
                                    std::to_string(a) + ")",
                                false});
    }

    // Affine AGL(d,p) whenever n = p^d.
    {
        std::uint32_t p = 0U;
        std::size_t e = 0U;
        if (as_prime_power(n, p, e)) {
            auto agl = agl_gens(e, p);
            if (agl.is_error())
                return fail<std::vector<MaximalCandidate>>(agl.error());
            raw.push_back(RawFamily{std::move(agl.value()),
                                    "AGL(" + std::to_string(e) + "," +
                                        std::to_string(p) + ") affine",
                                    false});
        }
    }

    // Projective PGL(d,p) on P^{d−1}(F_p) whenever n = (p^d−1)/(p−1), d ≥ 3.
    for (std::uint32_t p = 2U; p <= static_cast<std::uint32_t>(n); ++p) {
        bool p_prime = p >= 2U;
        for (std::uint32_t dd = 2U; dd * dd <= p && p_prime; ++dd) {
            if (p % dd == 0U) p_prime = false;
        }
        if (!p_prime) continue;
        std::size_t points = 1U + p;  // d = 2 term; extend below
        for (std::size_t d = 3U; points <= n; ++d) {
            std::size_t pw = 1U;
            for (std::size_t i = 0U; i < d; ++i) pw *= p;
            points = (pw - 1U) / (p - 1U);
            if (points != n) continue;
            auto proj = projective_gl_gens(d, p);
            if (proj.is_error())
                return fail<std::vector<MaximalCandidate>>(proj.error());
            raw.push_back(RawFamily{std::move(proj.value()),
                                    "PGL(" + std::to_string(d) + "," +
                                        std::to_string(p) + ") on P^" +
                                        std::to_string(d - 1U) + "(F_" +
                                        std::to_string(p) + ")",
                                    false});
        }
    }

    // Möbius/semilinear families on P¹(F_q), q = n−1 a prime power.
    {
        std::uint32_t p = 0U;
        std::size_t e = 0U;
        if (as_prime_power(n - 1U, p, e)) {
            auto gf = build_gf(p, e);
            if (gf.is_error())
                return fail<std::vector<MaximalCandidate>>(gf.error());
            moebius_families(gf.value(), raw);
        }
    }

    // Subset actions of S_m and A_m whenever n = C(m,k), 2 ≤ k ≤ m/2 (the
    // complementary action on (m−k)-subsets is the same group up to the
    // complementation relabelling, so only the canonical half is built).
    for (std::size_t m = 5U; m <= 20U; ++m) {
        for (std::size_t k = 2U; 2U * k <= m; ++k) {
            if (binomial_u64(m, k) != static_cast<std::uint64_t>(n)) continue;
            auto sm = ksubset_action_gens(m, k, symmetric_gens(m));
            if (sm.is_error())
                return fail<std::vector<MaximalCandidate>>(sm.error());
            const std::string suffix = " on " + std::to_string(k) +
                                       "-subsets of " + std::to_string(m) +
                                       " points";
            // A_m-image is normal (index 2) in the S_m-image; if the latter
            // is odd, neither class splits in A_n.
            bool sm_odd = false;
            for (const Perm& g : sm.value()) sm_odd = sm_odd || is_odd(g);
            auto am = ksubset_action_gens(m, k, alternating_gens(m));
            if (am.is_error())
                return fail<std::vector<MaximalCandidate>>(am.error());
            raw.push_back(
                RawFamily{std::move(sm.value()), "S_" + std::to_string(m) + suffix, false});
            raw.push_back(RawFamily{std::move(am.value()),
                                    "A_" + std::to_string(m) + suffix,
                                    sm_odd});
        }
    }

    // ---- verify structurally and split by parity side ----
    const std::uint64_t full_order = factorial_u64(n);
    struct Verified {
        BsgsGroup group;
        std::string provenance;
        bool no_split_certificate;
    };
    std::vector<Verified> s_side;
    std::vector<Verified> a_side;

    auto classify = [&](std::vector<Perm> gens, std::string provenance,
                        bool cert) -> Result<void> {
        auto built = BsgsGroup::build(n, std::move(gens));
        if (built.is_error()) return fail<void>(built.error());
        BsgsGroup g = std::move(built.value());
        // A candidate must be able to contain a transitive G_f.
        if (!g.is_transitive()) return ok();
        const std::uint64_t o = g.order();
        if (o == full_order) return ok();  // the whole S_n: not a descent step
        if (!g.has_odd_element()) {
            if (o == full_order / 2U) return ok();  // A_n itself: the ambient
            a_side.push_back(
                Verified{std::move(g), std::move(provenance), cert});
        } else {
            s_side.push_back(
                Verified{std::move(g), std::move(provenance), cert});
        }
        return ok();
    };

    for (RawFamily& fam : raw) {
        auto r = classify(std::move(fam.gens), std::move(fam.provenance),
                          fam.no_split_certificate);
        if (r.is_error()) return fail<std::vector<MaximalCandidate>>(r.error());
    }

    // Even parts of the odd candidates feed the A_n side (Schreier); the odd
    // parent normalises its index-2 even part, so no twin is ever needed.
    const std::size_t s_count = s_side.size();
    for (std::size_t i = 0U; i < s_count; ++i) {
        auto ev = even_part_gens(n, s_side[i].group.generators());
        if (ev.is_error())
            return fail<std::vector<MaximalCandidate>>(ev.error());
        auto r = classify(std::move(ev.value()),
                          "even part of " + s_side[i].provenance, true);
        if (r.is_error()) return fail<std::vector<MaximalCandidate>>(r.error());
    }

    // Twins: an even candidate whose class may split in A_n gets its
    // (0 1)-conjugate as a second representative (skipped when the conjugate
    // is the same subgroup).
    {
        const Perm t01 = [&] {
            Perm t = identity(n);
            t[0] = 1U;
            t[1] = 0U;
            return t;
        }();
        const std::size_t a_count = a_side.size();
        for (std::size_t i = 0U; i < a_count; ++i) {
            if (a_side[i].no_split_certificate) continue;
            std::vector<Perm> tg = conj_all(t01, a_side[i].group.generators());
            bool same = true;
            for (const Perm& g : tg) same = same && a_side[i].group.contains(g);
            if (same) continue;
            auto built = BsgsGroup::build(n, std::move(tg));
            if (built.is_error())
                return fail<std::vector<MaximalCandidate>>(built.error());
            a_side.push_back(Verified{std::move(built.value()),
                                      a_side[i].provenance + " ((0 1)-twin)",
                                      true});
        }
    }

    // Prune direct containments within the requested side (same point
    // indexing): a candidate lying inside another is reachable through it
    // later in the descent, so listing it here only duplicates work.
    std::vector<Verified>& side =
        (ambient == AmbientGroup::Symmetric) ? s_side : a_side;
    std::vector<bool> drop(side.size(), false);
    for (std::size_t i = 0U; i < side.size(); ++i) {
        for (std::size_t j = 0U; j < side.size() && !drop[i]; ++j) {
            if (i == j || drop[j]) continue;
            const std::uint64_t oi = side[i].group.order();
            const std::uint64_t oj = side[j].group.order();
            if (oi > oj || (oi == oj && i < j)) continue;
            bool contained = true;
            for (const Perm& g : side[i].group.generators()) {
                if (!side[j].group.contains(g)) {
                    contained = false;
                    break;
                }
            }
            if (contained) drop[i] = true;
        }
    }

    // ---- final assembly (A_n first for the symmetric ambient) ----
    std::vector<MaximalCandidate> out;
    if (ambient == AmbientGroup::Symmetric) {
        auto an = BsgsGroup::build(n, alternating_gens(n));
        if (an.is_error())
            return fail<std::vector<MaximalCandidate>>(an.error());
        out.push_back(MaximalCandidate{std::move(an.value()),
                                       "A_" + std::to_string(n), true});
    }
    for (std::size_t i = 0U; i < side.size(); ++i) {
        if (drop[i]) continue;
        auto prim =
            is_primitive(n, side[i].group.generators());
        if (prim.is_error())
            return fail<std::vector<MaximalCandidate>>(prim.error());
        out.push_back(MaximalCandidate{std::move(side[i].group),
                                       std::move(side[i].provenance),
                                       prim.value()});
    }
    return ok(std::move(out));
}

}  // namespace cas::algebra::permgrp
