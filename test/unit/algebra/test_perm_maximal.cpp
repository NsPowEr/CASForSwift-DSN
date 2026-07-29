// A6 — Maximal-transitive candidate tests, Brick 2 of the deg ≥ 8 Stauduhar
// closure.
//
// The decisive check is the n ≤ 7 CROSS-CHECK against the exhaustively
// generated subgroup lattice (which is itself table-free): every candidate
// must be a genuine proper transitive class, every transitive class must lie
// under some candidate (up to ambient conjugacy — with σ restricted to A_n
// on the alternating side, exercising the twin mechanism), and no transitive
// class may sit strictly between a candidate and the ambient (maximality).
// For 8 ≤ n ≤ 10 the lattice is out of reach by design, so the tests verify
// the structural facts (orders vs classical formulas, parity, primitivity,
// twin distinctness) that the assembly contract promises.

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "../../../src/algebra/perm_bsgs_internal.hpp"
#include "../../../src/algebra/perm_group_internal.hpp"
#include "../../../src/algebra/perm_maximal_internal.hpp"
#include "cas/error.hpp"

using namespace cas;
using namespace cas::algebra::permgrp;

namespace {

[[nodiscard]] PermGroup dense_of(std::size_t n, const std::vector<Perm>& gens) {
    auto d = PermGroup::closure(n, gens, factorial_u64(n));
    EXPECT_TRUE(d.is_ok());
    return std::move(d.value());
}

[[nodiscard]] std::vector<Perm> conj_gens(const Perm& t,
                                          const std::vector<Perm>& gens) {
    std::vector<Perm> out;
    out.reserve(gens.size());
    for (const Perm& g : gens) out.push_back(compose(t, compose(g, inverse(t))));
    return out;
}

// ∃ σ (even σ only, if requested) with T^σ ⊆ H — exhaustive over S_n,
// exact. T is given by generators, H as a dense group.
[[nodiscard]] bool conj_contained(std::size_t n,
                                  const std::vector<Perm>& t_gens,
                                  std::uint64_t t_order, const PermGroup& h,
                                  bool even_sigma_only) {
    if (t_order == 0U || h.order() % t_order != 0U) return false;  // Lagrange
    const std::uint64_t full = factorial_u64(n);
    for (std::uint64_t sr = 0U; sr < full; ++sr) {
        const Perm sigma = lehmer_unrank(static_cast<std::uint32_t>(sr), n);
        if (even_sigma_only && is_odd(sigma)) continue;
        const Perm sigma_inv = inverse(sigma);
        bool all_in = true;
        for (const Perm& g : t_gens) {
            if (!h.contains(compose(sigma, compose(g, sigma_inv)))) {
                all_in = false;
                break;
            }
        }
        if (all_in) return true;
    }
    return false;
}

void cross_check(std::size_t n, AmbientGroup ambient, std::uint64_t max_ops) {
    auto cands_r = maximal_transitive_candidates(ambient, n);
    ASSERT_TRUE(cands_r.is_ok()) << cands_r.error().message;
    const std::vector<MaximalCandidate>& cands = cands_r.value();
    ASSERT_FALSE(cands.empty());

    auto lattice_r = transitive_subgroup_classes(
        n, LatticeBudget{.max_degree = 7U, .max_ops = max_ops});
    ASSERT_TRUE(lattice_r.is_ok()) << lattice_r.error().message;
    const std::vector<PermGroup>& lattice = lattice_r.value();

    const std::uint64_t full = factorial_u64(n);
    const std::uint64_t amb_order =
        (ambient == AmbientGroup::Symmetric) ? full : full / 2U;

    // ── soundness: each candidate is a proper transitive class of the
    //    lattice, on the right parity side ─────────────────────────────────
    std::vector<PermGroup> dense_cands;
    for (const MaximalCandidate& c : cands) {
        PermGroup d = dense_of(n, c.group.generators());
        EXPECT_TRUE(d.is_transitive()) << c.provenance;
        EXPECT_LT(d.order(), amb_order + (ambient == AmbientGroup::Symmetric
                                              ? 0U
                                              : 1U))
            << c.provenance;
        EXPECT_NE(d.order(), full) << c.provenance;
        if (ambient == AmbientGroup::Alternating) {
            EXPECT_FALSE(d.has_odd_element()) << c.provenance;
            EXPECT_LT(d.order(), amb_order) << c.provenance;
        }
        bool in_lattice = false;
        for (const PermGroup& cls : lattice) {
            if (conjugate_in_sn(d, cls)) {
                in_lattice = true;
                break;
            }
        }
        EXPECT_TRUE(in_lattice) << c.provenance << " not a lattice class";
        dense_cands.push_back(std::move(d));
    }

    // ── coverage: every transitive class below the ambient lies under a
    //    candidate (alternating side: under A_n-conjugacy, both A_n-class
    //    representatives) ────────────────────────────────────────────────────
    const Perm t01 = [n] {
        Perm t = identity(n);
        t[0] = 1U;
        t[1] = 0U;
        return t;
    }();
    for (const PermGroup& cls : lattice) {
        if (cls.order() == full) continue;  // S_n itself
        if (ambient == AmbientGroup::Symmetric) {
            bool covered = false;
            for (const PermGroup& h : dense_cands) {
                if (conj_contained(n, cls.generators(), cls.order(), h,
                                   false)) {
                    covered = true;
                    break;
                }
            }
            EXPECT_TRUE(covered)
                << "S-side coverage failed for class of order " << cls.order();
        } else {
            if (cls.has_odd_element()) continue;   // not below A_n… it is,
                                                   // but not on this side's
                                                   // descent (odd ⇒ ⊄ A_n)
            if (cls.order() == amb_order) continue;  // A_n itself
            const std::vector<std::vector<Perm>> reps{
                cls.generators(), conj_gens(t01, cls.generators())};
            for (std::size_t ri = 0U; ri < reps.size(); ++ri) {
                bool covered = false;
                for (const PermGroup& h : dense_cands) {
                    if (conj_contained(n, reps[ri], cls.order(), h, true)) {
                        covered = true;
                        break;
                    }
                }
                EXPECT_TRUE(covered)
                    << "A-side coverage failed for class of order "
                    << cls.order() << " (A_n-class rep " << ri << ")";
            }
        }
    }

    // ── near-maximality tripwire (exact at n ≤ 7): a listed candidate MAY
    //    be non-maximal (the contract allows benign redundancy — e.g. F_21 =
    //    (AGL(1,7) ∩ A_7) sits inside PSL(3,2), the Sylow-7 normaliser), but
    //    then the intermediate class must itself be conjugate to a listed
    //    candidate. An intermediate that is NOT in the list would flag a
    //    family gap even though coverage still holds through it. ───────────
    for (std::size_t ci = 0U; ci < cands.size(); ++ci) {
        const std::uint64_t ho = dense_cands[ci].order();
        if (ho == full / 2U && ambient == AmbientGroup::Symmetric) continue;
        for (const PermGroup& k : lattice) {
            if (k.order() <= ho || k.order() >= amb_order) continue;
            if (ambient == AmbientGroup::Alternating && k.has_odd_element())
                continue;
            const PermGroup dense_k = dense_of(n, k.generators());
            if (!conj_contained(n, cands[ci].group.generators(), ho, dense_k,
                                false)) {
                continue;  // no class strictly between: candidate is maximal
            }
            bool intermediate_is_listed = false;
            for (const PermGroup& h : dense_cands) {
                if (conjugate_in_sn(h, dense_k)) {
                    intermediate_is_listed = true;
                    break;
                }
            }
            EXPECT_TRUE(intermediate_is_listed)
                << cands[ci].provenance
                << " lies strictly inside an UNLISTED class of order "
                << k.order();
        }
    }
}

[[nodiscard]] std::vector<std::uint64_t> sorted_orders(
    const std::vector<MaximalCandidate>& cands) {
    std::vector<std::uint64_t> o;
    o.reserve(cands.size());
    for (const auto& c : cands) o.push_back(c.group.order());
    std::sort(o.begin(), o.end());
    return o;
}

}  // namespace

// ── exhaustive cross-check against the generated lattice (n ≤ 7) ────────────

TEST(PermMaximalTest, CrossCheckDegree5BothAmbients) {
    cross_check(5U, AmbientGroup::Symmetric, 400'000'000ULL);
    cross_check(5U, AmbientGroup::Alternating, 400'000'000ULL);
}

TEST(PermMaximalTest, CrossCheckDegree6BothAmbients) {
    cross_check(6U, AmbientGroup::Symmetric, 4'000'000'000ULL);
    cross_check(6U, AmbientGroup::Alternating, 4'000'000'000ULL);
}

TEST(PermMaximalTest, CrossCheckDegree7BothAmbients) {
    cross_check(7U, AmbientGroup::Symmetric, 40'000'000'000ULL);
    cross_check(7U, AmbientGroup::Alternating, 40'000'000'000ULL);
}

// ── structural verification for 8 ≤ n ≤ 10 (lattice out of reach) ───────────

TEST(PermMaximalTest, Degree8Structure) {
    auto sym = maximal_transitive_candidates(AmbientGroup::Symmetric, 8U);
    ASSERT_TRUE(sym.is_ok()) << sym.error().message;
    // A_8, S_2≀S_4 (384), S_4≀S_2 (1152), PGL(2,7) (336).
    EXPECT_EQ(sorted_orders(sym.value()),
              (std::vector<std::uint64_t>{336U, 384U, 1152U, 20160U}));
    EXPECT_EQ(sym.value()[0].provenance, "A_8");

    auto alt = maximal_transitive_candidates(AmbientGroup::Alternating, 8U);
    ASSERT_TRUE(alt.is_ok()) << alt.error().message;
    // PSL(2,7) (168), (S_4≀S_2)∩A_8 = 2^4:(S_3×S_3)-shaped (576), AGL(3,2)
    // + its (0 1)-twin (1344 ×2 — the S_8-class splits into two
    // A_8-classes). (S_2≀S_4)∩A_8 (192) is correctly pruned: with the
    // canonical encodings the pair blocks {2j, 2j+1} are the cosets of the
    // subspace {0, e_0} ⊂ F_2³ and the even pair-preserving group is exactly
    // the affine stabiliser of that partition, 8·24 = 192 ⊂ AGL(3,2).
    EXPECT_EQ(sorted_orders(alt.value()),
              (std::vector<std::uint64_t>{168U, 576U, 1344U, 1344U}));
    // The two order-1344 candidates are genuinely different subgroups.
    std::vector<const MaximalCandidate*> affine;
    for (const auto& c : alt.value()) {
        if (c.group.order() == 1344U) affine.push_back(&c);
        if (c.group.order() == 1344U || c.group.order() == 168U) {
            EXPECT_TRUE(c.primitive) << c.provenance;
        }
        if (c.group.order() == 192U || c.group.order() == 576U) {
            EXPECT_FALSE(c.primitive) << c.provenance;
        }
    }
    ASSERT_EQ(affine.size(), 2U);
    bool identical = true;
    for (const Perm& g : affine[0]->group.generators()) {
        identical = identical && affine[1]->group.contains(g);
    }
    EXPECT_FALSE(identical);
}

TEST(PermMaximalTest, Degree9Structure) {
    auto sym = maximal_transitive_candidates(AmbientGroup::Symmetric, 9U);
    ASSERT_TRUE(sym.is_ok()) << sym.error().message;
    // A_9, S_3≀S_3 (1296), AGL(2,3) (432).
    EXPECT_EQ(sorted_orders(sym.value()),
              (std::vector<std::uint64_t>{432U, 1296U, 181440U}));

    auto alt = maximal_transitive_candidates(AmbientGroup::Alternating, 9U);
    ASSERT_TRUE(alt.is_ok()) << alt.error().message;
    // even parts (216, 648), PΓL(2,8) + twin (1512 ×2); PSL(2,8) pruned
    // (it lies inside PΓL(2,8) with the same point indexing).
    EXPECT_EQ(sorted_orders(alt.value()),
              (std::vector<std::uint64_t>{216U, 648U, 1512U, 1512U}));
}

TEST(PermMaximalTest, Degree10Structure) {
    auto sym = maximal_transitive_candidates(AmbientGroup::Symmetric, 10U);
    ASSERT_TRUE(sym.is_ok()) << sym.error().message;
    // A_10, S_2≀S_5 (3840), S_5≀S_2 (28800), PΓL(2,9) (1440), S_5 on pairs
    // (120); PGL(2,9) and PSL(2,9).⟨φ⟩ pruned inside PΓL(2,9).
    EXPECT_EQ(sorted_orders(sym.value()),
              (std::vector<std::uint64_t>{120U, 1440U, 3840U, 28800U,
                                          1814400U}));

    auto alt = maximal_transitive_candidates(AmbientGroup::Alternating, 10U);
    ASSERT_TRUE(alt.is_ok()) << alt.error().message;
    // A_5 on pairs (60), PSL(2,9).⟨γφ⟩ (720, even), even wreath parts
    // (1920, 14400); PSL(2,9) (360) pruned inside the 720.
    EXPECT_EQ(sorted_orders(alt.value()),
              (std::vector<std::uint64_t>{60U, 720U, 1920U, 14400U}));
    for (const auto& c : alt.value()) {
        if (c.group.order() == 720U || c.group.order() == 60U) {
            EXPECT_TRUE(c.primitive) << c.provenance;
        } else {
            EXPECT_FALSE(c.primitive) << c.provenance;
        }
    }
}

TEST(PermMaximalTest, DegreesOutsideRangeAreStructuredErrors) {
    auto low = maximal_transitive_candidates(AmbientGroup::Symmetric, 4U);
    ASSERT_FALSE(low.is_ok());
    EXPECT_EQ(low.error().kind, CASErrorKind::Unimplemented);
    auto high = maximal_transitive_candidates(AmbientGroup::Alternating, 11U);
    ASSERT_FALSE(high.is_ok());
    EXPECT_EQ(high.error().kind, CASErrorKind::Unimplemented);
}
