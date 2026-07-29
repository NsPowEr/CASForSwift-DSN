// A6 Brick 3b — certified relative invariants + coset transversals.
// Every expectation is a mathematical identity on GENERATED groups (Brick
// 1/2 machinery), never a transcribed value:
//   • [G:H] = |G|/|H| coset representatives, images pairwise distinct;
//   • exactly one representative (the identity) fixes F — Stab_G(F) = H
//     verified by counting;
//   • the wreath candidate is separated by a low-degree subset invariant,
//     while H = A_n (k-homogeneous for every k) provably exhausts tier 1
//     and must take the guaranteed Galois-monomial fallback;
//   • structured failures: H ⊄ G, H = G, ops budget.

#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <set>
#include <utility>
#include <vector>

#include "../../../src/algebra/galois_invariant_internal.hpp"
#include "../../../src/algebra/perm_bsgs_internal.hpp"
#include "../../../src/algebra/perm_group_internal.hpp"
#include "../../../src/algebra/perm_maximal_internal.hpp"
#include "cas/error.hpp"

using namespace cas;
using namespace cas::algebra;
using namespace cas::algebra::galois_invariant;
using namespace cas::algebra::permgrp;

namespace {

constexpr std::uint64_t kOps = 1ULL << 26U;

[[nodiscard]] Perm cycle(std::size_t n,
                         std::initializer_list<std::size_t> cyc) {
    Perm p = identity(n);
    if (cyc.size() >= 2U) {
        auto it = cyc.begin();
        const std::size_t first = *it;
        std::size_t prev = first;
        ++it;
        for (; it != cyc.end(); ++it) {
            p[prev] = static_cast<std::uint8_t>(*it);
            prev = *it;
        }
        p[prev] = static_cast<std::uint8_t>(first);
    }
    return p;
}

[[nodiscard]] BsgsGroup symmetric(std::size_t n) {
    std::vector<std::size_t> full(n);
    for (std::size_t i = 0U; i < n; ++i) full[i] = i;
    Perm ncyc = identity(n);
    for (std::size_t i = 0U; i < n; ++i) {
        ncyc[i] = static_cast<std::uint8_t>((i + 1U) % n);
    }
    auto g = BsgsGroup::build(n, {cycle(n, {0U, 1U}), ncyc});
    EXPECT_TRUE(g.is_ok());
    return g.value();
}

// Checks the full exactness contract of a computed invariant.
void expect_exact(const BsgsGroup& G, const RelativeInvariant& inv,
                  std::uint64_t index) {
    ASSERT_EQ(inv.coset_reps.size(), index);
    std::set<std::vector<Monomial>> images;
    std::size_t fixing = 0U;
    for (const auto& s : inv.coset_reps) {
        EXPECT_TRUE(G.contains(s));
        auto img = apply_perm_to_invariant(s, inv.monomials);
        if (img == inv.monomials) ++fixing;
        images.insert(std::move(img));
    }
    EXPECT_EQ(images.size(), index) << "coset images must be distinct";
    EXPECT_EQ(fixing, 1U) << "only the identity coset fixes F";
    EXPECT_EQ(inv.coset_reps[0], identity(G.degree()));
}

TEST(GaloisInvariant, WreathInS8LowDegreeInvariant) {
    const BsgsGroup s8 = symmetric(8U);
    auto cands = maximal_transitive_candidates(AmbientGroup::Symmetric, 8U);
    ASSERT_TRUE(cands.is_ok());
    bool found = false;
    for (const auto& c : cands.value()) {
        if (c.group.order() != 1152U) continue;  // S_4 wr S_2
        found = true;
        auto inv = relative_invariant(s8, c.group, kOps, nullptr);
        ASSERT_TRUE(inv.is_ok()) << inv.error().message;
        // Within-block pairs: the block partition is recoverable from the
        // pair graph, so a degree-2 subset invariant must suffice.
        EXPECT_EQ(inv.value().total_degree, 2U);
        expect_exact(s8, inv.value(), 40320U / 1152U);
    }
    EXPECT_TRUE(found);
}

TEST(GaloisInvariant, AlternatingForcesGaloisMonomialFallback) {
    // A_6 is k-homogeneous for every k, so every tier-1 subset orbit is the
    // full (symmetric) one and the guaranteed Galois monomial must fire:
    // degree n(n−1)/2 = 15, index 2, and Σ_{h∈A_6} h·m* has |A_6| = 360
    // distinct monomials.
    const BsgsGroup s6 = symmetric(6U);
    auto a6 = BsgsGroup::build(
        6U, {cycle(6U, {0U, 1U, 2U}), cycle(6U, {1U, 2U, 3U, 4U, 5U})});
    ASSERT_TRUE(a6.is_ok());
    ASSERT_EQ(a6.value().order(), 360U);
    auto inv = relative_invariant(s6, a6.value(), kOps, nullptr);
    ASSERT_TRUE(inv.is_ok()) << inv.error().message;
    EXPECT_EQ(inv.value().total_degree, 15U);
    EXPECT_EQ(inv.value().monomials.size(), 360U);
    expect_exact(s6, inv.value(), 2U);
}

TEST(GaloisInvariant, CyclicSubgroupExactStabiliser) {
    // H = C_6 ≤ S_6 is far from maximal (index 120): the certificate must
    // still pin Stab_G(F) = C_6 exactly — never an intermediate D_6.
    const BsgsGroup s6 = symmetric(6U);
    auto c6 =
        BsgsGroup::build(6U, {cycle(6U, {0U, 1U, 2U, 3U, 4U, 5U})});
    ASSERT_TRUE(c6.is_ok());
    ASSERT_EQ(c6.value().order(), 6U);
    auto inv = relative_invariant(s6, c6.value(), kOps, nullptr);
    ASSERT_TRUE(inv.is_ok()) << inv.error().message;
    expect_exact(s6, inv.value(), 120U);
}

TEST(GaloisInvariant, TwinAglClassesInA8) {
    // Both A_8-classes of AGL(3,2) (the twin mechanism of Brick 2) must
    // receive an exact index-15 invariant below A_8.
    auto cands =
        maximal_transitive_candidates(AmbientGroup::Alternating, 8U);
    ASSERT_TRUE(cands.is_ok());
    auto a8 = maximal_transitive_candidates(AmbientGroup::Symmetric, 8U);
    ASSERT_TRUE(a8.is_ok());
    ASSERT_FALSE(a8.value().empty());
    const BsgsGroup& a8g = a8.value().front().group;  // A_8 listed first
    ASSERT_EQ(a8g.order(), 20160U);
    std::size_t agl_seen = 0U;
    for (const auto& c : cands.value()) {
        if (c.group.order() != 1344U) continue;  // AGL(3,2)
        ++agl_seen;
        auto inv = relative_invariant(a8g, c.group, kOps, nullptr);
        ASSERT_TRUE(inv.is_ok()) << inv.error().message;
        expect_exact(a8g, inv.value(), 15U);
    }
    EXPECT_EQ(agl_seen, 2U);  // the twin pair
}

TEST(GaloisInvariant, StructuredFailures) {
    const BsgsGroup s6 = symmetric(6U);
    auto a6 = BsgsGroup::build(
        6U, {cycle(6U, {0U, 1U, 2U}), cycle(6U, {1U, 2U, 3U, 4U, 5U})});
    ASSERT_TRUE(a6.is_ok());
    // H ⊄ G: a transposition is not in A_6.
    auto t = BsgsGroup::build(6U, {cycle(6U, {0U, 1U})});
    ASSERT_TRUE(t.is_ok());
    auto bad = relative_invariant(a6.value(), t.value(), kOps, nullptr);
    ASSERT_TRUE(bad.is_error());
    EXPECT_EQ(bad.error().kind, CASErrorKind::InvalidArgument);
    // H = G: index 1.
    auto same = relative_invariant(s6, s6, kOps, nullptr);
    ASSERT_TRUE(same.is_error());
    EXPECT_EQ(same.error().kind, CASErrorKind::InvalidArgument);
    // Ops budget: structured Unimplemented, never a hang.
    auto tiny = relative_invariant(s6, a6.value(), 4U, nullptr);
    ASSERT_TRUE(tiny.is_error());
    EXPECT_EQ(tiny.error().kind, CASErrorKind::Unimplemented);
}

}  // namespace
