// A6 / CAS-L3-18 — Exact permutation-group engine tests.
//
// Every group here is *generated* from first-principles constructions
// (cycles, affine maps over F_p, Möbius maps over P¹(F_5)) and validated
// against theorems (orders, transitivity, orbit counts), never against
// copied tables.

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <vector>

#include "../../../src/algebra/perm_group_internal.hpp"
#include "cas/error.hpp"

using namespace cas;
using namespace cas::algebra::permgrp;

namespace {

// Build a permutation of degree n from a single cycle (a0 a1 … ak).
[[nodiscard]] Perm cycle(std::size_t n, std::initializer_list<std::size_t> cyc) {
    Perm p = identity(n);
    if (cyc.size() >= 2U) {
        auto it = cyc.begin();
        std::size_t first = *it;
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

TEST(PermGroupTest, ComposeInverseIdentity) {
    const Perm s = cycle(5U, {0U, 1U});           // transposition (01)
    const Perm c = cycle(5U, {0U, 1U, 2U, 3U, 4U});  // 5-cycle
    EXPECT_EQ(compose(s, inverse(s)), identity(5U));
    EXPECT_EQ(compose(inverse(c), c), identity(5U));
    // (01)∘(01234): apply cycle first → 0→1→0, so image of 0 is 0.
    const Perm sc = compose(s, c);
    EXPECT_EQ(sc[0], 0U);
    EXPECT_EQ(sc[1], 2U);
}

TEST(PermGroupTest, ParityAndCycleType) {
    EXPECT_TRUE(is_odd(cycle(5U, {0U, 1U})));               // transposition
    EXPECT_FALSE(is_odd(cycle(5U, {0U, 1U, 2U})));          // 3-cycle even
    EXPECT_FALSE(is_odd(cycle(5U, {0U, 1U, 2U, 3U, 4U})));  // 5-cycle even
    EXPECT_TRUE(is_odd(cycle(6U, {0U, 1U, 2U, 3U, 4U, 5U})));  // 6-cycle odd
    const auto ct = cycle_type(compose(cycle(6U, {0U, 1U}), cycle(6U, {2U, 3U, 4U})));
    ASSERT_EQ(ct.size(), 3U);
    EXPECT_EQ(ct[0], 3U);
    EXPECT_EQ(ct[1], 2U);
    EXPECT_EQ(ct[2], 1U);
}

TEST(PermGroupTest, LehmerRankRoundTripS4) {
    // Exhaustive round-trip over all 24 elements of S_4, and bijectivity.
    std::vector<bool> hit(24U, false);
    for (std::uint32_t r = 0U; r < 24U; ++r) {
        const Perm p = lehmer_unrank(r, 4U);
        EXPECT_TRUE(is_valid_perm(p));
        EXPECT_EQ(lehmer_rank(p), r);
        hit[r] = true;
    }
    EXPECT_TRUE(std::all_of(hit.begin(), hit.end(), [](bool b) { return b; }));
    EXPECT_EQ(lehmer_rank(identity(4U)), 0U);
}

TEST(PermGroupTest, ClosureS5FullSymmetric) {
    auto g = PermGroup::closure(
        5U, {cycle(5U, {0U, 1U}), cycle(5U, {0U, 1U, 2U, 3U, 4U})},
        factorial_u64(5U));
    ASSERT_TRUE(g.is_ok());
    EXPECT_EQ(g.value().order(), 120U);
    EXPECT_TRUE(g.value().is_transitive());
    EXPECT_TRUE(g.value().has_odd_element());
}

TEST(PermGroupTest, ClosureA5Alternating) {
    // ⟨(012), (01234)⟩ — both even → A_5, order 60.
    auto g = PermGroup::closure(
        5U, {cycle(5U, {0U, 1U, 2U}), cycle(5U, {0U, 1U, 2U, 3U, 4U})},
        factorial_u64(5U));
    ASSERT_TRUE(g.is_ok());
    EXPECT_EQ(g.value().order(), 60U);
    EXPECT_TRUE(g.value().is_transitive());
    EXPECT_FALSE(g.value().has_odd_element());
}

TEST(PermGroupTest, ClosureAffineGroupF7) {
    // AGL(1,7) = ⟨x→x+1, x→3x⟩ over F_7 (3 is a primitive root mod 7).
    // Order = 7·6 = 42, sharply 2-transitive.
    Perm add1(7U), mul3(7U);
    for (std::size_t i = 0U; i < 7U; ++i) {
        add1[i] = static_cast<std::uint8_t>((i + 1U) % 7U);
        mul3[i] = static_cast<std::uint8_t>((3U * i) % 7U);
    }
    auto g = PermGroup::closure(7U, {add1, mul3}, factorial_u64(7U));
    ASSERT_TRUE(g.is_ok());
    EXPECT_EQ(g.value().order(), 42U);
    EXPECT_TRUE(g.value().is_transitive());
    // 2-transitivity ⇒ single orbit on 2-subsets (C(7,2)=21).
    const auto pair_orbits = g.value().orbit_lengths_on_ksubsets(2U);
    ASSERT_EQ(pair_orbits.size(), 1U);
    EXPECT_EQ(pair_orbits[0], 21U);
}

TEST(PermGroupTest, ClosurePGL25OnProjectiveLine) {
    // PGL(2,5) acting on P¹(F_5) = {0,1,2,3,4,∞(=5)} via Möbius maps:
    //   z→z+1, z→2z, z→1/z.  |PGL(2,5)| = 120; the exotic transitive
    //   embedding of S_5 into S_6 (action of S_5 on its six Sylow-5s).
    Perm add1 = identity(6U), mul2 = identity(6U), inv = identity(6U);
    for (std::size_t z = 0U; z < 5U; ++z) {
        add1[z] = static_cast<std::uint8_t>((z + 1U) % 5U);
        mul2[z] = static_cast<std::uint8_t>((2U * z) % 5U);
    }
    add1[5] = 5U;  // ∞ fixed by translation
    mul2[5] = 5U;  // ∞ fixed by scaling
    // z→1/z: 0↔∞, 1→1, 2→3 (2·3=6≡1), 3→2, 4→4 (4·4=16≡1).
    inv[0] = 5U;
    inv[5] = 0U;
    inv[2] = 3U;
    inv[3] = 2U;
    auto g = PermGroup::closure(6U, {add1, mul2, inv}, factorial_u64(6U));
    ASSERT_TRUE(g.is_ok());
    EXPECT_EQ(g.value().order(), 120U);
    EXPECT_TRUE(g.value().is_transitive());
    // PGL(2,5) is 3-transitive on P¹ ⇒ single orbit on 2-subsets (15).
    const auto pair_orbits = g.value().orbit_lengths_on_ksubsets(2U);
    ASSERT_EQ(pair_orbits.size(), 1U);
    EXPECT_EQ(pair_orbits[0], 15U);
}

TEST(PermGroupTest, OrbitLengthsCyclicC6OnPairs) {
    // C_6 = ⟨(012345)⟩ on 2-subsets: distance classes d=1,2 give orbits of
    // length 6, the diameter class d=3 gives length 3. Total 15 = C(6,2).
    auto g = PermGroup::closure(6U, {cycle(6U, {0U, 1U, 2U, 3U, 4U, 5U})},
                                factorial_u64(6U));
    ASSERT_TRUE(g.is_ok());
    EXPECT_EQ(g.value().order(), 6U);
    const auto orb = g.value().orbit_lengths_on_ksubsets(2U);
    ASSERT_EQ(orb.size(), 3U);
    EXPECT_EQ(orb[0], 6U);
    EXPECT_EQ(orb[1], 6U);
    EXPECT_EQ(orb[2], 3U);
}

TEST(PermGroupTest, ClosureBudgetExceededIsStructuredError) {
    // S_5 has order 120 > 50 → structured Unimplemented, no silent truncation.
    auto g = PermGroup::closure(
        5U, {cycle(5U, {0U, 1U}), cycle(5U, {0U, 1U, 2U, 3U, 4U})}, 50U);
    ASSERT_TRUE(g.is_error());
    EXPECT_EQ(g.error().kind, CASErrorKind::Unimplemented);
}

TEST(PermGroupTest, CycleTypeDistributionS3) {
    auto g = PermGroup::closure(
        3U, {cycle(3U, {0U, 1U}), cycle(3U, {0U, 1U, 2U})}, factorial_u64(3U));
    ASSERT_TRUE(g.is_ok());
    const auto dist = g.value().cycle_type_distribution();
    // S_3: 1×identity (1,1,1), 3×transpositions (2,1), 2×3-cycles (3).
    ASSERT_EQ(dist.size(), 3U);
    std::size_t total = 0U;
    for (const auto& [ct, cnt] : dist) total += cnt;
    EXPECT_EQ(total, 6U);
    for (const auto& [ct, cnt] : dist) {
        if (ct == std::vector<std::size_t>{3U}) EXPECT_EQ(cnt, 2U);
        if (ct == std::vector<std::size_t>{2U, 1U}) EXPECT_EQ(cnt, 3U);
        if (ct == std::vector<std::size_t>{1U, 1U, 1U}) EXPECT_EQ(cnt, 1U);
    }
}

TEST(PermGroupTest, ContainsAndSubgroupRelations) {
    auto s5 = PermGroup::closure(
        5U, {cycle(5U, {0U, 1U}), cycle(5U, {0U, 1U, 2U, 3U, 4U})},
        factorial_u64(5U));
    auto a5 = PermGroup::closure(
        5U, {cycle(5U, {0U, 1U, 2U}), cycle(5U, {0U, 1U, 2U, 3U, 4U})},
        factorial_u64(5U));
    ASSERT_TRUE(s5.is_ok());
    ASSERT_TRUE(a5.is_ok());
    EXPECT_TRUE(s5.value().contains_group(a5.value()));
    EXPECT_FALSE(a5.value().contains_group(s5.value()));
    EXPECT_FALSE(a5.value().contains(cycle(5U, {0U, 1U})));
    EXPECT_TRUE(a5.value().contains(cycle(5U, {0U, 1U, 2U})));
    EXPECT_FALSE(a5.value().equals(s5.value()));
}

}  // namespace
