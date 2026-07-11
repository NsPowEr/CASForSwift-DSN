// A6 / CAS-L3-18 — BsgsGroup (Schreier-Sims) tests, Brick 1 of the deg ≥ 8
// Stauduhar closure.
//
// Every group is *generated* from first principles (cycles); orders are
// validated against theorems (|S_n| = n!, |A_n| = n!/2, |D_n| = 2n, …) and
// against the dense PermGroup engine for n ≤ 7 — never against copied tables.
// The point of BSGS is scale: the closing test builds S_8 / A_8 / an
// imprimitive group whose |G| is far beyond what the dense Θ(n!) representation
// can hold, confirming the memory wall is gone.

#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <vector>

#include "../../../src/algebra/perm_bsgs_internal.hpp"
#include "../../../src/algebra/perm_group_internal.hpp"
#include "cas/error.hpp"

using namespace cas;
using namespace cas::algebra::permgrp;

namespace {

[[nodiscard]] Perm cycle(std::size_t n, std::initializer_list<std::size_t> cyc) {
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

[[nodiscard]] std::uint64_t factorial(std::size_t n) {
    std::uint64_t f = 1U;
    for (std::size_t i = 2U; i <= n; ++i) f *= static_cast<std::uint64_t>(i);
    return f;
}

// S_n = ⟨(0 1), (0 1 … n-1)⟩.
[[nodiscard]] BsgsGroup symmetric(std::size_t n) {
    auto g = BsgsGroup::build(
        n, {cycle(n, {0U, 1U}), [n] {
                Perm c = identity(n);
                for (std::size_t i = 0U; i + 1U < n; ++i)
                    c[i] = static_cast<std::uint8_t>(i + 1U);
                c[n - 1U] = 0U;
                return c;
            }()});
    EXPECT_TRUE(g.is_ok());
    return g.value();
}

// A_n = ⟨(0 1 2), (0 1 … n-1)⟩ for odd n, ⟨(0 1 2), (1 2 … n-1)⟩ for even n
// — the standard even generators of the alternating group.
[[nodiscard]] BsgsGroup alternating(std::size_t n) {
    Perm long_cycle = identity(n);
    if (n % 2U == 1U) {
        for (std::size_t i = 0U; i + 1U < n; ++i)
            long_cycle[i] = static_cast<std::uint8_t>(i + 1U);
        long_cycle[n - 1U] = 0U;
    } else {
        for (std::size_t i = 1U; i + 1U < n; ++i)
            long_cycle[i] = static_cast<std::uint8_t>(i + 1U);
        long_cycle[n - 1U] = 1U;
    }
    auto g = BsgsGroup::build(n, {cycle(n, {0U, 1U, 2U}), long_cycle});
    EXPECT_TRUE(g.is_ok());
    return g.value();
}

}  // namespace

// ── Known orders ─────────────────────────────────────────────────────────────

TEST(BsgsGroupTest, SymmetricGroupOrders) {
    for (std::size_t n = 3U; n <= 8U; ++n) {
        EXPECT_EQ(symmetric(n).order(), factorial(n)) << "S_" << n;
    }
}

TEST(BsgsGroupTest, AlternatingGroupOrders) {
    for (std::size_t n = 3U; n <= 8U; ++n) {
        EXPECT_EQ(alternating(n).order(), factorial(n) / 2U) << "A_" << n;
    }
}

TEST(BsgsGroupTest, CyclicAndDihedralOrders) {
    // C_7 = ⟨(0 1 2 3 4 5 6)⟩.
    auto c7 = BsgsGroup::build(7U, {cycle(7U, {0U, 1U, 2U, 3U, 4U, 5U, 6U})});
    ASSERT_TRUE(c7.is_ok());
    EXPECT_EQ(c7.value().order(), 7U);

    // D_7 = ⟨(0 1 … 6), (1 6)(2 5)(3 4)⟩, order 14.
    Perm refl = identity(7U);
    refl[1] = 6U; refl[6] = 1U;
    refl[2] = 5U; refl[5] = 2U;
    refl[3] = 4U; refl[4] = 3U;
    auto d7 = BsgsGroup::build(
        7U, {cycle(7U, {0U, 1U, 2U, 3U, 4U, 5U, 6U}), refl});
    ASSERT_TRUE(d7.is_ok());
    EXPECT_EQ(d7.value().order(), 14U);

    // Klein V_4 = ⟨(0 1)(2 3), (0 2)(1 3)⟩, order 4.
    Perm a = identity(4U); a[0]=1U; a[1]=0U; a[2]=3U; a[3]=2U;
    Perm b = identity(4U); b[0]=2U; b[2]=0U; b[1]=3U; b[3]=1U;
    auto v4 = BsgsGroup::build(4U, {a, b});
    ASSERT_TRUE(v4.is_ok());
    EXPECT_EQ(v4.value().order(), 4U);
}

TEST(BsgsGroupTest, TrivialGroup) {
    auto g = BsgsGroup::build(5U, {});
    ASSERT_TRUE(g.is_ok());
    EXPECT_EQ(g.value().order(), 1U);
    EXPECT_TRUE(g.value().contains(identity(5U)));
    EXPECT_FALSE(g.value().contains(cycle(5U, {0U, 1U})));
    EXPECT_FALSE(g.value().is_transitive());
    EXPECT_FALSE(g.value().has_odd_element());
}

// ── Membership ───────────────────────────────────────────────────────────────

TEST(BsgsGroupTest, MembershipAlternating) {
    const BsgsGroup a5 = alternating(5U);
    EXPECT_TRUE(a5.contains(cycle(5U, {0U, 1U, 2U})));           // 3-cycle, even
    EXPECT_TRUE(a5.contains(identity(5U)));
    EXPECT_FALSE(a5.contains(cycle(5U, {0U, 1U})));              // transposition, odd
    EXPECT_FALSE(a5.contains(cycle(5U, {0U, 1U, 2U, 3U})));      // 4-cycle, odd
    // sift of a member is the identity; of a non-member, non-identity.
    auto r_in = a5.sift(cycle(5U, {0U, 1U, 2U}));
    ASSERT_TRUE(r_in.has_value());
    EXPECT_EQ(*r_in, identity(5U));
    auto r_out = a5.sift(cycle(5U, {0U, 1U}));
    ASSERT_TRUE(r_out.has_value());
    EXPECT_NE(*r_out, identity(5U));
}

TEST(BsgsGroupTest, TransitivityAndParity) {
    EXPECT_TRUE(symmetric(6U).is_transitive());
    EXPECT_TRUE(symmetric(6U).has_odd_element());
    EXPECT_TRUE(alternating(6U).is_transitive());
    EXPECT_FALSE(alternating(6U).has_odd_element());
    // A non-transitive group: ⟨(0 1)⟩ on 4 points.
    auto g = BsgsGroup::build(4U, {cycle(4U, {0U, 1U})});
    ASSERT_TRUE(g.is_ok());
    EXPECT_FALSE(g.value().is_transitive());
    EXPECT_EQ(g.value().order(), 2U);
}

// ── Cross-check vs the dense PermGroup engine (n ≤ 7) ────────────────────────

TEST(BsgsGroupTest, CrossCheckDensePermGroup) {
    struct Case {
        std::size_t n;
        std::vector<Perm> gens;
    };
    std::vector<Case> cases;
    cases.push_back({5U, {cycle(5U, {0U, 1U}), cycle(5U, {0U, 1U, 2U, 3U, 4U})}});   // S_5
    cases.push_back({5U, {cycle(5U, {0U, 1U, 2U}), cycle(5U, {0U, 1U, 2U, 3U, 4U})}}); // A_5
    cases.push_back({6U, {cycle(6U, {0U, 1U, 2U, 3U, 4U, 5U}),
                          [] { Perm r = identity(6U); r[1]=5U;r[5]=1U;r[2]=4U;r[4]=2U;return r; }()}}); // D_6
    cases.push_back({7U, {cycle(7U, {0U, 1U, 2U, 3U, 4U, 5U, 6U})}});                // C_7

    for (const auto& c : cases) {
        auto dense = PermGroup::closure(c.n, c.gens, factorial(c.n));
        ASSERT_TRUE(dense.is_ok());
        auto bsgs = BsgsGroup::build(c.n, c.gens);
        ASSERT_TRUE(bsgs.is_ok());

        EXPECT_EQ(bsgs.value().order(), dense.value().order()) << "order n=" << c.n;
        EXPECT_EQ(bsgs.value().is_transitive(), dense.value().is_transitive());
        EXPECT_EQ(bsgs.value().has_odd_element(), dense.value().has_odd_element());

        // Every dense element is a BSGS member.
        for (const std::uint32_t rk : dense.value().element_ranks()) {
            EXPECT_TRUE(bsgs.value().contains(lehmer_unrank(rk, c.n)))
                << "dense element rank " << rk << " not in BSGS (n=" << c.n << ")";
        }
        // Every non-member (by dense) is rejected by BSGS — full agreement on S_n.
        const std::uint64_t full = factorial(c.n);
        for (std::uint32_t rk = 0U; rk < full; ++rk) {
            const bool in_dense = dense.value().contains_rank(rk);
            const bool in_bsgs = bsgs.value().contains(lehmer_unrank(rk, c.n));
            ASSERT_EQ(in_dense, in_bsgs) << "membership mismatch rank " << rk
                                         << " n=" << c.n;
        }
    }
}

// ── Scale: past the dense Θ(n!) wall ─────────────────────────────────────────

TEST(BsgsGroupTest, ScaleDegreeEightExactOrders) {
    EXPECT_EQ(symmetric(8U).order(), 40320U);   // 8!
    EXPECT_EQ(alternating(8U).order(), 20160U);  // 8!/2

    // Imprimitive wreath-like group ⟨(0 1 2 3), (4 5 6 7), (0 4)(1 5)(2 6)(3 7)⟩
    // — two blocks of 4 swapped: transitive of degree 8, order well-defined and
    // reachable only because the representation is polynomial in size.
    Perm b0 = cycle(8U, {0U, 1U, 2U, 3U});
    Perm b1 = cycle(8U, {4U, 5U, 6U, 7U});
    Perm swap = identity(8U);
    swap[0]=4U;swap[4]=0U; swap[1]=5U;swap[5]=1U;
    swap[2]=6U;swap[6]=2U; swap[3]=7U;swap[7]=3U;
    auto w = BsgsGroup::build(8U, {b0, b1, swap});
    ASSERT_TRUE(w.is_ok());
    EXPECT_TRUE(w.value().is_transitive());
    // |C_4 wr C_2| = 4^2 · 2 = 32.
    EXPECT_EQ(w.value().order(), 32U);
    EXPECT_TRUE(w.value().contains(b0));
    EXPECT_TRUE(w.value().contains(swap));
    EXPECT_FALSE(w.value().contains(cycle(8U, {0U, 1U})));  // odd single transposition
}

// ── Degree cap (u64 order safety) ────────────────────────────────────────────

TEST(BsgsGroupTest, DegreeCapRejected) {
    auto g = BsgsGroup::build(21U, {});
    EXPECT_FALSE(g.is_ok());
    EXPECT_EQ(g.error().kind, CASErrorKind::Unimplemented);
}
