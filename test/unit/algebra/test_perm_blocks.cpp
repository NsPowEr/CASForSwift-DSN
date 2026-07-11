// A6 — Block-system tests (Atkinson), Brick 2 of the deg ≥ 8 Stauduhar
// closure. Groups are generated from first principles; expected block
// structures are classical facts (D_4/C_4 diagonal pairs, V_4's three
// pairings, primitivity of S_n/A_n in the natural action).

#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <vector>

#include "../../../src/algebra/perm_blocks_internal.hpp"
#include "../../../src/algebra/perm_bsgs_internal.hpp"
#include "../../../src/algebra/perm_construct_internal.hpp"
#include "../../../src/algebra/perm_group_internal.hpp"
#include "cas/error.hpp"

using namespace cas;
using namespace cas::algebra::permgrp;

namespace {

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

}  // namespace

TEST(PermBlocksTest, CyclicC4HasTheDiagonalPairSystem) {
    auto r = minimal_block_systems(4U, {cycle(4U, {0U, 1U, 2U, 3U})});
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    ASSERT_EQ(r.value().size(), 1U);
    const BlockSystem& sys = r.value()[0];
    EXPECT_EQ(sys.num_blocks, 2U);
    // {0,2} and {1,3}: the only C_4-congruence.
    EXPECT_EQ(sys.block_of,
              (std::vector<std::size_t>{0U, 1U, 0U, 1U}));
}

TEST(PermBlocksTest, KleinV4HasThreePairSystems) {
    Perm a = identity(4U);
    a[0] = 1U; a[1] = 0U; a[2] = 3U; a[3] = 2U;  // (0 1)(2 3)
    Perm b = identity(4U);
    b[0] = 2U; b[2] = 0U; b[1] = 3U; b[3] = 1U;  // (0 2)(1 3)
    auto r = minimal_block_systems(4U, {a, b});
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value().size(), 3U);
    for (const BlockSystem& sys : r.value()) EXPECT_EQ(sys.num_blocks, 2U);
}

TEST(PermBlocksTest, DihedralD4HasOnlyTheDiagonalSystem) {
    // D_4 = ⟨(0 1 2 3), (1 3)⟩: the diagonals {0,2}/{1,3} are the unique
    // non-trivial congruence.
    auto r = minimal_block_systems(
        4U, {cycle(4U, {0U, 1U, 2U, 3U}), cycle(4U, {1U, 3U})});
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    ASSERT_EQ(r.value().size(), 1U);
    EXPECT_EQ(r.value()[0].block_of,
              (std::vector<std::size_t>{0U, 1U, 0U, 1U}));
}

TEST(PermBlocksTest, NaturalSymmetricAndAlternatingArePrimitive) {
    for (std::size_t n = 4U; n <= 9U; ++n) {
        auto ps = is_primitive(n, symmetric_gens(n));
        ASSERT_TRUE(ps.is_ok());
        EXPECT_TRUE(ps.value()) << "S_" << n;
        auto pa = is_primitive(n, alternating_gens(n));
        ASSERT_TRUE(pa.is_ok());
        EXPECT_TRUE(pa.value()) << "A_" << n;
    }
}

TEST(PermBlocksTest, WreathBlockStructureAndBlockAction) {
    // S_2 ≀ S_4 on 8 points: unique minimal system = the four pairs; the
    // induced action on the blocks is the full S_4.
    auto w24 = wreath_gens(2U, 4U);
    ASSERT_TRUE(w24.is_ok());
    auto sys24 = minimal_block_systems(8U, w24.value());
    ASSERT_TRUE(sys24.is_ok()) << sys24.error().message;
    ASSERT_EQ(sys24.value().size(), 1U);
    EXPECT_EQ(sys24.value()[0].num_blocks, 4U);
    EXPECT_EQ(sys24.value()[0].block_of,
              (std::vector<std::size_t>{0U, 0U, 1U, 1U, 2U, 2U, 3U, 3U}));

    auto act = block_action_gens(8U, w24.value(), sys24.value()[0]);
    ASSERT_TRUE(act.is_ok()) << act.error().message;
    auto top = BsgsGroup::build(4U, act.value());
    ASSERT_TRUE(top.is_ok());
    EXPECT_EQ(top.value().order(), 24U);  // S_4 on the blocks
    EXPECT_TRUE(top.value().is_transitive());

    // S_4 ≀ S_2: unique minimal system = the two blocks of four.
    auto w42 = wreath_gens(4U, 2U);
    ASSERT_TRUE(w42.is_ok());
    auto sys42 = minimal_block_systems(8U, w42.value());
    ASSERT_TRUE(sys42.is_ok());
    ASSERT_EQ(sys42.value().size(), 1U);
    EXPECT_EQ(sys42.value()[0].num_blocks, 2U);

    auto prim = is_primitive(8U, w24.value());
    ASSERT_TRUE(prim.is_ok());
    EXPECT_FALSE(prim.value());
}

TEST(PermBlocksTest, IntransitiveIsRejected) {
    auto r = minimal_block_systems(4U, {cycle(4U, {0U, 1U})});
    ASSERT_FALSE(r.is_ok());
    EXPECT_EQ(r.error().kind, CASErrorKind::InvalidArgument);
}

TEST(PermBlocksTest, NonInvariantPartitionIsRejected) {
    // {{0,1},{2,3}} is not a block system of C_4.
    BlockSystem bad;
    bad.block_of = {0U, 0U, 1U, 1U};
    bad.num_blocks = 2U;
    auto r = block_action_gens(4U, {cycle(4U, {0U, 1U, 2U, 3U})}, bad);
    ASSERT_FALSE(r.is_ok());
    EXPECT_EQ(r.error().kind, CASErrorKind::InvalidArgument);
}
