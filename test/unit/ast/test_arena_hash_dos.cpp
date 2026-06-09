// F7.0-A3.7 — Hash-DoS defence sanity tests.
//
// The collision-chain detector is preventive: normal CAS workloads stay
// far below MAX_COLLISION_CHAIN (128) due to std::unordered_set's
// load-factor-driven rehash. These tests verify the API surface and the
// reset behaviour, not adversarial saturation (which would require
// thousands of crafted same-hash nodes).

#include <gtest/gtest.h>

#include "cas/ast.hpp"

using namespace cas;

namespace {

TEST(AstArenaHashDos, DefaultFlagClear) {
    AstArena arena;
    EXPECT_FALSE(arena.hash_dos_detected());
}

TEST(AstArenaHashDos, ManyAllocationsDontTripFalsePositive) {
    AstArena arena;
    for (int i = 0; i < 1000; ++i) {
        (void)arena.make<IntegerLit>(BigInt(i));
        (void)arena.make<Symbol>("s_" + std::to_string(i));
    }
    EXPECT_FALSE(arena.hash_dos_detected());
}

TEST(AstArenaHashDos, ClearResets) {
    AstArena arena;
    arena.clear_hash_dos_flag();
    EXPECT_FALSE(arena.hash_dos_detected());
    arena.reset();
    EXPECT_FALSE(arena.hash_dos_detected());
}

TEST(AstArenaHashDos, MaxCollisionChainExposedConstant) {
    EXPECT_EQ(AstArena::MAX_COLLISION_CHAIN, 128U);
}

}  // namespace
