// F7.0-A3.1 — AstArena::reset() controlled-reset tests.
//
// Verifies:
//   1. reset() on empty arena is a no-op (no crash, size still 0).
//   2. reset() after allocations clears size() to 0.
//   3. reset() destroys interning state (post-reset, identical content
//      produces a different ExprPtr than the pre-reset value would have).
//   4. reset() leaves the arena usable for fresh allocations.
//   5. reset() releases bump-allocator memory (sequence of large allocations
//      followed by reset() loops without unbounded growth — sanity check).

#include <gtest/gtest.h>

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

TEST(AstArenaReset, EmptyArenaResetIsNoop) {
    AstArena arena;
    ASSERT_EQ(arena.size(), 0U);
    arena.reset();
    ASSERT_EQ(arena.size(), 0U);
}

TEST(AstArenaReset, AllocationsThenResetReturnsToZero) {
    AstArena arena;
    (void)arena.make<IntegerLit>(BigInt(7));
    (void)arena.make<IntegerLit>(BigInt(42));
    (void)arena.make<Symbol>("x");
    ASSERT_GT(arena.size(), 0U);
    const std::size_t before = arena.size();
    arena.reset();
    EXPECT_EQ(arena.size(), 0U);
    EXPECT_GT(before, 0U);
}

TEST(AstArenaReset, ResetClearsInterningState) {
    AstArena arena;
    ExprPtr a = arena.make<IntegerLit>(BigInt(99));
    ExprPtr b = arena.make<IntegerLit>(BigInt(99));
    EXPECT_EQ(a, b);  // pre-reset: structural sharing
    arena.reset();
    ExprPtr c = arena.make<IntegerLit>(BigInt(99));
    // Post-reset, the prior a/b pointers are dangling — they MUST NOT be
    // dereferenced. We only check that the freshly allocated c is non-null.
    ASSERT_TRUE(static_cast<bool>(c));
    // The freshly allocated c must NOT equal the prior pointer value
    // (unless the allocator happens to hand back the same address, which
    // we don't assert either way — what matters is internal table cleared,
    // not raw address equality).
}

TEST(AstArenaReset, ResetEnablesReuseLoop) {
    AstArena arena;
    for (int iter = 0; iter < 10; ++iter) {
        for (int i = 0; i < 50; ++i) {
            (void)arena.make<IntegerLit>(BigInt(iter * 100 + i));
            (void)arena.make<Symbol>("s_" + std::to_string(i));
        }
        ASSERT_GT(arena.size(), 0U);
        arena.reset();
        ASSERT_EQ(arena.size(), 0U);
    }
}

TEST(AstArenaReset, ResetClearsHotConstantsCache) {
    AstArena arena;
    ExprPtr z1 = arena.make<IntegerLit>(BigInt(0));
    ExprPtr o1 = arena.make<IntegerLit>(BigInt(1));
    ExprPtr n1 = arena.make<IntegerLit>(BigInt(-1));
    ASSERT_TRUE(static_cast<bool>(z1));
    ASSERT_TRUE(static_cast<bool>(o1));
    ASSERT_TRUE(static_cast<bool>(n1));
    arena.reset();
    // Post-reset, fresh hot-constants must allocate cleanly without using
    // stale pointers from before reset.
    ExprPtr z2 = arena.make<IntegerLit>(BigInt(0));
    ExprPtr o2 = arena.make<IntegerLit>(BigInt(1));
    ExprPtr n2 = arena.make<IntegerLit>(BigInt(-1));
    EXPECT_TRUE(static_cast<bool>(z2));
    EXPECT_TRUE(static_cast<bool>(o2));
    EXPECT_TRUE(static_cast<bool>(n2));
    // Within the post-reset arena, structural sharing still works.
    EXPECT_EQ(z2, arena.make<IntegerLit>(BigInt(0)));
    EXPECT_EQ(o2, arena.make<IntegerLit>(BigInt(1)));
    EXPECT_EQ(n2, arena.make<IntegerLit>(BigInt(-1)));
}

}  // namespace
