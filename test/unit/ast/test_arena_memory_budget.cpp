// F7.0-A3.5 — AstArena memory budget guard tests.
//
// Prevents the OS OOM-killer from terminating the host process: the arena
// refuses heap allocation once the configured budget would be exceeded.

#include <gtest/gtest.h>

#include "cas/ast.hpp"

using namespace cas;

namespace {

TEST(AstArenaMemoryBudget, DefaultUnlimited) {
    AstArena arena;
    EXPECT_EQ(arena.max_memory_budget_bytes(), 0U);
    EXPECT_FALSE(arena.budget_exhausted());
}

TEST(AstArenaMemoryBudget, SetGetBudget) {
    AstArena arena;
    arena.set_max_memory_budget_bytes(1024U * 1024U);
    EXPECT_EQ(arena.max_memory_budget_bytes(), 1024U * 1024U);
}

TEST(AstArenaMemoryBudget, AllocationsTracked) {
    AstArena arena;
    EXPECT_EQ(arena.bytes_allocated(), 0U);
    (void)arena.make<IntegerLit>(BigInt(42));
    EXPECT_GT(arena.bytes_allocated(), 0U);
}

TEST(AstArenaMemoryBudget, BudgetBlocksAllocationAfterCap) {
    AstArena arena;
    // Tiny budget — first DEFAULT_BLOCK_BYTES=64KB block requests 64KB.
    arena.set_max_memory_budget_bytes(100U);
    ExprPtr p = arena.make<IntegerLit>(BigInt(7));
    EXPECT_FALSE(static_cast<bool>(p));   // allocation refused → null
    EXPECT_TRUE(arena.budget_exhausted());
}

TEST(AstArenaMemoryBudget, ResetClearsBudgetExhaustedAndAllocBytes) {
    AstArena arena;
    arena.set_max_memory_budget_bytes(100U);
    (void)arena.make<IntegerLit>(BigInt(7));
    EXPECT_TRUE(arena.budget_exhausted());

    arena.reset();
    EXPECT_FALSE(arena.budget_exhausted());
    EXPECT_EQ(arena.bytes_allocated(), 0U);
    EXPECT_EQ(arena.max_memory_budget_bytes(), 100U);  // budget cap preserved
}

TEST(AstArenaMemoryBudget, RaisingBudgetUnblocks) {
    AstArena arena;
    arena.set_max_memory_budget_bytes(100U);
    ExprPtr a = arena.make<IntegerLit>(BigInt(1));
    EXPECT_FALSE(static_cast<bool>(a));

    arena.set_max_memory_budget_bytes(1024U * 1024U);  // 1MB
    arena.reset();  // restart bookkeeping
    ExprPtr b = arena.make<IntegerLit>(BigInt(1));
    EXPECT_TRUE(static_cast<bool>(b));
}

}  // namespace
