// F1.3-NEW — Arena shard-lock + make_fresh_symbol thread-safety tests.
//
// Verifies:
//   1. 4 threads × 1000 random intern calls → no races, no duplicate nodes.
//   2. make_fresh_symbol with same prefix → 100 distinct names from one context.
//   3. make_fresh_symbol never collides with previously define()-d names.

#include <gtest/gtest.h>

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"

#include <atomic>
#include <set>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

using namespace cas;
using namespace cas::symbolic;

namespace {

// -------------------------------------------------------------------------
// Test 1: concurrent intern — 4 threads each intern 1000 symbols.
// All interned Symbol("x"), Symbol("y"), ... must be unique ExprPtr values
// (since the arena guarantees one canonical pointer per structurally-equal node).
// -------------------------------------------------------------------------
TEST(ArenaShardLocksTest, ConcurrentInternNoDuplicatePointers) {
    AstArena arena;

    // 10 distinct symbol names shared across all threads.
    const std::vector<std::string> names = {
        "a", "b", "c", "d", "e", "f", "g", "h", "i_sym", "j"
    };

    constexpr int N_THREADS = 4;
    constexpr int N_ITER    = 1000;

    // Each thread collects the ExprPtrs it gets back.
    std::vector<std::vector<ExprPtr>> thread_results(N_THREADS);
    for (auto& r : thread_results) r.reserve(N_ITER);

    std::vector<std::thread> threads;
    threads.reserve(N_THREADS);

    for (int t = 0; t < N_THREADS; ++t) {
        threads.emplace_back([&, t]() {
            for (int i = 0; i < N_ITER; ++i) {
                const std::string& name = names[static_cast<std::size_t>(i) % names.size()];
                ExprPtr ptr = arena.make<Symbol>(name);
                thread_results[static_cast<std::size_t>(t)].push_back(ptr);
            }
        });
    }

    for (auto& th : threads) th.join();

    // Verify: for the same symbol name all threads must have gotten the same pointer.
    for (const auto& name : names) {
        ExprPtr canonical;
        bool first = true;
        for (int t = 0; t < N_THREADS; ++t) {
            for (int i = 0; i < N_ITER; ++i) {
                const std::string& n = names[static_cast<std::size_t>(i) % names.size()];
                if (n != name) continue;
                ExprPtr ptr = thread_results[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)];
                if (first) {
                    canonical = ptr;
                    first = false;
                } else {
                    EXPECT_EQ(ptr, canonical)
                        << "Symbol(\"" << name << "\") returned different pointers from different threads";
                }
            }
        }
    }
}

// -------------------------------------------------------------------------
// Test 2: concurrent intern of IntegerLit values — hot-cache atoms are safe.
// -------------------------------------------------------------------------
TEST(ArenaShardLocksTest, ConcurrentInternIntegerLitHotCache) {
    AstArena arena;

    constexpr int N_THREADS = 4;
    constexpr int N_ITER    = 500;

    std::vector<std::thread> threads;
    threads.reserve(N_THREADS);
    std::atomic<int> failures{0};

    for (int t = 0; t < N_THREADS; ++t) {
        threads.emplace_back([&]() {
            for (int i = 0; i < N_ITER; ++i) {
                ExprPtr z1 = arena.make<IntegerLit>(BigInt(0));
                ExprPtr z2 = arena.make<IntegerLit>(BigInt(0));
                if (z1 != z2) failures.fetch_add(1, std::memory_order_relaxed);

                ExprPtr o1 = arena.make<IntegerLit>(BigInt(1));
                ExprPtr o2 = arena.make<IntegerLit>(BigInt(1));
                if (o1 != o2) failures.fetch_add(1, std::memory_order_relaxed);
            }
        });
    }

    for (auto& th : threads) th.join();

    EXPECT_EQ(failures.load(), 0)
        << "IntegerLit hot-cache returned different pointers under concurrent access";
}

// -------------------------------------------------------------------------
// Test 3: make_fresh_symbol — 100 calls with same prefix → all distinct.
// -------------------------------------------------------------------------
TEST(MakeFreshSymbolTest, HundredCallsDistinct) {
    CASContext ctx;

    constexpr int N = 100;
    std::set<std::string> seen;

    for (int i = 0; i < N; ++i) {
        Symbol s = ctx.make_fresh_symbol("C");
        EXPECT_TRUE(seen.insert(s.name).second)
            << "make_fresh_symbol(\"C\") produced duplicate name: " << s.name;
    }

    EXPECT_EQ(static_cast<int>(seen.size()), N);
}

// -------------------------------------------------------------------------
// Test 4: make_fresh_symbol never collides with user-defined variables.
// -------------------------------------------------------------------------
TEST(MakeFreshSymbolTest, NoCollisionWithDefinedVariables) {
    CASContext ctx;

    // Pre-define variables that look like fresh-symbol candidates.
    AstArena& arena = ctx.arena();
    ExprPtr val = arena.make<IntegerLit>(BigInt(42));
    ctx.define(Symbol("C_1"), val);
    ctx.define(Symbol("C_2"), val);
    ctx.define(Symbol("C_3"), val);

    // Now generate fresh symbols — none should be C_1, C_2, or C_3.
    for (int i = 0; i < 20; ++i) {
        Symbol s = ctx.make_fresh_symbol("C");
        EXPECT_NE(s.name, "C_1") << "Fresh symbol collided with defined C_1";
        EXPECT_NE(s.name, "C_2") << "Fresh symbol collided with defined C_2";
        EXPECT_NE(s.name, "C_3") << "Fresh symbol collided with defined C_3";
    }
}

}  // namespace
