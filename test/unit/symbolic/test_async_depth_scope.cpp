// F7.0-A3.2 — AsyncDepthScope propagation helper tests.
//
// Verifies the async-aware depth propagation primitive that allows a worker
// thread to inherit the parent thread's simplification_depth (preventing
// thread_local bypass of the recursion budget).

#include <gtest/gtest.h>

// Internal header — includes detail::current_simplify_depth + AsyncDepthScope.
#include "../../../src/symbolic/simplify_impl.hpp"

#include <future>
#include <thread>

using namespace cas::symbolic::detail;

namespace {

TEST(AsyncDepthScope, FreshThreadStartsAtZero) {
    EXPECT_EQ(current_simplify_depth(), 0);
}

TEST(AsyncDepthScope, ScopeSetsAndRestores) {
    EXPECT_EQ(current_simplify_depth(), 0);
    {
        AsyncDepthScope scope(42);
        EXPECT_EQ(current_simplify_depth(), 42);
    }
    EXPECT_EQ(current_simplify_depth(), 0);
}

TEST(AsyncDepthScope, NestedScopesNestCorrectly) {
    {
        AsyncDepthScope outer(10);
        EXPECT_EQ(current_simplify_depth(), 10);
        {
            AsyncDepthScope inner(20);
            EXPECT_EQ(current_simplify_depth(), 20);
        }
        EXPECT_EQ(current_simplify_depth(), 10);
    }
    EXPECT_EQ(current_simplify_depth(), 0);
}

TEST(AsyncDepthScope, WorkerThreadInheritsParentDepth) {
    // Without AsyncDepthScope, the worker would see 0.
    // With it, the worker observes the inherited value.
    constexpr int kParentDepth = 137;
    auto fut = std::async(std::launch::async, [&] {
        EXPECT_EQ(current_simplify_depth(), 0);  // worker starts at 0
        AsyncDepthScope scope(kParentDepth);
        return current_simplify_depth();
    });
    EXPECT_EQ(fut.get(), kParentDepth);
}

TEST(AsyncDepthScope, RestoresAfterExceptionInScope) {
    EXPECT_EQ(current_simplify_depth(), 0);
    try {
        AsyncDepthScope scope(99);
        EXPECT_EQ(current_simplify_depth(), 99);
        throw std::runtime_error("unwind test");
    } catch (const std::runtime_error&) {
        // Scope destructor must have run, restoring depth.
    }
    EXPECT_EQ(current_simplify_depth(), 0);
}

}  // namespace
