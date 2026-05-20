// CAS-L2-21 — Branch-cut policy global toggle test.
//
// Verifies that strict_branch_cuts() flag exists, defaults false,
// roundtrips via setter. Future integration: simplifier rules should
// consult this flag to refuse positivity-dependent identities under
// strict mode.

#include <gtest/gtest.h>

#include "cas/symbolic.hpp"

using namespace cas;

namespace {

TEST(BranchCutsGlobalTest, DefaultIsLenient) {
    symbolic::CASContext ctx;
    EXPECT_FALSE(ctx.strict_branch_cuts());
}

TEST(BranchCutsGlobalTest, SetStrictRoundtrips) {
    symbolic::CASContext ctx;
    ctx.set_strict_branch_cuts(true);
    EXPECT_TRUE(ctx.strict_branch_cuts());
    ctx.set_strict_branch_cuts(false);
    EXPECT_FALSE(ctx.strict_branch_cuts());
}

}  // namespace
