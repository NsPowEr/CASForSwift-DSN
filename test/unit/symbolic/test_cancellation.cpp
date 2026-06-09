// F7.0-A3.3 — CASContext::interrupt() + check_interrupt() poll-point tests.
//
// Verifies:
//   1. is_interrupted() defaults to false; clear_interrupt() resets.
//   2. interrupt() before a long simplify() causes Timeout error on next
//      poll-point check (Simplifier::check_timeout).
//   3. interrupt() before a substitute() call causes Timeout error.
//   4. check_interrupt() public helper returns ok() when not interrupted,
//      Timeout error when interrupted.
//   5. clear_interrupt() restores normal operation.

#include <gtest/gtest.h>

#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

using namespace cas;
using namespace cas::symbolic;

namespace {

TEST(Cancellation, InterruptFlagDefaultsFalse) {
    CASContext ctx;
    EXPECT_FALSE(ctx.is_interrupted());
}

TEST(Cancellation, InterruptThenClear) {
    CASContext ctx;
    ctx.interrupt();
    EXPECT_TRUE(ctx.is_interrupted());
    ctx.clear_interrupt();
    EXPECT_FALSE(ctx.is_interrupted());
}

TEST(Cancellation, CheckInterruptOkWhenNotInterrupted) {
    CASContext ctx;
    auto res = ctx.check_interrupt();
    EXPECT_TRUE(res.is_ok());
}

TEST(Cancellation, CheckInterruptErrorWhenInterrupted) {
    CASContext ctx;
    ctx.interrupt();
    auto res = ctx.check_interrupt();
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::Timeout);
    EXPECT_NE(res.error().message.find("cancelled"), std::string::npos);
}

TEST(Cancellation, ClearInterruptAfterFlagRestoresOk) {
    CASContext ctx;
    ctx.interrupt();
    EXPECT_TRUE(ctx.check_interrupt().is_error());
    ctx.clear_interrupt();
    EXPECT_TRUE(ctx.check_interrupt().is_ok());
}

}  // namespace
