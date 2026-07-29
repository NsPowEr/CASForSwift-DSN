#include <gtest/gtest.h>
#include "cas/error.hpp"
#include "cas/error_helpers.hpp"
#include "cas/result.hpp"
#include "cas/unimplemented_info.hpp"

using namespace cas;

TEST(UnimplementedInfoTest, InfoStructCreationAndFields) {
    UnimplementedInfo info{
        .module      = "test_mod",
        .function    = "test_fn",
        .input_shape = "test_shape",
        .reason      = "TEST_REASON",
        .suggestion  = "test_sugg",
        .ticket      = "T-100"
    };

    EXPECT_EQ(info.module, "test_mod");
    EXPECT_EQ(info.function, "test_fn");
    EXPECT_EQ(info.input_shape, "test_shape");
    EXPECT_EQ(info.reason, "TEST_REASON");
    EXPECT_EQ(info.suggestion, "test_sugg");
    EXPECT_EQ(info.ticket, "T-100");
    EXPECT_EQ(info.reason_code(), "TEST_REASON");
    EXPECT_EQ(info.ticket_id(), "T-100");
}

TEST(UnimplementedInfoTest, MakeUnimplementedErrorWithInfo) {
    UnimplementedInfo info{
        .module      = "calculus",
        .function    = "solve_ode",
        .input_shape = "non_linear",
        .reason      = "ODE_UNSUPPORTED",
        .suggestion  = "Use numeric solver",
        .ticket      = "A23"
    };

    auto err = make_unimplemented_error(info);
    ASSERT_EQ(err.kind, CASErrorKind::Unimplemented);
    ASSERT_TRUE(err.payload.has_value());
    EXPECT_EQ(err.payload->module, "calculus");
    EXPECT_EQ(err.payload->function, "solve_ode");
    EXPECT_EQ(err.payload->input_shape, "non_linear");
    EXPECT_EQ(err.payload->reason, "ODE_UNSUPPORTED");
    EXPECT_EQ(err.payload->suggestion, "Use numeric solver");
    EXPECT_EQ(err.payload->ticket, "A23");

    std::string user_msg = err.format_user_message();
    EXPECT_NE(user_msg.find("[Unimplemented] module=calculus function=solve_ode"), std::string::npos);
    EXPECT_NE(user_msg.find("Input shape: non_linear"), std::string::npos);
    EXPECT_NE(user_msg.find("Reason: ODE_UNSUPPORTED"), std::string::npos);
    EXPECT_NE(user_msg.find("Suggestion: Use numeric solver"), std::string::npos);
    EXPECT_NE(user_msg.find("Ticket: A23"), std::string::npos);
}

TEST(UnimplementedInfoTest, ResultUnimplementedStaticMethod) {
    UnimplementedInfo info{
        .module      = "algebra",
        .function    = "factor_poly",
        .input_shape = "deg_100",
        .reason      = "DEGREE_TOO_LARGE",
        .suggestion  = "Use modular factorizer",
        .ticket      = "A23-2"
    };

    auto res = Result<int>::unimplemented(info);
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
    ASSERT_TRUE(res.error().payload.has_value());
    EXPECT_EQ(res.error().payload->module, "algebra");
}

static Result<int> inner_propagate_test() {
    return Result<int>::unimplemented(UnimplementedInfo{
        .module      = "core",
        .function    = "inner",
        .input_shape = "shape",
        .reason      = "REASON",
        .suggestion  = "SUGG",
        .ticket      = "T-PROP"
    });
}

static Result<int> outer_propagate_test() {
    auto r = inner_propagate_test();
    if (r.is_error()) return fail<int>(r.error());
    return ok(42);
}

TEST(UnimplementedInfoTest, PropagationPreservesPayload) {
    auto res = outer_propagate_test();
    ASSERT_TRUE(res.is_error());
    ASSERT_TRUE(res.error().payload.has_value());
    EXPECT_EQ(res.error().payload->ticket, "T-PROP");
}
