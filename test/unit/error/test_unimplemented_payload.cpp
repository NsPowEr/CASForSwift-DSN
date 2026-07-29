#include <gtest/gtest.h>
#include "cas/error.hpp"
#include "cas/error_helpers.hpp"

using namespace cas;

TEST(UnimplementedPayloadTest, PayloadIsCorrectlyFormed) {
    auto err = make_unimplemented_error(
        "test_mod", "test_fn", "test_shape",
        error::reason_codes::GENERIC,
        "test_sugg", "T-123", "Brief msg");
    
    ASSERT_TRUE(err.payload.has_value());
    EXPECT_EQ(err.payload->module, "test_mod");
    EXPECT_EQ(err.payload->function, "test_fn");
    EXPECT_EQ(err.payload->reason, error::reason_codes::GENERIC);
    EXPECT_EQ(err.payload->ticket, "T-123");
    EXPECT_EQ(err.message, "Brief msg");
}

TEST(UnimplementedPayloadTest, ResultHelperWorks) {
    auto res = make_unimplemented<int>(
        "mod", "fn", "shape", "CODE", "sugg", "T001");
    
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
    EXPECT_TRUE(res.error().payload.has_value());
}
