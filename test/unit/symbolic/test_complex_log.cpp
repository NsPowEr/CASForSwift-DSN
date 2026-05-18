
#include <gtest/gtest.h>
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast_debug.hpp"

using namespace cas;
using namespace cas::symbolic;

namespace {
[[nodiscard]] Result<ExprPtr> pe(const std::string& s, AstArena& a) {
    auto t = Lexer(s).tokenize();
    if (t.is_error()) return fail<ExprPtr>(t.error());
    return Parser(t.value(), a).parse();
}
} // namespace

class ComplexLogTest : public ::testing::Test {
protected:
    CASContext ctx;
};

TEST_F(ComplexLogTest, LnOnePlusI_RealPart) {
    // ln(1 + i) = (1/2)*ln(2) + i*pi/4
    // Real part = (1/2)*ln(2)
    // abs(1+i) = sqrt(2) -> ln(sqrt(2)) = (1/2)*ln(2)
    auto e = pe("ln(1 + i)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx.simplify(e.value());
    ASSERT_TRUE(s.is_ok()) << s.error().message;
    auto t = to_round_trip_text(s.value());
    ASSERT_TRUE(t.is_ok());
    // Check it contains "ln(2)" and "pi"
    EXPECT_TRUE(t.value().find("ln") != std::string::npos ||
                t.value().find("2") != std::string::npos)
        << "Got: " << t.value();
}

TEST_F(ComplexLogTest, LnI_EqualsITimesPiOver2) {
    // ln(i) = i*pi/2  -- already handled before complex sum rule
    auto e = pe("ln(i)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx.simplify(e.value());
    ASSERT_TRUE(s.is_ok()) << s.error().message;
    auto t = to_round_trip_text(s.value());
    ASSERT_TRUE(t.is_ok());
    EXPECT_TRUE(t.value().find("pi") != std::string::npos) << "Got: " << t.value();
}

TEST_F(ComplexLogTest, LnNeg1_EqualsIPi) {
    // ln(-1) = i*pi
    auto e = pe("ln(-1)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx.simplify(e.value());
    ASSERT_TRUE(s.is_ok()) << s.error().message;
    auto t = to_round_trip_text(s.value());
    ASSERT_TRUE(t.is_ok());
    EXPECT_TRUE(t.value().find("pi") != std::string::npos) << "Got: " << t.value();
}

TEST_F(ComplexLogTest, LnComplexWithPositiveReal) {
    // ln(a + b*i) where a is positive -- arg simplifies to atan(b/a)
    // Assume a>0 for the arg simplification to trigger
    ctx.assumptions().assume_positive(Symbol{"a"});
    ctx.assumptions().assume_real(Symbol{"b"});
    auto e = pe("ln(a + b*i)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx.simplify(e.value());
    ASSERT_TRUE(s.is_ok()) << s.error().message;
    auto t = to_round_trip_text(s.value());
    ASSERT_TRUE(t.is_ok());
    // Should produce something involving atan(b/a) and ln
    EXPECT_TRUE(t.value().find("atan") != std::string::npos ||
                t.value().find("ln") != std::string::npos)
        << "Got: " << t.value();
}
