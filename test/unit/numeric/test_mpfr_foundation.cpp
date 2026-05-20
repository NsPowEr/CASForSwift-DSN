// L3-01 STEP 6 — MPFR arbitrary-precision foundation.
// Verifies eval_mpfr() at arbitrary decimal precision for constants
// and symbolic expressions. Foundation already implemented (commit
// 0f6832c BigFloat/MPFR). This test certifies coverage.

#include <gtest/gtest.h>

#include "cas/lexer.hpp"
#include "cas/numeric.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::numeric;

namespace {

class MpfrFoundationTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
};

TEST_F(MpfrFoundationTest, PiAt50Digits) {
    auto e = parse("pi");
    auto r = eval_mpfr(e, 50);
    ASSERT_TRUE(r.is_ok()) << (r.is_error() ? r.error().message : "");
    // First 15 digits of pi: 3.14159265358979
    EXPECT_TRUE(r.value().starts_with("3.14159265358979"))
        << "got: " << r.value();
}

TEST_F(MpfrFoundationTest, PiAt100Digits) {
    auto e = parse("pi");
    auto r = eval_mpfr(e, 100);
    ASSERT_TRUE(r.is_ok());
    // First ~20 digits should match.
    EXPECT_TRUE(r.value().starts_with("3.14159265358979323846"))
        << "got: " << r.value();
    // Length should reflect ~100 significant digits.
    EXPECT_GE(r.value().size(), 95U);
}

TEST_F(MpfrFoundationTest, ESqrtTwoExpOne) {
    auto e = parse("exp(1)");
    auto r = eval_mpfr(e, 30);
    ASSERT_TRUE(r.is_ok());
    // exp(1) = e ≈ 2.71828182845904523536
    EXPECT_TRUE(r.value().starts_with("2.71828182"))
        << "got: " << r.value();
}

TEST_F(MpfrFoundationTest, SqrtTwoAt50Digits) {
    auto e = parse("sqrt(2)");
    auto r = eval_mpfr(e, 50);
    ASSERT_TRUE(r.is_ok());
    // sqrt(2) ≈ 1.41421356237309504880168872420969807856967187537694
    EXPECT_TRUE(r.value().starts_with("1.41421356237309504880"))
        << "got: " << r.value();
}

TEST_F(MpfrFoundationTest, RationalExpression) {
    // 22/7 ≈ pi (Archimedes approx) — verify rational eval
    auto e = parse("22/7");
    auto r = eval_mpfr(e, 20);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(r.value().starts_with("3.142857142857"))
        << "got: " << r.value();
}

TEST_F(MpfrFoundationTest, SineAtPiOverFour) {
    // sin(pi/4) = sqrt(2)/2 ≈ 0.707106781186547524...
    auto e = parse("sin(pi/4)");
    auto r = eval_mpfr(e, 30);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(r.value().starts_with("0.70710678"))
        << "got: " << r.value();
}

}  // namespace
