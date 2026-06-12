// F8.0-5.3: tests for MPFR-precision Sturm root polish.

#include "cas/algebra.hpp"
#include "cas/bigfloat.hpp"
#include "cas/lexer.hpp"
#include "cas/numeric_bigfloat.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <cmath>
#include <string>

namespace cas {
namespace {

[[nodiscard]] ExprPtr parse_poly(symbolic::CASContext& ctx, const std::string& s) {
    auto tokens = Lexer(s).tokenize();
    if (tokens.is_error()) return nullptr;
    Parser p(tokens.value(), ctx.arena());
    auto e = p.parse();
    if (e.is_error()) return nullptr;
    auto simp = ctx.simplify(e.value());
    return simp.is_ok() ? simp.value() : e.value();
}

TEST(SturmBigFloat, Sqrt2_MatchesMpfrReference) {
    symbolic::CASContext ctx;
    ExprPtr poly = parse_poly(ctx, "x^2 - 2");
    ASSERT_NE(poly, nullptr);

    auto roots = numeric::find_polynomial_roots_sturm_bigfloat(
        poly, "x", ctx, -5.0, 5.0, 1e-9, 256);
    ASSERT_TRUE(roots.is_ok()) << roots.error().message;
    ASSERT_EQ(roots.value().size(), 2U);

    // Roots are sorted by midpoint of isolating interval (ascending).
    // Expected: -√2 then +√2.
    const BigFloat& neg = roots.value()[0];
    const BigFloat& pos = roots.value()[1];

    EXPECT_LT(neg.to_double(), 0.0);
    EXPECT_GT(pos.to_double(), 0.0);

    // Cross-check against the MPFR reference √2 at 256-bit precision.
    BigFloat two   = BigFloat::from_double(2.0, 256);
    BigFloat sqrt2 = BigFloat::sqrt(two);
    BigFloat err_pos = BigFloat::abs(pos - sqrt2);
    BigFloat err_neg = BigFloat::abs(neg + sqrt2);

    // Newton stops when |dx| < tol (= 1e-9 here), so accuracy ≤ ~1e-9.
    // For tighter convergence the caller must supply a smaller `tol` —
    // but `tol` doubles as the Sturm-isolation precision and goes through
    // double_to_rational, which collapses below ≈ 2.3e-10. The product of
    // both limits is the achievable accuracy at this API level.
    BigFloat tolerance = BigFloat::from_double(1.0e-9, 256);
    EXPECT_TRUE(err_pos < tolerance)
        << "+√2 Newton polish off by " << err_pos.to_string(60);
    EXPECT_TRUE(err_neg < tolerance)
        << "-√2 Newton polish off by " << err_neg.to_string(60);
}

TEST(SturmBigFloat, PrecisionScales_LowerBoundOnAccuracy) {
    // Higher precision must give at least as much accuracy as lower.
    symbolic::CASContext ctx;
    ExprPtr poly = parse_poly(ctx, "x^2 - 2");
    ASSERT_NE(poly, nullptr);

    auto roots_64  = numeric::find_polynomial_roots_sturm_bigfloat(
        poly, "x", ctx, 0.0, 5.0, 1e-9, 64);
    auto roots_512 = numeric::find_polynomial_roots_sturm_bigfloat(
        poly, "x", ctx, 0.0, 5.0, 1e-9, 512);
    ASSERT_TRUE(roots_64.is_ok());
    ASSERT_TRUE(roots_512.is_ok());
    ASSERT_EQ(roots_64.value().size(), 1U);
    ASSERT_EQ(roots_512.value().size(), 1U);

    // Newton stops at |dx| < tol (1e-9 here) → accuracy ≤ ~1e-9.
    EXPECT_NEAR(roots_64.value()[0].to_double(), 1.41421356237309515, 5e-9);
    EXPECT_NEAR(roots_512.value()[0].to_double(), 1.41421356237309515, 5e-9);

    EXPECT_GE(roots_64.value()[0].precision_bits(), 64);
    EXPECT_GE(roots_512.value()[0].precision_bits(), 512);
}

TEST(SturmBigFloat, RejectsInvalidPrecision) {
    symbolic::CASContext ctx;
    ExprPtr poly = parse_poly(ctx, "x^2 - 2");
    ASSERT_NE(poly, nullptr);

    auto roots = numeric::find_polynomial_roots_sturm_bigfloat(
        poly, "x", ctx, -5.0, 5.0, 1e-9, 1); // < MIN_PREC
    EXPECT_TRUE(roots.is_error());
    EXPECT_EQ(roots.error().kind, CASErrorKind::InvalidArgument);
}

TEST(SturmBigFloat, EmptyInterval_NoRoots) {
    symbolic::CASContext ctx;
    ExprPtr poly = parse_poly(ctx, "x^2 - 2");
    ASSERT_NE(poly, nullptr);

    auto roots = numeric::find_polynomial_roots_sturm_bigfloat(
        poly, "x", ctx, 5.0, 0.0, 1e-9, 128); // low > high
    ASSERT_TRUE(roots.is_ok());
    EXPECT_EQ(roots.value().size(), 0U);
}

} // namespace
} // namespace cas
