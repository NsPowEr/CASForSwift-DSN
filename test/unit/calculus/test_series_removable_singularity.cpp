// test_series_removable_singularity.cpp — A37.
//
// taylor_series used to fail with "ComplexRational: division by zero" on every
// quotient whose denominator vanishes at the expansion point, because it built
// coefficients by substituting the point into f and its derivatives (a 0/0 form
// there). These functions are analytic at the point — the singularity is
// removable — so the expansion must exist.
//
// Validation is mathematical, never textual: each computed polynomial is
// compared to the closed-form expected series via mathematically_equal.

#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

#include <string>

namespace cas::calculus {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

// Expands `input` around 0 to `order` and asserts the polynomial equals
// `expected_series` mathematically.
void expect_series(const std::string& input,
                   unsigned int order,
                   const std::string& expected_series) {
    symbolic::CASContext ctx;
    auto f = parse_expr(input, ctx.arena());
    ASSERT_TRUE(f.is_ok()) << "parse failed: " << input;
    auto zero = ctx.arena().make<IntegerLit>(BigInt(0));

    auto series = taylor_series(f.value(), Symbol("x"), zero, order, ctx);
    ASSERT_TRUE(series.is_ok())
        << input << " -> " << series.error().message;

    auto expected = parse_expr(expected_series, ctx.arena());
    ASSERT_TRUE(expected.is_ok());

    auto equal = symbolic::mathematically_equal(series.value().polynomial,
                                                expected.value(), ctx);
    ASSERT_TRUE(equal.is_ok()) << "equality inconclusive for " << input;
    EXPECT_TRUE(equal.value()) << input << ": series does not match expected\n"
        << "  got:      " << debug_print(series.value().polynomial) << "\n"
        << "  expected: " << debug_print(expected.value());
}

TEST(SeriesRemovableSingularityA37, ExpMinusOneOverX) {
    // (e^x - 1)/x = 1 + x/2 + x^2/6 + x^3/24 + x^4/120
    expect_series("(exp(x) - 1)/x", 4,
                  "1 + x/2 + x^2/6 + x^3/24 + x^4/120");
}

TEST(SeriesRemovableSingularityA37, SinhOverX) {
    // sinh(x)/x = 1 + x^2/6 + x^4/120
    expect_series("sinh(x)/x", 5, "1 + x^2/6 + x^4/120");
}

TEST(SeriesRemovableSingularityA37, LogOnePlusXOverX) {
    // log(1+x)/x = 1 - x/2 + x^2/3 - x^3/4 + x^4/5
    expect_series("log(1 + x)/x", 4,
                  "1 - x/2 + x^2/3 - x^3/4 + x^4/5");
}

TEST(SeriesRemovableSingularityA37, ArctanOverX) {
    // atan(x)/x = 1 - x^2/3 + x^4/5
    expect_series("atan(x)/x", 5, "1 - x^2/3 + x^4/5");
}

TEST(SeriesRemovableSingularityA37, XOverExpMinusOne) {
    // x/(e^x - 1) = 1 - x/2 + x^2/12 - x^4/720  (Bernoulli generating function)
    expect_series("x/(exp(x) - 1)", 4, "1 - x/2 + x^2/12 - x^4/720");
}

// A genuine pole is not a removable singularity: no Taylor expansion exists.
// The engine must say so explicitly instead of returning a wrong series or an
// opaque arithmetic error.
TEST(SeriesRemovableSingularityA37, GenuinePoleReportsUnimplemented) {
    symbolic::CASContext ctx;
    auto f = parse_expr("sin(x)/x^2", ctx.arena());
    ASSERT_TRUE(f.is_ok());
    auto zero = ctx.arena().make<IntegerLit>(BigInt(0));

    auto series = taylor_series(f.value(), Symbol("x"), zero, 4U, ctx);
    ASSERT_FALSE(series.is_ok()) << "a function with a pole has no Taylor series";
    EXPECT_EQ(series.error().kind, CASErrorKind::Unimplemented);
}

// The Laurent machinery does represent that same function.
TEST(SeriesRemovableSingularityA37, GenuinePoleHasLaurentExpansion) {
    symbolic::CASContext ctx;
    auto f = parse_expr("sin(x)/x^2", ctx.arena());
    ASSERT_TRUE(f.is_ok());
    auto zero = ctx.arena().make<IntegerLit>(BigInt(0));

    auto lau = laurent_series(f.value(), Symbol("x"), zero, 3U, ctx);
    ASSERT_TRUE(lau.is_ok()) << lau.error().message;
    // The expansion must record a pole (leading_order < 0); the exact bound
    // is the machinery's own, and need not be tight.
    EXPECT_LT(lau.value().leading_order, 0);
}

}  // namespace
}  // namespace cas::calculus
