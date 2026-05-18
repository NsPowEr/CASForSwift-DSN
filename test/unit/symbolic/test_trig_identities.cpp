#include <gtest/gtest.h>
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::symbolic;

namespace {
[[nodiscard]] Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

// Simplify, then compare text with expected string.
// Two Arenas: one for input parse, one for expected parse — both via the same CASContext.
[[nodiscard]] bool trig_simplifies_to(const char* input_expr, const char* expected_expr) {
    CASContext ctx;
    auto in_e = parse_expr(input_expr, ctx.arena());
    if (!in_e.is_ok()) return false;
    auto simp = ctx.simplify(in_e.value());
    if (!simp.is_ok()) return false;

    CASContext ctx_ex;
    auto ex_e = parse_expr(expected_expr, ctx_ex.arena());
    if (!ex_e.is_ok()) return false;
    auto simp_ex = ctx_ex.simplify(ex_e.value());
    if (!simp_ex.is_ok()) return false;

    // Compare text forms of both simplified results.
    auto t1 = to_round_trip_text(simp.value());
    auto t2 = to_round_trip_text(simp_ex.value());
    if (!t1.is_ok() || !t2.is_ok()) return false;
    return t1.value() == t2.value();
}
} // namespace

// ── L2-07: Addition formula (phase shifts) ───────────────────────────────────

TEST(TrigIdentitiesTest, SinXPlusHalfPiEqualsCosx) {
    // sin(x + pi/2) = cos(x)
    EXPECT_TRUE(trig_simplifies_to("sin(x + pi/2)", "cos(x)"));
}

TEST(TrigIdentitiesTest, CosXPlusHalfPiEqualsNegSinx) {
    // cos(x + pi/2) = -sin(x)
    EXPECT_TRUE(trig_simplifies_to("cos(x + pi/2)", "-sin(x)"));
}

TEST(TrigIdentitiesTest, SinXPlusPiEqualsNegSinx) {
    // sin(x + pi) = -sin(x)
    EXPECT_TRUE(trig_simplifies_to("sin(x + pi)", "-sin(x)"));
}

TEST(TrigIdentitiesTest, CosXPlusPiEqualsNegCosx) {
    // cos(x + pi) = -cos(x)
    EXPECT_TRUE(trig_simplifies_to("cos(x + pi)", "-cos(x)"));
}

TEST(TrigIdentitiesTest, SinXMinusHalfPiEqualsNegCosx) {
    // sin(x - pi/2) = -cos(x)
    EXPECT_TRUE(trig_simplifies_to("sin(x - pi/2)", "-cos(x)"));
}

TEST(TrigIdentitiesTest, CosXMinusHalfPiEqualsSinx) {
    // cos(x - pi/2) = sin(x)
    EXPECT_TRUE(trig_simplifies_to("cos(x - pi/2)", "sin(x)"));
}

TEST(TrigIdentitiesTest, SinXPlusThirdPiExpansion) {
    // sin(x + pi/3) = (1/2)*sin(x) + (sqrt(3)/2)*cos(x)
    // Verify via round-trip: simplify both sides and compare.
    EXPECT_TRUE(trig_simplifies_to("sin(x + pi/3)",
        "1/2 * sin(x) + sqrt(3)/2 * cos(x)"));
}

TEST(TrigIdentitiesTest, CosXPlusSixthPiExpansion) {
    // cos(x + pi/6) = (sqrt(3)/2)*cos(x) - (1/2)*sin(x)
    EXPECT_TRUE(trig_simplifies_to("cos(x + pi/6)",
        "sqrt(3)/2 * cos(x) - 1/2 * sin(x)"));
}

// ── L2-07: Double-angle compaction ───────────────────────────────────────────

TEST(TrigIdentitiesTest, TwoSinCosEqualsDoubleAngle) {
    // 2*sin(x)*cos(x) = sin(2*x)
    EXPECT_TRUE(trig_simplifies_to("2*sin(x)*cos(x)", "sin(2*x)"));
}

TEST(TrigIdentitiesTest, SinCosProductIsHalfDoubleAngle) {
    // sin(x)*cos(x) = sin(2*x)/2
    EXPECT_TRUE(trig_simplifies_to("sin(x)*cos(x)", "sin(2*x)/2"));
}

// ── L2-07: Anti-hardcode — Gauss constructible angle ─────────────────────────

TEST(TrigIdentitiesTest, SinXPlusFifthPiExpansion) {
    // sin(x + pi/5) = sin(pi/5)*cos(x) + cos(pi/5)*sin(x)
    // Uses Gauss-constructible sin(pi/5) and cos(pi/5).
    // Verify consistency: simplify(sin(x+pi/5)) == simplify(sin(pi/5)*cos(x) + cos(pi/5)*sin(x))
    EXPECT_TRUE(trig_simplifies_to("sin(x + pi/5)",
        "sin(pi/5)*cos(x) + cos(pi/5)*sin(x)"));
}

// ── Power reduction (already in simplify_power) — regression guard ───────────

TEST(TrigIdentitiesTest, SinSquaredPowerReduction) {
    // sin^2(x) = 1/2 - 1/2*cos(2*x)
    EXPECT_TRUE(trig_simplifies_to("sin(x)^2",
        "1/2 - 1/2 * cos(2*x)"));
}

TEST(TrigIdentitiesTest, CosSquaredPowerReduction) {
    // cos^2(x) = 1/2 + 1/2*cos(2*x)
    EXPECT_TRUE(trig_simplifies_to("cos(x)^2",
        "1/2 + 1/2 * cos(2*x)"));
}
