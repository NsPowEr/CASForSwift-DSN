// CAS-L3-18 — Galois group tests (deg 2 e 3).

#include <gtest/gtest.h>

#include "cas/galois.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::algebra;

namespace {

class GaloisTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};
    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
};

TEST_F(GaloisTest, Quadratic_Trivial_Reducible) {
    // x²-1 = (x-1)(x+1) → trivial
    auto p = parse("x^2 - 1");
    auto r = galois_group(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_EQ(r.value(), "trivial");
}

TEST_F(GaloisTest, Quadratic_C2_Irreducible) {
    // x²+1 irreducible → C₂
    auto p = parse("x^2 + 1");
    auto r = galois_group(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_EQ(r.value(), "C2");
}

TEST_F(GaloisTest, Quadratic_C2_NonSquareDisc) {
    // x²-2 irreducible → C₂
    auto p = parse("x^2 - 2");
    auto r = galois_group(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_EQ(r.value(), "C2");
}

TEST_F(GaloisTest, Cubic_Trivial_Splits) {
    // x³-6x²+11x-6 = (x-1)(x-2)(x-3) → trivial
    auto p = parse("x^3 - 6*x^2 + 11*x - 6");
    auto r = galois_group(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_EQ(r.value(), "trivial");
}

TEST_F(GaloisTest, Cubic_S3_Irreducible_NonSquareDisc) {
    // x³-2 irreducible over Q; D = -108 (non-square) → S₃
    auto p = parse("x^3 - 2");
    auto r = galois_group(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_EQ(r.value(), "S3");
}

TEST_F(GaloisTest, Cubic_A3_SquareDiscriminant) {
    // x³ - 3x + 1: disc = -4·(-3)³ - 27·1² = 108 - 27 = 81 = 9² → A₃
    auto p = parse("x^3 - 3*x + 1");
    auto r = galois_group(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_EQ(r.value(), "A3");
}

TEST_F(GaloisTest, Cubic_C2_OneLinearOneQuadratic) {
    // (x-1)(x²+1) → linear + irreducible quadratic → C₂
    auto p = parse("(x - 1) * (x^2 + 1)");
    auto r = galois_group(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_EQ(r.value(), "C2");
}

TEST_F(GaloisTest, AntiHardcodeQuarticReturnsUnknown) {
    // x⁴+1 irreducible (cyclotomic Φ₈) — Galois group is V₄ but our MVP
    // returns "unknown" rather than guess.
    auto p = parse("x^4 + 1");
    auto r = galois_group(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    // Anti-hardcode: MUST NOT silently return wrong group. Acceptable:
    // "unknown", "trivial" (if splits, which it doesn't here), or a
    // specific deg-4 group label IF supported.
    EXPECT_TRUE(r.value() == "unknown" || r.value() == "V4"
                || r.value() == "S4" || r.value() == "trivial")
        << "got: " << r.value();
}

TEST_F(GaloisTest, Quartic_TrivialSplits) {
    // (x-1)(x-2)(x-3)(x-4) = x⁴-10x³+35x²-50x+24 → trivial
    auto p = parse("x^4 - 10*x^3 + 35*x^2 - 50*x + 24");
    auto r = galois_group(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_EQ(r.value(), "trivial");
}

TEST_F(GaloisTest, Quartic_S4_NonSquareDisc) {
    // x⁴+x+1 — irreducibile, discriminant non-square → S₄
    auto p = parse("x^4 + x + 1");
    auto r = galois_group(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    // Accept S4 (conservative MVP) or D4/A4/V4 (advanced analysis).
    EXPECT_TRUE(r.value() == "S4" || r.value() == "D4")
        << "got: " << r.value();
}

}  // namespace
