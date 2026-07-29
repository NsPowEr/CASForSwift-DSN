// HC-F8-PENDING-20 / Branch_Cut_Propagation.md §3.2 — direction-limit table
// at branch-cut edges.
//
// Contract pinned here: for a one-sided limit whose inner argument lands on
// the negative real axis (the sqrt/ln cut) with a var-dependent imaginary
// part, the limit takes the CUT-EDGE value of the approached side. Before
// this table the direction was silently ignored: substitution always gave
// the principal (top-edge) value, wrong for bottom-edge approaches.

#include <gtest/gtest.h>

#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class LimitBranchCutTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol t{"t"};

    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto tok = Lexer(s).tokenize();
        EXPECT_TRUE(tok.is_ok()) << s;
        Parser p(tok.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }

    [[nodiscard]] Result<ExprPtr> lim(const std::string& s, LimitDirection dir) {
        return calculus::limit(parse(s), t, parse("0"), dir, ctx);
    }

    void expect_equals(Result<ExprPtr> got, const std::string& expected_src) {
        ASSERT_TRUE(got.is_ok()) << got.error().message;
        auto expected = ctx.simplify(parse(expected_src));
        ASSERT_TRUE(expected.is_ok());
        auto eq = symbolic::mathematically_equal(got.value(), expected.value(), ctx);
        ASSERT_TRUE(eq.is_ok()) << eq.error().message;
        EXPECT_TRUE(eq.value())
            << "got " << debug_print(got.value())
            << " != expected " << expected_src;
    }
};

// ── Spec §3.2 canonical example: sqrt(x ± iε), ε → 0±, x < 0 ────────────────

TEST_F(LimitBranchCutTest, SqrtTopEdge_FromAbove) {
    // Im(w) = t → 0⁺ approaches from Im > 0: principal edge, +2i.
    expect_equals(lim("sqrt(-4 + i*t)", LimitDirection::Right), "2*i");
}

TEST_F(LimitBranchCutTest, SqrtBottomEdge_FromBelow) {
    // t → 0⁻ ⇒ Im(w) < 0: bottom edge, −2i. Plain substitution would give
    // +2i (principal) — the silent-wrong this table fixes.
    expect_equals(lim("sqrt(-4 + i*t)", LimitDirection::Left), "-2*i");
}

TEST_F(LimitBranchCutTest, SqrtFlippedImaginary_RightIsBottom) {
    // w = −4 − i·t: t → 0⁺ approaches from Im < 0 → bottom edge.
    expect_equals(lim("sqrt(-4 - i*t)", LimitDirection::Right), "-2*i");
}

TEST_F(LimitBranchCutTest, SqrtBothSides_DisagreeingEdges_Undefined) {
    auto r = lim("sqrt(-4 + i*t)", LimitDirection::Both);
    ASSERT_FALSE(r.is_ok());
    EXPECT_EQ(r.error().kind, CASErrorKind::Undefined);
}

TEST_F(LimitBranchCutTest, SqrtEvenValuation_BothSidesTopEdge) {
    // Im(w) = t² > 0 on both sides → both directions approach the top edge:
    // the two-sided limit exists and equals +2i.
    expect_equals(lim("sqrt(-4 + i*t^2)", LimitDirection::Both), "2*i");
    expect_equals(lim("sqrt(-4 + i*t^2)", LimitDirection::Left), "2*i");
}

// ── ln on the cut ───────────────────────────────────────────────────────────

// The simplifier decomposes ln(z) = ln|z| + i·arg(z) before the limit engine
// runs, so the ln cut is exercised through arg and the Sum decomposition.
// Expected values are written as 1/2·ln(16) (≡ ln 4): the equality helpers
// cannot yet fold ln(4) − ln(16)/2 → 0 (log-power extraction gap, noted in
// HC-F8-PENDING-20).

TEST_F(LimitBranchCutTest, LnTopEdge_FromAbove) {
    expect_equals(lim("ln(-4 + i*t)", LimitDirection::Right),
                  "1/2*ln(16) + i*pi");
}

TEST_F(LimitBranchCutTest, LnBottomEdge_FromBelow) {
    expect_equals(lim("ln(-4 + i*t)", LimitDirection::Left),
                  "1/2*ln(16) - i*pi");
}

// ── Composite: sound decomposition or explicit diagnostic ──────────────────

TEST_F(LimitBranchCutTest, CompositeSum_DecomposesSoundly) {
    // sqrt(−4+i·t) + 1: substituting blindly would return 1+2i for BOTH
    // directions; the Sum decomposition must pick the bottom edge on Left.
    expect_equals(lim("sqrt(-4 + i*t) + 1", LimitDirection::Left), "1 - 2*i");
}

TEST_F(LimitBranchCutTest, IndeterminateThroughCut_ExplicitUnimplemented) {
    // (sqrt(−4+i·t) − 2i)/t is 0/0 THROUGH the cut: no sound decomposition
    // exists at this level — the engine must refuse, not answer top-edge.
    auto r = lim("(sqrt(-4 + i*t) - 2*i)/t", LimitDirection::Left);
    ASSERT_FALSE(r.is_ok());
    EXPECT_EQ(r.error().kind, CASErrorKind::Unimplemented);
}

// ── Regressions: cut-free paths unaffected ──────────────────────────────────

TEST_F(LimitBranchCutTest, RealSqrtAwayFromCut_Unchanged) {
    expect_equals(lim("sqrt(4 + t)", LimitDirection::Left), "2");
}

TEST_F(LimitBranchCutTest, RealPathOnCut_PrincipalBranch_Unchanged) {
    // No imaginary dependence on t: the path runs along the cut itself and
    // the principal value stands (top edge by convention).
    expect_equals(lim("sqrt(-4 + t)", LimitDirection::Right), "2*i");
}

}  // namespace
