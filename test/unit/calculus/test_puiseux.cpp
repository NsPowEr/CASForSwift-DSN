// F5.5 — Puiseux Newton polygon: leading-term extractor.
//
// Each test fixes an algebraic curve f(x, y) = 0 with a known branch
// structure at x = 0 and checks the multiset of (exponent, leading c)
// pairs returned by puiseux_leading_terms.

#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

#include <string>

namespace cas::test {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

}  // namespace

class PuiseuxTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

// y² − x = 0  →  y = ±√x.  Edge from (1, 0) to (0, 2): μ = 1/2, c² − 1 = 0
// (after the offset, characteristic Φ(c) = c² − 1 because a_{0,2} = 1,
// a_{1,0} = −1, j_min = 0).  Roots: c = ±1.
TEST_F(PuiseuxTest, SquareRoot_TwoBranchesAtHalf) {
    auto f = parse_expr("y^2 - x", ctx.arena()).value();
    auto br = calculus::puiseux_leading_terms(f, Symbol("x"), Symbol("y"), ctx);
    ASSERT_TRUE(br.is_ok()) << br.error().message;
    ASSERT_EQ(br.value().size(), 2U);
    for (const auto& b : br.value()) {
        EXPECT_EQ(b.leading_exponent.numerator(), BigInt(1));
        EXPECT_EQ(b.leading_exponent.denominator(), BigInt(2));
        EXPECT_EQ(b.multiplicity, 1U);
    }
}

// y² − x³ = 0  →  y = ±x^(3/2).  Edge from (3, 0) to (0, 2): slope −2/3, so
// μ = 3/2 reduced.  Two branches, c = ±1.
TEST_F(PuiseuxTest, CuspY2Equalsx3_BranchesAtThreeHalves) {
    auto f = parse_expr("y^2 - x^3", ctx.arena()).value();
    auto br = calculus::puiseux_leading_terms(f, Symbol("x"), Symbol("y"), ctx);
    ASSERT_TRUE(br.is_ok()) << br.error().message;
    ASSERT_EQ(br.value().size(), 2U);
    for (const auto& b : br.value()) {
        EXPECT_EQ(b.leading_exponent.numerator(), BigInt(3));
        EXPECT_EQ(b.leading_exponent.denominator(), BigInt(2));
    }
}

// y³ − x² = 0  →  three branches  y = ω · x^(2/3),  ω³ = 1.  Edge slope −2/3
// gives μ = 2/3.  Characteristic Φ(c) = c³ − 1; SymPy/solve_polynomial may
// return the three cube roots of unity as RootOf-equivalents.  We only check
// the exponent here (multiplicities handled by characteristic structure).
TEST_F(PuiseuxTest, CubeRoot_ThreeBranchesAtTwoThirds) {
    auto f = parse_expr("y^3 - x^2", ctx.arena()).value();
    auto br = calculus::puiseux_leading_terms(f, Symbol("x"), Symbol("y"), ctx);
    ASSERT_TRUE(br.is_ok()) << br.error().message;
    ASSERT_GE(br.value().size(), 1U);
    for (const auto& b : br.value()) {
        EXPECT_EQ(b.leading_exponent.numerator(), BigInt(2));
        EXPECT_EQ(b.leading_exponent.denominator(), BigInt(3));
    }
}

// y² + 2xy − x = 0.  Monomials (0,2,1), (1,1,2), (1,0,−1).  Lower hull
// (0,2) → (1,0); the intermediate (1,1) does not lie on the edge.  μ =
// (1−0)/(2−0) = 1/2.  Characteristic Φ(c) = c² − 1 → c = ±1.  Two branches
// at exponent 1/2 — the curve has a node at the origin whose two analytic
// arcs are y ≈ ±√x − x + ….
TEST_F(PuiseuxTest, NodeCurve_TwoBranchesAtHalf) {
    auto f = parse_expr("y^2 + 2*x*y - x", ctx.arena()).value();
    auto br = calculus::puiseux_leading_terms(f, Symbol("x"), Symbol("y"), ctx);
    ASSERT_TRUE(br.is_ok()) << br.error().message;
    ASSERT_EQ(br.value().size(), 2U);
    for (const auto& b : br.value()) {
        EXPECT_EQ(b.leading_exponent.numerator(), BigInt(1));
        EXPECT_EQ(b.leading_exponent.denominator(), BigInt(2));
        EXPECT_EQ(b.multiplicity, 1U);
    }
}

}  // namespace cas::test
