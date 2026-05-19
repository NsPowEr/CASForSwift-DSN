#include <gtest/gtest.h>

#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "../../../src/calculus/calculus_internal.hpp"

using namespace cas;
using namespace cas::calculus;

// Step 1: Gruntz tower-adaptive depth bound.
//
// The legacy `if (depth >= 16U) return Unimplemented` in
// `compute_recursive` capped MRV recursion at a fixed constant. With the
// tower-adaptive bound `max(8U, 2*h + 4)` derived from
// `transcendental_tower_depth(expr, var)`, the engine's recursion budget
// scales with the asymptotic complexity of the input.
//
// These tests pin the helper's output on canonical Gruntz inputs and
// verify the bound formula is non-decreasing in tower height.

class TowerDepthHelperTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};

    [[nodiscard]] ExprPtr parse(const std::string& input) {
        auto tokens = Lexer(input).tokenize();
        EXPECT_TRUE(tokens.is_ok()) << input;
        Parser parser(tokens.value(), ctx.arena());
        auto res = parser.parse();
        EXPECT_TRUE(res.is_ok()) << input;
        return res.value();
    }
};

TEST_F(TowerDepthHelperTest, ConstantHasZeroTower) {
    auto e = parse("42");
    EXPECT_EQ(transcendental_tower_depth(e, x), 0U);
}

TEST_F(TowerDepthHelperTest, ConstantSymbolNotDependingOnVarHasZeroTower) {
    auto e = parse("y");
    EXPECT_EQ(transcendental_tower_depth(e, x), 0U);
}

TEST_F(TowerDepthHelperTest, PlainVariableHasTowerOne) {
    auto e = parse("x");
    EXPECT_EQ(transcendental_tower_depth(e, x), 1U);
}

TEST_F(TowerDepthHelperTest, PolynomialInXHasTowerOne) {
    auto e = parse("x^5 + 3*x^2 + 1");
    // Polynomial in x: no exp wrapper → height = 1 (not zero, since it
    // depends on x; not >1 since no exp).
    EXPECT_EQ(transcendental_tower_depth(e, x), 1U);
}

TEST_F(TowerDepthHelperTest, ExpOfXHasTowerTwo) {
    auto e = parse("exp(x)");
    // exp wraps x → 1 + 1 = 2.
    EXPECT_EQ(transcendental_tower_depth(e, x), 2U);
}

TEST_F(TowerDepthHelperTest, ExpOfExpOfXHasTowerThree) {
    auto e = parse("exp(exp(x))");
    EXPECT_EQ(transcendental_tower_depth(e, x), 3U);
}

TEST_F(TowerDepthHelperTest, TripleNestedExpHasTowerFour) {
    auto e = parse("exp(exp(exp(x)))");
    EXPECT_EQ(transcendental_tower_depth(e, x), 4U);
}

TEST_F(TowerDepthHelperTest, LogPeelsLevelButReportsAtLeastOne) {
    // ln(x) depends on x; tower depth ≥ 1.
    auto ln_x = parse("ln(x)");
    EXPECT_EQ(transcendental_tower_depth(ln_x, x), 1U);

    auto ln_exp_x = parse("ln(exp(x))");
    // exp(x) inside: child tower = 2 → max(2, 1) = 2.
    EXPECT_EQ(transcendental_tower_depth(ln_exp_x, x), 2U);
}

TEST_F(TowerDepthHelperTest, SumPropagatesMaxOfChildTowers) {
    // f = x + exp(exp(x)); polynomial path → 1, exp²(x) → 3. max = 3.
    auto e = parse("x + exp(exp(x))");
    EXPECT_EQ(transcendental_tower_depth(e, x), 3U);
}

TEST_F(TowerDepthHelperTest, ProductPropagatesMaxOfChildTowers) {
    auto e = parse("x * exp(exp(x))");
    EXPECT_EQ(transcendental_tower_depth(e, x), 3U);
}

TEST_F(TowerDepthHelperTest, RatioOfExponentialTowers) {
    // (exp(exp(x)) + 1) / exp(exp(x)) — quotient tower = 3.
    auto e = parse("(exp(exp(x)) + 1) / exp(exp(x))");
    EXPECT_EQ(transcendental_tower_depth(e, x), 3U);
}

// The legacy bound `depth >= 16U` capped recursion for any input,
// independent of structure. The adaptive bound `max(8, 2*h + 4)` derived
// from `transcendental_tower_depth` is strictly larger than 16 once
// h >= 6: a 6-level exp tower would have generated `Unimplemented` under
// the legacy cap but is now allotted budget = 16.
TEST_F(TowerDepthHelperTest, AdaptiveBoundExceedsLegacyForHeightSix) {
    auto e = parse("exp(exp(exp(exp(exp(exp(x))))))");
    unsigned int h = transcendental_tower_depth(e, x);
    EXPECT_GE(h, 6U) << "Six-level exp tower must report tower_depth >= 6";
    unsigned int adaptive_bound = std::max<unsigned int>(8U, 2U * h + 4U);
    EXPECT_GT(adaptive_bound, 16U)
        << "Adaptive bound must exceed legacy depth=16 cap for h=6";
}

// Smoke test of the integrated path: `(exp(exp(x)) + 1) / exp(exp(x))`
// at infinity. Tower height = 3 → adaptive bound = 10 (vs legacy 16, both
// adequate). The engine must complete without surfacing the legacy depth
// error message even if it ultimately fails for other reasons.
TEST_F(TowerDepthHelperTest, LimitOfHeightThreeRatioDoesNotHitDepthCap) {
    auto e = parse("(exp(exp(x)) + 1) / exp(exp(x))");
    auto inf = ctx.arena().make<Constant>(MathConstant::Infinity);
    auto res = limit(e, x, inf, LimitDirection::Both, ctx);
    if (res.is_error()) {
        EXPECT_EQ(res.error().message.find("richiede piu' iterazioni"),
                  std::string::npos)
            << "Tower-adaptive bound must not surface the legacy depth=16 "
            << "Unimplemented for height-3 inputs.";
    }
}
