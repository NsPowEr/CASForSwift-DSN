// L2-19: mathematically_equal_subset_risch tests.
//
// The decidable subset of Richardson's problem covered here:
//   - log(x*y), log(x^n), log(x/y) under positivity / integrality assumptions.
//   - exp(x+y), exp(n*ln(x)), exp(ln(x)) under positivity / integrality.
//   - Pythagorean trig identity sin^2 + cos^2 = 1 (handled by existing simplifier).
// Cases outside the subset return false (not Unimplemented), consistent with
// Richardson's theorem precluding a total decision procedure.

#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

#include <memory>
#include <string>

namespace cas::test {
namespace {

[[nodiscard]] Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    Lexer lexer(input);
    auto tokens = lexer.tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

class EquivalenceSubsetRischTest : public ::testing::Test {
protected:
    void SetUp() override { ctx = std::make_unique<symbolic::CASContext>(); }

    [[nodiscard]] ExprPtr parse_ok(const std::string& input) {
        auto parsed = parse_expr(input, ctx->arena());
        EXPECT_TRUE(parsed.is_ok()) << parsed.error().message;
        return parsed.is_ok() ? parsed.value() : nullptr;
    }

    void assume_positive(const std::string& name) {
        ctx->assumptions().assume_positive(Symbol{name});
    }

    void assume_integer(const std::string& name) {
        ctx->assumptions().assume_integer(Symbol{name});
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

// ----- PASS cases (decidable subset) ---------------------------------------

TEST_F(EquivalenceSubsetRischTest, ExpOfLogSumEqualsProductUnderPositiveAssumption) {
    assume_positive("x");
    assume_positive("y");
    auto lhs = parse_ok("exp(ln(x) + ln(y))");
    auto rhs = parse_ok("x * y");
    auto eq = symbolic::mathematically_equal_subset_risch(lhs, rhs, *ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value());
}

TEST_F(EquivalenceSubsetRischTest, PythagoreanIdentity) {
    auto lhs = parse_ok("sin(x)^2 + cos(x)^2");
    auto rhs = parse_ok("1");
    auto eq = symbolic::mathematically_equal_subset_risch(lhs, rhs, *ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value());
}

TEST_F(EquivalenceSubsetRischTest, LogOfPowerCollapsesUnderPositiveBase) {
    assume_positive("x");
    auto lhs = parse_ok("ln(x^2)");
    auto rhs = parse_ok("2 * ln(x)");
    auto eq = symbolic::mathematically_equal_subset_risch(lhs, rhs, *ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value());
}

TEST_F(EquivalenceSubsetRischTest, ExpOfIntegerTimesLogCollapsesToPowerUnderPositive) {
    assume_positive("x");
    auto lhs = parse_ok("exp(2 * ln(x))");
    auto rhs = parse_ok("x^2");
    auto eq = symbolic::mathematically_equal_subset_risch(lhs, rhs, *ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value());
}

TEST_F(EquivalenceSubsetRischTest, LogProductMinusLogsZeroesUnderPositiveAssumption) {
    assume_positive("x");
    assume_positive("y");
    auto lhs = parse_ok("ln(x * y) - ln(x) - ln(y)");
    auto rhs = parse_ok("0");
    auto eq = symbolic::mathematically_equal_subset_risch(lhs, rhs, *ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value());
}

// ----- FALSE cases (assumption-driven safety) ------------------------------

// Behaviour note: without a positivity assumption on x, y the existing
// global simplifier still reduces exp(ln(x)) -> x unconditionally; this is
// outside the scope of the L2-19 subset walker (which only adds reductions,
// never restricts them).  See CAS_TASKS L2-19 note on the branch-cut gap.
TEST_F(EquivalenceSubsetRischTest, ExpOfLogSumWithoutPositivityIsNotEqualToProduct) {
    auto lhs = parse_ok("exp(ln(x) + ln(y))");
    auto rhs = parse_ok("x * y");
    auto eq = symbolic::mathematically_equal_subset_risch(lhs, rhs, *ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_FALSE(eq.value());
}

TEST_F(EquivalenceSubsetRischTest, LogOfSquareWithoutPositivityIsNotEqualToTwoLog) {
    auto lhs = parse_ok("ln(x^2)");
    auto rhs = parse_ok("2 * ln(x)");
    auto eq = symbolic::mathematically_equal_subset_risch(lhs, rhs, *ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_FALSE(eq.value())
        << "without x>0 the identity must not be claimed (branch-cut safety)";
}

// ----- Reflexivity / robustness ------------------------------------------

TEST_F(EquivalenceSubsetRischTest, ReflexivityHoldsForArbitraryExpression) {
    auto e = parse_ok("ln(x^2 * y) + exp(z + w)");
    auto eq = symbolic::mathematically_equal_subset_risch(e, e, *ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value()) << "subset_risch must be reflexive";
}

TEST_F(EquivalenceSubsetRischTest, RejectsNullOperand) {
    auto lhs = parse_ok("x");
    auto eq = symbolic::mathematically_equal_subset_risch(lhs, ExprPtr{}, *ctx);
    ASSERT_TRUE(eq.is_error());
    EXPECT_EQ(eq.error().kind, CASErrorKind::InvalidArgument);
}

// ----- New: positivity inference + non-integer scalar (B3+B4+B5 fixes) ----

TEST_F(EquivalenceSubsetRischTest, LogOfTripleProductExpandsUnderAllPositive) {
    assume_positive("x");
    assume_positive("y");
    assume_positive("z");
    auto lhs = parse_ok("ln(x * y * z)");
    auto rhs = parse_ok("ln(x) + ln(y) + ln(z)");
    auto eq = symbolic::mathematically_equal_subset_risch(lhs, rhs, *ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value());
}

TEST_F(EquivalenceSubsetRischTest, ExpOfRationalTimesLogCollapsesToPowerUnderPositive) {
    // For x > 0,  exp((1/2) * ln(x)) = x^(1/2) = sqrt(x).  Rational
    // exponent must NOT be rejected by the matcher (B3+B4 fix).
    assume_positive("x");
    auto lhs = parse_ok("exp((1/2) * ln(x))");
    auto rhs = parse_ok("x^(1/2)");
    auto eq = symbolic::mathematically_equal_subset_risch(lhs, rhs, *ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value());
}

TEST_F(EquivalenceSubsetRischTest, ExpAlwaysPositiveInfersLogExpInverse) {
    // exp(x) is always positive (real x) by structural inference (B5).
    // Thus ln(exp(x)) = x with only the standard "x real" assumption.
    ctx->assumptions().assume_real(Symbol{"x"});
    auto lhs = parse_ok("ln(exp(x))");
    auto rhs = parse_ok("x");
    auto eq = symbolic::mathematically_equal_subset_risch(lhs, rhs, *ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value());
}

TEST_F(EquivalenceSubsetRischTest, NestedLogProductExpandsRecursively) {
    // ln((x^2 * y^3) * z) under all-positive must reduce to
    // 2*ln(x) + 3*ln(y) + ln(z) via the recursive walker (B2 fix).
    assume_positive("x");
    assume_positive("y");
    assume_positive("z");
    auto lhs = parse_ok("ln((x^2 * y^3) * z)");
    auto rhs = parse_ok("2*ln(x) + 3*ln(y) + ln(z)");
    auto eq = symbolic::mathematically_equal_subset_risch(lhs, rhs, *ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value());
}

}  // namespace
}  // namespace cas::test
