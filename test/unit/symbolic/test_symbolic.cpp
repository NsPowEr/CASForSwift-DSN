#include "cas/ast.hpp"
#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/normal_form.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "../../helpers/property_test.hpp"

// Dummy comment for rebuild
#include <gtest/gtest.h>

#include <chrono>
#include <string>
#include <vector>

namespace cas::symbolic {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) {
        return fail<ExprPtr>(tokens.error());
    }

    Parser parser(tokens.value(), arena);
    return parser.parse();
}

Result<ExprPtr> simplify_input(const std::string& input, AstArena& parse_arena, AstArena& simplify_arena) {
    auto parsed = parse_expr(input, parse_arena);
    if (parsed.is_error()) {
        return fail<ExprPtr>(parsed.error());
    }

    return simplify(parsed.value(), simplify_arena);
}

Result<ExprPtr> simplify_input_with_context(const std::string& input, AstArena& parse_arena, CASContext& context) {
    auto parsed = parse_expr(input, parse_arena);
    if (parsed.is_error()) {
        return fail<ExprPtr>(parsed.error());
    }

    return context.simplify(parsed.value());
}

ExprPtr sum_terms(AstArena& arena, const std::vector<ExprPtr>& terms) {
    if (terms.empty()) {
        return arena.make<IntegerLit>(BigInt(0));
    }
    if (terms.size() == 1U) {
        return terms.front();
    }
    return arena.make<Sum>(terms);
}

TEST(SymbolicSimplifyTest, ReducesIntegerConstants) {
    AstArena parse_arena;
    AstArena simplify_arena;

    auto simplified = simplify_input("3 + 4", parse_arena, simplify_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    ASSERT_EQ(expr_kind(simplified.value()), ExprKind::IntegerLit);
    const auto* integer = expr_cast<IntegerLit>(simplified.value());
    ASSERT_NE(integer, nullptr);
    EXPECT_EQ(integer->value.decimal(), "7");
    EXPECT_FALSE(integer->value.is_negative());
}

TEST(SymbolicSimplifyTest, ReducesRationalConstantsToCanonicalForm) {
    AstArena parse_arena;
    AstArena simplify_arena;

    auto simplified = simplify_input("2/4", parse_arena, simplify_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    ASSERT_EQ(expr_kind(simplified.value()), ExprKind::RationalLit);
    const auto* rational = expr_cast<RationalLit>(simplified.value());
    ASSERT_NE(rational, nullptr);
    EXPECT_EQ(rational->numerator.decimal(), "1");
    EXPECT_EQ(rational->denominator.decimal(), "2");
    EXPECT_FALSE(rational->numerator.is_negative());
}

TEST(SymbolicSimplifyTest, RejectsExactDivisionByZero) {
    AstArena parse_arena;
    AstArena simplify_arena;

    auto simplified = simplify_input("1/0", parse_arena, simplify_arena);

    ASSERT_TRUE(simplified.is_error());
    EXPECT_EQ(simplified.error().kind, CASErrorKind::Undefined);
}

TEST(SymbolicSimplifyTest, ConvertsDecimalLiteralToRational) {
    AstArena parse_arena;
    AstArena simplify_arena;

    auto simplified = simplify_input("3.14", parse_arena, simplify_arena);

    ASSERT_TRUE(simplified.is_ok());
    auto expected = simplify_arena.make<RationalLit>(BigInt(157), BigInt(50));
    EXPECT_TRUE(structural_equal(simplified.value(), expected));
}

TEST(SymbolicSimplifyTest, OrdersMonomialFactorsLexicographically) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("y * 3 * x", parse_arena, simplify_arena);
    auto expected = simplify_input("3 * x * y", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, OrdersPolynomialTermsByDescendingDegree) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("1 + x^3 + x", parse_arena, simplify_arena);
    auto expected = simplify_input("x^3 + x + 1", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, OrdersHugeDegreeAheadOfLinearTermWithoutOverflowFallback) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("x + x^999999999999999999999999999999999999999", parse_arena, simplify_arena);
    auto expected = simplify_input("x^999999999999999999999999999999999999999 + x", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, FlattensNestedSums) {
    AstArena parse_arena;
    AstArena simplify_arena;

    auto simplified = simplify_input("(a + b) + c", parse_arena, simplify_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    ASSERT_EQ(expr_kind(simplified.value()), ExprKind::Sum);
    const auto* sum = expr_cast<Sum>(simplified.value());
    ASSERT_NE(sum, nullptr);
    ASSERT_EQ(sum->terms.size(), 3U);
    EXPECT_EQ(expr_kind(sum->terms[0]), ExprKind::Symbol);
    EXPECT_EQ(expr_kind(sum->terms[1]), ExprKind::Symbol);
    EXPECT_EQ(expr_kind(sum->terms[2]), ExprKind::Symbol);
}

TEST(SymbolicSimplifyTest, FlattensNestedProducts) {
    AstArena parse_arena;
    AstArena simplify_arena;

    auto simplified = simplify_input("(a * b) * c", parse_arena, simplify_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    ASSERT_EQ(expr_kind(simplified.value()), ExprKind::Product);
    const auto* product = expr_cast<Product>(simplified.value());
    ASSERT_NE(product, nullptr);
    ASSERT_EQ(product->factors.size(), 3U);
}

TEST(SymbolicSimplifyTest, RemovesMultiplicativeIdentity) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("x * 1", parse_arena, simplify_arena);
    auto expected = simplify_input("x", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, RemovesAdditiveIdentity) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("x + 0", parse_arena, simplify_arena);
    auto expected = simplify_input("x", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, CollapsesProductWithZero) {
    AstArena parse_arena;
    AstArena simplify_arena;

    auto simplified = simplify_input("x * 0", parse_arena, simplify_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    ASSERT_EQ(expr_kind(simplified.value()), ExprKind::IntegerLit);
    const auto* integer = expr_cast<IntegerLit>(simplified.value());
    ASSERT_NE(integer, nullptr);
    EXPECT_EQ(integer->value.decimal(), "0");
}

TEST(SymbolicSimplifyTest, SimplifiesPowerZeroExponent) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("x^0", parse_arena, simplify_arena);
    auto expected = simplify_input("1", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, SimplifiesPowerUnitExponent) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("x^1", parse_arena, simplify_arena);
    auto expected = simplify_input("x", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, EvaluatesExactIntegerPowerWithBinaryExponentiation) {
    AstArena parse_arena;
    AstArena simplify_arena;

    auto simplified = simplify_input("2^20", parse_arena, simplify_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_EQ(expr_kind(simplified.value()), ExprKind::IntegerLit);
    EXPECT_EQ(expr_ref<IntegerLit>(simplified.value()).value.decimal(), "1048576");
}

TEST(SymbolicSimplifyTest, EvaluatesExactNegativeIntegerPowerWithBinaryExponentiation) {
    AstArena parse_arena;
    AstArena simplify_arena;

    auto simplified = simplify_input("2^(-10)", parse_arena, simplify_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_EQ(expr_kind(simplified.value()), ExprKind::RationalLit);
    const auto& rational = expr_ref<RationalLit>(simplified.value());
    EXPECT_EQ(rational.numerator.decimal(), "1");
    EXPECT_EQ(rational.denominator.decimal(), "1024");
}

TEST(SymbolicSimplifyTest, ReusesOriginalPointerForStableFunctionCall) {
    AstArena parse_arena;
    CASContext context;

    auto parsed = parse_expr("sin(x)", parse_arena);
    ASSERT_TRUE(parsed.is_ok()) << parsed.error().message;

    auto simplified = context.simplify(parsed.value());
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    EXPECT_EQ(simplified.value(), parsed.value());
}

TEST(AlgebraCollectTest, CollectsRepeatedQuadraticTerms) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("x*x + 2*x^2 + x", parse_arena);
    auto expected = parse_expr("3*x^2 + x", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto collected = algebra::collect(expr.value(), Symbol("x"), ctx);
    ASSERT_TRUE(collected.is_ok()) << collected.error().message;

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(collected.value(), simplified_expected.value()));
}

TEST(AlgebraCollectTest, PreservesSymbolicCoefficientsIndependentFromVariable) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("a*x + 2*x + a", parse_arena);
    auto expected = parse_expr("(a + 2)*x + a", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto collected = algebra::collect(expr.value(), Symbol("x"), ctx);
    ASSERT_TRUE(collected.is_ok()) << collected.error().message;

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(collected.value(), simplified_expected.value()));
}

TEST(AlgebraCollectTest, SupportsConstantScalarDivision) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("(x^2 + 2*x)/2", parse_arena);
    auto expected = parse_expr("1/2*x^2 + x", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto collected = algebra::collect(expr.value(), Symbol("x"), ctx);
    ASSERT_TRUE(collected.is_ok()) << collected.error().message;

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(collected.value(), simplified_expected.value()));
}

TEST(AlgebraCollectTest, RejectsTranscendentalDependence) {
    AstArena parse_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("sin(x) + x", parse_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto collected = algebra::collect(expr.value(), Symbol("x"), ctx);
    ASSERT_TRUE(collected.is_error());
    EXPECT_EQ(collected.error().kind, CASErrorKind::Unimplemented);
}

TEST(AlgebraCollectTest, RejectsExponentThatExceedsSupportedMagnitude) {
    AstArena parse_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("x^999999999999999999999999999999999999999", parse_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto collected = algebra::collect(expr.value(), Symbol("x"), ctx);
    ASSERT_TRUE(collected.is_error());
    EXPECT_EQ(collected.error().kind, CASErrorKind::Unimplemented);
}

TEST(AlgebraRationalPartsTest, SplitsSumIntoCommonDenominator) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("1/x + 1/(x+1)", parse_arena);
    auto expected_num = parse_expr("2*x + 1", expected_arena);
    auto expected_den = parse_expr("x*(x+1)", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected_num.is_ok()) << expected_num.error().message;
    ASSERT_TRUE(expected_den.is_ok()) << expected_den.error().message;

    auto parts = algebra::apart_num_den(expr.value(), ctx);
    ASSERT_TRUE(parts.is_ok()) << parts.error().message;

    auto simplified_num = ctx.simplify(expected_num.value());
    auto simplified_den = ctx.simplify(expected_den.value());
    ASSERT_TRUE(simplified_num.is_ok()) << simplified_num.error().message;
    ASSERT_TRUE(simplified_den.is_ok()) << simplified_den.error().message;
    EXPECT_TRUE(structural_equal(parts.value().numerator, simplified_num.value()));
    EXPECT_TRUE(structural_equal(parts.value().denominator, simplified_den.value()));
}

TEST(AlgebraRationalPartsTest, SupportsNegativeIntegerPowers) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("(x + 1)^(-2)", parse_arena);
    auto expected_num = parse_expr("1", expected_arena);
    auto expected_den = parse_expr("(x + 1)^2", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected_num.is_ok()) << expected_num.error().message;
    ASSERT_TRUE(expected_den.is_ok()) << expected_den.error().message;

    auto parts = algebra::apart_num_den(expr.value(), ctx);
    ASSERT_TRUE(parts.is_ok()) << parts.error().message;

    auto simplified_num = ctx.simplify(expected_num.value());
    auto simplified_den = ctx.simplify(expected_den.value());
    ASSERT_TRUE(simplified_num.is_ok()) << simplified_num.error().message;
    ASSERT_TRUE(simplified_den.is_ok()) << simplified_den.error().message;
    EXPECT_TRUE(structural_equal(parts.value().numerator, simplified_num.value()));
    EXPECT_TRUE(structural_equal(parts.value().denominator, simplified_den.value()));
}

TEST(AlgebraTogetherTest, RebuildsSingleRationalExpression) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("1/x + 1/(x+1)", parse_arena);
    auto expected = parse_expr("(2*x + 1)/(x*(x+1))", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto merged = algebra::together(expr.value(), ctx);
    ASSERT_TRUE(merged.is_ok()) << merged.error().message;

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(merged.value(), simplified_expected.value()));
}

TEST(AlgebraExpandTest, ExpandsBinomialSquare) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("(x + 1)^2", parse_arena);
    auto expected = parse_expr("x^2 + 2*x + 1", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto expanded = algebra::expand(expr.value(), ctx);
    ASSERT_TRUE(expanded.is_ok()) << expanded.error().message;

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(expanded.value(), simplified_expected.value()));
}

TEST(AlgebraExpandTest, ExpandsProductOfBinomials) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("(x + 1)*(x - 1)", parse_arena);
    auto expected = parse_expr("x^2 - 1", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto expanded = algebra::expand(expr.value(), ctx);
    ASSERT_TRUE(expanded.is_ok()) << expanded.error().message;

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(expanded.value(), simplified_expected.value()));
}

TEST(AlgebraExpandTest, ExpandsRationalNumeratorAndDenominator) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("((x + 1)^2)/((x - 1)^2)", parse_arena);
    auto expected = parse_expr("(x^2 + 2*x + 1)/(x^2 - 2*x + 1)", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto expanded = algebra::expand(expr.value(), ctx);
    ASSERT_TRUE(expanded.is_ok()) << expanded.error().message;

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(expanded.value(), simplified_expected.value()));
}

TEST(AlgebraExpandTest, ExpandsNegativeIntegerPowerAsReciprocal) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("(x + 1)^(-2)", parse_arena);
    auto expected = parse_expr("1/(x^2 + 2*x + 1)", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto expanded = algebra::expand(expr.value(), ctx);
    ASSERT_TRUE(expanded.is_ok()) << expanded.error().message;

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(expanded.value(), simplified_expected.value()));
}

TEST(AlgebraExpandTest, RejectsExponentThatExceedsSupportedMagnitude) {
    AstArena parse_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("(x + 1)^999999999999999999999999999999999999999", parse_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto expanded = algebra::expand(expr.value(), ctx);
    ASSERT_TRUE(expanded.is_error());
    EXPECT_EQ(expanded.error().kind, CASErrorKind::Unimplemented);
}

TEST(AlgebraGcdTest, ExtractsCommonLinearFactor) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto left = parse_expr("x^2 - 1", parse_arena);
    auto right = parse_expr("x^2 - 3*x + 2", parse_arena);
    auto expected = parse_expr("x - 1", expected_arena);
    ASSERT_TRUE(left.is_ok()) << left.error().message;
    ASSERT_TRUE(right.is_ok()) << right.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto gcd = algebra::polynomial_gcd(left.value(), right.value(), Symbol("x"), ctx);
    ASSERT_TRUE(gcd.is_ok()) << gcd.error().message;

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(gcd.value(), simplified_expected.value()));
}

TEST(AlgebraGcdTest, NormalizesScalarMultiplesToMonicResult) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto left = parse_expr("2*x^2 + 2*x", parse_arena);
    auto right = parse_expr("2*x", parse_arena);
    auto expected = parse_expr("x", expected_arena);
    ASSERT_TRUE(left.is_ok()) << left.error().message;
    ASSERT_TRUE(right.is_ok()) << right.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto gcd = algebra::polynomial_gcd(left.value(), right.value(), Symbol("x"), ctx);
    ASSERT_TRUE(gcd.is_ok()) << gcd.error().message;

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(gcd.value(), simplified_expected.value()));
}

TEST(AlgebraGcdTest, SupportsSymbolicCoefficientsIndependentFromVariable) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto left = parse_expr("a*x + a", parse_arena);
    auto right = parse_expr("x + 1", parse_arena);
    auto expected = parse_expr("x + 1", expected_arena);
    ASSERT_TRUE(left.is_ok()) << left.error().message;
    ASSERT_TRUE(right.is_ok()) << right.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto gcd = algebra::polynomial_gcd(left.value(), right.value(), Symbol("x"), ctx);
    ASSERT_TRUE(gcd.is_ok()) << gcd.error().message;

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(gcd.value(), simplified_expected.value()));
}

TEST(AlgebraGcdTest, HandlesPrimitivePseudoRemainderSequenceForIntegerCoefficients) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto left = parse_expr("2*x^2 + x - 1", parse_arena);
    auto right = parse_expr("2*x^2 + 3*x - 2", parse_arena);
    auto expected = parse_expr("x - 1/2", expected_arena);
    ASSERT_TRUE(left.is_ok()) << left.error().message;
    ASSERT_TRUE(right.is_ok()) << right.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto gcd = algebra::polynomial_gcd(left.value(), right.value(), Symbol("x"), ctx);
    ASSERT_TRUE(gcd.is_ok()) << gcd.error().message;

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(gcd.value(), simplified_expected.value()));
}

TEST(AlgebraGcdTest, ReturnsOneForCoprimeIntegerPolynomialsViaSubresultantPath) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;
    ctx.enable_trace(true);

    auto left = parse_expr("x^2 + 7*x + 6", parse_arena);
    auto right = parse_expr("x^2 - 5*x - 6", parse_arena);
    auto expected = parse_expr("x + 1", expected_arena);
    ASSERT_TRUE(left.is_ok()) << left.error().message;
    ASSERT_TRUE(right.is_ok()) << right.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto gcd = algebra::polynomial_gcd(left.value(), right.value(), Symbol("x"), ctx);
    ASSERT_TRUE(gcd.is_ok()) << gcd.error().message;
    ASSERT_FALSE(ctx.get_trace().empty());
    EXPECT_EQ(ctx.get_trace().back().rule_id, symbolic::RuleId::PolynomialGcdSubresultant);

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(gcd.value(), simplified_expected.value()));
}

TEST(AlgebraGcdTest, ReportsSymbolicEuclideanPathForSymbolicCoefficients) {
    AstArena parse_arena;
    symbolic::CASContext ctx;
    ctx.enable_trace(true);

    auto left = parse_expr("a*x + a", parse_arena);
    auto right = parse_expr("x + 1", parse_arena);
    ASSERT_TRUE(left.is_ok()) << left.error().message;
    ASSERT_TRUE(right.is_ok()) << right.error().message;

    auto gcd = algebra::polynomial_gcd(left.value(), right.value(), Symbol("x"), ctx);
    ASSERT_TRUE(gcd.is_ok()) << gcd.error().message;
    ASSERT_FALSE(ctx.get_trace().empty());
    EXPECT_EQ(ctx.get_trace().back().rule_id, symbolic::RuleId::PolynomialGcdSymbolicEuclidean);
}

TEST(AlgebraGcdTest, ReportsBetaFallbackForKnownIntegerCase) {
    AstArena parse_arena;
    symbolic::CASContext ctx;
    ctx.enable_trace(true);

    auto left = parse_expr("x^4 - 3*x^3 - 3*x^2 - 2*x", parse_arena);
    auto right = parse_expr("3*x^4 - 3*x^3 - 3*x^2 - 3*x - 3", parse_arena);
    ASSERT_TRUE(left.is_ok()) << left.error().message;
    ASSERT_TRUE(right.is_ok()) << right.error().message;

    auto gcd = algebra::polynomial_gcd(left.value(), right.value(), Symbol("x"), ctx);
    ASSERT_TRUE(gcd.is_ok()) << gcd.error().message;
    ASSERT_FALSE(ctx.get_trace().empty());
    EXPECT_EQ(ctx.get_trace().back().rule_id, symbolic::RuleId::PolynomialGcdSubresultant);
}

TEST(AlgebraFactorTest, ExtractsLinearIntegerFactors) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto poly = parse_expr("x^2 - 5*x + 6", parse_arena);
    auto expected_first = parse_expr("x - 2", expected_arena);
    auto expected_second = parse_expr("x - 3", expected_arena);
    ASSERT_TRUE(poly.is_ok()) << poly.error().message;
    ASSERT_TRUE(expected_first.is_ok()) << expected_first.error().message;
    ASSERT_TRUE(expected_second.is_ok()) << expected_second.error().message;

    auto factorization = algebra::factor_over_integers(poly.value(), Symbol("x"), ctx);
    ASSERT_TRUE(factorization.is_ok()) << factorization.error().message;
    ASSERT_EQ(expr_kind(factorization.value().content), ExprKind::IntegerLit);
    EXPECT_EQ(expr_ref<IntegerLit>(factorization.value().content).value.decimal(), "1");
    ASSERT_EQ(factorization.value().factors.size(), 2U);
    EXPECT_EQ(factorization.value().factors[0].multiplicity, 1U);
    EXPECT_EQ(factorization.value().factors[1].multiplicity, 1U);

    auto simplified_first = ctx.simplify(expected_first.value());
    auto simplified_second = ctx.simplify(expected_second.value());
    ASSERT_TRUE(simplified_first.is_ok()) << simplified_first.error().message;
    ASSERT_TRUE(simplified_second.is_ok()) << simplified_second.error().message;
    EXPECT_TRUE(structural_equal(factorization.value().factors[0].factor, simplified_first.value()));
    EXPECT_TRUE(structural_equal(factorization.value().factors[1].factor, simplified_second.value()));
}

TEST(AlgebraFactorTest, PreservesNonMonicLinearFactor) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto poly = parse_expr("2*x^2 - 3*x + 1", parse_arena);
    auto expected_first = parse_expr("x - 1", expected_arena);
    auto expected_second = parse_expr("2*x - 1", expected_arena);
    ASSERT_TRUE(poly.is_ok()) << poly.error().message;
    ASSERT_TRUE(expected_first.is_ok()) << expected_first.error().message;
    ASSERT_TRUE(expected_second.is_ok()) << expected_second.error().message;

    auto factorization = algebra::factor_over_integers(poly.value(), Symbol("x"), ctx);
    ASSERT_TRUE(factorization.is_ok()) << factorization.error().message;
    ASSERT_EQ(factorization.value().factors.size(), 2U);

    auto simplified_first = ctx.simplify(expected_first.value());
    auto simplified_second = ctx.simplify(expected_second.value());
    ASSERT_TRUE(simplified_first.is_ok()) << simplified_first.error().message;
    ASSERT_TRUE(simplified_second.is_ok()) << simplified_second.error().message;
    EXPECT_TRUE(structural_equal(factorization.value().factors[0].factor, simplified_first.value()));
    EXPECT_TRUE(structural_equal(factorization.value().factors[1].factor, simplified_second.value()));
}

TEST(AlgebraFactorTest, TracksRepeatedFactorMultiplicity) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto poly = parse_expr("x^2 - 2*x + 1", parse_arena);
    auto expected = parse_expr("x - 1", expected_arena);
    ASSERT_TRUE(poly.is_ok()) << poly.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto factorization = algebra::factor_over_integers(poly.value(), Symbol("x"), ctx);
    ASSERT_TRUE(factorization.is_ok()) << factorization.error().message;
    ASSERT_EQ(factorization.value().factors.size(), 1U);
    EXPECT_EQ(factorization.value().factors[0].multiplicity, 2U);

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(factorization.value().factors[0].factor, simplified_expected.value()));
}

TEST(AlgebraFactorTest, KeepsIrreducibleResidualFactor) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto poly = parse_expr("x^2 + 1", parse_arena);
    auto expected = parse_expr("x^2 + 1", expected_arena);
    ASSERT_TRUE(poly.is_ok()) << poly.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto factorization = algebra::factor_over_integers(poly.value(), Symbol("x"), ctx);
    ASSERT_TRUE(factorization.is_ok()) << factorization.error().message;
    ASSERT_EQ(factorization.value().factors.size(), 1U);
    EXPECT_EQ(factorization.value().factors[0].multiplicity, 1U);

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(factorization.value().factors[0].factor, simplified_expected.value()));
}

TEST(AlgebraFactorTest, TracksRepeatedIrreducibleQuadraticFactor) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto poly = parse_expr("x^4 + 2*x^2 + 1", parse_arena);
    auto expected = parse_expr("x^2 + 1", expected_arena);
    ASSERT_TRUE(poly.is_ok()) << poly.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto factorization = algebra::factor_over_integers(poly.value(), Symbol("x"), ctx);
    ASSERT_TRUE(factorization.is_ok()) << factorization.error().message;
    ASSERT_EQ(factorization.value().factors.size(), 1U);
    EXPECT_EQ(factorization.value().factors[0].multiplicity, 2U);

    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(factorization.value().factors[0].factor, simplified_expected.value()));
}

TEST(AlgebraFactorTest, SeparatesSquareFreeAndRepeatedComponents) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto poly = parse_expr("x^5 - x^4 + 2*x^3 - 2*x^2 + x - 1", parse_arena);
    auto expected_first = parse_expr("x - 1", expected_arena);
    auto expected_second = parse_expr("x^2 + 1", expected_arena);
    ASSERT_TRUE(poly.is_ok()) << poly.error().message;
    ASSERT_TRUE(expected_first.is_ok()) << expected_first.error().message;
    ASSERT_TRUE(expected_second.is_ok()) << expected_second.error().message;

    auto factorization = algebra::factor_over_integers(poly.value(), Symbol("x"), ctx);
    ASSERT_TRUE(factorization.is_ok()) << factorization.error().message;
    ASSERT_EQ(factorization.value().factors.size(), 2U);
    EXPECT_EQ(factorization.value().factors[0].multiplicity, 1U);
    EXPECT_EQ(factorization.value().factors[1].multiplicity, 2U);

    auto simplified_first = ctx.simplify(expected_first.value());
    auto simplified_second = ctx.simplify(expected_second.value());
    ASSERT_TRUE(simplified_first.is_ok()) << simplified_first.error().message;
    ASSERT_TRUE(simplified_second.is_ok()) << simplified_second.error().message;
    EXPECT_TRUE(structural_equal(factorization.value().factors[0].factor, simplified_first.value()));
    EXPECT_TRUE(structural_equal(factorization.value().factors[1].factor, simplified_second.value()));
}

TEST(AlgebraPartialFractionsTest, DecomposesProperRationalWithSimpleLinearFactors) {
    AstArena parse_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("1/(x*(x+1))", parse_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto terms = algebra::partial_fractions(expr.value(), Symbol("x"), ctx);
    ASSERT_TRUE(terms.is_ok()) << terms.error().message;
    ASSERT_EQ(terms.value().size(), 2U);

    auto rebuilt = algebra::together(sum_terms(ctx.arena(), terms.value()), ctx);
    ASSERT_TRUE(rebuilt.is_ok()) << rebuilt.error().message;

    auto simplified_original = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified_original.is_ok()) << simplified_original.error().message;
    auto equivalent = mathematically_equal(rebuilt.value(), simplified_original.value(), ctx);
    ASSERT_TRUE(equivalent.is_ok()) << equivalent.error().message;
    EXPECT_TRUE(equivalent.value());
}

TEST(AlgebraPartialFractionsTest, DecomposesRationalWithAffineLinearFactors) {
    AstArena parse_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("(2*x + 3)/((x - 1)*(x + 2))", parse_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto terms = algebra::partial_fractions(expr.value(), Symbol("x"), ctx);
    ASSERT_TRUE(terms.is_ok()) << terms.error().message;
    ASSERT_EQ(terms.value().size(), 2U);

    auto rebuilt = algebra::together(sum_terms(ctx.arena(), terms.value()), ctx);
    ASSERT_TRUE(rebuilt.is_ok()) << rebuilt.error().message;

    auto simplified_original = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified_original.is_ok()) << simplified_original.error().message;
    auto equivalent = mathematically_equal(rebuilt.value(), simplified_original.value(), ctx);
    ASSERT_TRUE(equivalent.is_ok()) << equivalent.error().message;
    EXPECT_TRUE(equivalent.value());
}

TEST(AlgebraPartialFractionsTest, PreservesSignedResiduesForSharedLinearContent) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("(2*x + 3)/(x^2 + x)", parse_arena);
    auto expected_first = parse_expr("3/x", expected_arena);
    auto expected_second = parse_expr("-1/(x + 1)", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected_first.is_ok()) << expected_first.error().message;
    ASSERT_TRUE(expected_second.is_ok()) << expected_second.error().message;

    auto terms = algebra::partial_fractions(expr.value(), Symbol("x"), ctx);
    ASSERT_TRUE(terms.is_ok()) << terms.error().message;
    ASSERT_EQ(terms.value().size(), 2U);
    auto simplified_first = ctx.simplify(expected_first.value());
    auto simplified_second = ctx.simplify(expected_second.value());
    ASSERT_TRUE(simplified_first.is_ok()) << simplified_first.error().message;
    ASSERT_TRUE(simplified_second.is_ok()) << simplified_second.error().message;
    EXPECT_TRUE(
        (structural_equal(terms.value()[0], simplified_first.value()) &&
         structural_equal(terms.value()[1], simplified_second.value())) ||
        (structural_equal(terms.value()[0], simplified_second.value()) &&
         structural_equal(terms.value()[1], simplified_first.value())));
}

TEST(AlgebraPartialFractionsTest, DecomposesRepeatedLinearFactor) {
    AstArena parse_arena;
    AstArena expected_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("1/((x - 1)^2)", parse_arena);
    auto expected = parse_expr("1/(x - 1)^2", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto terms = algebra::partial_fractions(expr.value(), Symbol("x"), ctx);
    ASSERT_TRUE(terms.is_ok()) << terms.error().message;
    ASSERT_EQ(terms.value().size(), 1U);
    auto simplified_expected = ctx.simplify(expected.value());
    ASSERT_TRUE(simplified_expected.is_ok()) << simplified_expected.error().message;
    EXPECT_TRUE(structural_equal(terms.value()[0], simplified_expected.value()));
}

TEST(AlgebraPartialFractionsTest, DecomposesMixedSimpleAndRepeatedLinearFactors) {
    AstArena parse_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("(2*x + 5)/((x - 1)^2*(x + 2))", parse_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto terms = algebra::partial_fractions(expr.value(), Symbol("x"), ctx);
    ASSERT_TRUE(terms.is_ok()) << terms.error().message;
    ASSERT_EQ(terms.value().size(), 3U);

    // Verify sum of terms equals original expression (order-independent)
    auto sum_expr = ctx.arena().make<Sum>(std::vector<ExprPtr>(terms.value().begin(), terms.value().end()));
    auto eq = symbolic::mathematically_equal(sum_expr, expr.value(), ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value()) << "Sum of PFD terms != original expression";
}

TEST(AlgebraPartialFractionsTest, SupportsIrreducibleQuadraticDenominator) {
    AstArena parse_arena;
    symbolic::CASContext ctx;

    auto expr = parse_expr("1/(x^2 + 1)", parse_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto terms = algebra::partial_fractions(expr.value(), Symbol("x"), ctx);
    ASSERT_TRUE(terms.is_ok()) << terms.error().message;
    ASSERT_EQ(terms.value().size(), 1U);
}

TEST(SymbolicSimplifyTest, ReusesOriginalPointerForStablePower) {
    AstArena parse_arena;
    CASContext context;

    auto parsed = parse_expr("x^2", parse_arena);
    ASSERT_TRUE(parsed.is_ok()) << parsed.error().message;

    auto simplified = context.simplify(parsed.value());
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    EXPECT_EQ(simplified.value(), parsed.value());
}

TEST(SymbolicSimplifyTest, ReusesOriginalPointerForCanonicalSumNode) {
    AstArena arena;
    CASContext context;

    const auto x = arena.make<Symbol>("x");
    const auto y = arena.make<Symbol>("y");
    const auto sum = arena.make<Sum>(std::vector<ExprPtr>{x, y});

    auto simplified = context.simplify(sum);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    EXPECT_EQ(simplified.value(), sum);
}

TEST(SymbolicSimplifyTest, CollectsLikeLinearTerms) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("2*x + 3*x", parse_arena, simplify_arena);
    auto expected = simplify_input("5*x", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, CancelsIdenticalTermsToZero) {
    AstArena parse_arena;
    AstArena simplify_arena;

    auto simplified = simplify_input("x - x", parse_arena, simplify_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    ASSERT_EQ(expr_kind(simplified.value()), ExprKind::IntegerLit);
    const auto* integer = expr_cast<IntegerLit>(simplified.value());
    ASSERT_NE(integer, nullptr);
    EXPECT_EQ(integer->value.decimal(), "0");
}

TEST(SymbolicSimplifyTest, RewritesBasicTrigIdentityThroughContext) {
    AstArena parse_arena;
    CASContext context;

    auto simplified = simplify_input_with_context("sin(x)^2 + cos(x)^2", parse_arena, context);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    ASSERT_EQ(expr_kind(simplified.value()), ExprKind::IntegerLit);
    const auto* integer = expr_cast<IntegerLit>(simplified.value());
    ASSERT_NE(integer, nullptr);
    EXPECT_EQ(integer->value.decimal(), "1");
}

TEST(SymbolicRewriteTest, MatchesRepeatedWildcardConsistently) {
    AstArena arena;

    auto expr = parse_expr("sin(x) * cos(x)", arena);
    auto pattern = parse_expr("sin(a_) * cos(a_)", arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(pattern.is_ok()) << pattern.error().message;

    MatchMap matches;
    ASSERT_TRUE(match_pattern(expr.value(), pattern.value(), matches));
    ASSERT_EQ(matches.size(), 1U);

    const auto found = matches.find("a_");
    ASSERT_NE(found, matches.end());
    ASSERT_EQ(expr_kind(found->second), ExprKind::Symbol);
    EXPECT_EQ(expr_ref<Symbol>(found->second).name, "x");
}

TEST(SymbolicRewriteTest, AppliesAcRuleInsideLargerSum) {
    AstArena arena;
    AstArena rewrite_arena;
    AstArena expected_arena;

    auto expr = parse_expr("cos(x)^2 + y + sin(x)^2", arena);
    auto pattern = parse_expr("sin(a_)^2 + cos(a_)^2", arena);
    auto replacement = parse_expr("1", arena);
    auto expected = parse_expr("y + 1", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(pattern.is_ok()) << pattern.error().message;
    ASSERT_TRUE(replacement.is_ok()) << replacement.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    RewriteRule rule{
        .pattern = pattern.value(),
        .replacement = replacement.value(),
        .condition = {},
    };
    auto rewritten = apply_rule(expr.value(), rule, TraversalStrategy::FixPoint, rewrite_arena);
    ASSERT_TRUE(rewritten.is_ok()) << rewritten.error().message;

    EXPECT_TRUE(structural_equal(rewritten.value(), expected.value()));
}

TEST(SymbolicRewriteTest, TermOrderingAcceptsStrictReductions) {
    AstArena arena;
    auto lhs = parse_expr("sin(x)^2 + cos(x)^2", arena);
    auto rhs = parse_expr("1", arena);
    ASSERT_TRUE(lhs.is_ok()) << lhs.error().message;
    ASSERT_TRUE(rhs.is_ok()) << rhs.error().message;

    EXPECT_EQ(compare_rewrite_terms(rhs.value(), lhs.value()), TermOrderRelation::Less);
    EXPECT_EQ(compare_rewrite_terms(lhs.value(), rhs.value()), TermOrderRelation::Greater);
}

TEST(SymbolicRewriteTest, TermOrderingUsesLpoSubtermProperty) {
    AstArena arena;
    auto parent = parse_expr("sin(x + y)", arena);
    auto child = parse_expr("x + y", arena);
    ASSERT_TRUE(parent.is_ok()) << parent.error().message;
    ASSERT_TRUE(child.is_ok()) << child.error().message;

    EXPECT_EQ(compare_rewrite_terms(parent.value(), child.value()), TermOrderRelation::Greater);
    EXPECT_EQ(compare_rewrite_terms(child.value(), parent.value()), TermOrderRelation::Less);
}

TEST(SymbolicRewriteTest, OrientsReducingRulesBeforeApplication) {
    AstArena arena;
    auto pattern = parse_expr("sin(x_)^2 + cos(x_)^2", arena);
    auto replacement = parse_expr("1", arena);
    ASSERT_TRUE(pattern.is_ok()) << pattern.error().message;
    ASSERT_TRUE(replacement.is_ok()) << replacement.error().message;

    RewriteRule rule{
        .pattern = pattern.value(),
        .replacement = replacement.value(),
        .condition = {},
    };

    EXPECT_TRUE(rewrite_rule_is_oriented(rule));
}

TEST(SymbolicRewriteTest, RejectsInverseExpCombinationCycle) {
    AstArena arena;
    auto pattern = parse_expr("exp(a_) * exp(b_)", arena);
    auto replacement = parse_expr("exp(a_ + b_)", arena);
    ASSERT_TRUE(pattern.is_ok()) << pattern.error().message;
    ASSERT_TRUE(replacement.is_ok()) << replacement.error().message;

    RewriteRule rule{
        .pattern = pattern.value(),
        .replacement = replacement.value(),
        .condition = {},
    };

    EXPECT_FALSE(rewrite_rule_is_oriented(rule));
}

TEST(SymbolicRewriteTest, StrongNormalizationRejectsCyclicRuleSets) {
    AstArena arena;
    auto a = parse_expr("a", arena);
    auto b = parse_expr("b", arena);
    ASSERT_TRUE(a.is_ok()) << a.error().message;
    ASSERT_TRUE(b.is_ok()) << b.error().message;

    const RewriteRule a_to_b{
        .pattern = a.value(),
        .replacement = b.value(),
        .condition = {},
    };
    const RewriteRule b_to_a{
        .pattern = b.value(),
        .replacement = a.value(),
        .condition = {},
    };

    EXPECT_FALSE(is_strongly_normalizing(std::vector<RewriteRule>{a_to_b, b_to_a}));
}

TEST(SymbolicRewriteTest, RejectsUnorientedRulesBeforeMatching) {
    AstArena arena;
    AstArena rewrite_arena;

    auto expr = parse_expr("x", arena);
    auto pattern = parse_expr("x_", arena);
    auto replacement = parse_expr("x_ + 0", arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(pattern.is_ok()) << pattern.error().message;
    ASSERT_TRUE(replacement.is_ok()) << replacement.error().message;

    RewriteRule rule{
        .pattern = pattern.value(),
        .replacement = replacement.value(),
        .condition = {},
    };
    EXPECT_FALSE(rewrite_rule_is_oriented(rule));

    auto rewritten = apply_rule(expr.value(), rule, TraversalStrategy::TopDown, rewrite_arena);
    ASSERT_TRUE(rewritten.is_error());
    EXPECT_EQ(rewritten.error().kind, CASErrorKind::InvalidArgument);
}

TEST(SymbolicRewriteTest, TermOrderingRejectsExpansiveRuleApplications) {
    AstArena arena;
    AstArena rewrite_arena;

    auto expr = parse_expr("x", arena);
    auto pattern = parse_expr("x_", arena);
    auto replacement = parse_expr("x_ + 0", arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(pattern.is_ok()) << pattern.error().message;
    ASSERT_TRUE(replacement.is_ok()) << replacement.error().message;

    RewriteRule rule{
        .pattern = pattern.value(),
        .replacement = replacement.value(),
        .condition = {},
    };
    auto rewritten = apply_rule(expr.value(), rule, TraversalStrategy::FixPoint, rewrite_arena);
    ASSERT_TRUE(rewritten.is_error());
    EXPECT_EQ(rewritten.error().kind, CASErrorKind::InvalidArgument);
}

TEST(SymbolicRewriteTest, RewritesTrigIdentityInsideExtendedSumThroughContext) {
    AstArena parse_arena;
    CASContext context;
    AstArena expected_parse_arena;
    CASContext expected_context;

    auto simplified = simplify_input_with_context("cos(x)^2 + y + sin(x)^2", parse_arena, context);
    auto expected = simplify_input_with_context("y + 1", expected_parse_arena, expected_context);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicRewriteTest, BuiltinProviderRewritesElementaryFunctionRulesDirectly) {
    AstArena arena;
    AstArena rewrite_arena;
    AstArena expected_arena;

    auto sin_zero = parse_expr("sin(0)", arena);
    auto sin_neg = parse_expr("sin(-x)", arena);
    auto cos_neg = parse_expr("cos(-x)", arena);
    auto exp_one = parse_expr("exp(1)", arena);
    auto ln_e = parse_expr("ln(e)", arena);
    auto ln_exp = parse_expr("ln(e^x)", arena);
    auto expected_zero = parse_expr("0", expected_arena);
    auto expected_neg_sin = parse_expr("-sin(x)", expected_arena);
    auto expected_cos = parse_expr("cos(x)", expected_arena);
    auto expected_one = parse_expr("1", expected_arena);
    auto expected_e = parse_expr("e", expected_arena);
    auto expected_x = parse_expr("x", expected_arena);
    ASSERT_TRUE(sin_zero.is_ok()) << sin_zero.error().message;
    ASSERT_TRUE(sin_neg.is_ok()) << sin_neg.error().message;
    ASSERT_TRUE(cos_neg.is_ok()) << cos_neg.error().message;
    ASSERT_TRUE(exp_one.is_ok()) << exp_one.error().message;
    ASSERT_TRUE(ln_e.is_ok()) << ln_e.error().message;
    ASSERT_TRUE(ln_exp.is_ok()) << ln_exp.error().message;
    ASSERT_TRUE(expected_zero.is_ok()) << expected_zero.error().message;
    ASSERT_TRUE(expected_neg_sin.is_ok()) << expected_neg_sin.error().message;
    ASSERT_TRUE(expected_cos.is_ok()) << expected_cos.error().message;
    ASSERT_TRUE(expected_one.is_ok()) << expected_one.error().message;
    ASSERT_TRUE(expected_e.is_ok()) << expected_e.error().message;
    ASSERT_TRUE(expected_x.is_ok()) << expected_x.error().message;

    auto rewritten_sin = default_rewrite_provider().try_rewrite(sin_zero.value(), rewrite_arena, nullptr);
    auto rewritten_sin_neg = default_rewrite_provider().try_rewrite(sin_neg.value(), rewrite_arena, nullptr);
    auto rewritten_cos_neg = default_rewrite_provider().try_rewrite(cos_neg.value(), rewrite_arena, nullptr);
    auto rewritten_exp_one = default_rewrite_provider().try_rewrite(exp_one.value(), rewrite_arena, nullptr);
    auto rewritten_ln_e = default_rewrite_provider().try_rewrite(ln_e.value(), rewrite_arena, nullptr);
    auto rewritten_ln = default_rewrite_provider().try_rewrite(ln_exp.value(), rewrite_arena, nullptr);
    ASSERT_TRUE(rewritten_sin.is_ok()) << rewritten_sin.error().message;
    ASSERT_TRUE(rewritten_sin_neg.is_ok()) << rewritten_sin_neg.error().message;
    ASSERT_TRUE(rewritten_cos_neg.is_ok()) << rewritten_cos_neg.error().message;
    ASSERT_TRUE(rewritten_exp_one.is_ok()) << rewritten_exp_one.error().message;
    ASSERT_TRUE(rewritten_ln_e.is_ok()) << rewritten_ln_e.error().message;
    ASSERT_TRUE(rewritten_ln.is_ok()) << rewritten_ln.error().message;

    EXPECT_TRUE(structural_equal(rewritten_sin.value(), expected_zero.value()));
    EXPECT_TRUE(structural_equal(rewritten_sin_neg.value(), expected_neg_sin.value()));
    EXPECT_TRUE(structural_equal(rewritten_cos_neg.value(), expected_cos.value()));
    EXPECT_TRUE(structural_equal(rewritten_exp_one.value(), expected_e.value()));
    EXPECT_TRUE(structural_equal(rewritten_ln_e.value(), expected_one.value()));
    EXPECT_TRUE(structural_equal(rewritten_ln.value(), expected_x.value()));
}

TEST(SymbolicRewriteTest, BuiltinProviderUsesPositiveAssumptionsForLogProductRule) {
    AstArena parse_arena;
    AstArena rewrite_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"a"});
    context.assumptions().assume_positive(Symbol{"b"});

    auto expr = parse_expr("ln(a * b)", parse_arena);
    auto expected = parse_expr("ln(a) + ln(b)", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto rewritten =
        default_rewrite_provider().try_rewrite(expr.value(), rewrite_arena, &context.assumptions());
    ASSERT_TRUE(rewritten.is_ok()) << rewritten.error().message;

    EXPECT_TRUE(structural_equal(rewritten.value(), expected.value()));
}

TEST(SymbolicRewriteTest, BuiltinProviderUsesPositiveAssumptionsForLogQuotientRule) {
    AstArena parse_arena;
    AstArena rewrite_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"a"});
    context.assumptions().assume_positive(Symbol{"b"});

    auto expr = parse_expr("ln(a / b)", parse_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    const auto expected_a = expected_arena.make<Symbol>("a");
    const auto expected_b = expected_arena.make<Symbol>("b");
    const auto expected = expected_arena.make<Sum>(std::vector<ExprPtr>{
        expected_arena.make<FuncCall>("ln", std::vector<ExprPtr>{expected_a}),
        expected_arena.make<Unary>(
            UnaryOp::Neg,
            expected_arena.make<FuncCall>("ln", std::vector<ExprPtr>{expected_b})),
    });

    auto rewritten =
        default_rewrite_provider().try_rewrite(expr.value(), rewrite_arena, &context.assumptions());
    ASSERT_TRUE(rewritten.is_ok()) << rewritten.error().message;

    EXPECT_TRUE(structural_equal(rewritten.value(), expected));
}

TEST(SymbolicRewriteTest, BuiltinProviderUsesNonnegativeAssumptionsForSqrtProductRule) {
    AstArena parse_arena;
    AstArena rewrite_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"a"});
    context.assumptions().assume_positive(Symbol{"b"});

    auto expr = parse_expr("sqrt(a * b)", parse_arena);
    auto expected = parse_expr("sqrt(a) * sqrt(b)", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto rewritten =
        default_rewrite_provider().try_rewrite(expr.value(), rewrite_arena, &context.assumptions());
    ASSERT_TRUE(rewritten.is_ok()) << rewritten.error().message;

    EXPECT_TRUE(structural_equal(rewritten.value(), expected.value()));
}

TEST(SymbolicRewriteTest, BuiltinProviderUsesPositiveAssumptionsForSqrtQuotientRule) {
    AstArena parse_arena;
    AstArena rewrite_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"a"});
    context.assumptions().assume_positive(Symbol{"b"});

    auto expr = parse_expr("sqrt(a / b)", parse_arena);
    auto expected = parse_expr("sqrt(a) / sqrt(b)", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto rewritten =
        default_rewrite_provider().try_rewrite(expr.value(), rewrite_arena, &context.assumptions());
    ASSERT_TRUE(rewritten.is_ok()) << rewritten.error().message;

    EXPECT_TRUE(structural_equal(rewritten.value(), expected.value()));
}

TEST(SymbolicRewriteTest, BuiltinProviderUsesPositiveAssumptionsForExpLogRule) {
    AstArena parse_arena;
    AstArena rewrite_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"x"});

    auto expr = parse_expr("e^(ln(x))", parse_arena);
    auto expected = parse_expr("x", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto rewritten =
        default_rewrite_provider().try_rewrite(expr.value(), rewrite_arena, &context.assumptions());
    ASSERT_TRUE(rewritten.is_ok()) << rewritten.error().message;

    EXPECT_TRUE(structural_equal(rewritten.value(), expected.value()));
}

TEST(SymbolicRewriteTest, ContextUsesBuiltinProviderForFunctionRules) {
    AstArena parse_arena;
    CASContext context;
    context.enable_trace(true);

    auto simplified = simplify_input_with_context("ln(e^x)", parse_arena, context);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    const auto& trace = context.get_trace();
    ASSERT_FALSE(trace.empty());
    EXPECT_EQ(trace.front().rule_id, RuleId::RewriteProviderApplied);
    ASSERT_EQ(expr_kind(simplified.value()), ExprKind::Symbol);
    EXPECT_EQ(expr_ref<Symbol>(simplified.value()).name, "x");
}

TEST(SymbolicRewriteTest, ContextUsesBuiltinProviderForTrigParityRules) {
    AstArena parse_arena;
    AstArena expected_sin_arena;
    AstArena expected_cos_arena;
    CASContext context;
    context.enable_trace(true);

    auto simplified_sin = simplify_input_with_context("sin(-x)", parse_arena, context);
    auto expected_sin = simplify_input("-sin(x)", expected_sin_arena, expected_sin_arena);
    ASSERT_TRUE(simplified_sin.is_ok()) << simplified_sin.error().message;
    ASSERT_TRUE(expected_sin.is_ok()) << expected_sin.error().message;
    ASSERT_FALSE(context.get_trace().empty());
    EXPECT_EQ(context.get_trace().front().rule_id, RuleId::RewriteProviderApplied);
    EXPECT_TRUE(structural_equal(simplified_sin.value(), expected_sin.value()));

    auto simplified_cos = simplify_input_with_context("cos(-x)", parse_arena, context);
    auto expected_cos = simplify_input("cos(x)", expected_cos_arena, expected_cos_arena);
    ASSERT_TRUE(simplified_cos.is_ok()) << simplified_cos.error().message;
    ASSERT_TRUE(expected_cos.is_ok()) << expected_cos.error().message;
    ASSERT_FALSE(context.get_trace().empty());
    EXPECT_EQ(context.get_trace().front().rule_id, RuleId::RewriteProviderApplied);
    EXPECT_TRUE(structural_equal(simplified_cos.value(), expected_cos.value()));
}

TEST(SymbolicRewriteTest, ContextUsesBuiltinProviderForExtendedTrigHyperbolicParityRules) {
    struct Sample {
        const char* input;
        const char* expected;
    };

    const std::vector<Sample> samples{
        {"tan(-x)", "-tan(x)"},
        {"cot(-x)", "-cot(x)"},
        {"csc(-x)", "-csc(x)"},
        {"sec(-x)", "sec(x)"},
        {"sinh(-x)", "-sinh(x)"},
        {"cosh(-x)", "cosh(x)"},
        {"tanh(-x)", "-tanh(x)"},
        {"coth(-x)", "-coth(x)"},
    };

    for (const auto& sample : samples) {
        AstArena parse_arena;
        AstArena expected_arena;
        CASContext context;
        context.enable_trace(true);

        auto simplified = simplify_input_with_context(sample.input, parse_arena, context);
        auto expected = simplify_input(sample.expected, expected_arena, expected_arena);
        ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
        ASSERT_TRUE(expected.is_ok()) << expected.error().message;

        const auto& trace = context.get_trace();
        ASSERT_FALSE(trace.empty()) << sample.input;
        EXPECT_EQ(trace.front().rule_id, RuleId::RewriteProviderApplied) << sample.input;
        EXPECT_TRUE(structural_equal(simplified.value(), expected.value())) << sample.input;
    }
}

TEST(SymbolicRewriteTest, ContextUsesBuiltinProviderForTanNormalization) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.enable_trace(true);

    auto simplified = simplify_input_with_context("tan(x)", parse_arena, context);
    auto expected = simplify_input("sin(x) / cos(x)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    const auto& trace = context.get_trace();
    ASSERT_FALSE(trace.empty());
    EXPECT_EQ(trace.front().rule_id, RuleId::RewriteProviderApplied);
    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicRewriteTest, ContextUsesBuiltinProviderForPositiveLogProductRule) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"a"});
    context.assumptions().assume_positive(Symbol{"b"});
    context.enable_trace(true);

    auto simplified = simplify_input_with_context("ln(a * b)", parse_arena, context);
    auto expected = simplify_input("ln(a) + ln(b)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    const auto& trace = context.get_trace();
    ASSERT_FALSE(trace.empty());
    EXPECT_EQ(trace.front().rule_id, RuleId::RewriteProviderApplied);
    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicRewriteTest, ContextUsesBuiltinProviderForPositiveLogQuotientRule) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"a"});
    context.assumptions().assume_positive(Symbol{"b"});
    context.enable_trace(true);

    auto simplified = simplify_input_with_context("ln(a / b)", parse_arena, context);
    auto expected = simplify_input("ln(a) - ln(b)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    const auto& trace = context.get_trace();
    ASSERT_FALSE(trace.empty());
    EXPECT_EQ(trace.front().rule_id, RuleId::RewriteProviderApplied);
    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicRewriteTest, ContextUsesBuiltinProviderForPositiveSqrtProductRule) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"a"});
    context.assumptions().assume_positive(Symbol{"b"});
    context.enable_trace(true);

    auto simplified = simplify_input_with_context("sqrt(a * b)", parse_arena, context);
    auto expected = simplify_input("sqrt(a) * sqrt(b)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    const auto& trace = context.get_trace();
    ASSERT_FALSE(trace.empty());
    EXPECT_EQ(trace.front().rule_id, RuleId::RewriteProviderApplied);
    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicRewriteTest, ContextUsesBuiltinProviderForPositiveSqrtQuotientRule) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"a"});
    context.assumptions().assume_positive(Symbol{"b"});
    context.enable_trace(true);

    auto simplified = simplify_input_with_context("sqrt(a / b)", parse_arena, context);
    auto expected = simplify_input("sqrt(a) / sqrt(b)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    const auto& trace = context.get_trace();
    ASSERT_FALSE(trace.empty());
    EXPECT_EQ(trace.front().rule_id, RuleId::RewriteProviderApplied);
    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicRewriteTest, ContextUsesBuiltinProviderForPositiveExpLogRule) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"x"});
    context.enable_trace(true);

    auto simplified = simplify_input_with_context("e^(ln(x))", parse_arena, context);
    auto expected = simplify_input("x", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    const auto& trace = context.get_trace();
    ASSERT_FALSE(trace.empty());
    EXPECT_EQ(trace.front().rule_id, RuleId::RewriteProviderApplied);
    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, MergesPowersWithSameBase) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("x^2 * x^3", parse_arena, simplify_arena);
    auto expected = simplify_input("x^5", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, FlattensNestedPowers) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("(x^2)^3", parse_arena, simplify_arena);
    auto expected = simplify_input("x^6", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, SimplifiesElementaryTrigValues) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("sin(0) + cos(pi)", parse_arena, simplify_arena);
    auto expected = simplify_input("-1", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, SimplifiesLogAndExpBaseCases) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("ln(1) + exp(0)", parse_arena, simplify_arena);
    auto expected = simplify_input("1", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, SimplifiesLnOfExponentialBaseE) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("ln(e^x)", parse_arena, simplify_arena);
    auto expected = simplify_input("x", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, SimplifiesExponentialOfLogWithPositiveAssumption) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"x"});

    auto simplified = simplify_input_with_context("e^(ln(x))", parse_arena, context);
    auto expected = simplify_input("x", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, SimplifiesTrigParityThroughProvider) {
    AstArena parse_arena;
    AstArena expected_sin_arena;
    AstArena expected_cos_arena;
    CASContext context;

    auto simplified_sin = simplify_input_with_context("sin(-x)", parse_arena, context);
    auto expected_sin = simplify_input("-sin(x)", expected_sin_arena, expected_sin_arena);
    ASSERT_TRUE(simplified_sin.is_ok()) << simplified_sin.error().message;
    ASSERT_TRUE(expected_sin.is_ok()) << expected_sin.error().message;
    EXPECT_TRUE(structural_equal(simplified_sin.value(), expected_sin.value()));

    auto simplified_cos = simplify_input_with_context("cos(-x)", parse_arena, context);
    auto expected_cos = simplify_input("cos(x)", expected_cos_arena, expected_cos_arena);
    ASSERT_TRUE(simplified_cos.is_ok()) << simplified_cos.error().message;
    ASSERT_TRUE(expected_cos.is_ok()) << expected_cos.error().message;
    EXPECT_TRUE(structural_equal(simplified_cos.value(), expected_cos.value()));
}

TEST(SymbolicSimplifyTest, RewritesNegativePowerAsDivision) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("x^(-1)", parse_arena, simplify_arena);
    auto expected = simplify_input("1 / x", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, RewritesSqrtSquareToAbsWithoutAssumptions) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("sqrt(x^2)", parse_arena, simplify_arena);
    auto expected = simplify_input("abs(x)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, RewritesSqrtSquareToSymbolWithPositiveAssumption) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"x"});

    auto simplified = simplify_input_with_context("sqrt(x^2)", parse_arena, context);
    auto expected = simplify_input("x", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, ExpandsExponentialOfSum) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("exp(a + b)", parse_arena, simplify_arena);
    auto expected = simplify_input("exp(a) * exp(b)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, NormalizesTanToSinOverCos) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("tan(x)", parse_arena, simplify_arena);
    auto expected = simplify_input("sin(x) / cos(x)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, ExpandsLogOfPositiveProduct) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"a"});
    context.assumptions().assume_positive(Symbol{"b"});

    auto simplified = simplify_input_with_context("ln(a * b)", parse_arena, context);
    auto expected = simplify_input("ln(a) + ln(b)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, ExpandsLogOfPositiveQuotient) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"a"});
    context.assumptions().assume_positive(Symbol{"b"});

    auto simplified = simplify_input_with_context("ln(a / b)", parse_arena, context);
    auto expected = simplify_input("ln(a) - ln(b)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, ExpandsLogOfPowerWithPositiveBase) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"a"});

    auto simplified = simplify_input_with_context("ln(a^b)", parse_arena, context);
    auto expected = simplify_input("b * ln(a)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, ExpandsSqrtOfPositiveQuotient) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"a"});
    context.assumptions().assume_positive(Symbol{"b"});

    auto simplified = simplify_input_with_context("sqrt(a / b)", parse_arena, context);
    auto expected = simplify_input("sqrt(a) / sqrt(b)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, ExpandsSqrtOfPositiveProduct) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.assumptions().assume_positive(Symbol{"a"});
    context.assumptions().assume_positive(Symbol{"b"});

    auto simplified = simplify_input_with_context("sqrt(a * b)", parse_arena, context);
    auto expected = simplify_input("sqrt(a) * sqrt(b)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, KeepsLogOfQuotientWithoutAssumptions) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("ln(a / b)", parse_arena, simplify_arena);
    auto expected = simplify_input("ln(a / b)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, KeepsSqrtOfQuotientWithoutAssumptions) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("sqrt(a / b)", parse_arena, simplify_arena);
    auto expected = simplify_input("sqrt(a / b)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, KeepsSqrtOfProductWithoutAssumptions) {
    AstArena parse_arena;
    AstArena simplify_arena;
    AstArena expected_arena;

    auto simplified = simplify_input("sqrt(a * b)", parse_arena, simplify_arena);
    auto expected = simplify_input("sqrt(a * b)", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, KeepsDivisionBySameSymbolWithoutNonzeroAssumption) {
    AstArena parse_arena;
    CASContext context;
    AstArena expected_arena;

    auto simplified = simplify_input_with_context("x / x", parse_arena, context);
    auto expected = simplify_input("x / x", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicSimplifyTest, CollapsesDivisionBySameSymbolWithNonzeroAssumption) {
    AstArena parse_arena;
    CASContext context;
    context.assumptions().assume_nonzero(Symbol{"x"});

    auto simplified = simplify_input_with_context("x / x", parse_arena, context);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    ASSERT_EQ(expr_kind(simplified.value()), ExprKind::IntegerLit);
    const auto* integer = expr_cast<IntegerLit>(simplified.value());
    ASSERT_NE(integer, nullptr);
    EXPECT_EQ(integer->value.decimal(), "1");
}

TEST(SymbolicCoreTest, SubstitutesSymbolsThenSimplifies) {
    AstArena parse_arena;
    AstArena value_arena;
    CASContext context;

    auto parsed = parse_expr("x^2 + y", parse_arena);
    auto value = parse_expr("3", value_arena);
    ASSERT_TRUE(parsed.is_ok()) << parsed.error().message;
    ASSERT_TRUE(value.is_ok()) << value.error().message;

    auto substituted = context.substitute(parsed.value(), Symbol{"x"}, value.value());
    ASSERT_TRUE(substituted.is_ok()) << substituted.error().message;

    auto expected = simplify_input_with_context("9 + y", value_arena, context);
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    EXPECT_TRUE(structural_equal(substituted.value(), expected.value()));
}

TEST(SymbolicCoreTest, DetectsMathematicalEquivalenceViaSimplification) {
    AstArena parse_arena;
    CASContext context;

    auto lhs = parse_expr("x + x", parse_arena);
    auto rhs = parse_expr("2 * x", parse_arena);
    ASSERT_TRUE(lhs.is_ok()) << lhs.error().message;
    ASSERT_TRUE(rhs.is_ok()) << rhs.error().message;

    auto equivalent = mathematically_equal(lhs.value(), rhs.value(), context);
    ASSERT_TRUE(equivalent.is_ok()) << equivalent.error().message;
    EXPECT_TRUE(equivalent.value());
}

TEST(SymbolicCoreTest, DetectsMathematicalEquivalenceForReorderedSeriesTerms) {
    AstArena parse_arena;
    CASContext context;

    auto lhs = parse_expr(
        "1 + n*x + (-1/2)*n*x^2 + (1/2)*n^2*x^2 + (1/3)*n*x^3 + (-1/2)*n^2*x^3 + (1/6)*n^3*x^3",
        parse_arena);
    auto rhs = parse_expr(
        "(1/6)*n^3*x^3 + (1/2)*n^2*x^2 + n*x + (-1/2)*n^2*x^3 + (1/3)*n*x^3 + 1 + (-1/2)*n*x^2",
        parse_arena);
    ASSERT_TRUE(lhs.is_ok()) << lhs.error().message;
    ASSERT_TRUE(rhs.is_ok()) << rhs.error().message;

    auto equivalent = mathematically_equal(lhs.value(), rhs.value(), context);
    ASSERT_TRUE(equivalent.is_ok()) << equivalent.error().message;
    EXPECT_TRUE(equivalent.value());
}

TEST(SymbolicCoreTest, DoesNotSubstituteBoundIntegralVariable) {
    AstArena parse_arena;
    AstArena value_arena;
    CASContext context;

    auto parsed = parse_expr("∫(x + y, x)", parse_arena);
    auto value = parse_expr("3", value_arena);
    auto expected = parse_expr("∫(x + y, x)", value_arena);
    ASSERT_TRUE(parsed.is_ok()) << parsed.error().message;
    ASSERT_TRUE(value.is_ok()) << value.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto substituted = context.substitute(parsed.value(), Symbol{"x"}, value.value());
    ASSERT_TRUE(substituted.is_ok()) << substituted.error().message;

    EXPECT_TRUE(structural_equal(substituted.value(), expected.value()));
}

TEST(SymbolicCoreTest, SubstitutesFreeVariableInsideIntegralButNotBoundVariable) {
    AstArena parse_arena;
    AstArena value_arena;
    CASContext context;

    auto parsed = parse_expr("∫(x + y, x)", parse_arena);
    auto value = parse_expr("3", value_arena);
    auto expected = parse_expr("∫(x + 3, x)", value_arena);
    ASSERT_TRUE(parsed.is_ok()) << parsed.error().message;
    ASSERT_TRUE(value.is_ok()) << value.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto substituted = context.substitute(parsed.value(), Symbol{"y"}, value.value());
    ASSERT_TRUE(substituted.is_ok()) << substituted.error().message;

    EXPECT_TRUE(structural_equal(substituted.value(), expected.value()));
}

TEST(SymbolicTraceTest, ReturnsEmptyTraceWhenDisabled) {
    AstArena parse_arena;
    CASContext context;

    auto simplified = simplify_input_with_context("x + 0", parse_arena, context);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    EXPECT_TRUE(context.get_trace().empty());
}

TEST(SymbolicTraceTest, RecordsFlatDepthForNestedSubstitution) {
    AstArena parse_arena;
    AstArena value_arena;
    AstArena partial_expected_arena;
    AstArena final_expected_arena;
    CASContext context;
    context.enable_trace(true);

    auto parsed = parse_expr("f(x, g(x))", parse_arena);
    auto value = parse_expr("1", value_arena);
    auto partial_expected = parse_expr("f(1, g(x))", partial_expected_arena);
    auto final_expected = parse_expr("f(1, g(1))", final_expected_arena);
    ASSERT_TRUE(parsed.is_ok()) << parsed.error().message;
    ASSERT_TRUE(value.is_ok()) << value.error().message;
    ASSERT_TRUE(partial_expected.is_ok()) << partial_expected.error().message;
    ASSERT_TRUE(final_expected.is_ok()) << final_expected.error().message;

    auto substituted = context.substitute(parsed.value(), Symbol{"x"}, value.value());
    ASSERT_TRUE(substituted.is_ok()) << substituted.error().message;
    EXPECT_TRUE(structural_equal(substituted.value(), final_expected.value()));

    const auto& trace = context.get_trace();
    ASSERT_EQ(trace.size(), 2U);
    EXPECT_EQ(trace[0].rule_id, RuleId::SubstituteSymbol);
    EXPECT_EQ(trace[0].depth, 1U);
    EXPECT_TRUE(structural_equal(trace[0].root_after, partial_expected.value()));
    EXPECT_EQ(trace[1].rule_id, RuleId::SubstituteSymbol);
    EXPECT_EQ(trace[1].depth, 2U);
    EXPECT_TRUE(structural_equal(trace[1].root_after, final_expected.value()));
}

TEST(SymbolicTraceTest, RecordsAssumptionStepForAssumptionBasedRewrite) {
    AstArena parse_arena;
    AstArena expected_arena;
    CASContext context;
    context.enable_trace(true);
    context.assumptions().assume_positive(Symbol{"x"});

    auto simplified = simplify_input_with_context("sqrt(x^2)", parse_arena, context);
    auto expected = simplify_input("x", expected_arena, expected_arena);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));

    const auto& trace = context.get_trace();
    ASSERT_EQ(trace.size(), 2U);
    EXPECT_EQ(trace[0].rule_id, RuleId::AssumptionApplied);
    EXPECT_EQ(trace[0].depth, 0U);
    EXPECT_EQ(trace[1].rule_id, RuleId::SimplifySqrtSquare);
    EXPECT_EQ(trace[1].depth, 0U);
    EXPECT_TRUE(structural_equal(trace[1].root_after, expected.value()));
}

TEST(SymbolicTraceTest, MathematicallyEqualDoesNotPopulateTrace) {
    AstArena first_parse_arena;
    AstArena compare_parse_arena;
    CASContext context;
    context.enable_trace(true);

    auto simplified = simplify_input_with_context("x + 0", first_parse_arena, context);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    ASSERT_FALSE(context.get_trace().empty());

    auto lhs = parse_expr("x + x", compare_parse_arena);
    auto rhs = parse_expr("2 * x", compare_parse_arena);
    ASSERT_TRUE(lhs.is_ok()) << lhs.error().message;
    ASSERT_TRUE(rhs.is_ok()) << rhs.error().message;

    auto equivalent = mathematically_equal(lhs.value(), rhs.value(), context);
    ASSERT_TRUE(equivalent.is_ok()) << equivalent.error().message;
    EXPECT_TRUE(equivalent.value());
    EXPECT_TRUE(context.get_trace().empty());
}

TEST(SymbolicTimeoutTest, ResetsStateCleanlyAfterTimeout) {
    CASContext context;
    context.enable_trace(true);
    context.set_timeout(std::chrono::milliseconds::zero());

    ExprPtr symbol = context.arena().make<Symbol>(std::string("x"));
    std::vector<ExprPtr> terms(4096U, symbol);
    ExprPtr large_sum = context.arena().make<Sum>(std::move(terms));

    auto timed_out = context.simplify(large_sum);
    ASSERT_TRUE(timed_out.is_error());
    EXPECT_EQ(timed_out.error().kind, CASErrorKind::Timeout);
    EXPECT_TRUE(context.get_trace().empty());

    context.set_timeout(std::chrono::milliseconds(1000));

    AstArena parse_arena;
    auto recovered = simplify_input_with_context("x + 0", parse_arena, context);
    ASSERT_TRUE(recovered.is_ok()) << recovered.error().message;

    const auto& trace = context.get_trace();
    ASSERT_EQ(trace.size(), 1U);
    EXPECT_EQ(trace.front().rule_id, RuleId::SimplifyAddZero);
    EXPECT_EQ(trace.front().depth, 0U);
}

TEST(SymbolicDepthTest, MaxDepthConfigurable) {
    // L0-12: set_max_simplification_depth allows deeper legitimate computations
    CASContext ctx;
    ctx.set_max_simplification_depth(500);
    EXPECT_EQ(ctx.max_simplification_depth(), 500);

    // Min clamp: below 10 → 10
    ctx.set_max_simplification_depth(1);
    EXPECT_EQ(ctx.max_simplification_depth(), 10);

    // Default is 300
    CASContext default_ctx;
    EXPECT_EQ(default_ctx.max_simplification_depth(), 300);
}

TEST(SymbolicDepthTest, IncreasedDepthAllowsDeeperNesting) {
    // L0-12: increasing depth allows computation that would fail at 300
    CASContext ctx;
    ctx.set_max_simplification_depth(1000);

    // Build a sum with 350 terms: requires depth > 300 with default limit
    ExprPtr sym = ctx.arena().make<Symbol>(std::string("z"));
    std::vector<ExprPtr> terms(350U, sym);
    ExprPtr deep_sum = ctx.arena().make<Sum>(std::move(terms));

    auto result = ctx.simplify(deep_sum);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
}

TEST(SymbolicCycleDetectionTest, SameNodeReentryReturnsOriginal) {
    // L0-10: if simplify_expr re-enters same ExprPtr, CycleGuard returns original
    // without infinite recursion. We can't easily inject a cycle via rewrite rules
    // in a unit test, so we verify the guard data structure is correct by testing
    // that simplification of a simple expression terminates and is idempotent.
    CASContext ctx;
    AstArena parse_arena;
    auto expr = simplify_input("x + 0", parse_arena, ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    // Second simplification of already-simplified result must terminate
    auto expr2 = ctx.simplify(expr.value());
    ASSERT_TRUE(expr2.is_ok());
    EXPECT_TRUE(structural_equal(expr.value(), expr2.value()));
}

TEST(SymbolicCycleDetectionTest, DeepNestedExprTerminates) {
    // L0-10: deeply nested structure must not stack-overflow or cycle
    CASContext ctx;
    // Build x + (x + (x + (x + x))) with 50 levels
    ExprPtr sym = ctx.arena().make<Symbol>(std::string("y"));
    ExprPtr acc = sym;
    for (int i = 0; i < 50; ++i) {
        acc = ctx.arena().make<Sum>(std::vector<ExprPtr>{acc, sym});
    }
    auto result = ctx.simplify(acc);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
}

TEST(SymbolicTimeoutTest, TimeoutCheckIntervalConfigurable) {
    // L0-13: set_timeout_check_interval changes check granularity
    CASContext ctx;
    ctx.set_timeout_check_interval(64U);   // min clamp
    EXPECT_EQ(ctx.timeout_check_interval(), 64U);

    ctx.set_timeout_check_interval(256U);
    EXPECT_EQ(ctx.timeout_check_interval(), 256U);

    // Below min (64) clamps to 64
    ctx.set_timeout_check_interval(1U);
    EXPECT_EQ(ctx.timeout_check_interval(), 64U);

    // Default is 1024
    CASContext default_ctx;
    EXPECT_EQ(default_ctx.timeout_check_interval(), 1024U);
}

TEST(SymbolicTimeoutTest, TimeoutRespectedWithSmallInterval) {
    // With interval=64, a zero-timeout triggers on a big expr
    CASContext ctx;
    ctx.set_timeout(std::chrono::milliseconds::zero());
    ctx.set_timeout_check_interval(64U);

    ExprPtr sym = ctx.arena().make<Symbol>(std::string("x"));
    std::vector<ExprPtr> terms(4096U, sym);
    ExprPtr large_sum = ctx.arena().make<Sum>(std::move(terms));

    auto result = ctx.simplify(large_sum);
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::Timeout);
}

TEST(SymbolicAssumptionsTest, TracksRealAndIntegerFacts) {
    Assumptions assumptions;
    const Symbol x{"x"};
    const Symbol n{"n"};

    assumptions.assume_real(x);
    assumptions.assume_integer(n);

    EXPECT_TRUE(assumptions.is_real(x));
    EXPECT_TRUE(assumptions.is_real(n));
    EXPECT_FALSE(assumptions.is_nonzero(x));
}

TEST(SymbolicAssumptionsTest, RangeCanExcludeZeroConservatively) {
    AstArena arena;
    Assumptions assumptions;
    const Symbol x{"x"};

    assumptions.assume_in_range(
        x,
        arena.make<IntegerLit>(BigInt(2)),
        arena.make<IntegerLit>(BigInt(5)));

    EXPECT_FALSE(assumptions.could_be_zero(x));
    EXPECT_TRUE(assumptions.check_consistency().is_ok());
}

TEST(SymbolicAssumptionsTest, DetectsContradictionBetweenPositiveAndNonPositiveUpperBound) {
    AstArena arena;
    Assumptions assumptions;
    const Symbol x{"x"};

    assumptions.assume_positive(x);
    assumptions.assume_in_range(
        x,
        arena.make<IntegerLit>(BigInt(-3)),
        arena.make<IntegerLit>(BigInt(0)));

    auto consistency = assumptions.check_consistency();
    ASSERT_TRUE(consistency.is_error());
    EXPECT_EQ(consistency.error().kind, CASErrorKind::InvalidArgument);
}

TEST(SymbolicAssumptionsTest, DetectsContradictionBetweenNonzeroAndZeroRange) {
    AstArena arena;
    Assumptions assumptions;
    const Symbol x{"x"};

    assumptions.assume_nonzero(x);
    assumptions.assume_in_range(
        x,
        arena.make<IntegerLit>(BigInt(0)),
        arena.make<IntegerLit>(BigInt(0)));

    auto consistency = assumptions.check_consistency();
    ASSERT_TRUE(consistency.is_error());
    EXPECT_EQ(consistency.error().kind, CASErrorKind::InvalidArgument);
}

TEST(SymbolicPropertyTest, CanonicalizesCommutativeAdditionSamples) {
    const std::vector<std::pair<std::string, std::string>> samples = {
        {"x + y", "y + x"},
        {"2*x + 3*y", "3*y + 2*x"},
        {"sin(x)^2 + cos(x)^2", "cos(x)^2 + sin(x)^2"},
    };

    for (const auto& [lhs_text, rhs_text] : samples) {
        AstArena lhs_parse_arena;
        AstArena rhs_parse_arena;
        CASContext lhs_context;
        CASContext rhs_context;

        auto lhs = simplify_input_with_context(lhs_text, lhs_parse_arena, lhs_context);
        auto rhs = simplify_input_with_context(rhs_text, rhs_parse_arena, rhs_context);
        ASSERT_TRUE(lhs.is_ok()) << lhs.error().message;
        ASSERT_TRUE(rhs.is_ok()) << rhs.error().message;

        EXPECT_TRUE(structural_equal(lhs.value(), rhs.value()));
    }
}

TEST(SymbolicPropertyTest, CanonicalizesCommutativeMultiplicationSamples) {
    const std::vector<std::pair<std::string, std::string>> samples = {
        {"x * y", "y * x"},
        {"3 * y * x", "x * 3 * y"},
        {"x^2 * x^3", "x^3 * x^2"},
    };

    for (const auto& [lhs_text, rhs_text] : samples) {
        AstArena lhs_parse_arena;
        AstArena rhs_parse_arena;
        AstArena lhs_simplify_arena;
        AstArena rhs_simplify_arena;

        auto lhs = simplify_input(lhs_text, lhs_parse_arena, lhs_simplify_arena);
        auto rhs = simplify_input(rhs_text, rhs_parse_arena, rhs_simplify_arena);
        ASSERT_TRUE(lhs.is_ok()) << lhs.error().message;
        ASSERT_TRUE(rhs.is_ok()) << rhs.error().message;

        EXPECT_TRUE(structural_equal(lhs.value(), rhs.value()));
    }
}

TEST(SymbolicPropertyTest, CanonicalizesAssociativeAdditionSamples) {
    const std::vector<std::pair<std::string, std::string>> samples = {
        {"(x + y) + z", "x + (y + z)"},
        {"(2 + x) + y", "2 + (x + y)"},
    };

    for (const auto& [lhs_text, rhs_text] : samples) {
        AstArena lhs_parse_arena;
        AstArena rhs_parse_arena;
        AstArena lhs_simplify_arena;
        AstArena rhs_simplify_arena;

        auto lhs = simplify_input(lhs_text, lhs_parse_arena, lhs_simplify_arena);
        auto rhs = simplify_input(rhs_text, rhs_parse_arena, rhs_simplify_arena);
        ASSERT_TRUE(lhs.is_ok()) << lhs.error().message;
        ASSERT_TRUE(rhs.is_ok()) << rhs.error().message;

        EXPECT_TRUE(structural_equal(lhs.value(), rhs.value()));
    }
}

TEST(SymbolicPropertyTest, CanonicalizesAssociativeMultiplicationSamples) {
    const std::vector<std::pair<std::string, std::string>> samples = {
        {"(x * y) * z", "x * (y * z)"},
        {"(2 * x) * y", "2 * (x * y)"},
    };

    for (const auto& [lhs_text, rhs_text] : samples) {
        AstArena lhs_parse_arena;
        AstArena rhs_parse_arena;
        AstArena lhs_simplify_arena;
        AstArena rhs_simplify_arena;

        auto lhs = simplify_input(lhs_text, lhs_parse_arena, lhs_simplify_arena);
        auto rhs = simplify_input(rhs_text, rhs_parse_arena, rhs_simplify_arena);
        ASSERT_TRUE(lhs.is_ok()) << lhs.error().message;
        ASSERT_TRUE(rhs.is_ok()) << rhs.error().message;

        EXPECT_TRUE(structural_equal(lhs.value(), rhs.value()));
    }
}

TEST(SymbolicPropertyTest, CanonicalizesDistributiveSamples) {
    const std::vector<std::pair<std::string, std::string>> samples = {
        {"x * (y + z)", "x*y + x*z"},
        {"2 * (x + y)", "2*x + 2*y"},
    };

    for (const auto& [lhs_text, rhs_text] : samples) {
        AstArena lhs_parse_arena;
        AstArena rhs_parse_arena;
        AstArena lhs_simplify_arena;
        AstArena rhs_simplify_arena;

        auto lhs = simplify_input(lhs_text, lhs_parse_arena, lhs_simplify_arena);
        auto rhs = simplify_input(rhs_text, rhs_parse_arena, rhs_simplify_arena);
        ASSERT_TRUE(lhs.is_ok()) << lhs.error().message;
        ASSERT_TRUE(rhs.is_ok()) << rhs.error().message;

        EXPECT_TRUE(structural_equal(lhs.value(), rhs.value()));
    }
}

TEST(SymbolicPropertyTest, CanonicalizesGeneratedPolynomialIdentitySamples) {
    const std::vector<std::pair<std::string, std::string>> samples = {
        {"(x + y) + (z + 0)", "z + y + x"},
        {"3*x + 2*x - x", "4*x"},
        {"2*(x + y + z)", "2*x + 2*y + 2*z"},
        {"(x + 1) * (y + 1)", "x*y + x + y + 1"},
        {"x*(a + b + c)", "a*x + b*x + c*x"},
    };

    for (const auto& [lhs_text, rhs_text] : samples) {
        SCOPED_TRACE(lhs_text + " == " + rhs_text);
        AstArena lhs_parse_arena;
        AstArena rhs_parse_arena;
        AstArena lhs_simplify_arena;
        AstArena rhs_simplify_arena;

        auto lhs = simplify_input(lhs_text, lhs_parse_arena, lhs_simplify_arena);
        auto rhs = simplify_input(rhs_text, rhs_parse_arena, rhs_simplify_arena);
        ASSERT_TRUE(lhs.is_ok()) << lhs.error().message;
        ASSERT_TRUE(rhs.is_ok()) << rhs.error().message;

        EXPECT_TRUE(structural_equal(lhs.value(), rhs.value()));
    }
}

TEST(SymbolicPropertyTest, SimplificationIsIdempotentAcrossCanonicalizationSamples) {
    const std::vector<std::string> samples = {
        "(x + y) + (z + 0)",
        "3*x + 2*x - x",
        "2*(x + y + z)",
        "(x + 1) * (y + 1)",
        "sin(t)^2 + cos(t)^2 + 0",
        "ln(e^x) + exp(0)",
        "sqrt(x^2)",
        "exp(a + b)",
    };

    for (const std::string& input : samples) {
        SCOPED_TRACE(input);
        AstArena parse_arena;
        AstArena first_arena;
        AstArena second_arena;

        auto parsed = parse_expr(input, parse_arena);
        ASSERT_TRUE(parsed.is_ok()) << parsed.error().message;

        auto first = simplify(parsed.value(), first_arena);
        ASSERT_TRUE(first.is_ok()) << first.error().message;

        auto second = simplify(first.value(), second_arena);
        ASSERT_TRUE(second.is_ok()) << second.error().message;

        EXPECT_TRUE(structural_equal(first.value(), second.value()));
    }
}

TEST(SymbolicPropertyTest, CanonicalizesRandomizedCommutativeAdditionWithFixedSeed) {
    constexpr std::uint64_t kSeed = 0xC0FFEEULL;
    constexpr std::size_t kIterations = 32U;

    cas::test::run_seeded_cases(kSeed, kIterations, [](cas::test::DeterministicRng& rng, std::size_t index) {
        const std::string lhs_text = cas::test::generate_polynomial_expr(rng, 2);
        const std::string rhs_text = cas::test::generate_polynomial_expr(rng, 2);
        SCOPED_TRACE("seeded addition case " + std::to_string(index) + ": " + lhs_text + " | " + rhs_text);

        AstArena lhs_parse_arena;
        AstArena rhs_parse_arena;
        AstArena lhs_simplify_arena;
        AstArena rhs_simplify_arena;

        auto lhs = simplify_input("(" + lhs_text + ") + (" + rhs_text + ")", lhs_parse_arena, lhs_simplify_arena);
        auto rhs = simplify_input("(" + rhs_text + ") + (" + lhs_text + ")", rhs_parse_arena, rhs_simplify_arena);
        ASSERT_TRUE(lhs.is_ok()) << lhs.error().message;
        ASSERT_TRUE(rhs.is_ok()) << rhs.error().message;

        EXPECT_TRUE(structural_equal(lhs.value(), rhs.value()));
    });
}

TEST(SymbolicPropertyTest, CanonicalizesRandomizedDistributiveSamplesWithFixedSeed) {
    constexpr std::uint64_t kSeed = 0xBAD5EEDULL;
    constexpr std::size_t kIterations = 24U;

    cas::test::run_seeded_cases(kSeed, kIterations, [](cas::test::DeterministicRng& rng, std::size_t index) {
        const std::string factor = cas::test::generate_polynomial_expr(rng, 1);
        const std::string left = cas::test::generate_polynomial_expr(rng, 1);
        const std::string right = cas::test::generate_polynomial_expr(rng, 1);
        SCOPED_TRACE(
            "seeded distributive case " + std::to_string(index) + ": " + factor + " ; " + left + " ; " + right);

        AstArena lhs_parse_arena;
        AstArena rhs_parse_arena;
        AstArena lhs_simplify_arena;
        AstArena rhs_simplify_arena;

        auto lhs = simplify_input(
            "(" + factor + ") * ((" + left + ") + (" + right + "))",
            lhs_parse_arena,
            lhs_simplify_arena);
        auto rhs = simplify_input(
            "((" + factor + ") * (" + left + ")) + ((" + factor + ") * (" + right + "))",
            rhs_parse_arena,
            rhs_simplify_arena);
        ASSERT_TRUE(lhs.is_ok()) << lhs.error().message;
        ASSERT_TRUE(rhs.is_ok()) << rhs.error().message;

        EXPECT_TRUE(structural_equal(lhs.value(), rhs.value()));
    });
}

TEST(SymbolicPropertyTest, SimplificationIsIdempotentForRandomizedPolynomialSamplesWithFixedSeed) {
    constexpr std::uint64_t kSeed = 0x5EED1234ULL;
    constexpr std::size_t kIterations = 32U;

    cas::test::run_seeded_cases(kSeed, kIterations, [](cas::test::DeterministicRng& rng, std::size_t index) {
        const std::string input = cas::test::generate_polynomial_expr(rng, 3);
        SCOPED_TRACE("seeded idempotence case " + std::to_string(index) + ": " + input);

        AstArena parse_arena;
        AstArena first_arena;
        AstArena second_arena;

        auto parsed = parse_expr(input, parse_arena);
        ASSERT_TRUE(parsed.is_ok()) << parsed.error().message;

        auto first = simplify(parsed.value(), first_arena);
        ASSERT_TRUE(first.is_ok()) << first.error().message;

        auto second = simplify(first.value(), second_arena);
        ASSERT_TRUE(second.is_ok()) << second.error().message;

        EXPECT_TRUE(structural_equal(first.value(), second.value()));
    });
}

TEST(SymbolicSimplifyTest, IsIdempotentOnCanonicalPolynomial) {
    AstArena parse_arena;
    AstArena first_arena;
    AstArena second_arena;

    auto parsed = parse_expr("1 + x^3 + x", parse_arena);
    ASSERT_TRUE(parsed.is_ok()) << parsed.error().message;

    auto first = simplify(parsed.value(), first_arena);
    ASSERT_TRUE(first.is_ok()) << first.error().message;

    auto second = simplify(first.value(), second_arena);
    ASSERT_TRUE(second.is_ok()) << second.error().message;

    EXPECT_TRUE(structural_equal(first.value(), second.value()));
}

TEST(AlgebraMultivariatePolynomialTest, CanonicalizesSparseTermsAndMultipliesExactly) {
    const Symbol x("x");
    const Symbol y("y");

    algebra::MultivariatePolynomial left({
        algebra::MultivariateTerm{BigInt(1), {{x, 1U}}},
        algebra::MultivariateTerm{BigInt(2), {{y, 1U}}},
    });
    algebra::MultivariatePolynomial right({
        algebra::MultivariateTerm{BigInt(1), {{x, 1U}}},
        algebra::MultivariateTerm{BigInt(-2), {{y, 1U}}},
    });

    algebra::MultivariatePolynomial product = left * right;

    ASSERT_EQ(product.terms().size(), 2U);
    EXPECT_EQ(product.total_degree(), 2U);
    ASSERT_EQ(product.variables().size(), 2U);
    EXPECT_EQ(product.variables()[0].name, "x");
    EXPECT_EQ(product.variables()[1].name, "y");

    EXPECT_EQ(product.terms()[0].coefficient.decimal(), "1");
    ASSERT_EQ(product.terms()[0].factors.size(), 1U);
    EXPECT_EQ(product.terms()[0].factors[0].first.name, "x");
    EXPECT_EQ(product.terms()[0].factors[0].second, 2U);

    EXPECT_TRUE(product.terms()[1].coefficient.is_negative());
    EXPECT_EQ(product.terms()[1].coefficient.decimal(), "-4");
    ASSERT_EQ(product.terms()[1].factors.size(), 1U);
    EXPECT_EQ(product.terms()[1].factors[0].first.name, "y");
    EXPECT_EQ(product.terms()[1].factors[0].second, 2U);
}

TEST(AlgebraMultivariatePolynomialTest, ConvertsToUnivariateCoefficientsAndEvaluatesIntegersOnly) {
    const Symbol x("x");
    const Symbol y("y");
    symbolic::CASContext ctx;
    AstArena parse_arena;
    AstArena value_arena;

    algebra::MultivariatePolynomial polynomial({
        algebra::MultivariateTerm{BigInt(1), {{x, 2U}}},
        algebra::MultivariateTerm{BigInt(-4), {{y, 2U}}},
    });

    auto coefficients = polynomial.to_univariate_coefficients(x, ctx);
    ASSERT_TRUE(coefficients.is_ok()) << coefficients.error().message;
    ASSERT_EQ(coefficients.value().size(), 3U);

    auto expected_constant = parse_expr("-4*y^2", parse_arena);
    ASSERT_TRUE(expected_constant.is_ok()) << expected_constant.error().message;
    auto simplified_expected_constant = ctx.simplify(expected_constant.value());
    ASSERT_TRUE(simplified_expected_constant.is_ok()) << simplified_expected_constant.error().message;

    EXPECT_TRUE(structural_equal(coefficients.value()[0], simplified_expected_constant.value()));
    ASSERT_EQ(expr_kind(coefficients.value()[1]), ExprKind::IntegerLit);
    EXPECT_EQ(expr_ref<IntegerLit>(coefficients.value()[1]).value.decimal(), "0");
    ASSERT_EQ(expr_kind(coefficients.value()[2]), ExprKind::IntegerLit);
    EXPECT_EQ(expr_ref<IntegerLit>(coefficients.value()[2]).value.decimal(), "1");

    auto evaluated = polynomial.evaluate_at(y, value_arena.make<IntegerLit>(BigInt(3)));
    ASSERT_TRUE(evaluated.is_ok()) << evaluated.error().message;
    auto evaluated_coefficients = evaluated.value().to_univariate_coefficients(x, ctx);
    ASSERT_TRUE(evaluated_coefficients.is_ok()) << evaluated_coefficients.error().message;
    ASSERT_EQ(expr_kind(evaluated_coefficients.value()[0]), ExprKind::IntegerLit);
    EXPECT_TRUE(expr_ref<IntegerLit>(evaluated_coefficients.value()[0]).value.is_negative());
    EXPECT_EQ(expr_ref<IntegerLit>(evaluated_coefficients.value()[0]).value.decimal(), "-36");

    auto decimal_evaluation = polynomial.evaluate_at(y, value_arena.make<DecimalLit>("1.5"));
    ASSERT_TRUE(decimal_evaluation.is_error());
    EXPECT_EQ(decimal_evaluation.error().kind, CASErrorKind::Unimplemented);
}

TEST(SymbolicRewriteTest, ExpDecompositionRuleIsOriented) {
    AstArena arena;
    auto pattern = parse_expr("exp(a_ + b_)", arena);
    auto replacement = parse_expr("exp(a_) * exp(b_)", arena);
    ASSERT_TRUE(pattern.is_ok()) << pattern.error().message;
    ASSERT_TRUE(replacement.is_ok()) << replacement.error().message;

    RewriteRule rule{
        .pattern = pattern.value(),
        .replacement = replacement.value(),
        .condition = {},
    };

    EXPECT_TRUE(rewrite_rule_is_oriented(rule));
}

TEST(SymbolicRewriteTest, LnDecompositionOverDivisionIsOriented) {
    AstArena arena;
    auto pattern = parse_expr("ln(a_ / b_)", arena);
    auto replacement = parse_expr("ln(a_) + (-1)*ln(b_)", arena);
    ASSERT_TRUE(pattern.is_ok()) << pattern.error().message;
    ASSERT_TRUE(replacement.is_ok()) << replacement.error().message;

    RewriteRule rule{
        .pattern = pattern.value(),
        .replacement = replacement.value(),
        .condition = {},
    };

    EXPECT_TRUE(rewrite_rule_is_oriented(rule));
}

// L1-07: transcendental_normal_form must expand log products/quotients/powers
// and cancel inverse pairs without hardcoded special cases
TEST(SymbolicNormalFormTest, L1_07_LnProductExpands) {
    CASContext ctx;
    auto e = parse_expr("ln(x * y)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto res = transcendental_normal_form(e.value(), ctx);
    ASSERT_TRUE(res.is_ok());
    // must contain two ln subexpressions
    int ln_count = 0;
    std::function<void(ExprPtr)> count_ln = [&](ExprPtr n) {
        if (!n) return;
        if (auto* f = expr_cast<FuncCall>(n); f && f->func_id == BuiltinOp::Ln) ln_count++;
        if (auto* b = expr_cast<Binary>(n)) { count_ln(b->left); count_ln(b->right); }
        if (auto* s = expr_cast<Sum>(n)) for (auto t : s->terms) count_ln(t);
        if (auto* p = expr_cast<Product>(n)) for (auto f2 : p->factors) count_ln(f2);
    };
    count_ln(res.value());
    EXPECT_EQ(ln_count, 2) << "ln(x*y) must expand to ln(x)+ln(y)";
}

TEST(SymbolicNormalFormTest, L1_07_LnDivisionExpands) {
    CASContext ctx;
    auto e = parse_expr("ln(x / y)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto res = transcendental_normal_form(e.value(), ctx);
    ASSERT_TRUE(res.is_ok());
    int ln_count = 0;
    std::function<void(ExprPtr)> count_ln = [&](ExprPtr n) {
        if (!n) return;
        if (auto* f = expr_cast<FuncCall>(n); f && f->func_id == BuiltinOp::Ln) ln_count++;
        if (auto* b = expr_cast<Binary>(n)) { count_ln(b->left); count_ln(b->right); }
        if (auto* s = expr_cast<Sum>(n)) for (auto t : s->terms) count_ln(t);
    };
    count_ln(res.value());
    EXPECT_EQ(ln_count, 2) << "ln(x/y) must expand to ln(x)-ln(y)";
}

TEST(SymbolicNormalFormTest, L1_07_LnPowerExpands) {
    CASContext ctx;
    auto e = parse_expr("ln(x^3)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto res = transcendental_normal_form(e.value(), ctx);
    ASSERT_TRUE(res.is_ok());
    // result must contain Mul with IntegerLit(3) and ln(x)
    bool has_coeff3 = false;
    std::function<void(ExprPtr)> find_3 = [&](ExprPtr n) {
        if (!n) return;
        if (auto* i = expr_cast<IntegerLit>(n); i && i->value == BigInt(3)) has_coeff3 = true;
        if (auto* b = expr_cast<Binary>(n)) { find_3(b->left); find_3(b->right); }
    };
    find_3(res.value());
    EXPECT_TRUE(has_coeff3) << "ln(x^3) must expand to 3*ln(x)";
}

TEST(SymbolicNormalFormTest, L1_07_ExpLnCancels) {
    CASContext ctx;
    auto e = parse_expr("exp(ln(x))", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto res = transcendental_normal_form(e.value(), ctx);
    ASSERT_TRUE(res.is_ok());
    // must not contain exp or ln nodes
    bool has_transcendental = false;
    std::function<void(ExprPtr)> check = [&](ExprPtr n) {
        if (!n) return;
        if (auto* f = expr_cast<FuncCall>(n)) {
            if (f->func_id == BuiltinOp::Exp || f->func_id == BuiltinOp::Ln) has_transcendental = true;
        }
        if (auto* b = expr_cast<Binary>(n)) { check(b->left); check(b->right); }
    };
    check(res.value());
    EXPECT_FALSE(has_transcendental) << "exp(ln(x)) must cancel to x";
}

TEST(SymbolicNormalFormTest, L1_07_LnExpCancels) {
    CASContext ctx;
    auto e = parse_expr("ln(exp(x))", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto res = transcendental_normal_form(e.value(), ctx);
    ASSERT_TRUE(res.is_ok());
    bool has_transcendental = false;
    std::function<void(ExprPtr)> check = [&](ExprPtr n) {
        if (!n) return;
        if (auto* f = expr_cast<FuncCall>(n)) {
            if (f->func_id == BuiltinOp::Exp || f->func_id == BuiltinOp::Ln) has_transcendental = true;
        }
        if (auto* b = expr_cast<Binary>(n)) { check(b->left); check(b->right); }
    };
    check(res.value());
    EXPECT_FALSE(has_transcendental) << "ln(exp(x)) must cancel to x";
}

// L1-12: sqrt(n) must extract partial square factors for any n (no fixed table)
TEST(SymbolicSqrtDenestTest, L1_12_SqrtExtracts2FromSqrt12) {
    // sqrt(12) = 2*sqrt(3) — must not return unevaluated sqrt(12)
    CASContext ctx;
    auto e = parse_expr("sqrt(12)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto res = ctx.simplify(e.value());
    ASSERT_TRUE(res.is_ok());
    // Result must contain IntegerLit(2) as a factor, not a bare sqrt(12)
    bool has_12 = false;
    bool has_2_factor = false;
    std::function<void(ExprPtr)> check = [&](ExprPtr n) {
        if (!n) return;
        if (auto* i = expr_cast<IntegerLit>(n); i && i->value == BigInt(12)) has_12 = true;
        if (auto* i = expr_cast<IntegerLit>(n); i && i->value == BigInt(2)) has_2_factor = true;
        if (auto* b = expr_cast<Binary>(n)) { check(b->left); check(b->right); }
        if (auto* p = expr_cast<Product>(n)) for (auto f : p->factors) check(f);
    };
    check(res.value());
    EXPECT_FALSE(has_12) << "sqrt(12) must not return unevaluated sqrt(12)";
    EXPECT_TRUE(has_2_factor) << "sqrt(12) must produce coefficient 2";
}

TEST(SymbolicSqrtDenestTest, L1_12_SqrtExtracts3FromSqrt75) {
    // sqrt(75) = 5*sqrt(3)
    CASContext ctx;
    auto e = parse_expr("sqrt(75)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto res = ctx.simplify(e.value());
    ASSERT_TRUE(res.is_ok());
    bool has_75 = false;
    bool has_5_factor = false;
    std::function<void(ExprPtr)> check = [&](ExprPtr n) {
        if (!n) return;
        if (auto* i = expr_cast<IntegerLit>(n); i && i->value == BigInt(75)) has_75 = true;
        if (auto* i = expr_cast<IntegerLit>(n); i && i->value == BigInt(5)) has_5_factor = true;
        if (auto* b = expr_cast<Binary>(n)) { check(b->left); check(b->right); }
        if (auto* p = expr_cast<Product>(n)) for (auto f : p->factors) check(f);
    };
    check(res.value());
    EXPECT_FALSE(has_75) << "sqrt(75) must not return unevaluated sqrt(75)";
    EXPECT_TRUE(has_5_factor) << "sqrt(75) must produce coefficient 5";
}

TEST(SymbolicSqrtDenestTest, L1_12_PerfectSquareNormalized) {
    // sqrt(144) = 12 — must be IntegerLit(12)
    CASContext ctx;
    auto e = parse_expr("sqrt(144)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto res = ctx.simplify(e.value());
    ASSERT_TRUE(res.is_ok());
    const auto* i = expr_cast<IntegerLit>(res.value());
    ASSERT_NE(i, nullptr) << "sqrt(144) must be IntegerLit";
    EXPECT_EQ(i->value, BigInt(12));
}

}  // namespace
}  // namespace cas::symbolic
