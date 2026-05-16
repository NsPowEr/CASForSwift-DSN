#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include <gtest/gtest.h>

namespace cas::calculus {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

void expect_positive_infinity(ExprPtr expr) {
    const auto* constant = expr_cast<Constant>(expr);
    ASSERT_NE(constant, nullptr) << "Expected +Infinity, got: " << debug_print(expr);
    EXPECT_EQ(constant->value, MathConstant::Infinity);
}

void expect_equivalent(ExprPtr actual, ExprPtr expected) {
    symbolic::CASContext compare_context;
    auto normalized_actual = compare_context.simplify(actual);
    ASSERT_TRUE(normalized_actual.is_ok()) << normalized_actual.error().message;
    auto normalized_expected = compare_context.simplify(expected);
    ASSERT_TRUE(normalized_expected.is_ok()) << normalized_expected.error().message;

    auto equivalent = mathematically_equal(
        normalized_actual.value(),
        normalized_expected.value(),
        compare_context);
    ASSERT_TRUE(equivalent.is_ok()) << equivalent.error().message;
    EXPECT_TRUE(equivalent.value())
        << "actual=" << debug_print(normalized_actual.value())
        << " expected=" << debug_print(normalized_expected.value());
}

TEST(GruntzTest, ExponentialTowerReferenceCase) {
    symbolic::CASContext context;
    
    // lim(x->inf) exp(exp(x)) / exp(x^1000) = inf
    auto expr = parse_expr("exp(exp(x)) / exp(x^1000)", context.arena());
    ASSERT_TRUE(expr.is_ok());
    
    auto point = context.arena().make<Constant>(MathConstant::Infinity);
    auto res = limit(expr.value(), Symbol("x"), point, LimitDirection::Both, context);
    
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_positive_infinity(res.value());
}

TEST(GruntzTest, AntiHardcodeDifferentVariableReferenceCase) {
    symbolic::CASContext context;
    auto expr = parse_expr("exp(exp(k)) / exp(k^1000)", context.arena());
    ASSERT_TRUE(expr.is_ok());

    auto point = context.arena().make<Constant>(MathConstant::Infinity);
    auto res = limit(expr.value(), Symbol("k"), point, LimitDirection::Both, context);

    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_positive_infinity(res.value());
}

TEST(GruntzTest, AntiHardcodeEquivalentSyntaxReferenceCase) {
    symbolic::CASContext context;
    auto expr = parse_expr("exp(exp(kappa))*exp(-(kappa^1000))", context.arena());
    ASSERT_TRUE(expr.is_ok());

    auto point = context.arena().make<Constant>(MathConstant::Infinity);
    auto res = limit(expr.value(), Symbol("kappa"), point, LimitDirection::Both, context);

    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_positive_infinity(res.value());
}

TEST(GruntzTest, NegativeInfinityNormalizationForReferenceCase) {
    symbolic::CASContext context;
    auto expr = parse_expr("exp(exp(-u)) / exp((-u)^1000)", context.arena());
    ASSERT_TRUE(expr.is_ok());

    auto point = context.arena().make<Unary>(
        UnaryOp::Neg,
        context.arena().make<Constant>(MathConstant::Infinity));
    auto res = limit(expr.value(), Symbol("u"), point, LimitDirection::Both, context);

    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_positive_infinity(res.value());
}

TEST(GruntzTest, MrvSeriesExtractionForMonomialClass) {
    symbolic::CASContext context;
    auto expr = parse_expr("(exp(exp(t))^2) / exp(exp(t))", context.arena());
    ASSERT_TRUE(expr.is_ok());

    auto point = context.arena().make<Constant>(MathConstant::Infinity);
    auto res = limit(expr.value(), Symbol("t"), point, LimitDirection::Both, context);

    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_positive_infinity(res.value());
}

TEST(GruntzTest, ResolvesSameOrderMrvTermsWithExactNonzeroLeadingCoefficient) {
    symbolic::CASContext context;
    AstArena expected_arena;
    auto expr = parse_expr("(3*exp(y)+5)/(2*exp(y)-7)", context.arena());
    auto expected = parse_expr("3/2", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto point = context.arena().make<Constant>(MathConstant::Infinity);
    auto res = limit(expr.value(), Symbol("y"), point, LimitDirection::Both, context);

    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_equivalent(res.value(), expected.value());
}

TEST(GruntzTest, KeepsCancelledSameOrderMrvTermsConservative) {
    symbolic::CASContext context;
    AstArena expected_arena;
    auto expr = parse_expr("(exp(z)-exp(z)+1)/exp(z)", context.arena());
    auto expected = parse_expr("0", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto point = context.arena().make<Constant>(MathConstant::Infinity);
    auto res = limit(expr.value(), Symbol("z"), point, LimitDirection::Both, context);

    if (res.is_ok()) {
        expect_equivalent(res.value(), expected.value());
    } else {
        EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
    }
}

// L1-01 dynamic-rank anti-hardcode: triple exponential tower must dominate
// the double tower even when the inner argument is a high-degree polynomial.
// Previously rank was capped at 3 (Cat 10), making all exp(...) indistinguishable.
TEST(GruntzTest, TripleExponentialTowerDominatesDoubleTower_AntiHardcode) {
    symbolic::CASContext context;
    auto expr = parse_expr("exp(exp(exp(x))) / exp(exp(x^100))", context.arena());
    ASSERT_TRUE(expr.is_ok());
    auto point = context.arena().make<Constant>(MathConstant::Infinity);
    auto res = limit(expr.value(), Symbol("x"), point, LimitDirection::Both, context);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_positive_infinity(res.value());
}

TEST(GruntzTest, DoubleTowerDominatedByTripleTowerReciprocal_AntiHardcode) {
    symbolic::CASContext context;
    // The reciprocal must go to 0.
    auto expr = parse_expr("exp(exp(x^100)) / exp(exp(exp(x)))", context.arena());
    ASSERT_TRUE(expr.is_ok());
    auto point = context.arena().make<Constant>(MathConstant::Infinity);
    auto res = limit(expr.value(), Symbol("x"), point, LimitDirection::Both, context);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto zero = context.arena().make<IntegerLit>(BigInt(0));
    expect_equivalent(res.value(), zero);
}

}
}
