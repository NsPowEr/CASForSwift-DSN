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

void expect_equivalent(ExprPtr actual, ExprPtr expected) {
    symbolic::CASContext compare_context;
    auto normalized_actual = compare_context.simplify(actual);
    ASSERT_TRUE(normalized_actual.is_ok()) << normalized_actual.error().message;
    auto normalized_expected = compare_context.simplify(expected);
    ASSERT_TRUE(normalized_expected.is_ok()) << normalized_expected.error().message;
    auto equivalent = mathematically_equal(normalized_actual.value(), normalized_expected.value(), compare_context);
    ASSERT_TRUE(equivalent.is_ok()) << equivalent.error().message;
    EXPECT_TRUE(equivalent.value()) << "actual=" << debug_print(normalized_actual.value()) << " expected=" << debug_print(normalized_expected.value());
}

void expect_integration_oracle(const std::string& integrand_text, const std::string& variable) {
    symbolic::CASContext integrate_context;
    AstArena parse_arena;
    auto expr = parse_expr(integrand_text, parse_arena);
    ASSERT_TRUE(expr.is_ok());

    auto primitive = integrate(expr.value(), Symbol(variable), integrate_context);
    ASSERT_TRUE(primitive.is_ok()) << primitive.error().message;

    symbolic::CASContext differentiate_context;
    auto recovered = diff(primitive.value(), Symbol(variable), 1U, differentiate_context);
    ASSERT_TRUE(recovered.is_ok()) << recovered.error().message;

    expect_equivalent(recovered.value(), expr.value());
}

TEST(IntegrationAdvancedTest, TranscendentalExponentialPolynomial) {
    expect_integration_oracle("x*exp(x)", "x");
    expect_integration_oracle("(x^2 + 2*x + 1)*exp(x)", "x");
}

TEST(IntegrationAdvancedTest, TranscendentalLogarithmicPolynomial) {
    expect_integration_oracle("ln(x)", "x");
}

// New 2026-05-15: Risch logarithmic polynomial part (general k & coefficients).
// These exercise the descending recursion B_k = ∫ [a_k - (k+1)*B_{k+1}*u'/u] dx.
TEST(IntegrationAdvancedTest, LogarithmicPolynomialDegree1WithNonTrivialCoefficient) {
    // ∫ x * ln(x) dx  =  x^2/2 * ln(x) - x^2/4
    expect_integration_oracle("x*ln(x)", "x");
}

TEST(IntegrationAdvancedTest, LogarithmicPolynomialDegree2) {
    // ∫ ln(x)^2 dx  =  x*ln(x)^2 - 2*x*ln(x) + 2*x
    expect_integration_oracle("ln(x)^2", "x");
}

TEST(IntegrationAdvancedTest, LogarithmicPolynomialMixed) {
    // ∫ (2x + 3) * ln(x) dx  via Risch log poly part
    expect_integration_oracle("(2*x + 3)*ln(x)", "x");
}

TEST(IntegrationAdvancedTest, LogarithmicPolynomialQuadraticCoefficient) {
    // ∫ x^2 * ln(x) dx  =  x^3/3 * ln(x) - x^3/9
    expect_integration_oracle("x^2*ln(x)", "x");
}

} // namespace
} // namespace cas::calculus
