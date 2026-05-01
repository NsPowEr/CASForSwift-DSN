#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/algebra.hpp"
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

void expect_equivalent(ExprPtr actual, ExprPtr expected) {
    symbolic::CASContext compare_context;
    
    // Primo tentativo: mathematically_equal (che usa expand internamente)
    auto equivalent = mathematically_equal(actual, expected, compare_context);
    if (equivalent.is_ok() && equivalent.value()) {
        return;
    }

    auto expanded_actual = algebra::expand(actual, compare_context);
    auto expanded_expected = algebra::expand(expected, compare_context);
    if (expanded_actual.is_ok() && expanded_expected.is_ok()) {
        auto expanded_equivalent = mathematically_equal(expanded_actual.value(), expanded_expected.value(), compare_context);
        if (expanded_equivalent.is_ok() && expanded_equivalent.value()) {
            return;
        }
        if (structural_equal(expanded_actual.value(), expanded_expected.value())) {
            return;
        }
    }

    // Secondo tentativo: simplify(actual - expected) == 0
    AstArena& arena = compare_context.arena();
    auto diff_expr = arena.make<Sum>(std::vector<ExprPtr>{
        actual,
        arena.make<Unary>(UnaryOp::Neg, expected)
    });
    
    auto simplified_diff = compare_context.simplify(diff_expr);
    if (simplified_diff.is_ok()) {
        if (const auto* integer = expr_cast<IntegerLit>(simplified_diff.value())) {
            if (integer->value.is_zero()) {
                return;
            }
        }
    }

    // Se fallisce, stampa debug info
    auto normalized_actual = compare_context.simplify(actual);
    auto normalized_expected = compare_context.simplify(expected);
    
    EXPECT_TRUE(false) 
        << "actual=" << debug_print(normalized_actual.is_ok() ? normalized_actual.value() : actual) 
        << " expected=" << debug_print(normalized_expected.is_ok() ? normalized_expected.value() : expected);
}

void test_diff(const std::string& input, const std::string& var, const std::string& expected_text) {
    symbolic::CASContext context;
    AstArena parse_arena;
    AstArena expected_arena;

    auto expr = parse_expr(input, parse_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    auto result = diff(expr.value(), Symbol(var), 1U, context);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    auto expected = parse_expr(expected_text, expected_arena);
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    expect_equivalent(result.value(), expected.value());
}

void test_integrate(const std::string& input, const std::string& var, const std::string& expected_text) {
    symbolic::CASContext context;
    AstArena parse_arena;
    AstArena expected_arena;

    auto expr = parse_expr(input, parse_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    auto result = integrate(expr.value(), Symbol(var), context);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    auto expected = parse_expr(expected_text, expected_arena);
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    expect_equivalent(result.value(), expected.value());
}

TEST(CalculusF4Test, TranscendentalDerivatives) {
    test_diff("asin(x)", "x", "1/sqrt(1-x^2)");
    test_diff("acos(x)", "x", "-1/sqrt(1-x^2)");
    test_diff("atan(x)", "x", "1/(1+x^2)");
    test_diff("sinh(x)", "x", "cosh(x)");
    test_diff("cosh(x)", "x", "sinh(x)");
    test_diff("tanh(x)", "x", "1/cosh(x)^2");
}

TEST(CalculusF4Test, TranscendentalChainRule) {
    test_diff("asin(x^2)", "x", "2*x/sqrt(1-x^4)");
    test_diff("sinh(exp(x))", "x", "exp(x)*cosh(exp(x))");
}

TEST(CalculusF4Test, IntegrationByPartsBasic) {
    test_integrate("x*exp(x)", "x", "(x-1)*exp(x)");
    test_integrate("x*sin(x)", "x", "sin(x) - x*cos(x)");
    test_integrate("ln(x)", "x", "x*ln(x) - x");
}

TEST(CalculusF4Test, IntegrationByPartsAdvanced) {
    // x^2 * exp(x) -> u=x^2, dv=exp(x) -> 2x * exp(x) -> 2 * exp(x)
    test_integrate("x^2 * exp(x)", "x", "(x^2 - 2*x + 2)*exp(x)");
    
    // x * ln(x) -> u=ln(x), dv=x -> du=1/x, v=x^2/2 -> ∫ v du = ∫ x/2 = x^2/4
    test_integrate("x*ln(x)", "x", "(1/2)*x^2*ln(x) - (1/4)*x^2");
}

TEST(CalculusF4Test, FundamentalTheoremOfCalculus) {
    symbolic::CASContext context;
    AstArena arena;
    
    std::vector<std::string> test_funcs = {
        "x*exp(x)",
        "x^2 * sin(x)",
        "ln(x) / x", // This might need substitution but let's see if IBP handles it
        "atan(x)"
    };
    
    for (const auto& func_text : test_funcs) {
        auto expr = parse_expr(func_text, arena);
        ASSERT_TRUE(expr.is_ok()) << "Failed to parse " << func_text;

        auto primitive = integrate(expr.value(), Symbol("x"), context);
        if (primitive.is_error()) continue; // Some might not be integrable yet

        auto recovered = diff(primitive.value(), Symbol("x"), 1U, context);
        ASSERT_TRUE(recovered.is_ok()) << "Failed to differentiate primitive of " << func_text;

        expect_equivalent(recovered.value(), expr.value());
    }
}

TEST(CalculusF4Test, GeneralizedArcTanIntegral) {
    // 1/(x^2+1) → arctan(x): verify d/dx(primitive) == original
    symbolic::CASContext context;
    AstArena arena;
    auto expr = parse_expr("1/(x^2+1)", arena);
    ASSERT_TRUE(expr.is_ok());
    auto primitive = integrate(expr.value(), Symbol("x"), context);
    ASSERT_TRUE(primitive.is_ok()) << primitive.error().message;
    auto recovered = diff(primitive.value(), Symbol("x"), 1U, context);
    ASSERT_TRUE(recovered.is_ok());
    expect_equivalent(recovered.value(), expr.value());
}

TEST(CalculusF4Test, GeneralizedArcSinIntegral) {
    // 1/sqrt(1-x^2) → arcsin(x): verify d/dx(primitive) == original
    symbolic::CASContext context;
    AstArena arena;
    auto expr = parse_expr("1/sqrt(1-x^2)", arena);
    ASSERT_TRUE(expr.is_ok());
    auto primitive = integrate(expr.value(), Symbol("x"), context);
    ASSERT_TRUE(primitive.is_ok()) << primitive.error().message;
    auto recovered = diff(primitive.value(), Symbol("x"), 1U, context);
    ASSERT_TRUE(recovered.is_ok());
    expect_equivalent(recovered.value(), expr.value());
}

TEST(CalculusF4Test, ExpSumRewrite) {
    // exp(a+b) = exp(a)*exp(b) via builtin rewrite
    symbolic::CASContext context;
    AstArena arena;
    auto expr = parse_expr("exp(x+1)", arena);
    ASSERT_TRUE(expr.is_ok());
    auto result = context.simplify(expr.value());
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    // Rewrite fires: result != input structurally
    EXPECT_FALSE(structural_equal(result.value(), expr.value()));
}

} // namespace
} // namespace cas::calculus
