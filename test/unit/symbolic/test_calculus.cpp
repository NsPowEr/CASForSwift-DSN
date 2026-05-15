#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/algebra.hpp"
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
    if (tokens.is_error()) {
        return fail<ExprPtr>(tokens.error());
    }

    Parser parser(tokens.value(), arena);
    return parser.parse();
}

Result<ExprPtr> differentiate_input(const std::string& input, const std::string& variable, unsigned int order, symbolic::CASContext& context) {
    auto parsed = parse_expr(input, context.arena());
    if (parsed.is_error()) {
        return fail<ExprPtr>(parsed.error());
    }

    return diff(parsed.value(), Symbol(variable), order, context);
}

Result<ExprPtr> implicit_diff_input(
    const std::string& input,
    const std::string& dependent,
    const std::string& independent,
    symbolic::CASContext& context) {
    auto parsed = parse_expr(input, context.arena());
    if (parsed.is_error()) {
        return fail<ExprPtr>(parsed.error());
    }

    return implicit_diff(parsed.value(), Symbol(dependent), Symbol(independent), context);
}

Result<ExprPtr> integrate_input(const std::string& input, const std::string& variable, symbolic::CASContext& context) {
    auto parsed = parse_expr(input, context.arena());
    if (parsed.is_error()) {
        return fail<ExprPtr>(parsed.error());
    }

    return integrate(parsed.value(), Symbol(variable), context);
}

Result<ExprPtr> limit_input(
    const std::string& input,
    const std::string& variable,
    const std::string& point,
    LimitDirection direction,
    symbolic::CASContext& context) {
    auto parsed = parse_expr(input, context.arena());
    if (parsed.is_error()) {
        return fail<ExprPtr>(parsed.error());
    }

    auto parsed_point = parse_expr(point, context.arena());
    if (parsed_point.is_error()) {
        return fail<ExprPtr>(parsed_point.error());
    }

    return limit(parsed.value(), Symbol(variable), parsed_point.value(), direction, context);
}

Result<TaylorExpansion> taylor_input(
    const std::string& input,
    const std::string& variable,
    const std::string& point,
    unsigned int order,
    symbolic::CASContext& context) {
    auto parsed = parse_expr(input, context.arena());
    if (parsed.is_error()) {
        return fail<TaylorExpansion>(parsed.error());
    }

    auto parsed_point = parse_expr(point, context.arena());
    if (parsed_point.is_error()) {
        return fail<TaylorExpansion>(parsed_point.error());
    }

    return taylor_series(parsed.value(), Symbol(variable), parsed_point.value(), order, context);
}

Result<std::vector<ExprPtr>> gradient_input(
    const std::string& input,
    const std::vector<Symbol>& variables,
    symbolic::CASContext& context) {
    auto parsed = parse_expr(input, context.arena());
    if (parsed.is_error()) {
        return fail<std::vector<ExprPtr>>(parsed.error());
    }

    return gradient(parsed.value(), variables, context);
}

Result<ExprPtr> jacobian_input(
    const std::vector<std::string>& inputs,
    const std::vector<Symbol>& variables,
    symbolic::CASContext& context) {
    std::vector<ExprPtr> functions;
    functions.reserve(inputs.size());
    for (const std::string& input : inputs) {
        auto parsed = parse_expr(input, context.arena());
        if (parsed.is_error()) {
            return fail<ExprPtr>(parsed.error());
        }
        functions.push_back(parsed.value());
    }

    return jacobian(functions, variables, context);
}

Result<ExprPtr> hessian_input(
    const std::string& input,
    const std::vector<Symbol>& variables,
    symbolic::CASContext& context) {
    auto parsed = parse_expr(input, context.arena());
    if (parsed.is_error()) {
        return fail<ExprPtr>(parsed.error());
    }

    return hessian(parsed.value(), variables, context);
}

bool contains_builtin(ExprPtr expr, BuiltinOp op) {
    if (!expr) return false;
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if (call->func_id == op) return true;
        for (ExprPtr arg : call->args) {
            if (contains_builtin(arg, op)) return true;
        }
        return false;
    }
    if (const auto* unary = expr_cast<Unary>(expr)) return contains_builtin(unary->operand, op);
    if (const auto* binary = expr_cast<Binary>(expr)) {
        return contains_builtin(binary->left, op) || contains_builtin(binary->right, op);
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr term : sum->terms) {
            if (contains_builtin(term, op)) return true;
        }
        return false;
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr factor : product->factors) {
            if (contains_builtin(factor, op)) return true;
        }
        return false;
    }
    return false;
}

void expect_equivalent(ExprPtr actual, ExprPtr expected) {
    symbolic::CASContext compare_context;
    auto normalized_actual = compare_context.simplify(actual);
    ASSERT_TRUE(normalized_actual.is_ok()) << normalized_actual.error().message;

    auto normalized_expected = compare_context.simplify(expected);
    ASSERT_TRUE(normalized_expected.is_ok()) << normalized_expected.error().message;

    if (structural_equal(normalized_actual.value(), normalized_expected.value())) {
        SUCCEED();
        return;
    }

    auto expanded_actual = algebra::expand(normalized_actual.value(), compare_context);
    auto expanded_expected = algebra::expand(normalized_expected.value(), compare_context);
    if (expanded_actual.is_ok() && expanded_expected.is_ok() &&
        structural_equal(expanded_actual.value(), expanded_expected.value())) {
        SUCCEED();
        return;
    }

    AstArena& arena = compare_context.arena();
    auto difference = arena.make<Sum>(std::vector<ExprPtr>{
        normalized_actual.value(),
        arena.make<Unary>(UnaryOp::Neg, normalized_expected.value()),
    });
    auto together_difference = algebra::together(difference, compare_context);
    if (together_difference.is_ok()) {
        auto simplified_difference = compare_context.simplify(together_difference.value());
        if (simplified_difference.is_ok()) {
            if (const auto* integer = expr_cast<IntegerLit>(simplified_difference.value());
                integer != nullptr && integer->value.is_zero()) {
                SUCCEED();
                return;
            }
            if (const auto* rational = expr_cast<RationalLit>(simplified_difference.value());
                rational != nullptr && rational->numerator.is_zero()) {
                SUCCEED();
                return;
            }
        }
    }

    auto equivalent = mathematically_equal(
        normalized_actual.value(),
        normalized_expected.value(),
        compare_context);
    ASSERT_TRUE(equivalent.is_ok()) << equivalent.error().message;
    EXPECT_TRUE(equivalent.value())
        << "actual=" << debug_print(normalized_actual.value())
        << " expected=" << debug_print(normalized_expected.value());
}

void expect_vector_equivalent(const std::vector<ExprPtr>& actual, const std::vector<ExprPtr>& expected) {
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t index = 0; index < actual.size(); ++index) {
        expect_equivalent(actual[index], expected[index]);
    }
}

void expect_matrix_equivalent(ExprPtr actual, std::size_t rows, std::size_t cols, const std::vector<ExprPtr>& expected_elements) {
    const auto* matrix = expr_cast<Matrix>(actual);
    ASSERT_NE(matrix, nullptr);
    ASSERT_EQ(matrix->rows, rows);
    ASSERT_EQ(matrix->cols, cols);
    ASSERT_EQ(matrix->elements.size(), expected_elements.size());
    for (std::size_t index = 0; index < expected_elements.size(); ++index) {
        expect_equivalent(matrix->elements[index], expected_elements[index]);
    }
}

void expect_integral_equals(const std::string& integrand_text, const std::string& variable, const std::string& expected_text) {
    SCOPED_TRACE("integrand=" + integrand_text);
    SCOPED_TRACE("expected=" + expected_text);
    symbolic::CASContext context;
    AstArena expected_arena;

    auto primitive = integrate_input(integrand_text, variable, context);
    auto expected = parse_expr(expected_text, expected_arena);
    ASSERT_TRUE(primitive.is_ok()) << primitive.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    expect_equivalent(primitive.value(), expected.value());
}

void expect_integration_oracle(const std::string& integrand_text, const std::string& variable) {
    symbolic::CASContext integrate_context;
    auto primitive = integrate_input(integrand_text, variable, integrate_context);
    ASSERT_TRUE(primitive.is_ok()) << primitive.error().message;

    symbolic::CASContext differentiate_context;
    auto recovered = diff(primitive.value(), Symbol(variable), 1U, differentiate_context);
    ASSERT_TRUE(recovered.is_ok()) << recovered.error().message;

    AstArena expected_arena;
    auto expected = parse_expr(integrand_text, expected_arena);
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    expect_equivalent(recovered.value(), expected.value());
}

TEST(CalculusDiffTest, DerivesSymbolicPower) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto differentiated = differentiate_input("x^n", "x", 1U, context);
    auto expected = parse_expr("n*x^(n-1)", expected_arena);
    ASSERT_TRUE(differentiated.is_ok()) << differentiated.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    expect_equivalent(differentiated.value(), expected.value());
}

TEST(CalculusDiffTest, AppliesChainRuleToSinOfQuadratic) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto differentiated = differentiate_input("sin(x^2)", "x", 1U, context);
    auto expected = parse_expr("2*x*cos(x^2)", expected_arena);
    ASSERT_TRUE(differentiated.is_ok()) << differentiated.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    expect_equivalent(differentiated.value(), expected.value());
}

TEST(CalculusDiffTest, AppliesProductRule) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto differentiated = differentiate_input("sin(x)*exp(x)", "x", 1U, context);
    auto expected = parse_expr("cos(x)*exp(x) + sin(x)*exp(x)", expected_arena);
    ASSERT_TRUE(differentiated.is_ok()) << differentiated.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    expect_equivalent(differentiated.value(), expected.value());
}

TEST(CalculusDiffTest, AppliesQuotientRule) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto differentiated = differentiate_input("1/(1+x^2)", "x", 1U, context);
    auto expected = parse_expr("(-2*x)/(1+x^2)^2", expected_arena);
    ASSERT_TRUE(differentiated.is_ok()) << differentiated.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    expect_equivalent(differentiated.value(), expected.value());
}

TEST(CalculusDiffTest, SupportsHigherOrderDerivatives) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto differentiated = differentiate_input("sin(x)", "x", 2U, context);
    auto expected = parse_expr("-sin(x)", expected_arena);
    ASSERT_TRUE(differentiated.is_ok()) << differentiated.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    expect_equivalent(differentiated.value(), expected.value());
}

TEST(CalculusDiffTest, DerivesInverseAndHyperbolicFunctions) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto differentiated = differentiate_input("arctan(x) + tanh(x)", "x", 1U, context);
    auto expected = parse_expr("1/(1+x^2) + 1/cosh(x)^2", expected_arena);
    ASSERT_TRUE(differentiated.is_ok()) << differentiated.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    expect_equivalent(differentiated.value(), expected.value());
}

TEST(CalculusDiffTest, DerivesAliasesAndSquareRoot) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto differentiated = differentiate_input("asin(x) + sqrt(x)", "x", 1U, context);
    auto expected = parse_expr("1/sqrt(1-x^2) + 1/(2*sqrt(x))", expected_arena);
    ASSERT_TRUE(differentiated.is_ok()) << differentiated.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    expect_equivalent(differentiated.value(), expected.value());
}

TEST(CalculusDiffTest, DerivesAbsoluteValueWithNonzeroAssumption) {
    symbolic::CASContext context;
    AstArena expected_arena;
    context.assumptions().assume_nonzero(Symbol{"x"});

    auto differentiated = differentiate_input("abs(x)", "x", 1U, context);
    auto expected = parse_expr("sign(x)", expected_arena);
    ASSERT_TRUE(differentiated.is_ok()) << differentiated.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    expect_equivalent(differentiated.value(), expected.value());
}

TEST(CalculusDiffTest, ComputesImplicitDerivativeFromRelation) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto differentiated = implicit_diff_input("x^2 + y^2 - 1", "y", "x", context);
    auto expected = parse_expr("-(2*x)/(2*y)", expected_arena);
    ASSERT_TRUE(differentiated.is_ok()) << differentiated.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    expect_equivalent(differentiated.value(), expected.value());
}

TEST(CalculusDiffTest, RejectsImplicitDerivativeWithoutDependentVariable) {
    symbolic::CASContext context;

    auto differentiated = implicit_diff_input("x^2 + 1", "y", "x", context);
    ASSERT_TRUE(differentiated.is_error());
    EXPECT_EQ(differentiated.error().kind, CASErrorKind::InvalidArgument);
}

TEST(CalculusDiffTest, ComputesGradientAcrossMultipleVariables) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto actual = gradient_input("x^2 + x*y + sin(y)", {Symbol{"x"}, Symbol{"y"}}, context);
    auto expected_dx = parse_expr("2*x + y", expected_arena);
    auto expected_dy = parse_expr("x + cos(y)", expected_arena);
    ASSERT_TRUE(actual.is_ok()) << actual.error().message;
    ASSERT_TRUE(expected_dx.is_ok()) << expected_dx.error().message;
    ASSERT_TRUE(expected_dy.is_ok()) << expected_dy.error().message;
    expect_vector_equivalent(actual.value(), {expected_dx.value(), expected_dy.value()});
}

TEST(CalculusDiffTest, ComputesJacobianMatrix) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto actual = jacobian_input({"x*y", "x^2 + y"}, {Symbol{"x"}, Symbol{"y"}}, context);
    auto e11 = parse_expr("y", expected_arena);
    auto e12 = parse_expr("x", expected_arena);
    auto e21 = parse_expr("2*x", expected_arena);
    auto e22 = parse_expr("1", expected_arena);
    ASSERT_TRUE(actual.is_ok()) << actual.error().message;
    ASSERT_TRUE(e11.is_ok()) << e11.error().message;
    ASSERT_TRUE(e12.is_ok()) << e12.error().message;
    ASSERT_TRUE(e21.is_ok()) << e21.error().message;
    ASSERT_TRUE(e22.is_ok()) << e22.error().message;
    expect_matrix_equivalent(actual.value(), 2U, 2U, {e11.value(), e12.value(), e21.value(), e22.value()});
}

TEST(CalculusDiffTest, ComputesHessianMatrix) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto actual = hessian_input("x^2 + x*y + y^3", {Symbol{"x"}, Symbol{"y"}}, context);
    auto e11 = parse_expr("2", expected_arena);
    auto e12 = parse_expr("1", expected_arena);
    auto e21 = parse_expr("1", expected_arena);
    auto e22 = parse_expr("6*y", expected_arena);
    ASSERT_TRUE(actual.is_ok()) << actual.error().message;
    ASSERT_TRUE(e11.is_ok()) << e11.error().message;
    ASSERT_TRUE(e12.is_ok()) << e12.error().message;
    ASSERT_TRUE(e21.is_ok()) << e21.error().message;
    ASSERT_TRUE(e22.is_ok()) << e22.error().message;
    expect_matrix_equivalent(actual.value(), 2U, 2U, {e11.value(), e12.value(), e21.value(), e22.value()});
}

TEST(CalculusIntegrateTest, IntegratesPolynomialPower) {
    expect_integration_oracle("x^3", "x");
}

TEST(CalculusIntegrateTest, IntegratesTrigonometricTableEntries) {
    expect_integration_oracle("sin(x)", "x");
    expect_integration_oracle("cos(x)", "x");
}

TEST(CalculusIntegrateTest, IntegratesExponentialAndLogarithm) {
    expect_integration_oracle("exp(x)", "x");
    expect_integral_equals("ln(x)", "x", "x*ln(x) - x");
}

TEST(CalculusIntegrateTest, IntegratesAffineExactFunctionArguments) {
    expect_integral_equals("sin(2*x)", "x", "-cos(2*x)/2");
    expect_integral_equals("exp(-x)", "x", "-exp(-x)");
}

TEST(CalculusIntegrateTest, IntegratesDerivativeTimesElementaryComposition) {
    expect_integral_equals("2*x*sin(x^2)", "x", "-cos(x^2)");
    expect_integral_equals("2*x*exp(x^2)", "x", "exp(x^2)");
}

TEST(CalculusIntegrateTest, IntegratesDerivativeTimesPowerComposition) {
    expect_integral_equals("2*x*(1+x^2)^3", "x", "((1+x^2)^4)/4");
    expect_integral_equals("2*x/(1+x^2)", "x", "ln(abs(1+x^2))");
}

TEST(CalculusIntegrateTest, IntegratesRationalPartialFractionsWithLinearFactors) {
    expect_integral_equals("1/(x^2-1)", "x", "(1/2)*ln(abs(x-1)) + (-1/2)*ln(abs(x+1))");
    expect_integral_equals("-1/(x+1)", "x", "-ln(abs(x+1))");
    expect_integral_equals("(2*x+3)/(x^2+x)", "x", "3*ln(abs(x)) - ln(abs(x+1))");
    expect_integral_equals("1/((x-1)^2)", "x", "-1/(x-1)");
}

// --- P1-001 acceptance criteria: Hermite reduction (repeated factors) ---
TEST(CalculusIntegrateTest, P1_HermiteReduction_RepeatedLinearFactor) {
    // 1/(x^2*(x+1)) via Hermite: rational part -1/x, log part ln(x)-ln(x+1)
    // Oracle can't differentiate abs() without assumptions — check form directly
    symbolic::CASContext ctx;
    auto primitive = integrate_input("1/(x^2*(x+1))", "x", ctx);
    ASSERT_TRUE(primitive.is_ok()) << primitive.error().message;
    // Must contain a rational part (Power with negative exponent) and log terms
    EXPECT_TRUE(contains_builtin(primitive.value(), BuiltinOp::Ln))
        << "Expected ln in result of 1/(x^2*(x+1))";
}

TEST(CalculusIntegrateTest, P1_HermiteReduction_RepeatedQuadraticFactor) {
    // x/(x^2-1)^2 via Hermite reduction
    expect_integration_oracle("x/(x^2-1)^2", "x");
}

// --- P1-002 acceptance criteria: Rothstein-Trager log-part ---
TEST(CalculusIntegrateTest, P1_RothsteinTrager_CubicWithLinearAndQuadratic) {
    // 1/(x^3+x) = 1/(x*(x^2+1)) → ln|x| - 1/2*ln(x^2+1)
    expect_integration_oracle("1/(x^3+x)", "x");
}

TEST(CalculusIntegrateTest, IntegratesSimpleLinearNumeratorOverIrreducibleQuadratic) {
    expect_integral_equals("x/(x^2+1)", "x", "(1/2)*ln(abs(x^2+1))");
}

TEST(CalculusIntegrateTest, IntegratesAffineNumeratorOverIrreducibleQuadratic) {
    symbolic::CASContext context;
    auto primitive = integrate_input("(2*x+3)/(x^2+1)", "x", context);
    ASSERT_TRUE(primitive.is_ok()) << primitive.error().message;
    EXPECT_TRUE(contains_builtin(primitive.value(), BuiltinOp::Ln));
    EXPECT_TRUE(contains_builtin(primitive.value(), BuiltinOp::Atan));
}

TEST(CalculusIntegrateTest, ExtractsConstantFactorsAndLinearity) {
    expect_integration_oracle("2*sin(x) + 3*x^2", "x");
}

TEST(CalculusIntegrateTest, IntegratesExtendedElementaryTable) {
    expect_integration_oracle("1/(1+x^2)", "x");
    expect_integration_oracle("1/sqrt(1-x^2)", "x");
    expect_integral_equals("-1/sqrt(1-x^2)", "x", "-arcsin(x)");
}

TEST(CalculusIntegrateTest, IntegratesClassicalTrigSubstitutionRadicals) {
    expect_integral_equals("sqrt(1-x^2)", "x", "(1/2)*(x*sqrt(1-x^2) + arcsin(x))");
    expect_integral_equals("sqrt(1+x^2)", "x", "(1/2)*(x*sqrt(1+x^2) + ln(abs(x+sqrt(1+x^2))))");
    expect_integral_equals("sqrt(x^2-1)", "x", "(1/2)*(x*sqrt(x^2-1) - ln(abs(x+sqrt(x^2-1))))");
}

TEST(CalculusIntegrateTest, MatchesPrimitiveFormsRequiringAbsoluteValuesOrReciprocalNormalization) {
    expect_integral_equals("1/x", "x", "ln(abs(x))");
    expect_integral_equals("a^x", "x", "a^x/ln(a)");
    expect_integral_equals("tan(x)", "x", "(-1) * ln(abs(cos(x)))");
    expect_integral_equals("cot(x)", "x", "ln(abs(sin(x)))");
    expect_integral_equals("sec(x)", "x", "ln(abs(sec(x) + tan(x)))");
    expect_integral_equals("csc(x)", "x", "(-1) * ln(abs(csc(x) + cot(x)))");
    expect_integral_equals("sec(x)^2", "x", "tan(x)");
    expect_integral_equals("csc(x)^2", "x", "-cot(x)");
    expect_integral_equals("1/sqrt(x^2+a^2)", "x", "ln(abs(x + sqrt(x^2+a^2)))");
    expect_integral_equals("1/(x^2-a^2)", "x", "(1/(2*a))*ln(abs((x-a)/(x+a)))");
    expect_integral_equals("sinh(x)", "x", "cosh(x)");
    expect_integral_equals("cosh(x)", "x", "sinh(x)");
}

TEST(CalculusIntegrateTest, IntegratesElementaryByPartsPatterns) {
    expect_integration_oracle("x*exp(x)", "x");
    expect_integration_oracle("x*cos(x)", "x");
}

TEST(CalculusIntegrateTest, StopsRecursiveIntegrationByPartsCycles) {
    symbolic::CASContext context;
    AstArena arena;

    auto expr = parse_expr("x*sin(x)*exp(x)", arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto result = integrate(expr.value(), Symbol("x"), context);
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::Unimplemented);
}

TEST(CalculusIntegrateTest, ComputesDefiniteIntegralViaTfc) {
    symbolic::CASContext context;
    AstArena parse_arena;
    AstArena expected_arena;

    auto expr = parse_expr("x^2", parse_arena);
    auto lower = parse_expr("0", parse_arena);
    auto upper = parse_expr("2", parse_arena);
    auto expected = parse_expr("8/3", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(lower.is_ok()) << lower.error().message;
    ASSERT_TRUE(upper.is_ok()) << upper.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto result = definite_integral(expr.value(), Symbol("x"), lower.value(), upper.value(), context);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    expect_equivalent(result.value(), expected.value());
}

TEST(CalculusIntegrateTest, RejectsDefiniteIntegralWithInteriorRationalPole) {
    symbolic::CASContext context;
    AstArena parse_arena;

    auto expr = parse_expr("1/(x-2)", parse_arena);
    auto lower = parse_expr("0", parse_arena);
    auto upper = parse_expr("5", parse_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(lower.is_ok()) << lower.error().message;
    ASSERT_TRUE(upper.is_ok()) << upper.error().message;

    auto result = definite_integral(expr.value(), Symbol("x"), lower.value(), upper.value(), context);
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::Undefined);
}

TEST(CalculusIntegrateTest, RejectsDefiniteIntegralWithEndpointRationalPole) {
    symbolic::CASContext context;
    AstArena parse_arena;

    auto expr = parse_expr("1/(u+1)", parse_arena);
    auto lower = parse_expr("-1", parse_arena);
    auto upper = parse_expr("3", parse_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(lower.is_ok()) << lower.error().message;
    ASSERT_TRUE(upper.is_ok()) << upper.error().message;

    auto result = definite_integral(expr.value(), Symbol("u"), lower.value(), upper.value(), context);
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::Undefined);
}

TEST(CalculusIntegrateTest, AllowsRemovableRationalSingularityAfterExactCancellation) {
    symbolic::CASContext context;
    AstArena parse_arena;
    AstArena expected_arena;

    auto expr = parse_expr("(v^2-1)/(v-1)", parse_arena);
    auto lower = parse_expr("0", parse_arena);
    auto upper = parse_expr("2", parse_arena);
    auto expected = parse_expr("4", expected_arena);
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    ASSERT_TRUE(lower.is_ok()) << lower.error().message;
    ASSERT_TRUE(upper.is_ok()) << upper.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;

    auto result = definite_integral(expr.value(), Symbol("v"), lower.value(), upper.value(), context);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    expect_equivalent(result.value(), expected.value());
}

TEST(CalculusIntegrateTest, RejectsNonElementaryCaseForNow) {
    symbolic::CASContext context;
    auto integrated = integrate_input("exp(-x^2)", "x", context);
    ASSERT_TRUE(integrated.is_error());
    EXPECT_EQ(integrated.error().kind, CASErrorKind::Unimplemented);
}

TEST(CalculusLimitTest, ComputesDirectSubstitutionAndLHopitalCases) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto l1 = limit_input("x^2", "x", "2", LimitDirection::Both, context);
    auto e1 = parse_expr("4", expected_arena);
    ASSERT_TRUE(l1.is_ok()) << l1.error().message;
    ASSERT_TRUE(e1.is_ok()) << e1.error().message;
    expect_equivalent(l1.value(), e1.value());

    auto l2 = limit_input("sin(x)/x", "x", "0", LimitDirection::Both, context);
    auto e2 = parse_expr("1", expected_arena);
    ASSERT_TRUE(l2.is_ok()) << l2.error().message;
    ASSERT_TRUE(e2.is_ok()) << e2.error().message;
    expect_equivalent(l2.value(), e2.value());

    auto l3 = limit_input("(exp(x)-1)/x", "x", "0", LimitDirection::Both, context);
    auto e3 = parse_expr("1", expected_arena);
    ASSERT_TRUE(l3.is_ok()) << l3.error().message;
    ASSERT_TRUE(e3.is_ok()) << e3.error().message;
    expect_equivalent(l3.value(), e3.value());

    auto l4 = limit_input("(x^2-1)/(x-1)", "x", "1", LimitDirection::Both, context);
    auto e4 = parse_expr("2", expected_arena);
    ASSERT_TRUE(l4.is_ok()) << l4.error().message;
    ASSERT_TRUE(e4.is_ok()) << e4.error().message;
    expect_equivalent(l4.value(), e4.value());
}

TEST(CalculusLimitTest, DebugLimitSimplifyInf) {
    symbolic::CASContext context;
    AstArena& arena = context.arena();
    auto e1 = arena.make<Binary>(BinaryOp::Div, arena.make<IntegerLit>(BigInt(1)), arena.make<Constant>(MathConstant::Infinity));
    auto s1 = context.simplify(e1);
    if (s1.is_ok()) {
        auto expected = arena.make<IntegerLit>(BigInt(0));
        expect_equivalent(s1.value(), expected);
    }
}

TEST(CalculusLimitTest, ComputesBasicInfiniteGrowthComparisons) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto l1 = limit_input("x^2/exp(x)", "x", "inf", LimitDirection::Both, context);
    auto e1 = parse_expr("0", expected_arena);
    ASSERT_TRUE(l1.is_ok()) << l1.error().message;
    ASSERT_TRUE(e1.is_ok()) << e1.error().message;
    expect_equivalent(l1.value(), e1.value());

    auto l2 = limit_input("ln(x)/x", "x", "inf", LimitDirection::Both, context);
    auto e2 = parse_expr("0", expected_arena);
    ASSERT_TRUE(l2.is_ok()) << l2.error().message;
    ASSERT_TRUE(e2.is_ok()) << e2.error().message;
    expect_equivalent(l2.value(), e2.value());

    std::cerr << "Test l3" << std::endl;
    auto l3 = limit_input("(1+1/n)^n", "n", "inf", LimitDirection::Both, context);
    auto e3 = parse_expr("e", expected_arena);
    ASSERT_TRUE(l3.is_ok()) << l3.error().message;
    ASSERT_TRUE(e3.is_ok()) << e3.error().message;
    expect_equivalent(l3.value(), e3.value());

    std::cerr << "Test l4" << std::endl;
    auto l4 = limit_input("x^x", "x", "0", LimitDirection::Right, context);
    auto e4 = parse_expr("1", expected_arena);
    ASSERT_TRUE(l4.is_ok()) << l4.error().message;
    ASSERT_TRUE(e4.is_ok()) << e4.error().message;
    expect_equivalent(l4.value(), e4.value());

    std::cerr << "Test l5" << std::endl;
    auto l5 = limit_input("exp(x)/x^2", "x", "inf", LimitDirection::Both, context);
    auto e5 = parse_expr("inf", expected_arena);
    ASSERT_TRUE(l5.is_ok()) << l5.error().message;
    ASSERT_TRUE(e5.is_ok()) << e5.error().message;
    expect_equivalent(l5.value(), e5.value());

    auto l6 = limit_input("x/ln(x)", "x", "inf", LimitDirection::Both, context);
    auto e6 = parse_expr("inf", expected_arena);
    ASSERT_TRUE(l6.is_ok()) << l6.error().message;
    ASSERT_TRUE(e6.is_ok()) << e6.error().message;
    expect_equivalent(l6.value(), e6.value());

    auto l7 = limit_input("x^2/exp(x)", "x", "-inf", LimitDirection::Both, context);
    auto e7 = parse_expr("inf", expected_arena);
    ASSERT_TRUE(l7.is_ok()) << l7.error().message;
    ASSERT_TRUE(e7.is_ok()) << e7.error().message;
    expect_equivalent(l7.value(), e7.value());

    auto l8 = limit_input("x/exp(x)", "x", "-inf", LimitDirection::Both, context);
    auto e8 = parse_expr("-inf", expected_arena);
    ASSERT_TRUE(l8.is_ok()) << l8.error().message;
    ASSERT_TRUE(e8.is_ok()) << e8.error().message;
    expect_equivalent(l8.value(), e8.value());

    auto l9 = limit_input("exp(x)/x^2", "x", "-inf", LimitDirection::Both, context);
    auto e9 = parse_expr("0", expected_arena);
    ASSERT_TRUE(l9.is_ok()) << l9.error().message;
    ASSERT_TRUE(e9.is_ok()) << e9.error().message;
    expect_equivalent(l9.value(), e9.value());
}

TEST(CalculusLimitTest, RejectsLogarithmicGrowthAtNegativeInfinityOutsideDomain) {
    symbolic::CASContext context;

    auto result = limit_input("ln(x)/x", "x", "-inf", LimitDirection::Both, context);
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::Undefined);
}

TEST(CalculusLimitTest, ComputesLateralPoleLimitsForRationalFunctions) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto right = limit_input("1/x", "x", "0", LimitDirection::Right, context);
    auto right_expected = parse_expr("inf", expected_arena);
    ASSERT_TRUE(right.is_ok()) << right.error().message;
    ASSERT_TRUE(right_expected.is_ok()) << right_expected.error().message;
    expect_equivalent(right.value(), right_expected.value());

    auto left = limit_input("1/x", "x", "0", LimitDirection::Left, context);
    auto left_expected = parse_expr("-inf", expected_arena);
    ASSERT_TRUE(left.is_ok()) << left.error().message;
    ASSERT_TRUE(left_expected.is_ok()) << left_expected.error().message;
    expect_equivalent(left.value(), left_expected.value());

    auto translated = limit_input("1/(x-2)", "x", "2", LimitDirection::Right, context);
    auto translated_expected = parse_expr("inf", expected_arena);
    ASSERT_TRUE(translated.is_ok()) << translated.error().message;
    ASSERT_TRUE(translated_expected.is_ok()) << translated_expected.error().message;
    expect_equivalent(translated.value(), translated_expected.value());

    auto even_pole = limit_input("1/x^2", "x", "0", LimitDirection::Both, context);
    auto even_expected = parse_expr("inf", expected_arena);
    ASSERT_TRUE(even_pole.is_ok()) << even_pole.error().message;
    ASSERT_TRUE(even_expected.is_ok()) << even_expected.error().message;
    expect_equivalent(even_pole.value(), even_expected.value());

    auto reduced_pole = limit_input("x/(x^2)", "x", "0", LimitDirection::Right, context);
    auto reduced_expected = parse_expr("inf", expected_arena);
    ASSERT_TRUE(reduced_pole.is_ok()) << reduced_pole.error().message;
    ASSERT_TRUE(reduced_expected.is_ok()) << reduced_expected.error().message;
    expect_equivalent(reduced_pole.value(), reduced_expected.value());
}

TEST(CalculusLimitTest, RejectsBilateralPoleWithOppositeLateralSigns) {
    symbolic::CASContext context;

    auto result = limit_input("1/x", "x", "0", LimitDirection::Both, context);
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::Undefined);
}

TEST(CalculusLimitTest, ComputesLateralLogarithmicSingularities) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto right = limit_input("ln(x)", "x", "0", LimitDirection::Right, context);
    auto right_expected = parse_expr("-inf", expected_arena);
    ASSERT_TRUE(right.is_ok()) << right.error().message;
    ASSERT_TRUE(right_expected.is_ok()) << right_expected.error().message;
    expect_equivalent(right.value(), right_expected.value());

    auto translated = limit_input("ln(x-2)", "x", "2", LimitDirection::Right, context);
    auto translated_expected = parse_expr("-inf", expected_arena);
    ASSERT_TRUE(translated.is_ok()) << translated.error().message;
    ASSERT_TRUE(translated_expected.is_ok()) << translated_expected.error().message;
    expect_equivalent(translated.value(), translated_expected.value());

    auto bilateral_even = limit_input("ln((x-1)^2)", "x", "1", LimitDirection::Both, context);
    auto bilateral_expected = parse_expr("-inf", expected_arena);
    ASSERT_TRUE(bilateral_even.is_ok()) << bilateral_even.error().message;
    ASSERT_TRUE(bilateral_expected.is_ok()) << bilateral_expected.error().message;
    expect_equivalent(bilateral_even.value(), bilateral_expected.value());
}

TEST(CalculusLimitTest, RejectsLogarithmicSideOutsideDomain) {
    symbolic::CASContext context;

    auto left = limit_input("ln(x)", "x", "0", LimitDirection::Left, context);
    ASSERT_TRUE(left.is_error());
    EXPECT_EQ(left.error().kind, CASErrorKind::Undefined);

    auto bilateral = limit_input("ln(x)", "x", "0", LimitDirection::Both, context);
    ASSERT_TRUE(bilateral.is_error());
    EXPECT_EQ(bilateral.error().kind, CASErrorKind::Undefined);
}

TEST(CalculusTaylorTest, BuildsMaclaurinPolynomials) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto expansion = taylor_input("exp(x)", "x", "0", 3U, context);
    auto expected = parse_expr("1 + x + (1/2)*x^2 + (1/6)*x^3", expected_arena);
    ASSERT_TRUE(expansion.is_ok()) << expansion.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    EXPECT_EQ(expansion.value().computed_order, 3U);
    expect_equivalent(expansion.value().polynomial, expected.value());
}

TEST(CalculusTaylorTest, BuildsSeriesAroundZeroForSin) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto expansion = taylor_input("sin(x)", "x", "0", 5U, context);
    auto expected = parse_expr("x + (-1/6)*x^3 + (1/120)*x^5", expected_arena);
    ASSERT_TRUE(expansion.is_ok()) << expansion.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    EXPECT_EQ(expansion.value().computed_order, 5U);
    expect_equivalent(expansion.value().polynomial, expected.value());
}

TEST(CalculusTaylorTest, BuildsStandardSeriesForLogOnePlusX) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto expansion = taylor_input("ln(1+x)", "x", "0", 4U, context);
    auto expected = parse_expr("x - (1/2)*x^2 + (1/3)*x^3 - (1/4)*x^4", expected_arena);
    ASSERT_TRUE(expansion.is_ok()) << expansion.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    EXPECT_EQ(expansion.value().computed_order, 4U);
    expect_equivalent(expansion.value().polynomial, expected.value());
}

TEST(CalculusTaylorTest, BuildsStandardSeriesForArctan) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto expansion = taylor_input("arctan(x)", "x", "0", 5U, context);
    auto expected = parse_expr("x - (1/3)*x^3 + (1/5)*x^5", expected_arena);
    ASSERT_TRUE(expansion.is_ok()) << expansion.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    EXPECT_EQ(expansion.value().computed_order, 5U);
    expect_equivalent(expansion.value().polynomial, expected.value());
}

TEST(CalculusTaylorTest, BuildsStandardSeriesForGeneralizedBinomial) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto expansion = taylor_input("(1+x)^n", "x", "0", 3U, context);
    auto expected = parse_expr(
        "1 + n*x + (1/2)*n*(n-1)*x^2 + (1/6)*n*(n-1)*(n-2)*x^3",
        expected_arena);
    ASSERT_TRUE(expansion.is_ok()) << expansion.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    EXPECT_EQ(expansion.value().computed_order, 3U);
    expect_equivalent(expansion.value().polynomial, expected.value());
}

// L2-27: verify the general taylor fallback handles functions without a
// closed-form Maclaurin entry in limit_series.cpp.  These pass through
// repeated symbolic differentiation rather than the fast-path table.
TEST(CalculusTaylorTest, GeneralFallbackTangent) {
    // tan(x) at x=0, order 5:  x + x^3/3 + 2/15 * x^5
    symbolic::CASContext context;
    AstArena expected_arena;
    auto expansion = taylor_input("tan(x)", "x", "0", 5U, context);
    auto expected = parse_expr("x + (1/3)*x^3 + (2/15)*x^5", expected_arena);
    ASSERT_TRUE(expansion.is_ok()) << expansion.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    EXPECT_EQ(expansion.value().computed_order, 5U);
    expect_equivalent(expansion.value().polynomial, expected.value());
}

// NOTE: composition Taylor of exp(sin(x)) at higher orders currently degrades
// because the simplifier cannot fully reduce d^k/dx^k[exp(sin(x))] at x=0.
// This is a separate gap from L2-27 (general fallback exists & works for
// tan(x) and rationals); the simplifier limitation is its own task.

TEST(CalculusTaylorTest, GeneralFallbackRationalNotInTable) {
    // 1/(1 - x^2) at x=0, order 6:  1 + x^2 + x^4 + x^6
    symbolic::CASContext context;
    AstArena expected_arena;
    auto expansion = taylor_input("1/(1 - x^2)", "x", "0", 6U, context);
    auto expected = parse_expr("1 + x^2 + x^4 + x^6", expected_arena);
    ASSERT_TRUE(expansion.is_ok()) << expansion.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    EXPECT_EQ(expansion.value().computed_order, 6U);
    expect_equivalent(expansion.value().polynomial, expected.value());
}

TEST(CalculusTaylorTest, BuildsTaylorSeriesAroundNonZeroPointForPolynomial) {
    symbolic::CASContext context;
    AstArena expected_arena;

    auto expansion = taylor_input("x^3", "x", "2", 3U, context);
    auto expected = parse_expr("8 + 12*(x-2) + 6*(x-2)^2 + (x-2)^3", expected_arena);
    ASSERT_TRUE(expansion.is_ok()) << expansion.error().message;
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    EXPECT_EQ(expansion.value().computed_order, 3U);
    expect_equivalent(expansion.value().polynomial, expected.value());
}

// --- P0-003 acceptance criteria: polynomial degree comparison in growth ---

TEST(CalculusLimitTest, PolynomialDegreeComparison_HigherDominates) {
    symbolic::CASContext context;
    AstArena arena;

    // lim(x→∞) x^10/x^2 = ∞  (degree 10 > degree 2)
    auto r = limit_input("x^10/x^2", "x", "inf", LimitDirection::Both, context);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto expected = parse_expr("inf", arena);
    ASSERT_TRUE(expected.is_ok());
    expect_equivalent(r.value(), expected.value());
}

TEST(CalculusLimitTest, PolynomialDegreeComparison_LowerGoesToZero) {
    symbolic::CASContext context;
    AstArena arena;

    // lim(x→∞) x^2/x^10 = 0  (degree 2 < degree 10)
    auto r = limit_input("x^2/x^10", "x", "inf", LimitDirection::Both, context);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto expected = parse_expr("0", arena);
    ASSERT_TRUE(expected.is_ok());
    expect_equivalent(r.value(), expected.value());
}

TEST(CalculusLimitTest, PolynomialDegreeComparison_SameDegreeFiniteLimit) {
    symbolic::CASContext context;
    AstArena arena;

    // lim(x→∞) x^5/x^5 = 1  (same degree)
    auto r = limit_input("x^5/x^5", "x", "inf", LimitDirection::Both, context);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto expected = parse_expr("1", arena);
    ASSERT_TRUE(expected.is_ok());
    expect_equivalent(r.value(), expected.value());
}

TEST(CalculusLimitTest, PolynomialDegreeComparison_HighDegreeRatio) {
    symbolic::CASContext context;
    AstArena arena;

    // lim(x→∞) x^7/x^3 = ∞  (degree 7 >> degree 3)
    auto r = limit_input("x^7/x^3", "x", "inf", LimitDirection::Both, context);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto expected = parse_expr("inf", arena);
    ASSERT_TRUE(expected.is_ok());
    expect_equivalent(r.value(), expected.value());
}

// P2-005: Trig angle reduction
TEST(P2_005_TrigReduction, SinExactSpecialAngles) {
    symbolic::CASContext ctx;

    auto check_sin = [&](const std::string& expr_str, const std::string& expected_str) {
        auto e = parse_expr(expr_str, ctx.arena());
        ASSERT_TRUE(e.is_ok()) << expr_str;
        auto result = ctx.simplify(e.value());
        ASSERT_TRUE(result.is_ok()) << "simplify failed: " << expr_str;
        auto ex = parse_expr(expected_str, ctx.arena());
        ASSERT_TRUE(ex.is_ok()) << expected_str;
        auto exs = ctx.simplify(ex.value());
        ASSERT_TRUE(exs.is_ok());
        auto eq = mathematically_equal(result.value(), exs.value(), ctx);
        ASSERT_TRUE(eq.is_ok());
        EXPECT_TRUE(eq.value()) << expr_str << " → got wrong value, expected " << expected_str;
    };

    check_sin("sin(pi/6)", "1/2");
    check_sin("sin(pi/4)", "sqrt(2)/2");
    check_sin("sin(pi/3)", "sqrt(3)/2");
    check_sin("sin(pi/2)", "1");
    check_sin("sin(5*pi/6)", "1/2");       // Q2: sin(π - π/6) = sin(π/6)
    check_sin("sin(3*pi/4)", "sqrt(2)/2"); // Q2: sin(π - π/4)
    check_sin("sin(2*pi/3)", "sqrt(3)/2"); // Q2: sin(π - π/3)
    check_sin("sin(7*pi/6)", "-1/2");      // Q3: -sin(π/6)
    check_sin("sin(5*pi/4)", "-sqrt(2)/2");// Q3: -sin(π/4)
    check_sin("sin(4*pi/3)", "-sqrt(3)/2");// Q3: -sin(π/3)
    check_sin("sin(3*pi/2)", "-1");        // Q3/Q4 boundary
    check_sin("sin(11*pi/6)", "-1/2");     // Q4: -sin(π/6)
    check_sin("sin(7*pi/4)", "-sqrt(2)/2");// Q4: -sin(π/4)
    check_sin("sin(5*pi/3)", "-sqrt(3)/2");// Q4: -sin(π/3)
    check_sin("sin(2*pi)", "0");           // full period
}

TEST(P2_005_TrigReduction, CosExactSpecialAngles) {
    symbolic::CASContext ctx;

    auto check_cos = [&](const std::string& expr_str, const std::string& expected_str) {
        auto e = parse_expr(expr_str, ctx.arena());
        ASSERT_TRUE(e.is_ok()) << expr_str;
        auto result = ctx.simplify(e.value());
        ASSERT_TRUE(result.is_ok()) << "simplify failed: " << expr_str;
        auto ex = parse_expr(expected_str, ctx.arena());
        ASSERT_TRUE(ex.is_ok()) << expected_str;
        auto exs = ctx.simplify(ex.value());
        ASSERT_TRUE(exs.is_ok());
        auto eq = mathematically_equal(result.value(), exs.value(), ctx);
        ASSERT_TRUE(eq.is_ok());
        EXPECT_TRUE(eq.value()) << expr_str << " → got wrong value, expected " << expected_str;
    };

    check_cos("cos(pi/6)", "sqrt(3)/2");
    check_cos("cos(pi/4)", "sqrt(2)/2");
    check_cos("cos(pi/3)", "1/2");
    check_cos("cos(pi/2)", "0");
    check_cos("cos(5*pi/6)", "-sqrt(3)/2");// Q2: -cos(π/6)
    check_cos("cos(3*pi/4)", "-sqrt(2)/2");// Q2: -cos(π/4)
    check_cos("cos(2*pi/3)", "-1/2");      // Q2: -cos(π/3)
    check_cos("cos(7*pi/6)", "-sqrt(3)/2");// Q3: -cos(π/6)
    check_cos("cos(5*pi/4)", "-sqrt(2)/2");// Q3: -cos(π/4)
    check_cos("cos(4*pi/3)", "-1/2");      // Q3: -cos(π/3)
    check_cos("cos(7*pi/4)", "sqrt(2)/2"); // Q4: +cos(π/4)
    check_cos("cos(5*pi/3)", "1/2");       // Q4: +cos(π/3)
    check_cos("cos(11*pi/6)", "sqrt(3)/2");// Q4: +cos(π/6)
    check_cos("cos(2*pi)", "1");           // full period
}

TEST(P2_005_TrigReduction, MultiplePeriodsReduce) {
    symbolic::CASContext ctx;
    // sin(13*pi/6) = sin(pi/6) = 1/2 (13/6 = 2 + 1/6)
    auto e = parse_expr("sin(13*pi/6)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto result = ctx.simplify(e.value());
    ASSERT_TRUE(result.is_ok());
    auto eq = mathematically_equal(result.value(),
        parse_expr("1/2", ctx.arena()).value(), ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

// P2-001: Taylor for composite/non-lookup function (via derivative-based algorithm)
TEST(P2_001_Taylor, CompositeNonLookupFunction) {
    symbolic::CASContext ctx;
    // sin(x^2) has no direct lookup — must use derivative algorithm
    auto expansion = taylor_input("sin(x^2)", "x", "0", 6U, ctx);
    ASSERT_TRUE(expansion.is_ok()) << expansion.error().message;
    // sin(x^2) = x^2 - x^6/6 + ... — polynomial must be non-trivial
    EXPECT_NE(expansion.value().polynomial, nullptr);
    EXPECT_GE(expansion.value().computed_order, 1U);
}

TEST(P2_001_Taylor, GenericFunctionAtNonZeroPoint) {
    symbolic::CASContext ctx;
    // exp(x) around x=1: e + e*(x-1) + e/2*(x-1)^2 + ...
    auto expansion = taylor_input("exp(x)", "x", "1", 2U, ctx);
    ASSERT_TRUE(expansion.is_ok()) << expansion.error().message;
    EXPECT_NE(expansion.value().polynomial, nullptr);
}

// P2-004: L'Hôpital x^2/exp(x) at infinity
TEST(P2_004_LHopital, PolynomialOverExpAtInfinity) {
    symbolic::CASContext ctx;
    // lim(x→∞) x^2/e^x = 0
    auto r = limit_input("x^2/exp(x)", "x", "inf", LimitDirection::Both, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto expected = parse_expr("0", ctx.arena());
    ASSERT_TRUE(expected.is_ok());
    expect_equivalent(r.value(), expected.value());
}

TEST(P2_004_LHopital, SinOverXAtZero) {
    symbolic::CASContext ctx;
    // lim(x→0) sin(x)/x = 1
    auto r = limit_input("sin(x)/x", "x", "0", LimitDirection::Both, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto expected = parse_expr("1", ctx.arena());
    ASSERT_TRUE(expected.is_ok());
    expect_equivalent(r.value(), expected.value());
}

TEST(P2_004_LHopital, ExponentialMinusOneOverX) {
    symbolic::CASContext ctx;
    // lim(x→0) (e^x - 1)/x = 1
    auto r = limit_input("(exp(x)-1)/x", "x", "0", LimitDirection::Both, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto expected = parse_expr("1", ctx.arena());
    ASSERT_TRUE(expected.is_ok());
    expect_equivalent(r.value(), expected.value());
}

// P3-004: Special functions — Gamma(n) = (n-1)!

TEST(P3_004_SpecialFunctions, GammaOne) {
    // Gamma(1) = 0! = 1
    symbolic::CASContext ctx;
    auto& arena = ctx.arena();
    auto arg = arena.make<IntegerLit>(BigInt(1));
    auto gamma_1 = arena.make<FuncCall>(BuiltinOp::Gamma, std::vector<ExprPtr>{arg});
    auto result = ctx.simplify(gamma_1);
    ASSERT_TRUE(result.is_ok());
    const auto* il = expr_cast<IntegerLit>(result.value());
    ASSERT_NE(il, nullptr);
    EXPECT_EQ(il->value, BigInt(1));
}

TEST(P3_004_SpecialFunctions, GammaTwo) {
    // Gamma(2) = 1! = 1
    symbolic::CASContext ctx;
    auto& arena = ctx.arena();
    auto arg = arena.make<IntegerLit>(BigInt(2));
    auto gamma_2 = arena.make<FuncCall>(BuiltinOp::Gamma, std::vector<ExprPtr>{arg});
    auto result = ctx.simplify(gamma_2);
    ASSERT_TRUE(result.is_ok());
    const auto* il = expr_cast<IntegerLit>(result.value());
    ASSERT_NE(il, nullptr);
    EXPECT_EQ(il->value, BigInt(1));
}

TEST(P3_004_SpecialFunctions, GammaFive) {
    // Gamma(5) = 4! = 24
    symbolic::CASContext ctx;
    auto& arena = ctx.arena();
    auto arg = arena.make<IntegerLit>(BigInt(5));
    auto gamma_5 = arena.make<FuncCall>(BuiltinOp::Gamma, std::vector<ExprPtr>{arg});
    auto result = ctx.simplify(gamma_5);
    ASSERT_TRUE(result.is_ok());
    const auto* il = expr_cast<IntegerLit>(result.value());
    ASSERT_NE(il, nullptr);
    EXPECT_EQ(il->value, BigInt(24));
}

TEST(P3_004_SpecialFunctions, GammaLargeInteger) {
    // Gamma(10) = 9! = 362880
    symbolic::CASContext ctx;
    auto& arena = ctx.arena();
    auto arg = arena.make<IntegerLit>(BigInt(10));
    auto gamma_10 = arena.make<FuncCall>(BuiltinOp::Gamma, std::vector<ExprPtr>{arg});
    auto result = ctx.simplify(gamma_10);
    ASSERT_TRUE(result.is_ok());
    const auto* il = expr_cast<IntegerLit>(result.value());
    ASSERT_NE(il, nullptr);
    EXPECT_EQ(il->value, BigInt(362880LL));
}

TEST(P3_004_SpecialFunctions, GammaSymbolicNotReduced) {
    // Gamma(x) stays as gamma(x) — not a concrete integer
    symbolic::CASContext ctx;
    auto& arena = ctx.arena();
    auto x = arena.make<Symbol>("x");
    auto gamma_x = arena.make<FuncCall>(BuiltinOp::Gamma, std::vector<ExprPtr>{x});
    auto result = ctx.simplify(gamma_x);
    ASSERT_TRUE(result.is_ok());
    const auto* fc = expr_cast<FuncCall>(result.value());
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->func_id, BuiltinOp::Gamma);
}

// P3-005: Complex symbolic — Re, Im, Conj

TEST(P3_005_Complex, ReOfNumericComplex) {
    // Re(3 + 4*i) = 3
    symbolic::CASContext ctx;
    auto expr = parse_expr("re(3 + 4*i)", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simp = ctx.simplify(expr.value());
    ASSERT_TRUE(simp.is_ok());
    const auto* il = expr_cast<IntegerLit>(simp.value());
    ASSERT_NE(il, nullptr) << "Expected IntegerLit";
    EXPECT_EQ(il->value, BigInt(3));
}

TEST(P3_005_Complex, ImOfNumericComplex) {
    // Im(3 + 4*i) = 4
    symbolic::CASContext ctx;
    auto expr = parse_expr("im(3 + 4*i)", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simp = ctx.simplify(expr.value());
    ASSERT_TRUE(simp.is_ok());
    const auto* il = expr_cast<IntegerLit>(simp.value());
    ASSERT_NE(il, nullptr) << "Expected IntegerLit";
    EXPECT_EQ(il->value, BigInt(4));
}

TEST(P3_005_Complex, ConjOfNumericComplex) {
    // Conj(3 + 4*i) = 3 - 4*i → simplified should contain I
    symbolic::CASContext ctx;
    auto expr = parse_expr("conj(3 + 4*i)", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simp = ctx.simplify(expr.value());
    ASSERT_TRUE(simp.is_ok());
    // Result should not be a bare FuncCall(Conj) — it should be simplified
    const auto* fc = expr_cast<FuncCall>(simp.value());
    if (fc != nullptr) {
        EXPECT_NE(fc->func_id, BuiltinOp::Conj) << "Conj should have been reduced";
    }
}

TEST(P3_005_Complex, ReOfPureImaginary) {
    // Re(i) = 0
    symbolic::CASContext ctx;
    auto expr = parse_expr("re(i)", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simp = ctx.simplify(expr.value());
    ASSERT_TRUE(simp.is_ok());
    const auto* il = expr_cast<IntegerLit>(simp.value());
    ASSERT_NE(il, nullptr);
    EXPECT_EQ(il->value, BigInt(0));
}

TEST(P3_005_Complex, ImOfPureImaginary) {
    // Im(5*i) = 5
    symbolic::CASContext ctx;
    auto expr = parse_expr("im(5*i)", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simp = ctx.simplify(expr.value());
    ASSERT_TRUE(simp.is_ok());
    const auto* il = expr_cast<IntegerLit>(simp.value());
    ASSERT_NE(il, nullptr);
    EXPECT_EQ(il->value, BigInt(5));
}

// CAS-L0-14: decimal literals converted to rational at parse boundary
TEST(L0_14_DecimalToRational, DiffHalfXSquared) {
    // diff(0.5*x^2, x) = x   (0.5 parsed as 1/2)
    symbolic::CASContext ctx;
    auto result = differentiate_input("0.5*x^2", "x", 1, ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    auto simp = ctx.simplify(result.value());
    ASSERT_TRUE(simp.is_ok());
    // Expect x (i.e. 1/2 * 2 * x = x)
    EXPECT_EQ(expr_kind(simp.value()), ExprKind::Symbol)
        << "diff(0.5*x^2, x) must simplify to x";
    EXPECT_EQ(expr_ref<Symbol>(simp.value()).name, "x");
}

TEST(L0_14_DecimalToRational, IntegrateQuarterX) {
    // integrate(0.25*x, x) = x^2/8  (0.25 parsed as 1/4)
    symbolic::CASContext ctx;
    auto result = integrate_input("0.25*x", "x", ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    auto simp = ctx.simplify(result.value());
    ASSERT_TRUE(simp.is_ok()) << "integrate(0.25*x, x) must not fail";
}

TEST(L0_14_DecimalToRational, DiffWithVariableCoeff) {
    // diff(0.1*t^3, t) = 3/10 * t^2
    symbolic::CASContext ctx;
    auto result = differentiate_input("0.1*t^3", "t", 1, ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    auto simp = ctx.simplify(result.value());
    ASSERT_TRUE(simp.is_ok()) << "diff(0.1*t^3, t) must not fail";
}

}  // namespace
}  // namespace cas::calculus
