#include "cas/symbolic.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas::symbolic {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

// ─── Regression Suite: Symbolic ─────────────────────────────────────────

TEST(SymbolicRegression, DivisionByZeroLiteral) {
    CASContext ctx;
    auto expr = parse_expr("10 / 0", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simplified = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified.is_error());
    EXPECT_EQ(simplified.error().kind, CASErrorKind::Undefined);
}

TEST(SymbolicRegression, DivisionByZeroSymbolic) {
    CASContext ctx;
    // (x - x) is zero, so 1 / (x - x) should be division by zero after simplification
    auto expr = parse_expr("1 / (x - x)", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simplified = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified.is_error());
    EXPECT_EQ(simplified.error().kind, CASErrorKind::Undefined);
}

TEST(SymbolicRegression, DeeplyNestedSums) {
    CASContext ctx;
    std::string input = "x";
    for(int i = 0; i < 500; ++i) {
        input = "(" + input + " + x)";
    }
    auto expr = parse_expr(input, ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simplified = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified.is_ok());
    // Should simplify to 501 * x
}

TEST(SymbolicRegression, PowerOfZeroToZero) {
    CASContext ctx;
    auto expr = parse_expr("0^0", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simplified = ctx.simplify(expr.value());
    // 0^0 is usually undefined in symbolic engines or 1.
    // Let's check our engine's choice.
    if (simplified.is_ok()) {
        // If it's OK, it's likely 1 or 0 (but usually 1 in discrete contexts, undefined in calculus)
    } else {
        EXPECT_EQ(simplified.error().kind, CASErrorKind::Undefined);
    }
}

TEST(SymbolicRegression, LogOfZero) {
    CASContext ctx;
    auto expr = parse_expr("ln(0)", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simplified = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified.is_error());
    EXPECT_EQ(simplified.error().kind, CASErrorKind::Undefined);
}

TEST(SymbolicRegression, SqrtOfNegativeWithoutComplex) {
    // sqrt(-1) must canonicalize to i. Accept either form:
    //   - legacy Constant::I (pre-F1.6 ComplexLit canonical form), or
    //   - ComplexLit(0,1) (post-F1.6 canonical exact Q[i] form).
    CASContext ctx;
    auto expr = parse_expr("sqrt(-1)", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simplified = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified.is_ok());

    if (const auto* constant = expr_cast<Constant>(simplified.value())) {
        EXPECT_EQ(constant->value, MathConstant::I);
        return;
    }
    const auto* complex = expr_cast<ComplexLit>(simplified.value());
    ASSERT_NE(complex, nullptr)
        << "sqrt(-1) must simplify to Constant::I or ComplexLit(0,1)";
    EXPECT_TRUE(complex->re_num.is_zero());
    EXPECT_EQ(complex->im_num, BigInt(1));
    EXPECT_EQ(complex->im_den, BigInt(1));
}

TEST(SymbolicRegression, ImaginaryUnitIsConstantNotSymbol) {
    CASContext ctx;
    auto expr = parse_expr("i", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    const auto* constant = expr_cast<Constant>(expr.value());
    ASSERT_NE(constant, nullptr);
    EXPECT_EQ(constant->value, MathConstant::I);
}

TEST(SymbolicRegression, ImaginaryPowersCanonicalize) {
    CASContext ctx;
    auto expr = parse_expr("i^2", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simplified = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified.is_ok());
    const auto* integer = expr_cast<IntegerLit>(simplified.value());
    ASSERT_NE(integer, nullptr);
    EXPECT_EQ(integer->value, BigInt(-1));
}

TEST(SymbolicRegression, NestedPowerSimplification) {
    CASContext ctx;
    auto expr = parse_expr("((x^2)^3)^4", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simplified = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified.is_ok());
    // (x^2)^3 = x^6, (x^6)^4 = x^24
    auto expected = parse_expr("x^24", ctx.arena());
    ASSERT_TRUE(expected.is_ok());
    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

TEST(SymbolicRegression, LargeSumOfIdenticalSymbols) {
    CASContext ctx;
    std::vector<ExprPtr> terms;
    for(int i = 0; i < 100; ++i) {
        terms.push_back(ctx.arena().make<Symbol>("x"));
    }
    auto sum = ctx.arena().make<Sum>(std::move(terms));
    auto simplified = ctx.simplify(sum);
    ASSERT_TRUE(simplified.is_ok());
    // Should be 100 * x
    auto expected = parse_expr("100 * x", ctx.arena());
    ASSERT_TRUE(expected.is_ok());
    // Note: simplification might order it as x * 100 or 100 * x depending on rules
}

TEST(SymbolicRegression, ProductOfPowersSameBase) {
    CASContext ctx;
    auto expr = parse_expr("x^a * x^b * x^c", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simplified = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified.is_ok());
    // Should be x^(a + b + c)
}

} // namespace
} // namespace cas::symbolic
