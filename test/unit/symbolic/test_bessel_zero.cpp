#include "cas/calculus.hpp"
#include "cas/formatter.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

namespace cas::test {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

}  // namespace

class BesselZeroTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

TEST_F(BesselZeroTest, ParseBuildsDedicatedBuiltin) {
    auto expr = parse_expr("bessel_zero(0, 1)", ctx.arena());
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;
    const auto* call = expr_cast<FuncCall>(expr.value());
    ASSERT_NE(call, nullptr);
    EXPECT_EQ(call->func_id, BuiltinOp::BesselZero);
    ASSERT_EQ(call->args.size(), 2U);
}

TEST_F(BesselZeroTest, DerivativeIsZero) {
    auto expr = parse_expr("bessel_zero(0, 1)", ctx.arena());
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto deriv = calculus::diff(expr.value(), Symbol("x"), 1U, ctx);
    ASSERT_TRUE(deriv.is_ok()) << deriv.error().message;

    auto simplified = ctx.simplify(deriv.value());
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;
    const auto* integer = expr_cast<IntegerLit>(simplified.value());
    ASSERT_NE(integer, nullptr);
    EXPECT_EQ(integer->value, BigInt(0));
}

TEST_F(BesselZeroTest, TextAndLatexFormattingRoundTrip) {
    auto expr = parse_expr("bessel_zero(0, 1)", ctx.arena());
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    formatter::TextFormatter text;
    EXPECT_EQ(text.format(expr.value()), "BesselZero(0, 1)");

    formatter::LaTeXFormatter latex;
    EXPECT_EQ(latex.format(expr.value()), "j_{0,1}");

    auto reparsed = parse_expr(text.format(expr.value()), ctx.arena());
    ASSERT_TRUE(reparsed.is_ok()) << reparsed.error().message;
    const auto* call = expr_cast<FuncCall>(reparsed.value());
    ASSERT_NE(call, nullptr);
    EXPECT_EQ(call->func_id, BuiltinOp::BesselZero);
}

TEST_F(BesselZeroTest, DistinctIndicesRemainStructurallyDistinct) {
    auto first = parse_expr("bessel_zero(nu, 1)", ctx.arena());
    auto second = parse_expr("bessel_zero(nu, 2)", ctx.arena());
    ASSERT_TRUE(first.is_ok()) << first.error().message;
    ASSERT_TRUE(second.is_ok()) << second.error().message;
    EXPECT_FALSE(structural_equal(first.value(), second.value()));
}

}  // namespace cas::test
