#include "cas/lexer.hpp"
#include "cas/parser.hpp"

#include <gtest/gtest.h>

namespace cas {
namespace {

Result<ExprPtr> parse_input(const std::string& input, AstArena& arena);

struct ParserRootCase {
    const char* input;
    ExprKind root_kind;
};

class ParserRootKindTest : public ::testing::TestWithParam<ParserRootCase> {};

TEST_P(ParserRootKindTest, ParsesExpectedRootKind) {
    AstArena arena;
    const auto result = parse_input(GetParam().input, arena);

    ASSERT_TRUE(result.is_ok()) << GetParam().input;
    EXPECT_EQ(expr_kind(result.value()), GetParam().root_kind) << GetParam().input;
}

INSTANTIATE_TEST_SUITE_P(
    RootKinds,
    ParserRootKindTest,
    ::testing::Values(
        ParserRootCase{"0", ExprKind::IntegerLit},
        ParserRootCase{"42", ExprKind::IntegerLit},
        ParserRootCase{"3/4", ExprKind::RationalLit},
        ParserRootCase{"3.14", ExprKind::RationalLit},
        ParserRootCase{"x", ExprKind::Symbol},
        ParserRootCase{"pi", ExprKind::Constant},
        ParserRootCase{"∞", ExprKind::Constant},
        ParserRootCase{"-x", ExprKind::Unary},
        ParserRootCase{"x+y", ExprKind::Sum},
        ParserRootCase{"x*y", ExprKind::Product},
        ParserRootCase{"x-y", ExprKind::Binary},
        ParserRootCase{"x/y", ExprKind::Binary},
        ParserRootCase{"x^y", ExprKind::Binary},
        ParserRootCase{"x%2", ExprKind::Binary},
        ParserRootCase{"x!", ExprKind::Unary},
        ParserRootCase{"sin(x)", ExprKind::FuncCall},
        ParserRootCase{"foo(x,y,z)", ExprKind::FuncCall},
        ParserRootCase{"√(x)", ExprKind::FuncCall},
        ParserRootCase{"∫(x^2, x)", ExprKind::Integral},
        ParserRootCase{"diff(x^2, x)", ExprKind::Derivative},
        ParserRootCase{"d/dx(x^2)", ExprKind::Derivative},
        ParserRootCase{"lim(x, x, 0)", ExprKind::Limit},
        ParserRootCase{"RootOf(x^2 - 2, x)", ExprKind::RootOf},
        ParserRootCase{"[[1,2],[3,4]]", ExprKind::Matrix},
        ParserRootCase{"(x)", ExprKind::Symbol},
        ParserRootCase{"(x+1)", ExprKind::Sum},
        ParserRootCase{"((x))", ExprKind::Symbol},
        ParserRootCase{"x+y+z", ExprKind::Sum},
        ParserRootCase{"x*y*z", ExprKind::Product},
        ParserRootCase{"(x+y)*z", ExprKind::Product}));

class ParserMalformedInputTest : public ::testing::TestWithParam<const char*> {};

TEST_P(ParserMalformedInputTest, RejectsMalformedInputWithoutCrashing) {
    AstArena arena;
    const auto result = parse_input(GetParam(), arena);

    ASSERT_TRUE(result.is_error()) << GetParam();
    EXPECT_EQ(result.error().kind, CASErrorKind::ParseError);
}

INSTANTIATE_TEST_SUITE_P(
    MalformedInputs,
    ParserMalformedInputTest,
    ::testing::Values(
        "",
        "(",
        ")",
        "x + ",
        "sin(",
        "sin(,)",
        "foo(",
        "foo(x,)",
        "diff(",
        "diff(x",
        "diff(x,)",
        "diff(x, x, )",
        "lim(",
        "lim(x,)",
        "lim(x, x)",
        "[",
        "[[1,2]",
        "[[1],[2,3]]",
        "RootOf(x^2, 1)",
        "RootOf(x^2, x, -1)",
        "RootOf(x^2, x, 1, 2)",
        "∫(",
        "∫(x^2)",
        "∫(x^2, 1)",
        "x + * y",
        "x ^ ",
        "!",
        ",",
        "√",
        "d/d()"));

Result<ExprPtr> parse_input(const std::string& input, AstArena& arena) {
    Lexer lexer(input);
    auto tokens = lexer.tokenize();
    if (tokens.is_error()) {
        return Result<ExprPtr>(tokens.error());
    }

    Parser parser(tokens.value(), arena);
    return parser.parse();
}

TEST(ParserTest, ParsesSingleSymbol) {
    AstArena arena;
    auto result = parse_input("x", arena);

    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::Symbol);
    EXPECT_EQ(expr_ref<Symbol>(result.value()).name, "x");
}

TEST(ParserTest, ParsesDecimalAsRationalLiteral) {
    AstArena arena;
    auto result = parse_input("3.14", arena);

    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::RationalLit);
    const auto& rat = expr_ref<RationalLit>(result.value());
    EXPECT_EQ(rat.numerator.decimal(), "157");
    EXPECT_EQ(rat.denominator.decimal(), "50");
}

TEST(ParserTest, ParsesHalfDecimalAsRational) {
    AstArena arena;
    auto result = parse_input("0.5", arena);
    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::RationalLit);
    const auto& rat = expr_ref<RationalLit>(result.value());
    EXPECT_EQ(rat.numerator.decimal(), "1");
    EXPECT_EQ(rat.denominator.decimal(), "2");
}

TEST(ParserTest, ParsesQuarterDecimalAsRational) {
    AstArena arena;
    auto result = parse_input("0.25", arena);
    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::RationalLit);
    const auto& rat = expr_ref<RationalLit>(result.value());
    EXPECT_EQ(rat.numerator.decimal(), "1");
    EXPECT_EQ(rat.denominator.decimal(), "4");
}

TEST(ParserTest, ParsesOnePointZeroAsRational) {
    AstArena arena;
    auto result = parse_input("1.0", arena);
    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::RationalLit);
    const auto& rat = expr_ref<RationalLit>(result.value());
    EXPECT_EQ(rat.numerator.decimal(), "1");
    EXPECT_EQ(rat.denominator.decimal(), "1");
}

TEST(ParserTest, ParsesRationalLiteral) {
    AstArena arena;
    auto result = parse_input("3/4", arena);

    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::RationalLit);
    const auto& rational = expr_ref<RationalLit>(result.value());
    EXPECT_EQ(rational.numerator.decimal(), "3");
    EXPECT_EQ(rational.denominator.decimal(), "4");
}

TEST(ParserTest, PreservesOperatorPrecedence) {
    AstArena arena;
    auto result = parse_input("x + y * z", arena);

    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::Sum);
    const auto& sum = expr_ref<Sum>(result.value());
    ASSERT_EQ(sum.terms.size(), 2U);
    EXPECT_EQ(expr_kind(sum.terms[0]), ExprKind::Symbol);
    ASSERT_EQ(expr_kind(sum.terms[1]), ExprKind::Product);
    EXPECT_EQ(expr_ref<Symbol>(sum.terms[0]).name, "x");
    const auto& product = expr_ref<Product>(sum.terms[1]);
    ASSERT_EQ(product.factors.size(), 2U);
    EXPECT_EQ(expr_ref<Symbol>(product.factors[0]).name, "y");
    EXPECT_EQ(expr_ref<Symbol>(product.factors[1]).name, "z");
}

TEST(ParserTest, ParsesPowerAsRightAssociative) {
    AstArena arena;
    auto result = parse_input("x^y^z", arena);

    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::Binary);
    const auto& outer = expr_ref<Binary>(result.value());
    EXPECT_EQ(outer.op, BinaryOp::Pow);
    ASSERT_EQ(expr_kind(outer.right), ExprKind::Binary);
    EXPECT_EQ(expr_ref<Binary>(outer.right).op, BinaryOp::Pow);
}

TEST(ParserTest, ParsesFunctionCall) {
    AstArena arena;
    auto result = parse_input("sin(x^2 + 1)", arena);

    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::FuncCall);
    const auto& call = expr_ref<FuncCall>(result.value());
    EXPECT_EQ(call.name, "sin");
    ASSERT_EQ(call.args.size(), 1U);
    ASSERT_EQ(expr_kind(call.args[0]), ExprKind::Sum);
    const auto& sum = expr_ref<Sum>(call.args[0]);
    ASSERT_EQ(sum.terms.size(), 2U);
    ASSERT_EQ(expr_kind(sum.terms[0]), ExprKind::Binary);
    const auto& power = expr_ref<Binary>(sum.terms[0]);
    EXPECT_EQ(power.op, BinaryOp::Pow);
    EXPECT_EQ(expr_ref<Symbol>(power.left).name, "x");
    EXPECT_EQ(expr_ref<IntegerLit>(power.right).value.decimal(), "2");
    EXPECT_EQ(expr_ref<IntegerLit>(sum.terms[1]).value.decimal(), "1");
}

TEST(ParserTest, ParsesUnaryMinus) {
    AstArena arena;
    auto result = parse_input("-x", arena);

    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::Unary);
    EXPECT_EQ(expr_ref<Unary>(result.value()).op, UnaryOp::Neg);
}

TEST(ParserTest, ParsesProductOfGroupedExpressions) {
    AstArena arena;
    auto result = parse_input("(x + 1) * (x - 1)", arena);

    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::Product);
    const auto& product = expr_ref<Product>(result.value());
    ASSERT_EQ(product.factors.size(), 2U);
    ASSERT_EQ(expr_kind(product.factors[0]), ExprKind::Sum);
    ASSERT_EQ(expr_kind(product.factors[1]), ExprKind::Binary);
    const auto& left_sum = expr_ref<Sum>(product.factors[0]);
    ASSERT_EQ(left_sum.terms.size(), 2U);
    EXPECT_EQ(expr_ref<Symbol>(left_sum.terms[0]).name, "x");
    EXPECT_EQ(expr_ref<IntegerLit>(left_sum.terms[1]).value.decimal(), "1");
    const auto& right_sub = expr_ref<Binary>(product.factors[1]);
    EXPECT_EQ(right_sub.op, BinaryOp::Sub);
    EXPECT_EQ(expr_ref<Symbol>(right_sub.left).name, "x");
    EXPECT_EQ(expr_ref<IntegerLit>(right_sub.right).value.decimal(), "1");
}

TEST(ParserTest, ParsesIntegralNode) {
    AstArena arena;
    auto result = parse_input("∫(x^2, x)", arena);

    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::Integral);
    EXPECT_EQ(expr_ref<Integral>(result.value()).variable.name, "x");
}

TEST(ParserTest, ParsesDerivativeNode) {
    AstArena arena;
    auto result = parse_input("d/dx(sin(x))", arena);

    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::Derivative);
    EXPECT_EQ(expr_ref<Derivative>(result.value()).variable.name, "x");
}

TEST(ParserTest, ParsesGenericUnknownFunctionAsFuncCall) {
    AstArena arena;
    auto result = parse_input("foo(x, y)", arena);

    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::FuncCall);
    EXPECT_EQ(expr_ref<FuncCall>(result.value()).name, "foo");
}

TEST(ParserTest, ParsesOneDivZeroAsSyntacticDivision) {
    AstArena arena;
    auto result = parse_input("1/0", arena);

    ASSERT_TRUE(result.is_ok());
    ASSERT_EQ(expr_kind(result.value()), ExprKind::Binary);
    EXPECT_EQ(expr_ref<Binary>(result.value()).op, BinaryOp::Div);
}

TEST(ParserTest, RejectsEmptyInput) {
    AstArena arena;
    auto result = parse_input("", arena);

    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::ParseError);
    EXPECT_NE(result.error().message.find("1:1"), std::string::npos);
}

TEST(ParserTest, RejectsUnexpectedEof) {
    AstArena arena;
    auto result = parse_input("x + ", arena);

    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::ParseError);
    EXPECT_NE(result.error().message.find("1:5"), std::string::npos);
}

TEST(ParserTest, RejectsUnmatchedParenthesis) {
    AstArena arena;
    auto result = parse_input("sin(", arena);

    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::ParseError);
}

TEST(ParserTest, RejectsEmptyTokenStreamWithoutCrashing) {
    AstArena arena;
    const std::vector<Token> tokens;
    Parser parser(tokens, arena);

    const auto result = parser.parse();

    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::ParseError);
    EXPECT_NE(result.error().message.find("1:1"), std::string::npos);
}

TEST(ParserTest, RejectsOutOfRangeDerivativeOrder) {
    AstArena arena;
    auto result = parse_input("diff(x, x, 999999999999999999999999999999)", arena);

    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::ParseError);
}

TEST(ParserTest, RejectsOutOfRangeRootOfIndex) {
    AstArena arena;
    auto result = parse_input("RootOf(x^2 - 2, x, 999999999999999999999999999999)", arena);

    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::ParseError);
}

}  // namespace
}  // namespace cas
