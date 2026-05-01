#include "cas/lexer.hpp"
#include "cas/parser.hpp"

#include <gtest/gtest.h>

namespace cas {
namespace {

Result<ExprPtr> parse_roundtrip_input(const std::string& input, AstArena& arena) {
    Lexer lexer(input);
    auto tokens = lexer.tokenize();
    if (tokens.is_error()) {
        return Result<ExprPtr>(tokens.error());
    }

    Parser parser(tokens.value(), arena);
    return parser.parse();
}

TEST(ParserRoundTripPrinterTest, RejectsNullExpression) {
    const auto result = to_round_trip_text(ExprPtr{});

    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::InvalidArgument);
}

TEST(ParserRoundTripPrinterTest, EmitsNormalizedProtectedGrouping) {
    AstArena arena;
    const auto parsed = parse_roundtrip_input("(x + y) * z", arena);

    ASSERT_TRUE(parsed.is_ok());
    const auto text = to_round_trip_text(parsed.value());

    ASSERT_TRUE(text.is_ok());
    EXPECT_EQ(text.value(), "(x+y)*z");
}

TEST(ParserRoundTripPrinterTest, EmitsDedicatedSyntaxForSpecialNodes) {
    AstArena arena;
    auto integral = parse_roundtrip_input("int(x^2, x, 0, 1)", arena);
    auto derivative = parse_roundtrip_input("diff(sin(x), x, 2)", arena);
    auto limit = parse_roundtrip_input("lim(sin(x)/x, x, 0, right)", arena);
    auto matrix = parse_roundtrip_input("[[1,2],[3,4]]", arena);

    ASSERT_TRUE(integral.is_ok());
    ASSERT_TRUE(derivative.is_ok());
    ASSERT_TRUE(limit.is_ok());
    ASSERT_TRUE(matrix.is_ok());

    EXPECT_EQ(to_round_trip_text(integral.value()).value(), "int(x^2,x,0,1)");
    EXPECT_EQ(to_round_trip_text(derivative.value()).value(), "diff(sin(x),x,2)");
    EXPECT_EQ(to_round_trip_text(limit.value()).value(), "lim(sin(x)/x,x,0,right)");
    EXPECT_EQ(to_round_trip_text(matrix.value()).value(), "[[1,2],[3,4]]");
}

struct RoundTripCase {
    const char* input;
};

class ParserRoundTripTest : public ::testing::TestWithParam<RoundTripCase> {};

TEST_P(ParserRoundTripTest, ParsePrintParsePreservesStructure) {
    AstArena first_arena;
    AstArena second_arena;
    const auto original = parse_roundtrip_input(GetParam().input, first_arena);

    ASSERT_TRUE(original.is_ok()) << GetParam().input;

    const auto rendered = to_round_trip_text(original.value());
    ASSERT_TRUE(rendered.is_ok()) << GetParam().input;

    const auto reparsed = parse_roundtrip_input(rendered.value(), second_arena);
    ASSERT_TRUE(reparsed.is_ok()) << rendered.value();
    EXPECT_TRUE(structural_equal(original.value(), reparsed.value()))
        << "input=" << GetParam().input << " rendered=" << rendered.value();
}

INSTANTIATE_TEST_SUITE_P(
    RoundTripCore,
    ParserRoundTripTest,
    ::testing::Values(
        RoundTripCase{"x"},
        RoundTripCase{"3/4"},
        RoundTripCase{"3.14"},
        RoundTripCase{"x+y*z"},
        RoundTripCase{"(x+y)*z"},
        RoundTripCase{"2x"},
        RoundTripCase{"x(y+1)"},
        RoundTripCase{"2(x+y)^2"},
        RoundTripCase{"x^y^z"},
        RoundTripCase{"-x^2"},
        RoundTripCase{"sin(x^2+1)/(x-3)"},
        RoundTripCase{"int(x^2,x)"},
        RoundTripCase{"d/dx(sin(x))"},
        RoundTripCase{"diff(sin(x),x,2)"},
        RoundTripCase{"lim(sin(x)/x,x,0,left)"},
        RoundTripCase{"RootOf(x^2-2,x,1)"},
        RoundTripCase{"[[1,2],[3,4]]"}));

}  // namespace
}  // namespace cas
