#include "cas/bigint.hpp"
#include "cas/error.hpp"
#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <limits>

namespace cas {
namespace {

[[nodiscard]] BigInt parse_bigint_or_fail(const char* decimal) {
    auto parsed = BigInt::parse(decimal);
    EXPECT_TRUE(parsed.is_ok()) << parsed.error().message;
    return parsed.is_ok() ? std::move(parsed.value()) : BigInt(0);
}

TEST(FoundationAstTest, ArenaExposesApiVersion) {
    EXPECT_EQ(AstArena::API_VERSION, 1U);
}

TEST(FoundationAstTest, ArenaInternsHotLiteralAndSymbolNodes) {
    AstArena arena;

    const ExprPtr zero_a = arena.make<IntegerLit>(BigInt(0));
    const ExprPtr zero_b = arena.make<IntegerLit>(BigInt(0));
    const ExprPtr one_a = arena.make<IntegerLit>(BigInt(1));
    const ExprPtr one_b = arena.make<IntegerLit>(BigInt(1));
    const ExprPtr symbol_a = arena.make<Symbol>("x");
    const ExprPtr symbol_b = arena.make<Symbol>(std::string("x"));
    const ExprPtr pi_a = arena.make<Constant>(MathConstant::Pi);
    const ExprPtr pi_b = arena.make<Constant>(MathConstant::Pi);

    EXPECT_EQ(zero_a, zero_b);
    EXPECT_EQ(one_a, one_b);
    EXPECT_EQ(symbol_a, symbol_b);
    EXPECT_EQ(pi_a, pi_b);
}

TEST(FoundationResultTest, ValueResultStoresSuccessState) {
    auto result = ok(42);
    ASSERT_TRUE(result.is_ok());
    EXPECT_EQ(result.value(), 42);
}

TEST(FoundationResultTest, ErrorResultStoresTypedFailure) {
    auto result = fail<int>(CASError{
        .kind = CASErrorKind::InvalidArgument,
        .message = "invalid input",
        .hint = "check caller contract",
    });

    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::InvalidArgument);
    EXPECT_EQ(result.error().message, "invalid input");
    ASSERT_TRUE(result.error().hint.has_value());
    EXPECT_EQ(*result.error().hint, "check caller contract");
}

TEST(FoundationResultVoidTest, VoidResultSupportsSuccessAndFailure) {
    auto success = ok();
    ASSERT_TRUE(success.is_ok());

    auto failure = Result<void>(CASError{
        .kind = CASErrorKind::Timeout,
        .message = "timed out",
        .hint = std::nullopt,
    });

    ASSERT_TRUE(failure.is_error());
    EXPECT_EQ(failure.error().kind, CASErrorKind::Timeout);
}

TEST(FoundationBigIntTest, DefaultConstructedValueIsZero) {
    BigInt value;
    EXPECT_EQ(value.decimal(), "0");
    EXPECT_FALSE(value.is_negative());
}

TEST(FoundationBigIntTest, SignedConstructionPreservesMagnitude) {
    BigInt positive(123456789);
    BigInt negative(-42);

    EXPECT_EQ(positive.decimal(), "123456789");
    EXPECT_FALSE(positive.is_negative());
    EXPECT_EQ(negative.decimal(), "-42");
    EXPECT_TRUE(negative.is_negative());
}

TEST(FoundationBigIntTest, FromU64BuildsNativeMagnitudeWithoutDecimalParsing) {
    const BigInt value = BigInt::from_u64(std::numeric_limits<std::uint64_t>::max());

    EXPECT_EQ(value.decimal(), "18446744073709551615");
    EXPECT_FALSE(value.is_negative());
}

TEST(FoundationBigIntTest, DecimalConstructionCanonicalizesLeadingZeroes) {
    BigInt value = parse_bigint_or_fail("-000123");

    EXPECT_EQ(value.decimal(), "-123");
    EXPECT_TRUE(value.is_negative());
}

TEST(FoundationBigIntTest, SupportsSignedArithmeticAndDivisionRemainder) {
    const BigInt a = parse_bigint_or_fail("-12345678901234567890");
    const BigInt b = parse_bigint_or_fail("9876543210");

    EXPECT_EQ((a + b).decimal(), "-12345678891358024680");
    EXPECT_TRUE((a + b).is_negative());
    EXPECT_EQ((a - b).decimal(), "-12345678911111111100");
    EXPECT_TRUE((a - b).is_negative());
    EXPECT_EQ((a * b).decimal(), "-121932631124828532111263526900");
    EXPECT_TRUE((a * b).is_negative());
    EXPECT_EQ((a / b).decimal(), "-1249999989");
    EXPECT_TRUE((a / b).is_negative());
    EXPECT_EQ((a % b).decimal(), "2623456800");
    EXPECT_FALSE((a % b).is_negative());
}

TEST(FoundationBigIntTest, GcdUsesCanonicalPositiveMagnitude) {
    const BigInt left = parse_bigint_or_fail("-48");
    const BigInt right = parse_bigint_or_fail("180");

    const BigInt result = gcd(left, right);

    EXPECT_EQ(result.decimal(), "12");
    EXPECT_FALSE(result.is_negative());
}

TEST(FoundationBigIntTest, ParseReturnsTypedErrorForMalformedDecimalInputs) {
    const auto empty = BigInt::parse("");
    const auto invalid = BigInt::parse("12a");

    ASSERT_TRUE(empty.is_error());
    EXPECT_EQ(empty.error().kind, CASErrorKind::InvalidArgument);
    ASSERT_TRUE(invalid.is_error());
    EXPECT_EQ(invalid.error().kind, CASErrorKind::InvalidArgument);
}

TEST(FoundationBigIntTest, CheckedDivisionAndModuloReturnTypedErrors) {
    const BigInt value = parse_bigint_or_fail("42");
    const BigInt zero = parse_bigint_or_fail("0");

    const auto division = checked_divide(value, zero);
    const auto modulo = checked_mod(value, zero);

    ASSERT_TRUE(division.is_error());
    EXPECT_EQ(division.error().kind, CASErrorKind::Undefined);
    ASSERT_TRUE(modulo.is_error());
    EXPECT_EQ(modulo.error().kind, CASErrorKind::Undefined);
}

TEST(FoundationBigIntTest, CheckedDivisionWithRemainderReturnsTypedErrorForZeroDivisor) {
    const BigInt value = parse_bigint_or_fail("42");
    const BigInt zero = parse_bigint_or_fail("0");

    const auto division = checked_divide_with_remainder(value, zero);

    ASSERT_TRUE(division.is_error());
    EXPECT_EQ(division.error().kind, CASErrorKind::Undefined);
}

TEST(FoundationBigIntTest, CheckedDivisionWithRemainderMatchesOperatorContractOnValidInput) {
    const BigInt dividend = parse_bigint_or_fail("-12345678901234567890");
    const BigInt divisor = parse_bigint_or_fail("9876543210");

    const auto division = checked_divide_with_remainder(dividend, divisor);

    ASSERT_TRUE(division.is_ok()) << division.error().message;
    EXPECT_EQ(division.value().first.decimal(), "-1249999989");
    EXPECT_TRUE(division.value().first.is_negative());
    EXPECT_EQ(division.value().second.decimal(), "2623456800");
    EXPECT_FALSE(division.value().second.is_negative());
}

TEST(FoundationRationalTest, CanonicalizesSignAndReducesFraction) {
    const Rational value(parse_bigint_or_fail("6"), parse_bigint_or_fail("-8"));

    EXPECT_EQ(value.numerator().decimal(), "-3");
    EXPECT_TRUE(value.numerator().is_negative());
    EXPECT_EQ(value.denominator().decimal(), "4");
    EXPECT_FALSE(value.denominator().is_negative());
}

TEST(FoundationRationalTest, SupportsExactArithmetic) {
    const Rational lhs(parse_bigint_or_fail("2"), parse_bigint_or_fail("3"));
    const Rational rhs(parse_bigint_or_fail("5"), parse_bigint_or_fail("7"));

    const Rational sum = lhs + rhs;
    const Rational difference = lhs - rhs;
    const Rational product = lhs * rhs;
    const Rational quotient = lhs / rhs;

    EXPECT_EQ(sum.numerator().decimal(), "29");
    EXPECT_EQ(sum.denominator().decimal(), "21");
    EXPECT_EQ(difference.numerator().decimal(), "-1");
    EXPECT_TRUE(difference.numerator().is_negative());
    EXPECT_EQ(difference.denominator().decimal(), "21");
    EXPECT_EQ(product.numerator().decimal(), "10");
    EXPECT_EQ(product.denominator().decimal(), "21");
    EXPECT_EQ(quotient.numerator().decimal(), "14");
    EXPECT_EQ(quotient.denominator().decimal(), "15");
}

TEST(FoundationRationalTest, ZeroNormalizesToUnitDenominator) {
    const Rational value(parse_bigint_or_fail("0"), parse_bigint_or_fail("-999"));

    EXPECT_EQ(value.numerator().decimal(), "0");
    EXPECT_FALSE(value.numerator().is_negative());
    EXPECT_EQ(value.denominator().decimal(), "1");
    EXPECT_TRUE(value.is_integer());
}

TEST(FoundationRationalTest, MakeReturnsTypedErrorForZeroDenominator) {
    const auto result = Rational::make(parse_bigint_or_fail("1"), parse_bigint_or_fail("0"));

    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::Undefined);
}

TEST(FoundationRationalTest, CheckedDivisionReturnsTypedErrorForZeroRational) {
    const Rational lhs(parse_bigint_or_fail("3"), parse_bigint_or_fail("5"));
    const Rational zero(parse_bigint_or_fail("0"), parse_bigint_or_fail("7"));

    const auto result = checked_divide(lhs, zero);

    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::Undefined);
}

TEST(FoundationRationalTest, CheckedDoubleConversionReturnsTypedOverflow) {
    const Rational huge(parse_bigint_or_fail("1" "00000000000000000000000000000000000000000000000000"
                                             "00000000000000000000000000000000000000000000000000"
                                             "00000000000000000000000000000000000000000000000000"
                                             "00000000000000000000000000000000000000000000000000"
                                             "00000000000000000000000000000000000000000000000000"
                                             "00000000000000000000000000000000000000000000000000"
                                             "00000000000000000000000000000000000000000000000000"));

    const auto result = huge.to_double_checked();

    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::Overflow);
}

TEST(FoundationRationalTest, CheckedDoubleConversionHandlesLargeScaleRatio) {
    const Rational tiny(parse_bigint_or_fail("1"),
                        parse_bigint_or_fail("1" "00000000000000000000000000000000000000000000000000"
                                             "00000000000000000000000000000000000000000000000000"));

    const auto result = tiny.to_double_checked();

    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_GT(result.value(), 0.0);
    EXPECT_LT(result.value(), 1.0e-90);
}

}  // namespace
}  // namespace cas
