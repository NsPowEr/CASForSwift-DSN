#include <gtest/gtest.h>
#include "cas/bigfloat.hpp"
#include "cas/numeric.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::numeric;

// ── BigFloat unit tests ───────────────────────────────────────────────────────

TEST(BigFloatTest, ZeroConstruction) {
    BigFloat z(128);
    EXPECT_TRUE(z.is_zero());
    EXPECT_EQ(z.to_string(5), "0");
}

TEST(BigFloatTest, PiHighPrecision) {
    // Pi to 50 decimal digits — compare known leading digits
    BigFloat pi = BigFloat::pi(decimal_digits_to_bits(50));
    std::string s = pi.to_string(50);
    EXPECT_TRUE(s.find("3.14159265358979") != std::string::npos)
        << "Pi string: " << s;
}

TEST(BigFloatTest, EConstruction) {
    BigFloat e = BigFloat::e(decimal_digits_to_bits(30));
    std::string s = e.to_string(20);
    EXPECT_TRUE(s.find("2.718281828") != std::string::npos)
        << "e string: " << s;
}

TEST(BigFloatTest, ArithmeticAdd) {
    BigFloat a = BigFloat::from_double(1.0, 128);
    BigFloat b = BigFloat::from_double(2.0, 128);
    BigFloat c = a + b;
    EXPECT_NEAR(c.to_double(), 3.0, 1e-15);
}

TEST(BigFloatTest, ArithmeticMul) {
    BigFloat a = BigFloat::from_double(3.0, 128);
    BigFloat b = BigFloat::from_double(4.0, 128);
    EXPECT_NEAR((a * b).to_double(), 12.0, 1e-14);
}

TEST(BigFloatTest, ArithmeticDiv) {
    BigFloat a = BigFloat::from_double(1.0, 128);
    BigFloat b = BigFloat::from_double(3.0, 128);
    // 1/3 to 20 digits
    std::string s = (a / b).to_string(20);
    EXPECT_TRUE(s.find("0.333333333") != std::string::npos)
        << "1/3 string: " << s;
}

TEST(BigFloatTest, SqrtTwo) {
    BigFloat two = BigFloat::from_double(2.0, decimal_digits_to_bits(40));
    BigFloat r = BigFloat::sqrt(two);
    std::string s = r.to_string(30);
    EXPECT_TRUE(s.find("1.41421356237") != std::string::npos)
        << "sqrt(2) string: " << s;
}

TEST(BigFloatTest, SinPiOver2) {
    BigFloat pi = BigFloat::pi(decimal_digits_to_bits(30));
    BigFloat two = BigFloat::from_double(2.0, pi.precision_bits());
    BigFloat half_pi = pi / two;
    BigFloat s = BigFloat::sin(half_pi);
    EXPECT_NEAR(s.to_double(), 1.0, 1e-15);
}

TEST(BigFloatTest, ExpLnInverse) {
    BigFloat x = BigFloat::from_double(2.5, 128);
    BigFloat result = BigFloat::ln(BigFloat::exp(x));
    EXPECT_NEAR(result.to_double(), 2.5, 1e-14);
}

TEST(BigFloatTest, Comparisons) {
    BigFloat a = BigFloat::from_double(1.0, 128);
    BigFloat b = BigFloat::from_double(2.0, 128);
    EXPECT_TRUE(a < b);
    EXPECT_TRUE(b > a);
    EXPECT_FALSE(a == b);
    EXPECT_TRUE(a != b);
}

TEST(BigFloatTest, IsNegative) {
    BigFloat a = BigFloat::from_double(-3.0, 128);
    EXPECT_TRUE(a.is_negative());
    EXPECT_FALSE((-a).is_negative());
}

TEST(BigFloatTest, FromIntegerString) {
    // Large integer: 10^30
    std::string big = "1" + std::string(30, '0');
    BigFloat f = BigFloat::from_integer_string(big, decimal_digits_to_bits(35));
    std::string s = f.to_string(5);
    EXPECT_TRUE(s.find("1e+30") != std::string::npos || s.find("1.0e+30") != std::string::npos
        || s.find("1e30") != std::string::npos || s.find("1.e+30") != std::string::npos)
        << "10^30 string: " << s;
}

TEST(BigFloatTest, FromRationalParts) {
    // 1/7 — repeating decimal
    BigFloat f = BigFloat::from_rational_parts("1", "7", decimal_digits_to_bits(25));
    std::string s = f.to_string(20);
    EXPECT_TRUE(s.find("0.14285714285") != std::string::npos)
        << "1/7 string: " << s;
}

TEST(BigFloatTest, EulerGamma) {
    BigFloat gamma = BigFloat::euler_gamma(decimal_digits_to_bits(20));
    std::string s = gamma.to_string(15);
    EXPECT_TRUE(s.find("0.5772156649") != std::string::npos)
        << "Euler-Mascheroni string: " << s;
}

// ── eval_mpfr integration tests ───────────────────────────────────────────────

TEST(EvalMpfrTest, PiDigits50) {
    symbolic::CASContext ctx;
    auto pi_expr = ctx.arena().make<Constant>(MathConstant::Pi);
    auto res = eval_mpfr(pi_expr, 50);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    const std::string& s = res.value();
    // Verify first 20 digits of pi
    EXPECT_TRUE(s.find("3.14159265358979") != std::string::npos)
        << "pi@50: " << s;
    // String should be at least 20 chars long
    EXPECT_GE(s.size(), 20U);
}

TEST(EvalMpfrTest, SimpleArithmetic) {
    symbolic::CASContext ctx;
    // 1/3 + 2/3 = 1
    auto one_third = ctx.arena().make<RationalLit>(BigInt(1), BigInt(3));
    auto two_thirds = ctx.arena().make<RationalLit>(BigInt(2), BigInt(3));
    auto expr = ctx.arena().make<Sum>(std::vector<ExprPtr>{one_third, two_thirds});
    auto res = eval_mpfr(expr, 30);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(res.value().find("1") != std::string::npos)
        << "1/3+2/3: " << res.value();
}

TEST(EvalMpfrTest, SqrtTwo30Digits) {
    symbolic::CASContext ctx;
    auto two = ctx.arena().make<IntegerLit>(BigInt(2));
    auto expr = ctx.arena().make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{two});
    auto res = eval_mpfr(expr, 30);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(res.value().find("1.41421356237") != std::string::npos)
        << "sqrt(2)@30: " << res.value();
}

TEST(EvalMpfrTest, NBuiltinTwoArgs) {
    // Test N(pi, 50) via ctx.simplify
    symbolic::CASContext ctx;
    auto pi_expr = ctx.arena().make<Constant>(MathConstant::Pi);
    auto digits_expr = ctx.arena().make<IntegerLit>(BigInt(50));
    auto n_call = ctx.arena().make<FuncCall>(BuiltinOp::N,
        std::vector<ExprPtr>{pi_expr, digits_expr});
    auto res = ctx.simplify(n_call);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    const auto* dl = expr_cast<DecimalLit>(res.value());
    ASSERT_NE(dl, nullptr) << "Expected DecimalLit result";
    EXPECT_TRUE(dl->text.find("3.14159265358979") != std::string::npos)
        << "N(pi,50): " << dl->text;
}

TEST(EvalMpfrTest, NBuiltinLargeDigits) {
    // N(sqrt(2), 100) — 100 decimal digits
    symbolic::CASContext ctx;
    auto two = ctx.arena().make<IntegerLit>(BigInt(2));
    auto sqrt2 = ctx.arena().make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{two});
    auto digits = ctx.arena().make<IntegerLit>(BigInt(100));
    auto n_call = ctx.arena().make<FuncCall>(BuiltinOp::N,
        std::vector<ExprPtr>{sqrt2, digits});
    auto res = ctx.simplify(n_call);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    const auto* dl = expr_cast<DecimalLit>(res.value());
    ASSERT_NE(dl, nullptr);
    // Must start with 1.4142...
    EXPECT_TRUE(dl->text.find("1.4142135623730950488") != std::string::npos)
        << "N(sqrt(2),100): " << dl->text;
}
