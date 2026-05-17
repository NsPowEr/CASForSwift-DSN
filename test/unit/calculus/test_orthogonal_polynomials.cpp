#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

#include <string>

namespace cas::test {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

[[nodiscard]] bool is_integer_zero(ExprPtr expr) {
    const auto* lit = expr_cast<IntegerLit>(expr);
    return lit != nullptr && lit->value.is_zero();
}

}  // namespace

class OrthogonalPolynomialsTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    void expect_zero_or_match(ExprPtr result, const std::string& expected_src) {
        auto expected = parse_expr(expected_src, ctx.arena());
        ASSERT_TRUE(expected.is_ok());
        auto exp_s = ctx.simplify(expected.value());
        ASSERT_TRUE(exp_s.is_ok());
        auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, result, exp_s.value());
        auto delta_s = ctx.simplify(delta);
        ASSERT_TRUE(delta_s.is_ok());
        EXPECT_TRUE(is_integer_zero(delta_s.value()))
            << "result and expected differ; got expression that did not reduce to zero";
    }
};

// ─── Legendre P_n on [-1, 1] ────────────────────────────────────────────────

TEST_F(OrthogonalPolynomialsTest, LegendreCrossOrthogonality) {
    // ∫_{-1}^{1} P_1(x)·P_3(x) dx = 0.
    auto e = parse_expr("LegendreP(1, x) * LegendreP(3, x)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto lower = parse_expr("-1", ctx.arena()).value();
    auto upper = parse_expr("1", ctx.arena()).value();
    auto res = calculus::definite_integral(e.value(), Symbol("x"), lower, upper, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(is_integer_zero(res.value()));
}

TEST_F(OrthogonalPolynomialsTest, LegendreDiagonalNorm) {
    // ∫_{-1}^{1} P_2(x)² dx = 2/(2·2+1) = 2/5.
    auto e = parse_expr("LegendreP(2, x) * LegendreP(2, x)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto lower = parse_expr("-1", ctx.arena()).value();
    auto upper = parse_expr("1", ctx.arena()).value();
    auto res = calculus::definite_integral(e.value(), Symbol("x"), lower, upper, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_zero_or_match(res.value(), "2/5");
}

TEST_F(OrthogonalPolynomialsTest, LegendreAntiHardcodeHighIndex) {
    // ∫_{-1}^{1} P_7(x)² dx = 2/(2·7+1) = 2/15.  Anti-hardcode: index > 5.
    auto e = parse_expr("LegendreP(7, x) * LegendreP(7, x)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto lower = parse_expr("-1", ctx.arena()).value();
    auto upper = parse_expr("1", ctx.arena()).value();
    auto res = calculus::definite_integral(e.value(), Symbol("x"), lower, upper, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_zero_or_match(res.value(), "2/15");
}

// ─── Hermite physicist H_n on R, weight e^{-x²} ─────────────────────────────

TEST_F(OrthogonalPolynomialsTest, HermiteHCross) {
    // ∫_{-∞}^{∞} H_1(x)·H_2(x)·e^{-x²} dx = 0.
    auto e = parse_expr("HermiteH(1, x) * HermiteH(2, x) * exp(-(x^2))", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto lower = parse_expr("-inf", ctx.arena()).value();
    auto upper = parse_expr("inf", ctx.arena()).value();
    auto res = calculus::definite_integral(e.value(), Symbol("x"), lower, upper, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(is_integer_zero(res.value()));
}

TEST_F(OrthogonalPolynomialsTest, HermiteHDiagonal) {
    // ∫_{-∞}^{∞} H_3(x)²·e^{-x²} dx = 2³·3!·√π = 8·6·√π = 48·√π.
    auto e = parse_expr("HermiteH(3, x) * HermiteH(3, x) * exp(-(x^2))", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto lower = parse_expr("-inf", ctx.arena()).value();
    auto upper = parse_expr("inf", ctx.arena()).value();
    auto res = calculus::definite_integral(e.value(), Symbol("x"), lower, upper, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_zero_or_match(res.value(), "48 * sqrt(pi)");
}

// ─── Hermite probabilist He_n on R, weight e^{-x²/2} ────────────────────────

TEST_F(OrthogonalPolynomialsTest, HermiteHeCross) {
    // ∫_{-∞}^{∞} He_2(x)·He_4(x)·e^{-x²/2} dx = 0.
    auto e = parse_expr("HermiteHe(2, x) * HermiteHe(4, x) * exp(-(x^2)/2)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto lower = parse_expr("-inf", ctx.arena()).value();
    auto upper = parse_expr("inf", ctx.arena()).value();
    auto res = calculus::definite_integral(e.value(), Symbol("x"), lower, upper, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_TRUE(is_integer_zero(res.value()));
}

TEST_F(OrthogonalPolynomialsTest, HermiteHeDiagonalAntiHardcode) {
    // ∫_{-∞}^{∞} He_4(x)²·e^{-x²/2} dx = 4!·√(2π) = 24·√(2π).
    auto e = parse_expr("HermiteHe(4, x) * HermiteHe(4, x) * exp(-(x^2)/2)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto lower = parse_expr("-inf", ctx.arena()).value();
    auto upper = parse_expr("inf", ctx.arena()).value();
    auto res = calculus::definite_integral(e.value(), Symbol("x"), lower, upper, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_zero_or_match(res.value(), "24 * sqrt(2*pi)");
}

TEST_F(OrthogonalPolynomialsTest, LegendreWrongIntervalNotMatched) {
    // ∫_{0}^{1} P_2(x)² dx is NOT 2/5 — pattern must NOT fire on the wrong
    // interval.  Either succeeds with the correct restricted-domain value
    // (computed via the indefinite integral) or returns Unimplemented,
    // but must never produce the orthogonality answer 2/5.
    auto e = parse_expr("LegendreP(2, x) * LegendreP(2, x)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto lower = parse_expr("0", ctx.arena()).value();
    auto upper = parse_expr("1", ctx.arena()).value();
    auto res = calculus::definite_integral(e.value(), Symbol("x"), lower, upper, ctx);
    if (res.is_ok()) {
        // The correct value on [0,1] for P_2² is 1/5 — but mainly check it is NOT 2/5.
        auto two_over_five = parse_expr("2/5", ctx.arena()).value();
        auto two_over_five_s = ctx.simplify(two_over_five).value();
        auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, res.value(), two_over_five_s);
        auto delta_s = ctx.simplify(delta);
        ASSERT_TRUE(delta_s.is_ok());
        EXPECT_FALSE(is_integer_zero(delta_s.value()))
            << "orthogonality pattern fired on wrong interval";
    }
}

}  // namespace cas::test
