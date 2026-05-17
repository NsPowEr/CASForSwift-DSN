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

class BesselOrthogonalityTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

TEST_F(BesselOrthogonalityTest, CrossTermVanishesNu0) {
    // ∫_0^1 x · J_0(j_{0,1}·x) · J_0(j_{0,2}·x) dx = 0
    auto integrand = parse_expr(
        "x * bessel_j(0, bessel_zero(0, 1) * x) * bessel_j(0, bessel_zero(0, 2) * x)",
        ctx.arena());
    ASSERT_TRUE(integrand.is_ok()) << integrand.error().message;

    auto lower = parse_expr("0", ctx.arena());
    auto upper = parse_expr("1", ctx.arena());
    ASSERT_TRUE(lower.is_ok());
    ASSERT_TRUE(upper.is_ok());

    auto result = calculus::definite_integral(
        integrand.value(), Symbol("x"), lower.value(), upper.value(), ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_TRUE(is_integer_zero(result.value()));
}

TEST_F(BesselOrthogonalityTest, DiagonalNormSquaredNu0) {
    // ∫_0^1 x · J_0(j_{0,m}·x)^2 dx = (1/2) · J_1(j_{0,m})^2
    auto integrand = parse_expr(
        "x * bessel_j(0, bessel_zero(0, 3) * x) * bessel_j(0, bessel_zero(0, 3) * x)",
        ctx.arena());
    ASSERT_TRUE(integrand.is_ok()) << integrand.error().message;

    auto lower = parse_expr("0", ctx.arena());
    auto upper = parse_expr("1", ctx.arena());
    auto result = calculus::definite_integral(
        integrand.value(), Symbol("x"), lower.value(), upper.value(), ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    // Expected: (1/2) * J_1(j_{0,3})^2
    auto expected = parse_expr(
        "(1/2) * bessel_j(1, bessel_zero(0, 3))^2",
        ctx.arena());
    ASSERT_TRUE(expected.is_ok());
    auto expected_simpl = ctx.simplify(expected.value());
    ASSERT_TRUE(expected_simpl.is_ok());

    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, result.value(), expected_simpl.value());
    auto delta_simpl = ctx.simplify(delta);
    ASSERT_TRUE(delta_simpl.is_ok());
    EXPECT_TRUE(is_integer_zero(delta_simpl.value()))
        << "expected zero difference, got: "
        << "result and expected differ";
}

TEST_F(BesselOrthogonalityTest, DiagonalNuOne) {
    // Anti-hardcode: ν=1 must also be detected.
    // ∫_0^1 x · J_1(j_{1,2}·x)^2 dx = (1/2) · J_2(j_{1,2})^2
    auto integrand = parse_expr(
        "x * bessel_j(1, bessel_zero(1, 2) * x) * bessel_j(1, bessel_zero(1, 2) * x)",
        ctx.arena());
    ASSERT_TRUE(integrand.is_ok());

    auto lower = parse_expr("0", ctx.arena());
    auto upper = parse_expr("1", ctx.arena());
    auto result = calculus::definite_integral(
        integrand.value(), Symbol("x"), lower.value(), upper.value(), ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    auto expected = parse_expr(
        "(1/2) * bessel_j(2, bessel_zero(1, 2))^2",
        ctx.arena());
    ASSERT_TRUE(expected.is_ok());
    auto expected_simpl = ctx.simplify(expected.value());
    ASSERT_TRUE(expected_simpl.is_ok());

    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, result.value(), expected_simpl.value());
    auto delta_simpl = ctx.simplify(delta);
    ASSERT_TRUE(delta_simpl.is_ok());
    EXPECT_TRUE(is_integer_zero(delta_simpl.value()));
}

TEST_F(BesselOrthogonalityTest, SymbolicScaleFiveByFive) {
    // 5x5 orthogonality matrix on [0,1] for ν=0.
    auto lower = parse_expr("0", ctx.arena()).value();
    auto upper = parse_expr("1", ctx.arena()).value();
    for (int m = 1; m <= 5; ++m) {
        for (int k = 1; k <= 5; ++k) {
            const std::string s = "x * bessel_j(0, bessel_zero(0, " + std::to_string(m)
                + ") * x) * bessel_j(0, bessel_zero(0, " + std::to_string(k) + ") * x)";
            auto e = parse_expr(s, ctx.arena());
            ASSERT_TRUE(e.is_ok());
            auto res = calculus::definite_integral(e.value(), Symbol("x"), lower, upper, ctx);
            ASSERT_TRUE(res.is_ok()) << "m=" << m << " k=" << k << " err=" << res.error().message;
            if (m != k) {
                EXPECT_TRUE(is_integer_zero(res.value())) << "m=" << m << " k=" << k;
            }
        }
    }
}

TEST_F(BesselOrthogonalityTest, DoesNotMatchUnrelatedIntegrand) {
    // Anti-hardcode negative: x*J_0(α·x)*J_0(β·x) with α not equal to BesselZero/upper
    // must NOT trigger the pattern; falls back to the generic indefinite path.
    auto integrand = parse_expr(
        "x * bessel_j(0, 2 * x) * bessel_j(0, 3 * x)",
        ctx.arena());
    ASSERT_TRUE(integrand.is_ok());

    auto lower = parse_expr("0", ctx.arena()).value();
    auto upper = parse_expr("1", ctx.arena()).value();
    auto result = calculus::definite_integral(integrand.value(), Symbol("x"), lower, upper, ctx);
    // Either the indefinite fallback produces something (we don't assert correctness),
    // or it returns Unimplemented.  The only forbidden outcome is returning the
    // orthogonality answer (0 or J^2/2 expression) without justification.
    if (result.is_ok()) {
        EXPECT_FALSE(is_integer_zero(result.value()))
            << "Pattern incorrectly fired on non-orthogonality integrand";
    }
}

}  // namespace cas::test
