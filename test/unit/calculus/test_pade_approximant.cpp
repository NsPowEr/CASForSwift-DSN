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

[[nodiscard]] bool reduces_to_zero(ExprPtr actual, ExprPtr expected, symbolic::CASContext& ctx) {
    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, actual, expected);
    auto s = ctx.simplify(delta);
    if (s.is_error()) return false;
    const auto* lit = expr_cast<IntegerLit>(s.value());
    return lit != nullptr && lit->value.is_zero();
}

}  // namespace

class PadeApproximantTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

// ─── exp(x) at 0: well-known Padé table entries ─────────────────────────────

TEST_F(PadeApproximantTest, ExpZeroOnePade) {
    // [1/1] Pade of exp at 0:  (1 + x/2) / (1 − x/2).
    auto expr = parse_expr("exp(x)", ctx.arena()).value();
    auto center = parse_expr("0", ctx.arena()).value();
    auto pade = calculus::pade_approximant(expr, Symbol("x"), center, 1U, 1U, ctx);
    ASSERT_TRUE(pade.is_ok()) << pade.error().message;
    auto expected_num = parse_expr("1 + x/2", ctx.arena()).value();
    auto expected_den = parse_expr("1 - x/2", ctx.arena()).value();
    EXPECT_TRUE(reduces_to_zero(pade.value().numerator, ctx.simplify(expected_num).value(), ctx));
    EXPECT_TRUE(reduces_to_zero(pade.value().denominator, ctx.simplify(expected_den).value(), ctx));
}

TEST_F(PadeApproximantTest, ExpZeroTwoPadeAntiHardcode) {
    // [2/2] Pade of exp at 0: (1 + x/2 + x²/12) / (1 − x/2 + x²/12).
    auto expr = parse_expr("exp(x)", ctx.arena()).value();
    auto center = parse_expr("0", ctx.arena()).value();
    auto pade = calculus::pade_approximant(expr, Symbol("x"), center, 2U, 2U, ctx);
    ASSERT_TRUE(pade.is_ok()) << pade.error().message;
    // Verify via Taylor-truncation identity P − f·Q ≡ 0 mod x^5.
    auto pf = ctx.arena().make<Binary>(BinaryOp::Sub,
        pade.value().numerator,
        ctx.arena().make<Binary>(BinaryOp::Mul, expr, pade.value().denominator));
    auto pf_simp = ctx.simplify(pf);
    ASSERT_TRUE(pf_simp.is_ok());
    for (unsigned int k = 0; k <= 4U; ++k) {
        Result<ExprPtr> deriv = (k == 0U)
            ? ok(pf_simp.value())
            : calculus::diff(pf_simp.value(), Symbol("x"), k, ctx);
        ASSERT_TRUE(deriv.is_ok());
        auto val = ctx.substitute(deriv.value(), Symbol("x"), center);
        ASSERT_TRUE(val.is_ok());
        auto val_s = ctx.simplify(val.value());
        ASSERT_TRUE(val_s.is_ok());
        const auto* lit = expr_cast<IntegerLit>(val_s.value());
        EXPECT_TRUE(lit != nullptr && lit->value.is_zero())
            << "Pade [2/2] of exp failed Taylor-match at order " << k;
    }
}

TEST_F(PadeApproximantTest, GeometricSeriesReproducedExactly) {
    // 1/(1 − x) is itself a [0/1] rational; the [0/1] Pade must reproduce it.
    auto expr = parse_expr("1/(1 - x)", ctx.arena()).value();
    auto center = parse_expr("0", ctx.arena()).value();
    auto pade = calculus::pade_approximant(expr, Symbol("x"), center, 0U, 1U, ctx);
    ASSERT_TRUE(pade.is_ok()) << pade.error().message;
    auto expected_num = parse_expr("1", ctx.arena()).value();
    auto expected_den = parse_expr("1 - x", ctx.arena()).value();
    EXPECT_TRUE(reduces_to_zero(pade.value().numerator, ctx.simplify(expected_num).value(), ctx));
    EXPECT_TRUE(reduces_to_zero(pade.value().denominator, ctx.simplify(expected_den).value(), ctx));
}

TEST_F(PadeApproximantTest, ConsistencyAgainstTaylorTruncation) {
    // For an analytic f, P(x) − f(x)·Q(x) must vanish to order ≥ m + n + 1 at x = 0.
    // Verify by checking that the Taylor coefficients of P − f·Q up to order m+n
    // are all zero.  f = ln(1 + x), [2/2].
    auto expr = parse_expr("ln(1 + x)", ctx.arena()).value();
    auto center = parse_expr("0", ctx.arena()).value();
    auto pade = calculus::pade_approximant(expr, Symbol("x"), center, 2U, 2U, ctx);
    ASSERT_TRUE(pade.is_ok()) << pade.error().message;

    auto pf = ctx.arena().make<Binary>(BinaryOp::Sub,
        pade.value().numerator,
        ctx.arena().make<Binary>(BinaryOp::Mul, expr, pade.value().denominator));
    auto pf_simp = ctx.simplify(pf);
    ASSERT_TRUE(pf_simp.is_ok());
    for (unsigned int k = 0; k <= 4U; ++k) {
        Result<ExprPtr> deriv = (k == 0U)
            ? ok(pf_simp.value())
            : calculus::diff(pf_simp.value(), Symbol("x"), k, ctx);
        ASSERT_TRUE(deriv.is_ok());
        auto val = ctx.substitute(deriv.value(), Symbol("x"), center);
        ASSERT_TRUE(val.is_ok());
        auto val_s = ctx.simplify(val.value());
        ASSERT_TRUE(val_s.is_ok());
        const auto* lit = expr_cast<IntegerLit>(val_s.value());
        EXPECT_TRUE(lit != nullptr && lit->value.is_zero())
            << "Pade defect at order " << k;
    }
}

// ─── F5.5 / B5 — Non-Q symbolic coefficients ────────────────────────────────
//
// pade_approximant used to reject any Taylor coefficient that did not fold to
// IntegerLit/RationalLit, blocking expansion at algebraic / transcendental
// centres.  After the symbolic rewrite the solver works directly over the
// expression field; the tests below cover three orthogonal non-Q cases.

TEST_F(PadeApproximantTest, GeometricSeriesWithAlgebraicCoefficient) {
    // 1/(1 − √2 · x) at 0:  c_k = (√2)^k.  [0/1] Padé must reproduce 1/(1 − √2 x).
    auto expr = parse_expr("1/(1 - sqrt(2)*x)", ctx.arena()).value();
    auto center = parse_expr("0", ctx.arena()).value();
    auto pade = calculus::pade_approximant(expr, Symbol("x"), center, 0U, 1U, ctx);
    ASSERT_TRUE(pade.is_ok()) << pade.error().message;

    auto expected_num = parse_expr("1", ctx.arena()).value();
    auto expected_den = parse_expr("1 - sqrt(2)*x", ctx.arena()).value();
    EXPECT_TRUE(reduces_to_zero(pade.value().numerator, ctx.simplify(expected_num).value(), ctx));
    EXPECT_TRUE(reduces_to_zero(pade.value().denominator, ctx.simplify(expected_den).value(), ctx));
}

TEST_F(PadeApproximantTest, CosineAtAlgebraicCentreAntiHardcode) {
    // cos(x) at x = π/4.  Taylor coefficients carry √2/2; the legacy rational
    // path bailed out with Unimplemented.  Verify the Taylor-truncation
    // identity P − f·Q ≡ 0 mod (x − π/4)^{m+n+1} for [1/1].
    auto expr = parse_expr("cos(x)", ctx.arena()).value();
    auto center = parse_expr("pi/4", ctx.arena()).value();
    auto pade = calculus::pade_approximant(expr, Symbol("x"), center, 1U, 1U, ctx);
    ASSERT_TRUE(pade.is_ok()) << pade.error().message;

    auto pf = ctx.arena().make<Binary>(BinaryOp::Sub,
        pade.value().numerator,
        ctx.arena().make<Binary>(BinaryOp::Mul, expr, pade.value().denominator));
    auto pf_simp = ctx.simplify(pf);
    ASSERT_TRUE(pf_simp.is_ok());
    for (unsigned int k = 0; k <= 2U; ++k) {
        Result<ExprPtr> deriv = (k == 0U)
            ? ok(pf_simp.value())
            : calculus::diff(pf_simp.value(), Symbol("x"), k, ctx);
        ASSERT_TRUE(deriv.is_ok());
        auto val = ctx.substitute(deriv.value(), Symbol("x"), center);
        ASSERT_TRUE(val.is_ok());
        auto val_s = ctx.simplify(val.value());
        ASSERT_TRUE(val_s.is_ok());
        const auto* lit = expr_cast<IntegerLit>(val_s.value());
        EXPECT_TRUE(lit != nullptr && lit->value.is_zero())
            << "Pade [1/1] of cos at pi/4 defect at order " << k;
    }
}

TEST_F(PadeApproximantTest, ExpAtSymbolicCentreLnTwoAntiHardcode) {
    // exp(x) at x = ln(2).  c_k = 2 / k!  ∈ Q, but the substitution goes through
    // exp(ln(2)) which only folds to 2 via the simplifier — exercises the
    // symbolic path even though the final coefficients land in Q.
    auto expr = parse_expr("exp(x)", ctx.arena()).value();
    auto center = parse_expr("ln(2)", ctx.arena()).value();
    auto pade = calculus::pade_approximant(expr, Symbol("x"), center, 1U, 1U, ctx);
    ASSERT_TRUE(pade.is_ok()) << pade.error().message;
    // P(centre)/Q(centre) must equal exp(ln 2) = 2.
    auto num_at = ctx.substitute(pade.value().numerator, Symbol("x"), center);
    auto den_at = ctx.substitute(pade.value().denominator, Symbol("x"), center);
    ASSERT_TRUE(num_at.is_ok());
    ASSERT_TRUE(den_at.is_ok());
    auto ratio = ctx.arena().make<Binary>(BinaryOp::Div, num_at.value(), den_at.value());
    auto two = parse_expr("2", ctx.arena()).value();
    EXPECT_TRUE(reduces_to_zero(ratio, ctx.simplify(two).value(), ctx));
}

TEST_F(PadeApproximantTest, ShiftedCenterAntiHardcode) {
    // Centre at x = 1.  f(x) = 1/x = 1/(1 + (x−1)) → [1/1] Pade is exact.
    auto expr = parse_expr("1/x", ctx.arena()).value();
    auto center = parse_expr("1", ctx.arena()).value();
    auto pade = calculus::pade_approximant(expr, Symbol("x"), center, 0U, 1U, ctx);
    ASSERT_TRUE(pade.is_ok()) << pade.error().message;
    // Expected: P = 1, Q = 1 + (x − 1) = x  (Q is normalised at the centre, so
    // Q(centre) = 1 and Q(x) = 1 + (x − 1)).
    auto expected_num = parse_expr("1", ctx.arena()).value();
    EXPECT_TRUE(reduces_to_zero(pade.value().numerator, ctx.simplify(expected_num).value(), ctx));

    // P/Q at x = 2 must equal 1/2.
    auto eval_num = ctx.substitute(pade.value().numerator, Symbol("x"), parse_expr("2", ctx.arena()).value());
    auto eval_den = ctx.substitute(pade.value().denominator, Symbol("x"), parse_expr("2", ctx.arena()).value());
    ASSERT_TRUE(eval_num.is_ok());
    ASSERT_TRUE(eval_den.is_ok());
    auto ratio = ctx.arena().make<Binary>(BinaryOp::Div, eval_num.value(), eval_den.value());
    auto half = parse_expr("1/2", ctx.arena()).value();
    EXPECT_TRUE(reduces_to_zero(ratio, ctx.simplify(half).value(), ctx));
}

}  // namespace cas::test
