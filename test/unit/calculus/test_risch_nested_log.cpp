// A27 — log-derivative recognition for NESTED log towers.
// ∫ 1/(x·ln(x)·ln(ln(x))) dx = ln(ln(ln(x))): the integrand is D(g)/g with
// g = ln(ln(x)) a top generator whose argument is itself a log generator.
// Every accepted antiderivative is verified by differentiation round-trip
// (sound-by-verify — the assertions never trust the solver's own claim).

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <memory>
#include <string>

namespace cas {
namespace {

class RischNestedLogTest : public ::testing::Test {
protected:
    void SetUp() override {
        ctx = std::make_unique<symbolic::CASContext>();
    }

    [[nodiscard]] ExprPtr parse_ok(const std::string& input) {
        Lexer lexer(input);
        auto tokens = lexer.tokenize();
        EXPECT_TRUE(tokens.is_ok()) << tokens.error().message;
        if (tokens.is_error()) return nullptr;
        Parser parser(tokens.value(), ctx->arena());
        auto parsed = parser.parse();
        EXPECT_TRUE(parsed.is_ok()) << parsed.error().message;
        return parsed.is_ok() ? parsed.value() : nullptr;
    }

    // D(candidate) − integrand ≡ 0, checked structurally after together+simplify.
    void expect_antiderivative(ExprPtr candidate, ExprPtr integrand, const Symbol& x) {
        auto d = calculus::diff(candidate, x, 1U, *ctx);
        ASSERT_TRUE(d.is_ok()) << d.error().message;
        ExprPtr delta = ctx->arena().make<Binary>(BinaryOp::Sub, d.value(), integrand);
        auto tog = algebra::together(delta, *ctx);
        ExprPtr t = tog.is_ok() ? tog.value() : delta;
        auto simp = ctx->simplify(t);
        ASSERT_TRUE(simp.is_ok()) << simp.error().message;
        bool zero = false;
        if (const auto* il = expr_cast<IntegerLit>(simp.value())) zero = il->value.is_zero();
        if (const auto* rl = expr_cast<RationalLit>(simp.value())) zero = rl->numerator.is_zero();
        EXPECT_TRUE(zero) << "round-trip residual is not zero\n"
                          << "candidate: " << debug_print(candidate) << "\n"
                          << "residual:  " << debug_print(simp.value());
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

TEST_F(RischNestedLogTest, ThreeLevelChain_UnitConstant) {
    // ∫ 1/(x·ln(x)·ln(ln(x))) dx = ln(ln(ln(x)))
    Symbol x{"x"};
    ExprPtr integrand = parse_ok("1/(x*ln(x)*ln(ln(x)))");
    ASSERT_NE(integrand, nullptr);
    auto res = calculus::integrate(integrand, x, *ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_antiderivative(res.value(), integrand, x);
}

TEST_F(RischNestedLogTest, ThreeLevelChain_IntegerConstant) {
    // ∫ 3/(x·ln(x)·ln(ln(x))) dx = 3·ln(ln(ln(x))) — formal constant c=3 must
    // survive the nested-generator probe substitution (no closed constant set).
    Symbol x{"x"};
    ExprPtr integrand = parse_ok("3/(x*ln(x)*ln(ln(x)))");
    ASSERT_NE(integrand, nullptr);
    auto res = calculus::integrate(integrand, x, *ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_antiderivative(res.value(), integrand, x);
}

TEST_F(RischNestedLogTest, TwoLevelChain_RegressionGuard) {
    // ∫ 1/(x·ln(x)) dx = ln(ln(x)) — the already-working 2-level case must
    // keep working (HPP-007 closure regression guard).
    Symbol x{"x"};
    ExprPtr integrand = parse_ok("1/(x*ln(x))");
    ASSERT_NE(integrand, nullptr);
    auto res = calculus::integrate(integrand, x, *ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_antiderivative(res.value(), integrand, x);
}

}  // namespace
}  // namespace cas
