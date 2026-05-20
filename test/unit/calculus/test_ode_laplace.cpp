// CAS-L3-10 — ODE solver via Laplace transform tests.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include "../../../src/calculus/calculus_internal.hpp"

using namespace cas;
using namespace cas::calculus;

namespace {

class OdeLaplaceTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol t{"t"};
    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto tk = Lexer(s).tokenize();
        EXPECT_TRUE(tk.is_ok()) << s;
        Parser p(tk.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
    [[nodiscard]] ExprPtr int_e(long long v) {
        return ctx.arena().make<IntegerLit>(BigInt(v));
    }
    [[nodiscard]] bool equiv(ExprPtr a, ExprPtr b) {
        auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, a, b);
        auto tog = algebra::together(delta, ctx);
        auto simp = ctx.simplify(tog.is_ok() ? tog.value() : delta);
        if (simp.is_error()) return false;
        if (auto* il = expr_cast<IntegerLit>(simp.value())) return il->value.is_zero();
        return false;
    }
};

TEST_F(OdeLaplaceTest, FirstOrderHomogeneous) {
    // y' - y = 0, y(0)=1 → y = exp(t)
    // coeffs[0]=-1, coeffs[1]=1
    std::vector<ExprPtr> coeffs = {int_e(-1), int_e(1)};
    ExprPtr f = int_e(0);
    std::vector<ExprPtr> ics = {int_e(1)};
    auto r = solve_ode_laplace(coeffs, f, ics, t, ctx);
    if (!r.is_ok()) {
        FAIL() << "solve_ode_laplace error: " << r.error().message;
    }
    auto expected = parse("exp(t)");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(OdeLaplaceTest, FirstOrderForced) {
    // y' + y = 1, y(0)=0 → y = 1 - exp(-t).
    // Y(s) = 1/(s(s+1)) gestito via PFD + Unary(Neg) bridge.
    std::vector<ExprPtr> coeffs = {int_e(1), int_e(1)};
    ExprPtr f = int_e(1);
    std::vector<ExprPtr> ics = {int_e(0)};
    auto r = solve_ode_laplace(coeffs, f, ics, t, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    // Verify y' + y = 1.
    auto d = calculus::diff(r.value(), t, 1U, ctx);
    ASSERT_TRUE(d.is_ok());
    ExprPtr lhs = ctx.arena().make<Binary>(BinaryOp::Add, d.value(), r.value());
    auto tog = algebra::together(lhs, ctx);
    auto simp = ctx.simplify(tog.is_ok() ? tog.value() : lhs);
    ASSERT_TRUE(simp.is_ok());
    auto* lit = expr_cast<IntegerLit>(simp.value());
    EXPECT_TRUE(lit != nullptr && lit->value == BigInt(1))
        << "y' + y should equal 1";
}

TEST_F(OdeLaplaceTest, SecondOrderHarmonicOscillator) {
    // y'' + y = 0, y(0)=1, y'(0)=0 → y = cos(t)
    std::vector<ExprPtr> coeffs = {int_e(1), int_e(0), int_e(1)};
    ExprPtr f = int_e(0);
    std::vector<ExprPtr> ics = {int_e(1), int_e(0)};
    auto r = solve_ode_laplace(coeffs, f, ics, t, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto expected = parse("cos(t)");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(OdeLaplaceTest, SecondOrderHarmonicSine) {
    // y'' + y = 0, y(0)=0, y'(0)=1 → y = sin(t)
    std::vector<ExprPtr> coeffs = {int_e(1), int_e(0), int_e(1)};
    ExprPtr f = int_e(0);
    std::vector<ExprPtr> ics = {int_e(0), int_e(1)};
    auto r = solve_ode_laplace(coeffs, f, ics, t, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto expected = parse("sin(t)");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(OdeLaplaceTest, AntiHardcodeFirstOrderDecay) {
    // y' + 3y = 0, y(0)=2 → y = 2·exp(-3t)
    std::vector<ExprPtr> coeffs = {int_e(3), int_e(1)};
    ExprPtr f = int_e(0);
    std::vector<ExprPtr> ics = {int_e(2)};
    auto r = solve_ode_laplace(coeffs, f, ics, t, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto expected = parse("2 * exp(-3 * t)");
    EXPECT_TRUE(equiv(r.value(), expected));
}

}  // namespace
