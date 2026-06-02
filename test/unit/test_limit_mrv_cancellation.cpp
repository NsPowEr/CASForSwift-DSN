#include <gtest/gtest.h>
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "cas/formatter.hpp"

using namespace cas;
using namespace cas::calculus;

class LimitMrvTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};

    ExprPtr parse(const std::string& input) {
        auto tokens = Lexer(input).tokenize();
        if (tokens.is_error()) {
            throw std::runtime_error("Lex error: " + tokens.error().message);
        }

        Parser parser(tokens.value(), ctx.arena());
        auto res = parser.parse();
        if (res.is_error()) {
            throw std::runtime_error("Parse error: " + res.error().message);
        }
        return res.value();
    }

    void verify_limit(const std::string& expr_str, const std::string& expected_str) {
        auto expr = parse(expr_str);
        auto expected = parse(expected_str);
        auto inf = ctx.arena().make<Constant>(MathConstant::Infinity);
        
        auto res = limit(expr, x, inf, LimitDirection::Both, ctx);
        ASSERT_TRUE(res.is_ok()) << "Limit failed for " << expr_str << ": " << res.error().message;
        
        auto simplified_res = ctx.simplify(res.value());
        ASSERT_TRUE(simplified_res.is_ok());
        
        auto simplified_expected = ctx.simplify(expected);
        ASSERT_TRUE(simplified_expected.is_ok());

        EXPECT_TRUE(structural_equal(simplified_res.value(), simplified_expected.value()))
            << "For " << expr_str << " limit(x->inf) expected " << expected_str 
            << " but got " << formatter::TextFormatter{}.format(simplified_res.value());
    }
};

TEST_F(LimitMrvTest, CancellationSimple) {
    // (x + 1) - x -> 1
    verify_limit("(x + 1) - x", "1");
}

TEST_F(LimitMrvTest, CancellationExponential) {
    // (exp(x) + 1) - exp(x) -> 1
    verify_limit("(exp(x) + 1) - exp(x)", "1");
}

TEST_F(LimitMrvTest, CancellationTower) {
    // exp(x + exp(-x)) - exp(x) -> 1
    verify_limit("exp(x + exp(-x)) - exp(x)", "1");
}

TEST_F(LimitMrvTest, ComplexCancellation) {
    // (exp(x) + x) / (exp(x) + 1) -> 1
    verify_limit("(exp(x) + x) / (exp(x) + 1)", "1");
}

// ── F5.2 / B3 — Dynamic Gruntz growth-rank tests ────────────────────────────
//
// compare_growth (limit_mrv.cpp) is now the single source of truth used by
// both try_infinite_limit (Sum / Binary dominance) and the MRV engine.
// These tests exercise the dynamic recursive comparison at arbitrary tower
// depth, which the legacy static rank could not represent.

TEST_F(LimitMrvTest, GruntzDoubleExponentialDominatesPolynomialExp) {
    // exp(exp(x)) - exp(x^5)  → +∞  (exp(exp(x)) lives strictly above exp(x^N)).
    verify_limit("exp(exp(x)) - exp(x^5)", "inf");
}

TEST_F(LimitMrvTest, GruntzDoubleExponentialPolynomialExpReverse) {
    // exp(x^5) - exp(exp(x))  → -∞.
    verify_limit("exp(x^5) - exp(exp(x))", "-inf");
}

TEST_F(LimitMrvTest, GruntzTripleExponentialDominatesDouble) {
    // exp(exp(exp(x))) - exp(exp(x^3))  → +∞.
    verify_limit("exp(exp(exp(x))) - exp(exp(x^3))", "inf");
}

TEST_F(LimitMrvTest, GruntzLogPolyVsExp) {
    // ln(x) + x^10 vs exp(x):  exp(x) dominates → result is -∞ (since exp(x)
    // appears with negative sign).
    verify_limit("ln(x) + x^10 - exp(x)", "-inf");
}

TEST_F(LimitMrvTest, GruntzNestedLogVsPolynomial) {
    // exp(x) - ln(ln(x))·x^100  → +∞.  exp(x) sits in rank 3, ln(ln(x))·x^100
    // sits in rank 2 (polynomial × bounded slow log).
    verify_limit("exp(x) - ln(ln(x))*x^100", "inf");
}
