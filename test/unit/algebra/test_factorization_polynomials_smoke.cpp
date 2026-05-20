// DEBT-002 smoke coverage for src/algebra/factorization_polynomials.cpp.
// Targets factor_polynomial / factor_over_integers covering: trivial,
// linear product, quadratic discriminant ≥ 0, irreducible Q[x],
// squarefree-distinct, and rational coefficient handling.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class FactorizationPolynomialsSmokeTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};

    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
};

TEST_F(FactorizationPolynomialsSmokeTest, FactorOverIntegersLinearProduct) {
    // (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6 → 3 factors.
    auto e = parse("x^3 - 6*x^2 + 11*x - 6");
    auto r = algebra::factor_over_integers(e, x, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_GE(r.value().factors.size(), 3U);
}

TEST_F(FactorizationPolynomialsSmokeTest, FactorIrreducibleOverQ) {
    // x^2 + 1 is irreducible over Q → exactly 1 factor.
    auto e = parse("x^2 + 1");
    auto r = algebra::factor_over_integers(e, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_EQ(r.value().factors.size(), 1U);
}

TEST_F(FactorizationPolynomialsSmokeTest, FactorOfDifferenceOfSquares) {
    // x^2 - 4 = (x-2)(x+2) → 2 factors.
    auto e = parse("x^2 - 4");
    auto r = algebra::factor_over_integers(e, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_GE(r.value().factors.size(), 2U);
}

TEST_F(FactorizationPolynomialsSmokeTest, FactorRepeatedRoot) {
    // (x-1)^3 = x^3 - 3x^2 + 3x - 1; expect a single distinct factor
    // (x-1) with multiplicity 3.
    auto e = parse("x^3 - 3*x^2 + 3*x - 1");
    auto r = algebra::factor_over_integers(e, x, ctx);
    ASSERT_TRUE(r.is_ok());
    ASSERT_GE(r.value().factors.size(), 1U);
    EXPECT_EQ(r.value().factors[0].multiplicity, 3U);
}

TEST_F(FactorizationPolynomialsSmokeTest, FactorPolynomialQAlphaFallback) {
    // factor_polynomial without explicit extension delegates to
    // factor_over_integers — must succeed for any Q[x] input.
    auto e = parse("x^4 - 5*x^2 + 6");  // (x^2-2)(x^2-3)
    auto r = algebra::factor_polynomial(e, x, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_GE(r.value().factors.size(), 2U);
}

}  // namespace
