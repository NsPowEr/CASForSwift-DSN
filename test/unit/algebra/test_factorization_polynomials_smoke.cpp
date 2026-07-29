// DEBT-002 smoke coverage for src/algebra/factorization_polynomials.cpp.
// Targets factor_polynomial / factor_over_integers covering: trivial,
// linear product, quadratic discriminant ≥ 0, irreducible Q[x],
// squarefree-distinct, and rational coefficient handling.

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "algebra/polynomial_internal.hpp"
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

TEST_F(FactorizationPolynomialsSmokeTest,
       RecombinationComplementKeepsLargeModularFactor) {
    // Regression (found via A6, Φ₇ 3-set resolvent): mod 109 this factors
    // as degrees {1, 5, 2}. The historic recombination DROPPED modular
    // factors of degree > n/2 from the pool, so the complement `right` no
    // longer satisfied left·right ≡ f (mod p): hensel_lift then produced
    // garbage for every subset and the polynomial was silently declared
    // irreducible. Expected: {2, 6} — Maxima-verified
    // (y²+y+2)·(y⁶+3y⁵+2y⁴−y³+4y²−2y+1).
    auto e = parse(
        "x^8 + 4*x^7 + 7*x^6 + 7*x^5 + 7*x^4 + 7*x^2 - 3*x + 2");
    auto r = algebra::factor_over_integers(e, x, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    std::vector<std::size_t> degs;
    for (const auto& pf : r.value().factors) {
        auto pp = algebra::parse_polynomial(pf.factor, x, ctx);
        ASSERT_TRUE(pp.is_ok());
        for (unsigned int m = 0U; m < pf.multiplicity; ++m) {
            degs.push_back(algebra::poly_degree(pp.value()));
        }
    }
    std::sort(degs.begin(), degs.end());
    EXPECT_EQ(degs, (std::vector<std::size_t>{2U, 6U}));
}

}  // namespace
