#include <gtest/gtest.h>
#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "algebra/polynomial_internal.hpp"

using namespace cas;
using namespace cas::algebra;

static ExprPtr parse_string(const std::string& input, symbolic::CASContext& ctx) {
    Lexer lexer(input);
    auto tokens = lexer.tokenize();
    if (tokens.is_error()) throw std::runtime_error("Lex error: " + tokens.error().message);
    Parser parser(tokens.value(), ctx.arena());
    auto result = parser.parse();
    if (result.is_error()) {
        throw std::runtime_error("Parse error: " + result.error().message);
    }
    return result.value();
}

TEST(AlgebraLLLTest, BasicReduction4x4) {
    auto R = [](long long n) { return Rational(BigInt(n)); };
    LatticeMatrix basis = {
        {R(1), R(1),  R(1), R(1)},
        {R(-1), R(0), R(2), R(1)},
        {R(0),  R(1), R(2), R(3)},
        {R(1),  R(2), R(3), R(4)}
    };

    lll_reduction(basis);

    // Check basis[0] has small norm squared
    double norm_sq = 0.0;
    for (const auto& x : basis[0]) norm_sq += x.to_double() * x.to_double();
    EXPECT_LT(norm_sq, 10.0);
}

TEST(AlgebraHenselTest, UnivariateLifting) {
    // f = x^2 - 1, g = x-1, h = x+1 mod 5
    // Lift to 5^2 = 25
    IntPoly f({BigInt(-1), BigInt(0), BigInt(1)});
    IntPoly g({BigInt(-1), BigInt(1)});
    IntPoly h({BigInt(1), BigInt(1)});
    BigInt p(5);
    
    auto result = hensel_lift(f, g, h, p, 2);
    ASSERT_TRUE(result.is_ok());
    
    auto [G, H] = result.value();
    // G*H should be f mod 25
    // Actually for x^2-1 it should stay x-1, x+1
    EXPECT_EQ(G.degree(), 1U);
    EXPECT_EQ(H.degree(), 1U);
}

TEST(AlgebraFactorizationTest, Degree10LargeCoeffs) {
    symbolic::CASContext ctx;
    Symbol x("x");
    
    // (x^5 + 3x + 1)(x^5 - x^2 + 7)
    // = x^10 - x^7 + 3x^6 + 8x^5 - 3x^3 - x^2 + 21x + 7
    std::string poly_str = "x^10 - x^7 + 3*x^6 + 8*x^5 - 3*x^3 - x^2 + 21*x + 7";
    ExprPtr poly = parse_string(poly_str, ctx);

    auto result = factor_over_integers(poly, x, ctx);
    ASSERT_TRUE(result.is_ok());
    
    // Should have 2 factors
    // Note: one might be the negative of the other if content is -1, but here it should be 2 factors of degree 5.
    EXPECT_EQ(result.value().factors.size(), 2U);
    
    // Verify each factor degree is 5
    for (const auto& f : result.value().factors) {
        auto parsed = parse_polynomial(f.factor, x, ctx);
        ASSERT_TRUE(parsed.is_ok());
        EXPECT_EQ(parsed.value().degree(), 5U);
    }
}

// CAS-L0-05: prime selection must be hash-based, not fixed p=13
// Polynomial with lc=13 forces all prior fallback candidates to be skipped
TEST(AlgebraFactorizationTest, L0_05_HashBasedPrimeSelection_Lc13) {
    symbolic::CASContext ctx;
    Symbol x("x");
    // lc=13: p=13 divides lc. The hash-based selector must pick a different prime.
    // 13*(x^2 - 2) = 13*x^2 - 26; irreducible over Z
    ExprPtr poly = parse_string("13*x^2 - 26", ctx);
    auto result = factor_over_integers(poly, x, ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    // Content factor 13, plus (x^2 - 2) irreducible
    EXPECT_GE(result.value().factors.size(), 1U);
}

TEST(AlgebraFactorizationTest, L0_05_HashBasedPrimeSelection_DifferentPolys) {
    symbolic::CASContext ctx;
    Symbol x("x");
    // Two polynomials that are factorizations of (x-1)(x+1)=x^2-1 but with different
    // leading coefficients: lc=1 vs lc=2*(x^2-1) would pick different starting primes
    ExprPtr p1 = parse_string("x^2 - 1", ctx);
    ExprPtr p2 = parse_string("x^4 - 1", ctx);
    auto r1 = factor_over_integers(p1, x, ctx);
    auto r2 = factor_over_integers(p2, x, ctx);
    ASSERT_TRUE(r1.is_ok()) << r1.error().message;
    ASSERT_TRUE(r2.is_ok()) << r2.error().message;
    // x^2-1 = (x-1)(x+1) -> 2 factors
    EXPECT_EQ(r1.value().factors.size(), 2U);
    // x^4-1 = (x-1)(x+1)(x^2+1) -> 3 factors
    EXPECT_EQ(r2.value().factors.size(), 3U);
}
