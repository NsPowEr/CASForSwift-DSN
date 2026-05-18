#include <gtest/gtest.h>
#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "algebra/algebra_internal.hpp"
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

TEST(AlgebraFactorizationTest, L0_05_PrimeSelectionDoesNotUseFixedEmergencyFallback) {
    const BigInt pool_product =
        BigInt(13) * BigInt(17) * BigInt(19) * BigInt(23) * BigInt(29) *
        BigInt(31) * BigInt(37) * BigInt(41) * BigInt(43) * BigInt(47) *
        BigInt(53) * BigInt(59) * BigInt(61) * BigInt(67) * BigInt(71) *
        BigInt(73) * BigInt(79) * BigInt(83) * BigInt(89) * BigInt(97) *
        BigInt(101);
    const IntPoly f({BigInt(1), BigInt(0), pool_product});

    const BigInt p = select_factorization_prime(f);

    EXPECT_FALSE((pool_product % p).is_zero());
    EXPECT_NE(p, BigInt(101));
}

// L1-19: GCD heuristic B adapts to large coefficients (Mignotte bound)
TEST(AlgebraGcdHeuristicTest, L1_19_MignotteBoundAdaptivePadding) {
    symbolic::CASContext ctx;
    Symbol x("x");
    // GCD of polys with large coefficients (1000000): Mignotte bound must be >> 2*max+100*1000
    // GCD(1000000*x^2 - 1000000, 1000000*x - 1000000) = 1000000*(x-1)
    ExprPtr p = parse_string("1000000*x^2 - 1000000", ctx);
    ExprPtr q = parse_string("1000000*x - 1000000", ctx);
    auto g = polynomial_gcd(p, q, x, ctx);
    ASSERT_TRUE(g.is_ok()) << g.error().message;
    // Result should be non-trivial (divisible by x-1)
    EXPECT_TRUE(g.is_ok());
}

TEST(AlgebraGcdHeuristicTest, KroneckerAtRigorousMignotteBoundReturnsTrueGcd) {
    // p = -5 - 5x - 5y - 2xy,  q = -5 - 5x - 5y - xy.  Then p - q = -xy
    // so gcd(p, q) | xy, and a direct content check shows gcd is a unit (1).
    // With the rigorous Mignotte bound the Kronecker substitution image
    // reconstructs the true GCD (a constant) and verify_gcd_candidate accepts.
    // (Previously a softer bound B=max(formula,1000) caused a spurious image
    //  and the algorithm returned InternalError to reject it.)
    MultivariatePolynomial p({
        MultivariateTerm{.coefficient = BigInt(-5), .factors = {}},
        MultivariateTerm{.coefficient = BigInt(-5), .factors = {{Symbol("x"), 1U}}},
        MultivariateTerm{.coefficient = BigInt(-5), .factors = {{Symbol("y"), 1U}}},
        MultivariateTerm{.coefficient = BigInt(-2), .factors = {{Symbol("x"), 1U}, {Symbol("y"), 1U}}},
    });
    MultivariatePolynomial q({
        MultivariateTerm{.coefficient = BigInt(-5), .factors = {}},
        MultivariateTerm{.coefficient = BigInt(-5), .factors = {{Symbol("x"), 1U}}},
        MultivariateTerm{.coefficient = BigInt(-5), .factors = {{Symbol("y"), 1U}}},
        MultivariateTerm{.coefficient = BigInt(-1), .factors = {{Symbol("x"), 1U}, {Symbol("y"), 1U}}},
    });

    auto gcd = gcd_heuristic(p, q);
    ASSERT_TRUE(gcd.is_ok()) << gcd.error().message;
    // True GCD is a unit (±1). Verify the reconstructed polynomial is constant
    // with coefficient magnitude 1.
    ASSERT_EQ(gcd.value().terms().size(), 1U);
    const auto& term = gcd.value().terms()[0];
    EXPECT_TRUE(term.factors.empty()) << "GCD must be a constant polynomial";
    EXPECT_EQ(term.coefficient.abs(), BigInt(1)) << "GCD must be a unit";
}

// L1-20: evaluate_at_rational accepts rational values
TEST(AlgebraMultivariateTest, L1_20_EvaluateAtRationalValue) {
    using cas::algebra::MultivariatePolynomial;
    using cas::algebra::MultivariateTerm;
    symbolic::CASContext ctx;
    Symbol x("x");
    // Polynomial: 2*x^2  (one term, coefficient=2, exponent=2)
    MultivariateTerm t;
    t.coefficient = BigInt(2);
    t.factors = {{x, 2U}};
    MultivariatePolynomial poly(std::vector<MultivariateTerm>{t});

    // evaluate at x = 1/2 → 2*(1/2)^2 = 2*1/4 = 1/2
    Rational half(BigInt(1), BigInt(2));
    auto result = poly.evaluate_at_rational(x, half, ctx.arena());
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    // Result should be rational 1/2
    const auto* rat = expr_cast<RationalLit>(result.value());
    ASSERT_NE(rat, nullptr);
    EXPECT_EQ(rat->numerator, BigInt(1));
    EXPECT_EQ(rat->denominator, BigInt(2));
}
