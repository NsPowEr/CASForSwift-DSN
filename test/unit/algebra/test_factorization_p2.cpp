#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "algebra/polynomial_internal.hpp"
#include <gtest/gtest.h>
#include <vector>

namespace cas::algebra {

TEST(FactorPolynomialP2, EDFTraceDegree3) {
    // x^6 + x^5 + x^4 + x^3 + x^2 + x + 1 = (x^3 + x + 1) * (x^3 + x^2 + 1) mod 2
    IntPoly f({BigInt(1), BigInt(1), BigInt(1), BigInt(1), BigInt(1), BigInt(1), BigInt(1)});
    BigInt p(2);
    
    auto factors = factor_polynomial_mod_p(f, p);
    ASSERT_TRUE(factors.is_ok());
    // Should find 2 factors of degree 3
    EXPECT_EQ(factors.value().size(), 2U);
    for (const auto& fact : factors.value()) {
        EXPECT_EQ(fact.degree(), 3U);
    }
}

TEST(FactorPolynomialP2, SimpleDegree1) {
    // x^2 + x = x(x+1) mod 2
    IntPoly f({BigInt(0), BigInt(1), BigInt(1)});
    BigInt p(2);
    
    auto factors = factor_polynomial_mod_p(f, p);
    ASSERT_TRUE(factors.is_ok());
    EXPECT_EQ(factors.value().size(), 2U);
    for (const auto& fact : factors.value()) {
        EXPECT_EQ(fact.degree(), 1U);
    }
}

TEST(FactorPolynomialP2, Degree4) {
    // Irreducible factors of degree 4 in F2:
    // x^4 + x + 1
    // x^4 + x^3 + 1
    // x^4 + x^3 + x^2 + x + 1
    
    // f = (x^4 + x + 1) * (x^4 + x^3 + 1)
    // f = x^8 + x^7 + x^4 + x^5 + x^4 + x + x^4 + x^3 + 1
    // f = x^8 + x^7 + x^5 + x^4 + x^3 + x + 1
    
    IntPoly f({BigInt(1), BigInt(1), BigInt(0), BigInt(1), BigInt(1), BigInt(1), BigInt(0), BigInt(1), BigInt(1)});
    BigInt p(2);
    
    auto factors = factor_polynomial_mod_p(f, p);
    ASSERT_TRUE(factors.is_ok());
    EXPECT_EQ(factors.value().size(), 2U);
    for (const auto& fact : factors.value()) {
        EXPECT_EQ(fact.degree(), 4U);
    }
}

} // namespace cas::algebra
