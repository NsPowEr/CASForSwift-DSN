// F3.3 — Cyclotomic polynomial via Möbius inversion: power-gain tests.
//
// Pre-fix: `compute_cyclotomic` used recursive identity with thread-
// unsafe static cache and recursion depth ≤ log₂(n).
// Post-fix: direct Möbius construction Φ_n(x) = ∏ (x^d − 1)^μ(n/d),
// no recursion, no cache. Verifies known closed forms and stress on
// large n where pre-fix would have hit the kCyclotomicMaxN cap.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "../../../src/algebra/polynomial_internal.hpp"

using namespace cas;
using namespace cas::algebra;

namespace {

[[nodiscard]] IntPoly poly(std::initializer_list<long long> coeffs) {
    std::vector<BigInt> bs;
    bs.reserve(coeffs.size());
    for (auto c : coeffs) bs.emplace_back(c);
    return IntPoly(std::move(bs));
}

[[nodiscard]] bool poly_eq(const IntPoly& a, const IntPoly& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (a[i] != b[i]) return false;
    return true;
}

}  // namespace

// Φ_1(x) = x - 1
TEST(CyclotomicMobiusTest, Phi1) {
    auto v = is_cyclotomic(poly({-1, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 1);
}

// Φ_2(x) = x + 1
TEST(CyclotomicMobiusTest, Phi2) {
    auto v = is_cyclotomic(poly({1, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 2);
}

// Φ_3(x) = x² + x + 1
TEST(CyclotomicMobiusTest, Phi3) {
    auto v = is_cyclotomic(poly({1, 1, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 3);
}

// Φ_4(x) = x² + 1
TEST(CyclotomicMobiusTest, Phi4) {
    auto v = is_cyclotomic(poly({1, 0, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 4);
}

// Φ_6(x) = x² - x + 1 (canonical "tricky" case: μ(6/d) involves
// μ(6) = +1, μ(2)/μ(3) = -1, μ(1) = +1).
TEST(CyclotomicMobiusTest, Phi6) {
    auto v = is_cyclotomic(poly({1, -1, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 6);
}

// Φ_8(x) = x^4 + 1
TEST(CyclotomicMobiusTest, Phi8) {
    auto v = is_cyclotomic(poly({1, 0, 0, 0, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 8);
}

// Φ_12(x) = x^4 - x² + 1
TEST(CyclotomicMobiusTest, Phi12) {
    auto v = is_cyclotomic(poly({1, 0, -1, 0, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 12);
}

// Φ_15(x) = x^8 − x^7 + x^5 − x^4 + x³ − x + 1
TEST(CyclotomicMobiusTest, Phi15) {
    auto v = is_cyclotomic(poly({1, -1, 0, 1, -1, 1, 0, -1, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 15);
}

// Power-gain: Φ_97 (prime, degree 96). Pre-fix this would have spent
// ~97 recursive calls each allocating polynomial divisors; post-fix
// is a single Möbius pass. The polynomial is x^96 + x^95 + ... + 1.
TEST(CyclotomicMobiusTest, Phi97PrimeStillIdentified) {
    // Φ_p(x) = (x^p - 1)/(x - 1) = sum_{i=0}^{p-1} x^i  (for p prime)
    std::vector<BigInt> c97(97, BigInt(1));
    IntPoly phi_97(std::move(c97));
    auto v = is_cyclotomic(phi_97);
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 97);
}

// Power-gain: non-cyclotomic input must return nullopt without
// running forever or hitting the recursion cap.
TEST(CyclotomicMobiusTest, NonCyclotomicReturnsNullopt) {
    // (x+2) is not a cyclotomic polynomial.
    auto v = is_cyclotomic(poly({2, 1}));
    EXPECT_FALSE(v.has_value());
}
