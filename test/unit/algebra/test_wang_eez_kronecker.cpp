#include <gtest/gtest.h>
#include "cas/bigint.hpp"
#include "../../../src/algebra/polynomial_internal.hpp"
#include "../../../src/algebra/algebra_internal.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::algebra;

namespace {

IntPoly make_poly(std::vector<long long> coeffs) {
    std::vector<BigInt> c;
    c.reserve(coeffs.size());
    for (long long v : coeffs) c.emplace_back(v);
    IntPoly p(std::move(c));
    normalize_integer_poly(p);
    return p;
}

bool polys_equal(const IntPoly& a, const IntPoly& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

bool check_factorization(const IntPoly& original, const std::vector<IntPoly>& factors) {
    IntPoly prod(std::vector<BigInt>{BigInt(1)});
    for (const auto& f : factors) {
        // Each factor must have degree >= 1
        if (f.degree() < 1) return false;
        // Multiply prod * f
        IntPoly next_prod;
        next_prod.resize(prod.size() + f.size() - 1, BigInt(0));
        for (std::size_t i = 0; i < prod.size(); ++i) {
            for (std::size_t j = 0; j < f.size(); ++j) {
                next_prod[i + j] += prod[i] * f[j];
            }
        }
        normalize_integer_poly(next_prod);
        prod = std::move(next_prod);
    }
    return polys_equal(original, prod);
}

} // namespace

TEST(WangEezKroneckerTest, ContextParamsExposed) {
    symbolic::CASContext ctx;
    
    // Check default values
    EXPECT_EQ(ctx.max_hensel_lift_attempts(), 8U);
    EXPECT_EQ(ctx.kronecker_max_degree(), 8U);
    
    // Check setting and getting
    ctx.set_max_hensel_lift_attempts(15U);
    EXPECT_EQ(ctx.max_hensel_lift_attempts(), 15U);
    
    ctx.set_kronecker_max_degree(12U);
    EXPECT_EQ(ctx.kronecker_max_degree(), 12U);
}

TEST(WangEezKroneckerTest, KroneckerFactorizationQuadratic) {
    symbolic::CASContext ctx;
    
    // (x - 2)(x - 3) = x^2 - 5x + 6
    IntPoly f = make_poly({6, -5, 1});
    auto res = factorize_kronecker(f, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_EQ(res.value().size(), 2U);
    EXPECT_TRUE(check_factorization(f, res.value()));
    
    // x^2 + x + 1 (irreducible)
    IntPoly f_irr = make_poly({1, 1, 1});
    auto res_irr = factorize_kronecker(f_irr, ctx);
    ASSERT_TRUE(res_irr.is_ok());
    EXPECT_EQ(res_irr.value().size(), 1U);
    EXPECT_TRUE(polys_equal(res_irr.value()[0], f_irr));
}

TEST(WangEezKroneckerTest, KroneckerFactorizationQuarticProductOfQuadratics) {
    symbolic::CASContext ctx;
    
    // (x^2 + x + 1)(x^2 + 2x + 2) = x^4 + 3x^3 + 5x^2 + 4x + 2
    IntPoly f = make_poly({2, 4, 5, 3, 1});
    auto res = factorize_kronecker(f, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_EQ(res.value().size(), 2U);
    EXPECT_TRUE(check_factorization(f, res.value()));
}

TEST(WangEezKroneckerTest, HenselDispatcherNormalFactorization) {
    symbolic::CASContext ctx;
    
    // (x^2 - 2)(x^2 - 3) = x^4 - 5x^2 + 6
    IntPoly f = make_poly({6, 0, -5, 0, 1});
    auto res = factorize_univariate_hensel_or_kronecker(f, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_TRUE(check_factorization(f, res.value()));
}

TEST(WangEezKroneckerTest, HenselDispatcherFallbackToKronecker) {
    symbolic::CASContext ctx;
    
    // Force Hensel lifting attempts to 0 so it immediately falls back to Kronecker
    ctx.set_max_hensel_lift_attempts(0U);
    
    // (x^2 + x + 1)(x^2 + 2x + 2) = x^4 + 3x^3 + 5x^2 + 4x + 2
    IntPoly f = make_poly({2, 4, 5, 3, 1});
    auto res = factorize_univariate_hensel_or_kronecker(f, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_TRUE(check_factorization(f, res.value()));
}

TEST(WangEezKroneckerTest, HenselDispatcherBailoutOnLargeDegree) {
    symbolic::CASContext ctx;
    
    // Force Hensel lifting attempts to 0, and limit Kronecker max degree to 2
    ctx.set_max_hensel_lift_attempts(0U);
    ctx.set_kronecker_max_degree(2U);
    
    // (x^2 + x + 1)(x^2 + 2x + 2) = x^4 + 3x^3 + 5x^2 + 4x + 2 (degree 4 > 2)
    IntPoly f = make_poly({2, 4, 5, 3, 1});
    auto res = factorize_univariate_hensel_or_kronecker(f, ctx);
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
}
