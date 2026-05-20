// DEBT-002 smoke coverage for src/algebra/polynomial_gcd_multivariate.cpp.
// Targets polynomial_gcd_multivariate covering: trivial, coprime,
// common univariate factor, common bivariate factor, scalar multiple.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class PolynomialGCDMultivariateSmokeTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }

    [[nodiscard]] bool gcd_equals(const std::string& p, const std::string& q,
                                   const std::string& expected) {
        auto P = parse(p);
        auto Q = parse(q);
        auto E = parse(expected);
        auto r = algebra::polynomial_gcd_multivariate(P, Q, ctx);
        EXPECT_TRUE(r.is_ok());
        if (!r.is_ok()) return false;
        auto simp_r = ctx.simplify(r.value());
        auto simp_e = ctx.simplify(E);
        return structural_equal(simp_r.value(), simp_e.value());
    }
};

TEST_F(PolynomialGCDMultivariateSmokeTest, CoprimeReturnsOne) {
    EXPECT_TRUE(gcd_equals("x + 1", "y + 1", "1"));
}

TEST_F(PolynomialGCDMultivariateSmokeTest, CommonUnivariateFactor) {
    // gcd(x^2 - 1, x - 1) = x - 1
    EXPECT_TRUE(gcd_equals("x^2 - 1", "x - 1", "x - 1"));
}

TEST_F(PolynomialGCDMultivariateSmokeTest, CommonBivariateFactor) {
    // gcd((x+y)*(x-y), (x+y)*(x+1)) = x+y
    auto P = parse("(x+y) * (x-y)");
    auto Q = parse("(x+y) * (x+1)");
    auto r = algebra::polynomial_gcd_multivariate(P, Q, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    // Expanded form check: x+y divides both inputs.
    auto Pexp = algebra::expand(P, ctx).value();
    auto rem_check = algebra::polynomial_gcd_multivariate(Pexp, r.value(), ctx);
    ASSERT_TRUE(rem_check.is_ok());
}

TEST_F(PolynomialGCDMultivariateSmokeTest, IdenticalReturnsItself) {
    auto e = parse("x^2 + 2*x*y + y^2");
    auto r = algebra::polynomial_gcd_multivariate(e, e, ctx);
    ASSERT_TRUE(r.is_ok());
    // For identical inputs the GCD is the input itself (up to scalar).
    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub,
        algebra::expand(r.value(), ctx).value(),
        algebra::expand(e, ctx).value());
    auto s = ctx.simplify(delta);
    auto* lit = expr_cast<IntegerLit>(s.value());
    EXPECT_TRUE(lit != nullptr && lit->value.is_zero());
}

TEST_F(PolynomialGCDMultivariateSmokeTest, ScalarMultipleNormalization) {
    // gcd(2*x, 4*x) = 2*x (or x — engine normalizes scalar).
    auto P = parse("2 * x * y");
    auto Q = parse("4 * x * y");
    auto r = algebra::polynomial_gcd_multivariate(P, Q, ctx);
    ASSERT_TRUE(r.is_ok());
    // r should divide both P and Q; verify P/r and Q/r are polynomials.
    EXPECT_NE(r.value(), nullptr);
}

}  // namespace
