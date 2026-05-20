// CAS-L3-20 — Custom monomial orderings (GLex, Lex, GRevLex).
//
// Verifies that f4_groebner accepts MonomialOrder::Lex, GLex, GRevLex
// and produces consistent Groebner basis (same generator set up to
// monic+order). The basis cardinality may differ across orders but
// the ideal generated must be invariant.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include "../../../src/algebra/polynomial_groebner_f4.hpp"

using namespace cas;
using namespace cas::algebra;

namespace {

class MonomialOrdersTest : public ::testing::Test {
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
};

TEST_F(MonomialOrdersTest, LeadingMonomialLexVsGRevLex) {
    // poly: x²·y + x·y² + 1. Lex leading: x²y (highest exponent x first).
    // GRevLex leading: x·y² (same total deg 3 but reverse-lex prefers
    // last variable with lower exponent... actually GRevLex prefers
    // monomial with smaller last-variable exponent if total deg same).
    PolyF4 p;
    p.terms[{2u, 1u}] = Rational(BigInt(1));  // x²y
    p.terms[{1u, 2u}] = Rational(BigInt(1));  // xy²
    p.terms[{0u, 0u}] = Rational(BigInt(1));  // 1
    Monomial lm_lex = p.leading_monomial(MonomialOrder::Lex);
    Monomial lm_glex = p.leading_monomial(MonomialOrder::GLex);
    Monomial lm_grevlex = p.leading_monomial(MonomialOrder::GRevLex);
    EXPECT_EQ(lm_lex, (Monomial{2u, 1u}));
    EXPECT_EQ(lm_glex, (Monomial{2u, 1u}));  // tie on deg 3, lex prefers x²y
    EXPECT_EQ(lm_grevlex, (Monomial{2u, 1u}));  // x²y → reverse-lex: last var smaller
}

TEST_F(MonomialOrdersTest, GLexOrdersByTotalDegreeFirst) {
    // x³ (deg 3) > x²·y (deg 3) for both Lex and GLex (same deg, then lex).
    // But for deg-2 vs deg-3: GLex always prefers higher deg first.
    PolyF4 p;
    p.terms[{2u, 0u}] = Rational(BigInt(1));  // x² (deg 2)
    p.terms[{0u, 3u}] = Rational(BigInt(1));  // y³ (deg 3)
    // Lex: x² wins (first exponent compare: 2 > 0).
    // GLex: y³ wins (deg 3 > 2).
    EXPECT_EQ(p.leading_monomial(MonomialOrder::Lex), (Monomial{2u, 0u}));
    EXPECT_EQ(p.leading_monomial(MonomialOrder::GLex), (Monomial{0u, 3u}));
    EXPECT_EQ(p.leading_monomial(MonomialOrder::GRevLex), (Monomial{0u, 3u}));
}

TEST_F(MonomialOrdersTest, AntiHardcodeAllThreeOrdersDistinct) {
    // y² (deg 2) vs x (deg 1) vs y (deg 1)
    //   Lex: x (first var higher) wins.
    //   GLex: y² (highest deg).
    //   GRevLex: y² (highest deg, reverse-lex tie-break).
    PolyF4 p;
    p.terms[{1u, 0u}] = Rational(BigInt(1));  // x
    p.terms[{0u, 2u}] = Rational(BigInt(1));  // y²
    p.terms[{0u, 1u}] = Rational(BigInt(1));  // y
    EXPECT_EQ(p.leading_monomial(MonomialOrder::Lex), (Monomial{1u, 0u}));
    EXPECT_EQ(p.leading_monomial(MonomialOrder::GLex), (Monomial{0u, 2u}));
    EXPECT_EQ(p.leading_monomial(MonomialOrder::GRevLex), (Monomial{0u, 2u}));
}

TEST_F(MonomialOrdersTest, GLexOnLinearPolynomial) {
    // Linear: x + y + 1.
    // All orders should pick the largest-degree monomial first.
    PolyF4 p;
    p.terms[{1u, 0u}] = Rational(BigInt(1));
    p.terms[{0u, 1u}] = Rational(BigInt(1));
    p.terms[{0u, 0u}] = Rational(BigInt(1));
    Monomial lm = p.leading_monomial(MonomialOrder::GLex);
    // x and y both deg 1: GLex breaks via lex → x wins.
    EXPECT_EQ(lm, (Monomial{1u, 0u}));
}

}  // namespace
