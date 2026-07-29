// A6 / CAS-L3-18 — Galois group identification for irreducible septics.
//
// Both expected groups are backed by a theorem, not by an external table:
//   • x⁷ − x − 1: Gal = S₇ by Osada's theorem (the trinomial xⁿ − x − 1 has
//     the full symmetric group for every n). Its discriminant 776887 is not a
//     rational square, consistent with S₇ ⊄ A₇.
//   • x⁷ + x⁶ − 12x⁵ − 7x⁴ + 28x³ + 14x² − 9x + 1 is the minimal polynomial of
//     the Gaussian period of length 4 in Q(ζ₂₉) — i.e. it defines the unique
//     degree-7 subfield of Q(ζ₂₉), which is cyclic because Gal(Q(ζ₂₉)/Q) ≅
//     (Z/29Z)* is cyclic of order 28 and 7 | 28. So Gal = C₇, and indeed its
//     discriminant is 414613², a perfect square (C₇ ⊆ A₇).
//
// The first test in this binary to reach the deg-7 driver pays the one-time
// construction of the transitive lattice of S₇ (memoized process-wide).

#include <gtest/gtest.h>

#include <string>

#include "cas/galois.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::algebra;

namespace {

class GaloisDeg7Test : public ::testing::Test {
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

TEST_F(GaloisDeg7Test, OsadaTrinomial_S7) {
    auto p = parse("x^7 - x - 1");
    auto r = galois_group(p, x, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value(), "S7");
}

TEST_F(GaloisDeg7Test, GaussianPeriodCyclotomic29_C7) {
    auto p = parse("x^7 + x^6 - 12*x^5 - 7*x^4 + 28*x^3 + 14*x^2 - 9*x + 1");
    auto r = galois_group(p, x, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value(), "C7");
}

}  // namespace
