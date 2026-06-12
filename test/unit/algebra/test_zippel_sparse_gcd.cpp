#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class ZippelSparseGcdTest : public ::testing::Test {
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

TEST_F(ZippelSparseGcdTest, SparseBivariate) {
    ctx.set_zippel_density_threshold(2.0);
    // gcd(x^5 + y^5, x^5 + y^5) = x^5 + y^5 (which is sparse compared to full bivariate deg 5 polynomial)
    EXPECT_TRUE(gcd_equals("x^5 + y^5", "x^5 + y^5", "x^5 + y^5"));
}

TEST_F(ZippelSparseGcdTest, DenseBivariateFallback) {
    ctx.set_zippel_density_threshold(100.0); // very high, bypassing Zippel
    // gcd(x^2 + x*y + y^2, x^2 + 2*x*y + y^2) = 1
    EXPECT_TRUE(gcd_equals("x^2 + x*y + y^2", "x^2 + 2*x*y + y^2", "1"));
}

}  // namespace
