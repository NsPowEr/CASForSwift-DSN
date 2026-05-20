// L0-14 STEP 2 — DecimalLit→RationalLit conversion at parser boundary.
// Verifies regola CLAUDE.md "DecimalLit ammessi solo per preservare
// input utente; conversione esatta avviene al parser, non nel core".
#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class DecimalToRationalAtParserTest : public ::testing::Test {
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

TEST_F(DecimalToRationalAtParserTest, SimpleHalf) {
    auto e = parse("0.5");
    auto* rat = expr_cast<RationalLit>(e);
    ASSERT_NE(rat, nullptr) << "0.5 should parse to RationalLit, not DecimalLit";
    EXPECT_EQ(rat->numerator, BigInt(1));
    EXPECT_EQ(rat->denominator, BigInt(2));
}

TEST_F(DecimalToRationalAtParserTest, OneQuarter) {
    auto e = parse("0.25");
    auto* rat = expr_cast<RationalLit>(e);
    ASSERT_NE(rat, nullptr);
    EXPECT_EQ(rat->numerator, BigInt(1));
    EXPECT_EQ(rat->denominator, BigInt(4));
}

TEST_F(DecimalToRationalAtParserTest, OneTenthExact) {
    // 0.1 is exact in rational form (1/10), not lossy as float.
    auto e = parse("0.1");
    auto* rat = expr_cast<RationalLit>(e);
    ASSERT_NE(rat, nullptr);
    EXPECT_EQ(rat->numerator, BigInt(1));
    EXPECT_EQ(rat->denominator, BigInt(10));
}

TEST_F(DecimalToRationalAtParserTest, NegativeDecimal) {
    auto e = parse("-1.5");
    // Could be Unary(Neg, RationalLit) or directly negative RationalLit.
    if (auto* rat = expr_cast<RationalLit>(e)) {
        EXPECT_EQ(rat->numerator, BigInt(-3));
        EXPECT_EQ(rat->denominator, BigInt(2));
    } else if (auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
        auto* inner = expr_cast<RationalLit>(un->operand);
        ASSERT_NE(inner, nullptr);
        EXPECT_EQ(inner->numerator, BigInt(3));
        EXPECT_EQ(inner->denominator, BigInt(2));
    } else {
        FAIL() << "expected RationalLit or Unary(Neg, RationalLit)";
    }
}

TEST_F(DecimalToRationalAtParserTest, IntegerLikeFloatReduces) {
    // 1.0 → 1/1 (or IntegerLit(1) after simplify) — either is acceptable.
    auto e = parse("1.0");
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok());
    bool ok = false;
    if (auto* il = expr_cast<IntegerLit>(simp.value())) {
        ok = (il->value == BigInt(1));
    } else if (auto* rl = expr_cast<RationalLit>(simp.value())) {
        ok = (rl->numerator == BigInt(1) && rl->denominator == BigInt(1));
    }
    EXPECT_TRUE(ok);
}

TEST_F(DecimalToRationalAtParserTest, DiffOnDecimalNowWorks) {
    // Pre-L0-14: 0.5*x^2 emitted DecimalLit → diff returns Unimplemented.
    // Post-L0-14: 0.5*x^2 emits RationalLit(1/2) → diff returns x.
    auto e = parse("0.5 * x^2");
    auto d = calculus::diff(e, x, 1U, ctx);
    ASSERT_TRUE(d.is_ok()) << "diff failed on parsed 0.5*x^2: " << (d.is_error() ? d.error().message : "");
    // Verify result equivalent to x (could be x or RationalLit(1)*x, etc.)
    auto x_expr = ctx.arena().make<Symbol>("x");
    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, d.value(), x_expr);
    auto delta_simp = ctx.simplify(delta);
    ASSERT_TRUE(delta_simp.is_ok());
    auto* lit = expr_cast<IntegerLit>(delta_simp.value());
    EXPECT_TRUE(lit != nullptr && lit->value.is_zero())
        << "D(0.5*x^2) should equal x";
}

TEST_F(DecimalToRationalAtParserTest, IntegrateOnDecimalNowWorks) {
    auto e = parse("0.25 * x");
    auto r = calculus::integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok()) << "integrate failed: " << (r.is_error() ? r.error().message : "");
    // Verify D(result) = integrand
    auto d_back = calculus::diff(r.value(), x, 1U, ctx);
    ASSERT_TRUE(d_back.is_ok());
    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, d_back.value(), e);
    auto t = algebra::together(delta, ctx);
    ASSERT_TRUE(t.is_ok());
    auto simp = ctx.simplify(t.value());
    auto* lit = expr_cast<IntegerLit>(simp.value());
    EXPECT_TRUE(lit != nullptr && lit->value.is_zero());
}

TEST_F(DecimalToRationalAtParserTest, AntiHardcodeRandomFiniteDecimals) {
    // Verify generality: any finite decimal d.ddd parses to exactly num/10^k.
    struct Case { const char* str; long long num; long long den; };
    Case cases[] = {
        {"3.14", 157, 50},      // 314/100 = 157/50
        {"0.001", 1, 1000},
        {"100.0", 100, 1},
        {"0.000625", 1, 1600},  // 625/1000000 = 1/1600
    };
    for (auto& c : cases) {
        auto e = parse(c.str);
        auto* rat = expr_cast<RationalLit>(e);
        if (!rat) {
            // Could be IntegerLit if denominator reduces to 1
            auto simp = ctx.simplify(e);
            if (auto* il = expr_cast<IntegerLit>(simp.value())) {
                EXPECT_EQ(il->value, BigInt(c.num));
                EXPECT_EQ(c.den, 1);
                continue;
            }
            rat = expr_cast<RationalLit>(simp.value());
        }
        ASSERT_NE(rat, nullptr) << c.str;
        EXPECT_EQ(rat->numerator, BigInt(c.num)) << c.str;
        EXPECT_EQ(rat->denominator, BigInt(c.den)) << c.str;
    }
}

}  // namespace
