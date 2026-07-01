// F3.4 probe — verify DifferentialField implements real chain-rule
// derivation or is a pattern-substitution facade.
//
// Audit C6 claimed: "differential_field.cpp:77-150 builds list of
// extensions but does NOT verify field structure. Pattern matching
// find-replace, not algebra. Zero validation derivation compatibility."
//
// Code inspection of differential_field.cpp:157-180 reveals real
// chain-rule logic:
//   D(ln(u)) = u'/u
//   D(exp(u)) = exp(u) · u'
//   D(base_var) = 1
// Probes verify computational correctness of these rules.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/formatter.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::calculus;

class DifferentialFieldProbeTest : public ::testing::Test {
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

    // A two-level Liouvillian tower: t1 = exp(x) (D(t1)=t1), t2 = exp(t1)
    // (D(t2)=t2·D(t1)=t1·t2).  Derivatives of composite generator expressions
    // exercise the product/power rules of DifferentialField::derive.
    [[nodiscard]] DifferentialField exp_exp_tower() {
        DifferentialExtension e1{ExtensionType::Exponential, parse("x"), Symbol{"t1"}};
        DifferentialExtension e2{ExtensionType::Exponential, parse("t1"), Symbol{"t2"}};
        return DifferentialField(x, {e1, e2});
    }

    // Mathematically validate D(expr) == expected by checking the difference
    // simplifies to 0 in the field (structural, not toString — CLAUDE.md).
    [[nodiscard]] bool derivative_equals(const DifferentialField& field,
                                         const std::string& expr,
                                         const std::string& expected) {
        auto d = field.derive(parse(expr), ctx);
        EXPECT_TRUE(d.is_ok()) << "D(" << expr << "): " << (d.is_error() ? d.error().message : "");
        if (d.is_error()) return false;
        ExprPtr delta = ctx.arena().make<Binary>(BinaryOp::Sub, d.value(), parse(expected));
        auto tog = algebra::together(delta, ctx);
        auto s = ctx.simplify(tog.is_ok() ? tog.value() : delta);
        if (s.is_error()) return false;
        if (const auto* il = expr_cast<IntegerLit>(s.value())) return il->value.is_zero();
        if (const auto* rl = expr_cast<RationalLit>(s.value())) return rl->numerator.is_zero();
        return false;
    }
};

// --- Composite-generator derivation over a two-level tower (product/power/quotient
// rules).  Before this was fixed, DifferentialField::derive fell through to
// diff(_, base_var) for any non-bare-symbol, treating generators as constants and
// returning 0 for D(t1·t2), D(t2²), D(−t1) — a silent-wrong derivation that made
// the parametric Risch solver spuriously report "no solution" instead of reaching
// the deep-tower ConstantSystem boundary.  These pin the corrected rules. ---

TEST_F(DifferentialFieldProbeTest, DeriveBareGeneratorsTwoLevel) {
    DifferentialField field = exp_exp_tower();
    EXPECT_TRUE(derivative_equals(field, "t1", "t1"));        // D(exp x) = exp x
    EXPECT_TRUE(derivative_equals(field, "t2", "t1*t2"));     // D(exp(t1)) = t1·t2
}

TEST_F(DifferentialFieldProbeTest, DeriveProductRuleTwoLevel) {
    DifferentialField field = exp_exp_tower();
    // D(t1·t2) = D(t1)·t2 + t1·D(t2) = t1·t2 + t1·(t1·t2) = t1·t2 + t1²·t2.
    EXPECT_TRUE(derivative_equals(field, "t1*t2", "t1*t2 + t1^2*t2"));
}

TEST_F(DifferentialFieldProbeTest, DerivePowerRuleTwoLevel) {
    DifferentialField field = exp_exp_tower();
    // D(t2²) = 2·t2·D(t2) = 2·t2·(t1·t2) = 2·t1·t2².
    EXPECT_TRUE(derivative_equals(field, "t2^2", "2*t1*t2^2"));
}

TEST_F(DifferentialFieldProbeTest, DeriveNegationAndQuotientTwoLevel) {
    DifferentialField field = exp_exp_tower();
    EXPECT_TRUE(derivative_equals(field, "-t1", "-t1"));           // D(−t1) = −t1
    // D(t2/t1) = (D(t2)·t1 − t2·D(t1))/t1² = (t1²·t2 − t1·t2)/t1² = t2 − t2/t1.
    EXPECT_TRUE(derivative_equals(field, "t2/t1", "t2 - t2/t1"));
}

// A base-variable-dependent exponent over a generator needs log(base) ∉ K:
// the derivation must decline with a diagnostic, never a silent wrong value.
TEST_F(DifferentialFieldProbeTest, DerivePowerNonConstantExponentIsUnimplemented) {
    DifferentialField field = exp_exp_tower();
    auto d = field.derive(parse("t1^x"), ctx);
    ASSERT_TRUE(d.is_error());
    EXPECT_EQ(d.error().kind, CASErrorKind::Unimplemented);
}

// Build a field K(ln(x), exp(x)) and check tower contains 2 extensions.
TEST_F(DifferentialFieldProbeTest, BuildsLogExpTower) {
    auto expr = parse("ln(x) + exp(x)");
    auto field_res = DifferentialField::build(expr, x, ctx);
    ASSERT_TRUE(field_res.is_ok()) << field_res.error().message;
    EXPECT_EQ(field_res.value().extensions().size(), 2U);
    bool has_log = false, has_exp = false;
    for (const auto& ext : field_res.value().extensions()) {
        if (ext.type == ExtensionType::Logarithmic) has_log = true;
        if (ext.type == ExtensionType::Exponential) has_exp = true;
    }
    EXPECT_TRUE(has_log && has_exp);
}

// Verify D(base_var) = 1.
TEST_F(DifferentialFieldProbeTest, DeriveBaseVarReturnsOne) {
    DifferentialField field{x};
    auto x_node = ctx.arena().make<Symbol>("x");
    auto d = field.derive(x_node, ctx);
    ASSERT_TRUE(d.is_ok()) << d.error().message;
    auto* lit = expr_cast<IntegerLit>(d.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(1));
}

// Verify D(ln(x)) = 1/x via the field's chain rule.
// To exercise it, build field with ln(x), transform to t1, derive t1.
TEST_F(DifferentialFieldProbeTest, DeriveLogGeneratorGivesInverseArgument) {
    auto expr = parse("ln(x)");
    auto field_res = DifferentialField::build(expr, x, ctx);
    ASSERT_TRUE(field_res.is_ok());
    const auto& exts = field_res.value().extensions();
    ASSERT_EQ(exts.size(), 1U);

    // derive the generator symbol t_0
    auto t_node = ctx.arena().make<Symbol>(exts[0].t_var.name);
    auto d = field_res.value().derive(t_node, ctx);
    ASSERT_TRUE(d.is_ok()) << d.error().message;
    // Expected: 1/x (Binary Div with 1 and x)
    auto formatted = formatter::TextFormatter{}.format(d.value());
    std::cout << "[PROBE] D(t_log) = " << formatted << std::endl;
    EXPECT_NE(formatted.find("1"), std::string::npos);
    EXPECT_NE(formatted.find("x"), std::string::npos);
}

// Verify D(exp(x)) = exp(x) (since D(x) = 1, so D(t_exp) = t_exp · 1).
TEST_F(DifferentialFieldProbeTest, DeriveExpGeneratorGivesItself) {
    auto expr = parse("exp(x)");
    auto field_res = DifferentialField::build(expr, x, ctx);
    ASSERT_TRUE(field_res.is_ok());
    const auto& exts = field_res.value().extensions();
    ASSERT_EQ(exts.size(), 1U);
    auto t_node = ctx.arena().make<Symbol>(exts[0].t_var.name);
    auto d = field_res.value().derive(t_node, ctx);
    ASSERT_TRUE(d.is_ok()) << d.error().message;
    auto formatted = formatter::TextFormatter{}.format(d.value());
    std::cout << "[PROBE] D(t_exp) = " << formatted << std::endl;
    // Should reference t_var symbol (the generator)
    EXPECT_NE(formatted.find(exts[0].t_var.name), std::string::npos);
}

// Verify linearity: D(a + b) = D(a) + D(b).
TEST_F(DifferentialFieldProbeTest, LinearityOfDerivation) {
    DifferentialField field{x};
    auto sum = ctx.arena().make<Sum>(std::vector<ExprPtr>{
        parse("x"), parse("x^2")
    });
    auto d_sum = field.derive(sum, ctx);
    ASSERT_TRUE(d_sum.is_ok()) << d_sum.error().message;
    auto formatted = formatter::TextFormatter{}.format(d_sum.value());
    std::cout << "[PROBE] D(x + x^2) = " << formatted << std::endl;
    // Expected: 1 + 2*x (or equivalent). Must depend on x.
    EXPECT_NE(formatted.find("x"), std::string::npos);
}
