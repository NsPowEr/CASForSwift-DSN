// test_integrate_exp_rational.cpp — rational functions of exp(x): the
// hyperexponential residue term (Bronstein §5.9).
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Symbolic_Integration_I.md
//   §5.9 (IntegrateHyperexponential): after the residue reduction
//   (Rothstein-Trager), h − Dg₂ is a nonzero element of k in general when the
//   topmost monomial t is hyperexponential (deg_t(Dv) = deg_t(v)), and must be
//   integrated recursively over k.  Dropping it produced the silently wrong
//   antiderivative −ln|eˣ+1| for ∫ 1/(1+eˣ) dx (correct: x − ln(eˣ+1));
//   found by the 2026-07-16 golden refresh (Maxima oracle content diff).
//
// Every accepted antiderivative is verified by differentiation:
// simplify(F' − f) ≡ 0, with multi-point numeric probes as fallback
// (single-point probes can false-negative; see project memory).

#include <gtest/gtest.h>
#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast_debug.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include <string>

using namespace cas;
using namespace cas::calculus;

namespace {

class IntegrateExpRationalTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] ExprPtr parse(const std::string& src) {
        auto tokens = Lexer(src).tokenize();
        if (tokens.is_error()) return nullptr;
        Parser parser(tokens.value(), ctx.arena());
        auto res = parser.parse();
        return res.is_ok() ? res.value() : nullptr;
    }

    [[nodiscard]] static bool literal_zero(ExprPtr e) {
        if (const auto* il = expr_cast<IntegerLit>(e)) return il->value.is_zero();
        if (const auto* rl = expr_cast<RationalLit>(e)) return rl->numerator.is_zero();
        return false;
    }

    // simplify(F' − f) ≡ 0, falling back to numeric probes at several points.
    void expect_derivative_matches(ExprPtr F, ExprPtr integrand,
                                   const std::string& label) {
        Symbol x("x");
        auto Fp = diff(F, x, 1U, ctx);
        ASSERT_TRUE(Fp.is_ok()) << label << ": diff failed: " << Fp.error().message;
        AstArena& a = ctx.arena();
        ExprPtr resid = a.make<Binary>(BinaryOp::Sub, Fp.value(), integrand);
        // together first: rational combinations over a common exp-denominator
        // (e.g. −1/a − eˣ/a + 1, a = eˣ+1) collapse to 0 only after clearing
        // denominators; plain simplify treats eˣ as an opaque atom.
        if (auto tog = algebra::together(resid, ctx); tog.is_ok())
            resid = tog.value();
        auto rs = ctx.simplify(resid);
        ASSERT_TRUE(rs.is_ok()) << label << ": " << rs.error().message;
        if (literal_zero(rs.value())) return;
        // Multi-point probes: a residual that vanishes at all of them and is
        // built from exp/rational atoms only is identically zero in practice.
        for (const char* probe : {"1/2", "1", "-3/2"}) {
            ExprPtr pt = parse(probe);
            ASSERT_NE(pt, nullptr);
            auto sub = ctx.substitute(rs.value(), x, pt);
            ASSERT_TRUE(sub.is_ok()) << label << ": " << sub.error().message;
            if (auto tog = algebra::together(sub.value(), ctx); tog.is_ok())
                sub = ok(tog.value());
            auto ns = ctx.simplify(sub.value());
            ASSERT_TRUE(ns.is_ok()) << label << ": " << ns.error().message;
            EXPECT_TRUE(literal_zero(ns.value()))
                << label << ": F' - integrand ≠ 0 at x=" << probe
                << ";  residual: " << debug_print(ns.value())
                << ";  F = " << debug_print(F);
            if (!literal_zero(ns.value())) return;  // one failure is enough
        }
    }

    // The integral MUST succeed and be correct.
    void expect_integral_correct(const std::string& integrand_src) {
        Symbol x("x");
        ExprPtr integrand = parse(integrand_src);
        ASSERT_NE(integrand, nullptr) << "parse failure: " << integrand_src;
        auto F = integrate(integrand, x, ctx);
        ASSERT_TRUE(F.is_ok())
            << "integrate failed for " << integrand_src
            << ":  " << F.error().message;
        expect_derivative_matches(F.value(), integrand, integrand_src);
    }

    // The integral may cleanly bail (Unimplemented) but MUST NOT be wrong.
    void expect_correct_or_clean_bail(const std::string& integrand_src) {
        Symbol x("x");
        ExprPtr integrand = parse(integrand_src);
        ASSERT_NE(integrand, nullptr) << "parse failure: " << integrand_src;
        auto F = integrate(integrand, x, ctx);
        if (F.is_error()) {
            EXPECT_EQ(F.error().kind, CASErrorKind::Unimplemented)
                << integrand_src << ": non-Unimplemented failure: "
                << F.error().message;
            return;
        }
        expect_derivative_matches(F.value(), integrand, integrand_src);
    }
};

// ─── The golden-refresh bug and its immediate family (T1) ───────────────────

// ∫ 1/(1+eˣ) dx = x − ln(1+eˣ).  Engine used to return −ln|eˣ+1| (wrong).
TEST_F(IntegrateExpRationalTest, OneOverOnePlusExpX) {
    expect_integral_correct("1/(1 + exp(x))");
}

// ∫ eˣ/(1+eˣ) dx = ln(1+eˣ) — pure log part, residue 0.  Regression guard.
TEST_F(IntegrateExpRationalTest, ExpOverOnePlusExpX) {
    expect_integral_correct("exp(x)/(1 + exp(x))");
}

// ─── Systematic sweep of the 1/(a+b·e^{cx}) family (T2) ─────────────────────

TEST_F(IntegrateExpRationalTest, OneOverTwoPlusThreeExpX) {
    expect_correct_or_clean_bail("1/(2 + 3*exp(x))");
}

TEST_F(IntegrateExpRationalTest, OneOverOneMinusExpX) {
    expect_correct_or_clean_bail("1/(1 - exp(x))");
}

TEST_F(IntegrateExpRationalTest, OneOverOnePlusExpNegX) {
    expect_correct_or_clean_bail("1/(1 + exp(-x))");
}

TEST_F(IntegrateExpRationalTest, OneOverOnePlusExpTwoX) {
    expect_correct_or_clean_bail("1/(1 + exp(2*x))");
}

// Proper rational with quadratic denominator in t = eˣ:
// ∫ 1/((1+eˣ)·(2+eˣ)) dx — partial fractions over k(t), both factors normal.
TEST_F(IntegrateExpRationalTest, OneOverProductOfTwoExpLinears) {
    expect_correct_or_clean_bail("1/((1 + exp(x))*(2 + exp(x)))");
}

// Improper in t (deg num = deg den): ∫ (1+2eˣ)/(1+eˣ) dx.
TEST_F(IntegrateExpRationalTest, ImproperRationalInExp) {
    expect_correct_or_clean_bail("(1 + 2*exp(x))/(1 + exp(x))");
}

} // namespace
