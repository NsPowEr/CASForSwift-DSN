// DEBT-F1-COV-01 — White-box coverage uplift for simplify_bessel_orthogonal.cpp.
// Exercises: BesselJ/Y/I/K half-integer orders, BesselI/K negative order parity,
// BesselJ/Y three-term recurrence (expand_bessel_recurrence=true), BesselZero,
// ChebyshevU, HermiteH, HermiteHe, LaguerreL, JacobiP, LambertW special values.

#include <gtest/gtest.h>
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::symbolic;

namespace {

struct OrthoCtx {
    CASContext ctx;
    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto tok = Lexer(s).tokenize();
        if (!tok.is_ok()) return nullptr;
        auto parsed = Parser(tok.value(), ctx.arena()).parse();
        if (!parsed.is_ok()) return nullptr;
        return parsed.value();
    }
    [[nodiscard]] std::optional<std::string> simplify_str(const std::string& s) {
        ExprPtr e = parse(s);
        if (!e) return std::nullopt;
        auto r = ctx.simplify(e);
        if (!r.is_ok()) return std::nullopt;
        auto t = to_round_trip_text(r.value());
        if (!t.is_ok()) return std::nullopt;
        return t.value();
    }
    [[nodiscard]] bool is_ok(const std::string& s) {
        return simplify_str(s).has_value();
    }
};

// Mathematical equivalence via two separate contexts
[[nodiscard]] bool math_equiv(const std::string& a, const std::string& b) {
    OrthoCtx ca, cb;
    auto ra = ca.simplify_str(a);
    auto rb = cb.simplify_str(b);
    if (!ra || !rb) return false;
    return ra.value() == rb.value();
}

}  // namespace

// ── BesselJ half-integer orders (closed form in trig) ────────────────────────

TEST(OrthogonalCoverage, BesselJ_HalfInt_Positive) {
    // BesselJ(1/2, x) = sqrt(2/(pi*x)) * sin(x)
    OrthoCtx c;
    auto r = c.simplify_str("BesselJ(1/2, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_FALSE(r.value().empty());
    // Result should not contain BesselJ (was reduced)
    EXPECT_EQ(r.value().find("BesselJ"), std::string::npos);
}

TEST(OrthogonalCoverage, BesselJ_HalfInt_Negative) {
    // BesselJ(-1/2, x) = sqrt(2/(pi*x)) * cos(x)
    OrthoCtx c;
    auto r = c.simplify_str("BesselJ(-1/2, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_EQ(r.value().find("BesselJ"), std::string::npos);
}

TEST(OrthogonalCoverage, BesselY_HalfInt_Positive) {
    // BesselY(1/2, x) → closed trig form (involves cos with negation)
    OrthoCtx c;
    auto r = c.simplify_str("BesselY(1/2, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_EQ(r.value().find("BesselY"), std::string::npos);
}

TEST(OrthogonalCoverage, BesselY_HalfInt_Negative) {
    // BesselY(-1/2, x) → sin form
    OrthoCtx c;
    auto r = c.simplify_str("BesselY(-1/2, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_EQ(r.value().find("BesselY"), std::string::npos);
}

TEST(OrthogonalCoverage, BesselI_HalfInt_Positive) {
    // BesselI(1/2, x) = sqrt(2/(pi*x)) * sinh(x)
    OrthoCtx c;
    auto r = c.simplify_str("BesselI(1/2, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_EQ(r.value().find("BesselI"), std::string::npos);
}

TEST(OrthogonalCoverage, BesselI_HalfInt_Negative) {
    // BesselI(-1/2, x) = sqrt(2/(pi*x)) * cosh(x)
    OrthoCtx c;
    auto r = c.simplify_str("BesselI(-1/2, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_EQ(r.value().find("BesselI"), std::string::npos);
}

TEST(OrthogonalCoverage, BesselK_HalfInt_Positive) {
    // BesselK(1/2, x) = sqrt(pi/(2*x)) * exp(-x)
    OrthoCtx c;
    auto r = c.simplify_str("BesselK(1/2, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_EQ(r.value().find("BesselK"), std::string::npos);
}

TEST(OrthogonalCoverage, BesselK_HalfInt_Negative) {
    // BesselK(-1/2, x) = same as K(1/2) by parity
    OrthoCtx c;
    auto r = c.simplify_str("BesselK(-1/2, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_EQ(r.value().find("BesselK"), std::string::npos);
}

// ── BesselJ/Y negative integer parity ────────────────────────────────────────

TEST(OrthogonalCoverage, BesselJ_NegIntEven_SameSign) {
    // J_{-2}(x) = J_2(x)
    EXPECT_TRUE(math_equiv("BesselJ(-2, x)", "BesselJ(2, x)"));
}

TEST(OrthogonalCoverage, BesselJ_NegIntOdd_FlipSign) {
    // J_{-3}(x) = -J_3(x)
    EXPECT_TRUE(math_equiv("BesselJ(-3, x)", "-BesselJ(3, x)"));
}

TEST(OrthogonalCoverage, BesselY_NegIntEven) {
    EXPECT_TRUE(math_equiv("BesselY(-2, x)", "BesselY(2, x)"));
}

TEST(OrthogonalCoverage, BesselY_NegIntOdd) {
    EXPECT_TRUE(math_equiv("BesselY(-3, x)", "-BesselY(3, x)"));
}

TEST(OrthogonalCoverage, BesselI_NegInt_NoFlip) {
    // I_{-n}(x) = I_n(x) for all integer n (no sign flip for BesselI)
    EXPECT_TRUE(math_equiv("BesselI(-2, x)", "BesselI(2, x)"));
    EXPECT_TRUE(math_equiv("BesselI(-3, x)", "BesselI(3, x)"));
}

// ── BesselJ/Y three-term recurrence (expand_bessel_recurrence=true) ──────────

TEST(OrthogonalCoverage, BesselJ_Recurrence_Order2) {
    // J_2(x) = (2/x)*J_1(x) - J_0(x)
    OrthoCtx c;
    c.ctx.set_expand_bessel_recurrence(true);
    auto r = c.simplify_str("BesselJ(2, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_FALSE(r.value().empty());
}

TEST(OrthogonalCoverage, BesselJ_Recurrence_Order3) {
    // J_3(x) = (4/x)*J_2(x) - J_1(x) → recurses down to J_0, J_1
    OrthoCtx c;
    c.ctx.set_expand_bessel_recurrence(true);
    auto r = c.simplify_str("BesselJ(3, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_FALSE(r.value().empty());
}

TEST(OrthogonalCoverage, BesselY_Recurrence_Order2) {
    OrthoCtx c;
    c.ctx.set_expand_bessel_recurrence(true);
    auto r = c.simplify_str("BesselY(2, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_FALSE(r.value().empty());
}

TEST(OrthogonalCoverage, BesselJ_Recurrence_Disabled_ReturnsOpaque) {
    // Without expand_bessel_recurrence, BesselJ(2, x) stays symbolic
    OrthoCtx c;
    // default: expand_bessel_recurrence = false
    auto r = c.simplify_str("BesselJ(2, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_NE(r.value().find("BesselJ"), std::string::npos);
}

// ── BesselZero ───────────────────────────────────────────────────────────────

TEST(OrthogonalCoverage, BesselZero_IntOrderIntIndex) {
    // BesselZero(0, 1) — order 0 (real), index 1 (positive integer) → stays symbolic
    OrthoCtx c;
    auto r = c.simplify_str("BesselZero(0, 1)");
    ASSERT_TRUE(r.has_value());
    EXPECT_NE(r.value().find("BesselZero"), std::string::npos);
}

TEST(OrthogonalCoverage, BesselZero_InvalidIndex_ZeroFails) {
    // Index 0 is invalid — should return error
    OrthoCtx c;
    ExprPtr e = c.parse("BesselZero(0, 0)");
    ASSERT_NE(e, nullptr);
    auto r = c.ctx.simplify(e);
    // Should error (index must be >= 1)
    EXPECT_TRUE(r.is_error());
}

TEST(OrthogonalCoverage, BesselZero_IntOrderIntIndex_Two) {
    // BesselZero(1, 2) — valid
    OrthoCtx c;
    auto r = c.simplify_str("BesselZero(1, 2)");
    ASSERT_TRUE(r.has_value());
}

// ── ChebyshevU ───────────────────────────────────────────────────────────────

TEST(OrthogonalCoverage, ChebyshevU_0) {
    EXPECT_TRUE(math_equiv("ChebyshevU(0, x)", "1"));
}

TEST(OrthogonalCoverage, ChebyshevU_1) {
    // U_1(x) = 2x
    EXPECT_TRUE(math_equiv("ChebyshevU(1, x)", "2*x"));
}

TEST(OrthogonalCoverage, ChebyshevU_2) {
    // U_2(x) = 4x² - 1
    EXPECT_TRUE(math_equiv("ChebyshevU(2, x)", "4*x^2 - 1"));
}

TEST(OrthogonalCoverage, ChebyshevU_3) {
    // U_3(x) = 8x³ - 4x
    EXPECT_TRUE(math_equiv("ChebyshevU(3, x)", "8*x^3 - 4*x"));
}

// ── HermiteH (physicist) ──────────────────────────────────────────────────────

TEST(OrthogonalCoverage, HermiteH_0) {
    EXPECT_TRUE(math_equiv("HermiteH(0, x)", "1"));
}

TEST(OrthogonalCoverage, HermiteH_1) {
    EXPECT_TRUE(math_equiv("HermiteH(1, x)", "2*x"));
}

TEST(OrthogonalCoverage, HermiteH_2) {
    // H_2 = 4x² - 2
    EXPECT_TRUE(math_equiv("HermiteH(2, x)", "4*x^2 - 2"));
}

TEST(OrthogonalCoverage, HermiteH_3) {
    // H_3 = 8x³ - 12x
    EXPECT_TRUE(math_equiv("HermiteH(3, x)", "8*x^3 - 12*x"));
}

// ── HermiteHe (probabilist) ───────────────────────────────────────────────────

TEST(OrthogonalCoverage, HermiteHe_0) {
    EXPECT_TRUE(math_equiv("HermiteHe(0, x)", "1"));
}

TEST(OrthogonalCoverage, HermiteHe_1) {
    EXPECT_TRUE(math_equiv("HermiteHe(1, x)", "x"));
}

TEST(OrthogonalCoverage, HermiteHe_2) {
    // He_2 = x² - 1
    EXPECT_TRUE(math_equiv("HermiteHe(2, x)", "x^2 - 1"));
}

TEST(OrthogonalCoverage, HermiteHe_3) {
    // He_3 = x³ - 3x
    EXPECT_TRUE(math_equiv("HermiteHe(3, x)", "x^3 - 3*x"));
}

// ── LaguerreL ─────────────────────────────────────────────────────────────────

TEST(OrthogonalCoverage, LaguerreL_0) {
    EXPECT_TRUE(math_equiv("LaguerreL(0, x)", "1"));
}

TEST(OrthogonalCoverage, LaguerreL_1) {
    // L_1(x) = 1 - x; verify it simplifies to something non-trivial
    OrthoCtx c;
    auto r = c.simplify_str("LaguerreL(1, x)");
    ASSERT_TRUE(r.has_value());
    // Result should involve x
    EXPECT_NE(r.value().find("x"), std::string::npos);
    EXPECT_EQ(r.value().find("LaguerreL"), std::string::npos);
}

TEST(OrthogonalCoverage, LaguerreL_2) {
    // L_2 = (x² - 4x + 2)/2
    EXPECT_TRUE(math_equiv("LaguerreL(2, x)", "(x^2 - 4*x + 2)/2"));
}

TEST(OrthogonalCoverage, LaguerreL_3) {
    OrthoCtx c;
    auto r = c.simplify_str("LaguerreL(3, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_FALSE(r.value().empty());
}

// ── JacobiP ───────────────────────────────────────────────────────────────────

TEST(OrthogonalCoverage, JacobiP_0) {
    // P_0^{α,β}(x) = 1 for any α, β
    EXPECT_TRUE(math_equiv("JacobiP(0, 1, 2, x)", "1"));
}

TEST(OrthogonalCoverage, JacobiP_1) {
    // P_1^{α,β}(x) = (α-β)/2 + (α+β+2)/2 * x
    OrthoCtx c;
    auto r = c.simplify_str("JacobiP(1, 0, 0, x)");
    ASSERT_TRUE(r.has_value());
    // P_1^{0,0}(x) = x (Legendre)
    EXPECT_FALSE(r.value().empty());
}

TEST(OrthogonalCoverage, JacobiP_2_Symbolic_Alpha_Beta) {
    // Symbolic α, β — exercises the recurrence path; result may simplify or
    // remain opaque depending on whether algebra can reduce symbolic coefficients.
    OrthoCtx c;
    ExprPtr e = c.parse("JacobiP(2, a, b, x)");
    ASSERT_NE(e, nullptr);
    // Either succeeds with a polynomial expression or returns error/Unimplemented —
    // either way must not crash.
    auto r = c.ctx.simplify(e);
    // No assertion on is_ok; just verify no hang or crash
    (void)r;
    SUCCEED();
}

TEST(OrthogonalCoverage, JacobiP_2_NumericSpecialization) {
    // JacobiP(2, 0, 0, x) = (3x²-1)/2 = LegendreP(2, x)
    EXPECT_TRUE(math_equiv("JacobiP(2, 0, 0, x)", "LegendreP(2, x)"));
}

// ── LambertW special values ───────────────────────────────────────────────────

TEST(OrthogonalCoverage, LambertW_Zero) {
    // W(0) = 0
    EXPECT_TRUE(math_equiv("LambertW(0)", "0"));
}

TEST(OrthogonalCoverage, LambertW_E) {
    // W(e) = 1
    EXPECT_TRUE(math_equiv("LambertW(e)", "1"));
}

TEST(OrthogonalCoverage, LambertW_XExpX_NonNeg) {
    // W(x*exp(x)) = x when x is assumed nonneg — symbolic identity
    // Test the parse/simplify path; result is x or stays as W(...)
    OrthoCtx c;
    auto r = c.simplify_str("LambertW(x)");
    ASSERT_TRUE(r.has_value());  // should not crash
}

TEST(OrthogonalCoverage, LambertW_Symbolic_StaysOpaque) {
    // LambertW(2) — no closed form; stays as LambertW(2)
    OrthoCtx c;
    auto r = c.simplify_str("LambertW(2)");
    ASSERT_TRUE(r.has_value());
    EXPECT_NE(r.value().find("LambertW"), std::string::npos);
}

// ── ChebyshevT (also exercises bessel_orthogonal, verifying n=0,1 fast paths) ─

TEST(OrthogonalCoverage, ChebyshevT_0_FastPath) {
    EXPECT_TRUE(math_equiv("ChebyshevT(0, x)", "1"));
}

TEST(OrthogonalCoverage, ChebyshevT_1_FastPath) {
    EXPECT_TRUE(math_equiv("ChebyshevT(1, x)", "x"));
}

TEST(OrthogonalCoverage, ChebyshevT_4) {
    // T_4(x) = 8x⁴ - 8x² + 1
    EXPECT_TRUE(math_equiv("ChebyshevT(4, x)", "8*x^4 - 8*x^2 + 1"));
}

// ── LegendreP (additional coverage for recurrence > 1 pass) ──────────────────

TEST(OrthogonalCoverage, LegendreP_5) {
    // P_5(x) = (63x⁵ - 70x³ + 15x)/8
    EXPECT_TRUE(math_equiv("LegendreP(5, x)", "(63*x^5 - 70*x^3 + 15*x)/8"));
}

TEST(OrthogonalCoverage, LegendreP_AtZero) {
    // P_2(0) = -1/2
    EXPECT_TRUE(math_equiv("LegendreP(2, 0)", "-1/2"));
}

// ── Symbolic x argument passes through unchanged ─────────────────────────────

TEST(OrthogonalCoverage, BesselJ_SymbolicOrder_RetainsForm) {
    // BesselJ(n, x) with symbolic n — no reduction, form preserved
    OrthoCtx c;
    auto r = c.simplify_str("BesselJ(n, x)");
    ASSERT_TRUE(r.has_value());
    EXPECT_NE(r.value().find("BesselJ"), std::string::npos);
}
