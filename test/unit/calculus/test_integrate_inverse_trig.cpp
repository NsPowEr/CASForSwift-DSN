// test_integrate_inverse_trig.cpp — F7.5.B1 inverse trig + inverse hyperbolic
// standalone integral acceptance tests.
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Risch_Inverse_Trig.md
// Closes the spec's "≥14 cases" acceptance criterion:
//   asin/acos/atan/asinh/acosh/atanh/acoth × {standalone, scaled affine
//   argument, product `x · f(x)` via integration-by-parts}.

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

class IntegrateInverseTrigTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] ExprPtr parse(const std::string& src) {
        auto tokens = Lexer(src).tokenize();
        if (tokens.is_error()) return nullptr;
        Parser parser(tokens.value(), ctx.arena());
        auto res = parser.parse();
        return res.is_ok() ? res.value() : nullptr;
    }

    // Run ∫expr dx and verify by symbolic differentiation: simplify((F)' - expr) ≡ 0
    // (literal-zero check; falls back to single-point numeric residual when
    // the simplifier leaves a structurally non-zero but algebraically zero
    // residual — typical for sqrt expressions).
    void expect_integral_correct(const std::string& integrand_src,
                                 const char* probe_at = "1/2") {
        Symbol x("x");
        ExprPtr integrand = parse(integrand_src);
        ASSERT_NE(integrand, nullptr) << "parse failure: " << integrand_src;
        auto F = integrate(integrand, x, ctx);
        ASSERT_TRUE(F.is_ok())
            << "integrate failed for " << integrand_src
            << ":  " << F.error().message;
        auto Fp = diff(F.value(), x, 1U, ctx);
        ASSERT_TRUE(Fp.is_ok()) << Fp.error().message;
        AstArena& a = ctx.arena();
        ExprPtr resid = a.make<Binary>(BinaryOp::Sub, Fp.value(), integrand);
        auto rs = ctx.simplify(resid);
        ASSERT_TRUE(rs.is_ok()) << rs.error().message;
        if (auto* il = expr_cast<IntegerLit>(rs.value()))
            if (il->value.is_zero()) return;
        if (auto* rl = expr_cast<RationalLit>(rs.value()))
            if (rl->numerator.is_zero()) return;
        // Fall back to numeric probe at x = probe_at.
        ExprPtr pt = parse(probe_at);
        ASSERT_NE(pt, nullptr);
        auto sub = ctx.substitute(rs.value(), x, pt);
        ASSERT_TRUE(sub.is_ok()) << sub.error().message;
        auto ns = ctx.simplify(sub.value());
        ASSERT_TRUE(ns.is_ok()) << ns.error().message;
        bool zero = false;
        if (auto* il = expr_cast<IntegerLit>(ns.value()))
            zero = il->value.is_zero();
        else if (auto* rl = expr_cast<RationalLit>(ns.value()))
            zero = rl->numerator.is_zero();
        EXPECT_TRUE(zero) << "F' - integrand ≠ 0 at probe;  residual: "
            << debug_print(ns.value());
    }
};

// ─── Standalone primitives ──────────────────────────────────────────────────

TEST_F(IntegrateInverseTrigTest, AsinX)  { expect_integral_correct("asin(x)"); }
TEST_F(IntegrateInverseTrigTest, AcosX)  { expect_integral_correct("acos(x)"); }
TEST_F(IntegrateInverseTrigTest, AtanX)  { expect_integral_correct("atan(x)"); }
TEST_F(IntegrateInverseTrigTest, AsinhX) { expect_integral_correct("asinh(x)"); }
TEST_F(IntegrateInverseTrigTest, AcoshX) {
    // acosh derivative has sqrt(x²-1) — probe at x = 2 (real domain).
    expect_integral_correct("acosh(x)", "2");
}
TEST_F(IntegrateInverseTrigTest, AtanhX) {
    // atanh domain |x|<1; probe at x=1/3.
    expect_integral_correct("atanh(x)", "1/3");
}
TEST_F(IntegrateInverseTrigTest, AcothX) {
    // acoth domain |x|>1; probe at x=2.
    expect_integral_correct("acoth(x)", "2");
}

// ─── Affine-argument scaling  ∫f(a·x + b) dx = (1/a)·F(a·x + b) ────────────

TEST_F(IntegrateInverseTrigTest, AsinTwoX) {
    expect_integral_correct("asin(2*x)", "1/4");
}
TEST_F(IntegrateInverseTrigTest, AtanThreeX) {
    expect_integral_correct("atan(3*x)", "1/6");
}
TEST_F(IntegrateInverseTrigTest, AsinhFiveX) {
    expect_integral_correct("asinh(5*x)", "1/10");
}

// ─── Product x · f(x) routed through integration-by-parts (ILATE) ──────────
// ILATE selects the inverse function as `u` (HC-IBP-INVHYP fix in
// integrate_parts.cpp: `BuiltinOp::Unknown` name-match for arc*).  The
// recursive sub-integration of `v·du` (rational in x and 1/(1±x²) etc.)
// must close to elementary form; for some integrands this requires the
// rational-integration path to recognise patterns the current Risch/PFD
// pipeline does not.  Soft-skip when sub-integration fails.

TEST_F(IntegrateInverseTrigTest, XTimesAtanhX) {
    // Verified ok currently via PFD on 1/(1−x²).
    expect_integral_correct("x * atanh(x)", "1/3");
}

// HC-IBP-VDU closure (2026-06-15): the rational-fraction gate +
// `integrate_xsq_over_sqrt_quadratic` helper close the sub-integration of
// `x² · 1/(1±x²)` and `x² · 1/√(1±x²)` shapes IBP produces as v·du.
TEST_F(IntegrateInverseTrigTest, XTimesAtanX) {
    expect_integral_correct("x * atan(x)");
}
TEST_F(IntegrateInverseTrigTest, XTimesAsinX) {
    expect_integral_correct("x * asin(x)", "1/3");
}
TEST_F(IntegrateInverseTrigTest, XTimesAsinhX) {
    expect_integral_correct("x * asinh(x)", "1/2");
}

} // anonymous namespace

// ─── Probe sub-integrand patterns (HC-IBP-VDU diagnostic) ───────────────────
// These exercise the rational/sqrt-of-quadratic paths directly.  Skipped
// when integrate fails so they don't gate the build.
TEST_F(IntegrateInverseTrigTest, Probe_X2_Over_1_Plus_X2) {
    Symbol x("x");
    ExprPtr e = parse("x^2/(1+x^2)");
    ASSERT_NE(e, nullptr);
    auto r = integrate(e, x, ctx);
    if (r.is_error()) GTEST_SKIP() << r.error().message;
    std::cout << "x^2/(1+x^2) → " << debug_print(r.value()) << std::endl;
}
TEST_F(IntegrateInverseTrigTest, Probe_X2_Over_Sqrt_1_Plus_X2) {
    Symbol x("x");
    ExprPtr e = parse("x^2/sqrt(1+x^2)");
    ASSERT_NE(e, nullptr);
    auto r = integrate(e, x, ctx);
    if (r.is_error()) GTEST_SKIP() << r.error().message;
    std::cout << "x^2/sqrt(1+x^2) → " << debug_print(r.value()) << std::endl;
}
TEST_F(IntegrateInverseTrigTest, Probe_X2_Over_Sqrt_1_Minus_X2) {
    Symbol x("x");
    ExprPtr e = parse("x^2/sqrt(1-x^2)");
    ASSERT_NE(e, nullptr);
    auto r = integrate(e, x, ctx);
    if (r.is_error()) GTEST_SKIP() << r.error().message;
    std::cout << "x^2/sqrt(1-x^2) → " << debug_print(r.value()) << std::endl;
}
TEST_F(IntegrateInverseTrigTest, Probe_Half_X2_Over_1_Plus_X2) {
    Symbol x("x");
    ExprPtr e = parse("(1/2)*x^2/(1+x^2)");
    ASSERT_NE(e, nullptr);
    auto r = integrate(e, x, ctx);
    if (r.is_error()) GTEST_SKIP() << r.error().message;
    std::cout << "(1/2)x^2/(1+x^2) → " << debug_print(r.value()) << std::endl;
}
TEST_F(IntegrateInverseTrigTest, Probe_X2_Times_PowMinus1) {
    Symbol x("x");
    AstArena& a = ctx.arena();
    ExprPtr xsy = a.make<Symbol>("x");
    ExprPtr x2 = a.make<Binary>(BinaryOp::Pow, xsy, a.make<IntegerLit>(BigInt(2)));
    ExprPtr one_plus_x2 = a.make<Sum>(std::vector<ExprPtr>{a.make<IntegerLit>(BigInt(1)), x2});
    ExprPtr inv = a.make<Binary>(BinaryOp::Pow, one_plus_x2, a.make<IntegerLit>(BigInt(-1)));
    ExprPtr prod = a.make<Product>(std::vector<ExprPtr>{a.make<RationalLit>(BigInt(1), BigInt(2)), x2, inv});
    auto r = integrate(prod, x, ctx);
    if (r.is_error()) GTEST_SKIP() << r.error().message;
    std::cout << "Product → " << debug_print(r.value()) << std::endl;
}
TEST_F(IntegrateInverseTrigTest, Probe_SimplifyProduct) {
    AstArena& a = ctx.arena();
    ExprPtr xsy = a.make<Symbol>("x");
    ExprPtr x2 = a.make<Binary>(BinaryOp::Pow, xsy, a.make<IntegerLit>(BigInt(2)));
    ExprPtr one_plus_x2 = a.make<Sum>(std::vector<ExprPtr>{a.make<IntegerLit>(BigInt(1)), x2});
    ExprPtr inv = a.make<Binary>(BinaryOp::Pow, one_plus_x2, a.make<IntegerLit>(BigInt(-1)));
    ExprPtr prod = a.make<Product>(std::vector<ExprPtr>{a.make<RationalLit>(BigInt(1), BigInt(2)), x2, inv});
    auto s = ctx.simplify(prod);
    if (s.is_ok()) std::cout << "SIMP → " << debug_print(s.value()) << std::endl;
    else std::cout << "SIMP FAIL\n";
}
// (DependsOn probe removed — internal namespace not exposed.)
TEST_F(IntegrateInverseTrigTest, Probe_VduAsinhShape) {
    AstArena& a = ctx.arena();
    Symbol x("x");
    ExprPtr xsy = a.make<Symbol>("x");
    ExprPtr x2 = a.make<Binary>(BinaryOp::Pow, xsy, a.make<IntegerLit>(BigInt(2)));
    ExprPtr one_plus_x2 = a.make<Sum>(std::vector<ExprPtr>{a.make<IntegerLit>(BigInt(1)), x2});
    ExprPtr sq = a.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{one_plus_x2});
    ExprPtr inv_sq = a.make<Binary>(BinaryOp::Pow, sq, a.make<IntegerLit>(BigInt(-1)));
    ExprPtr prod = a.make<Product>(std::vector<ExprPtr>{
        a.make<RationalLit>(BigInt(1), BigInt(2)), x2, inv_sq});
    auto r = integrate(prod, x, ctx);
    if (r.is_error()) GTEST_SKIP() << r.error().message;
    std::cout << "Product → " << debug_print(r.value()) << std::endl;
}
TEST_F(IntegrateInverseTrigTest, Probe_XExpMinusXsq) {
    Symbol x("x");
    ExprPtr e = parse("x * exp(-x^2)");
    ASSERT_NE(e, nullptr);
    auto r = integrate(e, x, ctx);
    if (r.is_error()) GTEST_SKIP() << r.error().message;
    std::cout << "x*exp(-x^2) → " << debug_print(r.value()) << std::endl;
    // Verify
    auto Fp = diff(r.value(), x, 1U, ctx);
    ASSERT_TRUE(Fp.is_ok());
    AstArena& a = ctx.arena();
    ExprPtr resid = a.make<Binary>(BinaryOp::Sub, Fp.value(), e);
    auto rs = ctx.simplify(resid);
    std::cout << "residual: " << debug_print(rs.value()) << std::endl;
}
TEST_F(IntegrateInverseTrigTest, Weierstrass_InnerProbe) {
    Symbol t("t");
    // ∫ 2/(t^2+2*t+1) dt should be -2/(t+1) + C. NO log.
    ExprPtr e = parse("2/(t^2+2*t+1)");
    ASSERT_NE(e, nullptr);
    auto r = integrate(e, t, ctx);
    if (r.is_error()) GTEST_SKIP() << r.error().message;
    std::cout << "2/(t+1)^2 → " << debug_print(r.value()) << std::endl;
    // Also probe (1+t)^2 form
    ExprPtr e2 = parse("2/(1+t)^2");
    auto r2 = integrate(e2, t, ctx);
    if (r2.is_ok()) std::cout << "2/(1+t)^2 → " << debug_print(r2.value()) << std::endl;
    // The trigger:
    ExprPtr e3 = parse("1/(1+sin(x))");
    Symbol x("x");
    auto r3 = integrate(e3, x, ctx);
    if (r3.is_ok()) std::cout << "1/(1+sin x) → " << debug_print(r3.value()) << std::endl;
}
namespace cas::symbolic {
extern bool weierstrass_zero_diff(ExprPtr, CASContext&);
}
TEST_F(IntegrateInverseTrigTest, SinCos2_Probe) {
    Symbol x("x");
    ExprPtr e = parse("sin(x) * cos(x)^2");
    auto r = integrate(e, x, ctx);
    if (r.is_error()) { std::cout << "err: " << r.error().message << std::endl; return; }
    auto rs = ctx.simplify(r.value());
    std::cout << "integrate result: " << debug_print(rs.value()) << std::endl;
}
TEST_F(IntegrateInverseTrigTest, SqrtX2M1_Probe) {
    Symbol x("x");
    ExprPtr cas_anti = parse("1/2 * (x * sqrt(x^2 - 1) - acosh(x))");
    ExprPtr max_anti = parse("1/2 * x * sqrt(x^2 - 1) - 1/2 * log(2*x + 2*sqrt(x^2 - 1))");
    auto d_cas = diff(cas_anti, x, 1U, ctx);
    auto d_max = diff(max_anti, x, 1U, ctx);
    auto cs = ctx.simplify(d_cas.value());
    auto ms = ctx.simplify(d_max.value());
    std::cout << "d_cas: " << debug_print(cs.value()) << std::endl;
    std::cout << "d_max: " << debug_print(ms.value()) << std::endl;
    auto eq = cas::symbolic::mathematically_equal(cs.value(), ms.value(), ctx);
    std::cout << "eq: " << (eq.is_ok() ? (eq.value() ? "TRUE" : "FALSE") : "ERR") << std::endl;
}
TEST_F(IntegrateInverseTrigTest, LRT_OneOverX2P1_Probe) {
    Symbol x("x");
    // Direct LRT on simple 1/(x^2+1) — should give arctan(x).
    ExprPtr P = parse("1");
    ExprPtr Q = parse("x^2 + 1");
    auto r = cas::algebra::integrate_rational_lrt(P, Q, x, ctx);
    if (r.is_error()) { std::cout << "err: " << r.error().message << std::endl; return; }
    auto rs = ctx.simplify(r.value());
    std::cout << "LRT result: " << debug_print(rs.value()) << std::endl;
    // Test 1/(2(x^2+1))
    ExprPtr P2 = parse("1");
    ExprPtr Q2 = parse("2*x^2 + 2");
    auto r2 = cas::algebra::integrate_rational_lrt(P2, Q2, x, ctx);
    if (r2.is_ok()) std::cout << "LRT P=1, Q=2x²+2: " << debug_print(r2.value()) << std::endl;
    // Test P = 1/2, Q = x²+1 (the actual call site)
    ExprPtr P3 = parse("1/2");
    ExprPtr Q3 = parse("x^2 + 1");
    auto r3 = cas::algebra::integrate_rational_lrt(P3, Q3, x, ctx);
    if (r3.is_ok()) std::cout << "LRT P=1/2, Q=x²+1: " << debug_print(r3.value()) << std::endl;
}
TEST_F(IntegrateInverseTrigTest, PFD_X2_X2P1Sq_Probe) {
    Symbol x("x");
    ExprPtr e = parse("x^2/(x^2+1)^2");
    auto terms = cas::algebra::partial_fractions(e, x, ctx);
    if (terms.is_error()) { std::cout << "PFD err: " << terms.error().message << std::endl; return; }
    std::cout << "PFD terms: " << terms.value().size() << std::endl;
    for (auto t : terms.value())
        std::cout << "  term: " << debug_print(t) << std::endl;
}
TEST_F(IntegrateInverseTrigTest, X2OverX2P1Sq_Probe) {
    Symbol x("x");
    ExprPtr e = parse("x^2/(x^2+1)^2");
    auto r = integrate(e, x, ctx);
    if (r.is_error()) { std::cout << "err: " << r.error().message << std::endl; GTEST_SKIP(); }
    auto rs = ctx.simplify(r.value());
    std::cout << "result: " << debug_print(rs.value()) << std::endl;
    // Verify derivative.
    auto dr = diff(rs.value(), x, 1U, ctx);
    if (dr.is_ok()) {
        AstArena& a = ctx.arena();
        ExprPtr resid = a.make<Binary>(BinaryOp::Sub, dr.value(), e);
        auto rsd = ctx.simplify(resid);
        std::cout << "residual: " << debug_print(rsd.value()) << std::endl;
    }
}
TEST_F(IntegrateInverseTrigTest, SqrtX2P1Probe) {
    Symbol x("x");
    ExprPtr cas_anti = parse("1/2 * (x * sqrt(x^2+1) + ln(abs(x+sqrt(x^2+1))))");
    ExprPtr max_anti = parse("1/2 * x * sqrt(x^2+1) + 1/2 * asinh(x)");
    auto eq = cas::symbolic::mathematically_equal(cas_anti, max_anti, ctx);
    std::cout << "direct eq: " << (eq.is_ok() ? (eq.value() ? "TRUE" : "FALSE") : "ERR") << std::endl;
    auto d_cas = diff(cas_anti, x, 1U, ctx);
    auto d_max = diff(max_anti, x, 1U, ctx);
    auto cs = ctx.simplify(d_cas.value());
    auto ms = ctx.simplify(d_max.value());
    auto eq2 = cas::symbolic::mathematically_equal(cs.value(), ms.value(), ctx);
    std::cout << "deriv eq: " << (eq2.is_ok() ? (eq2.value() ? "TRUE" : "FALSE") : "ERR") << std::endl;
    std::cout << "d_cas: " << debug_print(cs.value()) << std::endl;
    std::cout << "d_max: " << debug_print(ms.value()) << std::endl;
}
TEST_F(IntegrateInverseTrigTest, WeierstrassSimpleProbe) {
    // sin(x) - 2*sin(x/2)*cos(x/2) ≡ 0
    ExprPtr e = parse("sin(x) - 2*sin(x/2)*cos(x/2)");
    bool z = cas::symbolic::weierstrass_zero_diff(e, ctx);
    std::cout << "double_angle: " << (z ? "TRUE" : "FALSE") << std::endl;
    // cos(x) - (1 - 2*sin(x/2)^2) ≡ 0
    ExprPtr e2 = parse("cos(x) - 1 + 2*sin(x/2)^2");
    bool z2 = cas::symbolic::weierstrass_zero_diff(e2, ctx);
    std::cout << "half_angle_cos: " << (z2 ? "TRUE" : "FALSE") << std::endl;
}
TEST_F(IntegrateInverseTrigTest, WeierstrassDirectProbe) {
    Symbol x("x");
    ExprPtr cas_anti = parse("-2/(cos(x/2)^(-1)*sin(x/2)+1)");
    ExprPtr max_anti = parse("-2/(sin(x)/(cos(x)+1)+1)");
    auto d_cas = diff(cas_anti, x, 1U, ctx);
    auto d_max = diff(max_anti, x, 1U, ctx);
    auto cs = ctx.simplify(d_cas.value());
    auto ms = ctx.simplify(d_max.value());
    AstArena& a = ctx.arena();
    ExprPtr d = a.make<Binary>(BinaryOp::Sub, cs.value(), ms.value());
    bool zero = cas::symbolic::weierstrass_zero_diff(d, ctx);
    std::cout << "weierstrass_zero_diff: " << (zero ? "TRUE" : "FALSE") << std::endl;
}
TEST_F(IntegrateInverseTrigTest, AntiderivDiffProbe) {
    Symbol x("x");
    ExprPtr cas_anti = parse("-2/(cos(x/2)^(-1)*sin(x/2)+1)");
    ExprPtr max_anti = parse("-2/(sin(x)/(cos(x)+1)+1)");
    auto d_cas = diff(cas_anti, x, 1U, ctx);
    auto d_max = diff(max_anti, x, 1U, ctx);
    auto cs = ctx.simplify(d_cas.value());
    auto ms = ctx.simplify(d_max.value());
    std::cout << "d_cas: " << debug_print(cs.value()) << std::endl;
    std::cout << "d_max: " << debug_print(ms.value()) << std::endl;
    AstArena& a = ctx.arena();
    ExprPtr d = a.make<Binary>(BinaryOp::Sub, cs.value(), ms.value());
    auto ds = ctx.simplify(d);
    std::cout << "diff: " << debug_print(ds.value()) << std::endl;
}
TEST_F(IntegrateInverseTrigTest, TrigDiffProbe) {
    Symbol x("x");
    ExprPtr d_cas = parse("1/(1 + sin(x))");
    ExprPtr d_max = parse("2*(cos(x)+1)/(sin(x)+cos(x)+1)^2");
    AstArena& a = ctx.arena();
    ExprPtr diff = a.make<Binary>(BinaryOp::Sub, d_cas, d_max);
    auto s = ctx.simplify(diff);
    if (s.is_ok()) std::cout << "DIFF: " << debug_print(s.value()) << std::endl;
}
TEST_F(IntegrateInverseTrigTest, NegPow_Probe) {
    ExprPtr e1 = parse("-x^2");
    std::cout << "RAW -x^2: " << debug_print(e1) << std::endl;
    auto s1 = ctx.simplify(e1);
    std::cout << "SIMP -x^2: " << debug_print(s1.value()) << std::endl;
    ExprPtr e2 = parse("exp(-x^2)");
    std::cout << "RAW exp(-x^2): " << debug_print(e2) << std::endl;
    auto s2 = ctx.simplify(e2);
    std::cout << "SIMP exp(-x^2): " << debug_print(s2.value()) << std::endl;
    ExprPtr e3 = parse("(-x)^2");
    std::cout << "RAW (-x)^2: " << debug_print(e3) << std::endl;
    auto s3 = ctx.simplify(e3);
    std::cout << "SIMP (-x)^2: " << debug_print(s3.value()) << std::endl;
}
