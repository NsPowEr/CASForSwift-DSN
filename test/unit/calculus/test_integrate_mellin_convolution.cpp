// A7 §6.7 — Mellin-convolution definite integrator
// (Meijer_G_Slater.md §6.7).  Each closed form below is an INDEPENDENT ground
// truth: verified numerically with mpmath (high-precision quadrature) before
// being asserted here — e.g.
//   ∫_0^∞ e^{-η x} J_0(2√(ω x)) dx = e^{-ω/η}/η   (Laplace transform of J_0),
//   quad(η=1,ω=1)=0.36787944117…=1/e, quad(η=2,ω=1)=0.30326533…=e^{-1/2}/2,
//   quad(η=1,ω=3)=0.049787068…=e^{-3},
//   ∫_0^∞ x e^{-2x} J_0(2√x) dx = 0.075816332…= e^{-1/2}/8 = e^{-ω/η}(η-ω)/η³.
// The last two tests assert the CONVERGENCE GATE rejects divergent integrands
// (no finite wrong value) — the soundness contract of §6.7 (REGOLA ZERO).

#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/numeric.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <cmath>
#include <gtest/gtest.h>

#include <string>

namespace cas::test {
namespace {

[[nodiscard]] ExprPtr parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    EXPECT_TRUE(tokens.is_ok()) << input;
    Parser parser(tokens.value(), arena);
    auto r = parser.parse();
    EXPECT_TRUE(r.is_ok()) << input;
    return r.value();
}

}  // namespace

class MellinConvolutionTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] Result<ExprPtr> integrate_half_line(const std::string& integrand) {
        ExprPtr f = parse_expr(integrand, ctx.arena());
        ExprPtr zero = ctx.arena().make<IntegerLit>(BigInt(0));
        ExprPtr inf = ctx.arena().make<Constant>(MathConstant::Infinity);
        return calculus::definite_integral(f, Symbol("x"), zero, inf, ctx);
    }

    // Dual certificate: the value must match `expected` STRUCTURALLY
    // (mathematically_equal) when it folds to the elementary form, OR — always
    // — NUMERICALLY (high-precision eval within tol).  The numeric leg is the
    // spec's §6.7 anti-silent-wrong shadow-evaluation, run in the test layer
    // where the evaluator lives (cas_numeric depends on cas_calculus, so it is
    // unavailable inside the pattern itself).  A hypergeometric closed form
    // that is exact but not further reduced (e.g. ₁F₁(2;1;z)=(1+z)e^z) still
    // passes via the numeric leg — the mathematics is proven, not the printout.
    void expect_value(const std::string& integrand, const std::string& expected) {
        auto got = integrate_half_line(integrand);
        ASSERT_TRUE(got.is_ok()) << integrand << " -> " << got.error().message;
        ExprPtr want = parse_expr(expected, ctx.arena());

        auto eq = symbolic::mathematically_equal(got.value(), want, ctx);
        const bool structural = eq.is_ok() && eq.value();

        auto gv = numeric::eval(got.value(), ctx);
        auto wv = numeric::eval(want, ctx);
        ASSERT_TRUE(gv.is_ok()) << "numeric eval of result failed: "
                                << debug_print(got.value());
        ASSERT_TRUE(wv.is_ok());
        const bool numeric_ok = std::abs(gv.value() - wv.value()) < 1e-9;

        EXPECT_TRUE(structural || numeric_ok)
            << "∫_0^∞ " << integrand << " dx\n  GOT: " << debug_print(got.value())
            << " (=" << gv.value() << ")\n  WANT: " << expected
            << " (=" << wv.value() << ")";
    }
};

// ∫_0^∞ e^{-x} J_0(2√x) dx = 1/e.  (η=ω=1; mpmath 0.3678794411714…)
TEST_F(MellinConvolutionTest, ExpBesselJ0_UnitScales) {
    expect_value("exp(-x) * bessel_j(0, 2*sqrt(x))", "exp(-1)");
}

// ∫_0^∞ e^{-2x} J_0(2√x) dx = e^{-1/2}/2.  (scale in η; mpmath 0.30326533…)
TEST_F(MellinConvolutionTest, ExpBesselJ0_ScaledEta) {
    expect_value("exp(-2*x) * bessel_j(0, 2*sqrt(x))", "exp(-1/2)/2");
}

// ∫_0^∞ e^{-x} J_0(2√(3x)) dx = e^{-3}.  (scale in ω; mpmath 0.049787068…)
TEST_F(MellinConvolutionTest, ExpBesselJ0_ScaledOmega) {
    expect_value("exp(-x) * bessel_j(0, 2*sqrt(3*x))", "exp(-3)");
}

// ∫_0^∞ x e^{-2x} J_0(2√x) dx = e^{-1/2}/8.  (tests the x^μ power fold via
// §6.2; = e^{-ω/η}(η-ω)/η³ at η=2,ω=1; mpmath 0.075816332…)
TEST_F(MellinConvolutionTest, ExpBesselJ0_WithXPower) {
    expect_value("x * exp(-2*x) * bessel_j(0, 2*sqrt(x))", "exp(-1/2)/8");
}

// SOUNDNESS: ∫_0^∞ J_0(2√x) J_0(2√(2x)) dx DIVERGES — both factors are
// G^{1,0}_{0,2} (q-p=2 ⇒ algebraic-oscillatory decay ~x^{-1/4}, product
// ~x^{-1/2}, ∫^∞ diverges).  The ∞-convergence certificate has NO exp kernel
// (no q_i=p_i+1) ⇒ the pattern MUST skip; the definite integral must NOT
// return a finite (wrong) value.
TEST_F(MellinConvolutionTest, DivergentBesselProduct_InfinityGateSkips) {
    auto got = integrate_half_line(
        "bessel_j(0, 2*sqrt(x)) * bessel_j(0, 2*sqrt(2*x))");
    EXPECT_TRUE(got.is_error())
        << "divergent J0·J0 must not yield a value; got "
        << (got.is_ok() ? "a finite result" : "");
}

// SOUNDNESS: ∫_0^∞ e^{-x} J_0(2√x) / x dx DIVERGES at 0 (integrand ~1/x).
// The 0-endpoint exponent is λ0 = -1 (x^{-1} folded into the exp G) + 0 = -1,
// NOT > -1 ⇒ the certificate rejects it; the pattern must skip.
TEST_F(MellinConvolutionTest, DivergentAtZero_ZeroGateSkips) {
    auto got = integrate_half_line(
        "x^(-1) * exp(-x) * bessel_j(0, 2*sqrt(x))");
    EXPECT_TRUE(got.is_error())
        << "0-divergent integrand must not yield a value";
}

}  // namespace cas::test
