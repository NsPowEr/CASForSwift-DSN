// Regression gate for the Risch integration pipeline on the Bronstein
// ch.6-8 transcendental corpus (T-007).  Two guarantees:
//
//   1. NO SILENT WRONG ANSWERS (REGOLA ZERO, "mai output sbagliato").  For
//      every integrand, if integrate() returns a closed form F, then F must
//      satisfy F'(x) = f(x) at every sampled point.  A mismatch fails the
//      test.  This pins the fix for ∫ e^{±x}·ln(x) dx, which previously
//      returned the wrong elementary form e^{∓x}·(x·ln x − x) (both integrals
//      are in fact non-elementary, Ei-log) — now they return Unimplemented.
//
//   2. COVERAGE FLOOR.  At least 28 of the listed integrals must still be
//      solved and verified, so a future change cannot quietly trade
//      correctness for lost coverage (the elementary set currently solves 30).
//
// Non-elementary kernels (Ei/erf/Si/dilog) are expected to return
// Unimplemented; that is counted as neither wrong nor covered.
//
// NOTE: two integrands — log(x + sqrt(x^2 + 1)) and sin(x)/x — previously
// drove integrate() into a non-terminating recursion (S5 hang).  Both are now
// pinned by dedicated regressions below: sin(x)/x returns Unimplemented
// (non-elementary Si) via the Weierstrass guard; log(x+√(x²+1)) is now SOLVED
// in closed form via conjugate rationalisation (HC-IBP-RADSUM-RATIONALIZE +
// the parse_integer_exponent robustness fix).

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

#include "cas/calculus.hpp"
#include "cas/formatter.hpp"
#include "cas/lexer.hpp"
#include "cas/numeric.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

struct Entry { const char* input; const char* ref; };

const std::vector<Entry> kCorpus = {
    {"log(x)", "ch6_ex1"},
    {"log(x)^2", "ch6_ex2"},
    {"x*log(x)", "ch6_ex3"},
    {"log(x^2 + 1)", "ch6_ex4"},
    {"log(x + 1)/(x + 1)", "ch6_ex5"},
    {"1/(x*log(x))", "ch6_ex6"},
    {"log(x)/x^2", "ch6_ex7"},
    {"x^2*log(x)", "ch6_ex8"},
    {"(log(x))^2/x", "ch6_ex10"},
    {"exp(x)/x", "ch7_ex1_ei"},
    {"exp(x)*x", "ch7_ex2"},
    {"exp(x^2)", "ch7_ex3_erf"},
    {"exp(x)*sin(x)", "ch7_ex5"},
    {"exp(x)*cos(x)", "ch7_ex6"},
    {"exp(-x)*log(x)", "ch7_ex8_eilog"},     // regression: was silently wrong
    {"exp(2*x)*sin(3*x)", "ch7_ex9"},
    {"x*exp(-x^2)", "ch7_ex10"},
    {"sin(x)^3", "ch8_ex1"},
    {"cos(x)^3", "ch8_ex2"},
    {"sin(x)^4", "ch8_ex3"},
    {"sin(x)*cos(x)^2", "ch8_ex4"},
    {"sin(x)^2*cos(x)^2", "ch8_ex5"},
    {"x*log(x^2 + 1)", "ch6_ex12"},
    {"log(x)/(x^2 + 1)", "ch6_ex14_dilog"},
    {"x^2*log(x)^2", "ch6_ex15"},
    {"log(x)^3/x", "ch6_ex16"},
    {"log(sqrt(x))", "ch6_ex17"},
    {"log(x + 1)^2", "ch6_ex19"},
    {"(log(x))^2 / (x^2)", "ch6_ex20"},
    {"x^2*exp(2*x)", "ch7_ex11"},
    {"x^3*exp(-x)", "ch7_ex12"},
    {"exp(x)*sin(2*x)", "ch7_ex13"},
    {"exp(-x)*cos(x)", "ch7_ex14"},
    {"x*exp(-x^2/2)", "ch7_ex16"},
    {"exp(2*x)*cos(3*x)", "ch7_ex17"},
    {"exp(x)*log(x)", "ch7_ex18_eilog"},     // regression: was silently wrong
    {"exp(x)/(x^2)", "ch7_ex19_ei"},
};

class BronsteinCorpus : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};

    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        if (!t.is_ok()) return nullptr;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        return r.is_ok() ? r.value() : nullptr;
    }
};

TEST_F(BronsteinCorpus, NoSilentWrongAnswers_And_CoverageFloor) {
    const std::vector<double> pts = {0.37, 0.81, 1.43, 2.17, 3.05};
    int n_verified = 0;
    std::vector<std::string> wrong;
    formatter::TextFormatter fmt;

    for (const auto& e : kCorpus) {
        ExprPtr f = parse(e.input);
        ASSERT_NE(f, nullptr) << e.input;
        auto r = calculus::integrate(f, x, ctx);
        if (r.is_error()) continue;  // Unimplemented: acceptable (non-elementary)

        ExprPtr F = r.value();
        auto dF = calculus::diff(F, x, 1U, ctx);
        if (dF.is_error()) continue;

        bool mismatch = false, eval_ok_any = false;
        for (double xv : pts) {
            numeric::NumericEnv env{{"x", xv}};
            auto fv = numeric::eval(f, env);
            auto Fv = numeric::eval(dF.value(), env);
            if (fv.is_error() || Fv.is_error()) continue;
            eval_ok_any = true;
            double d = std::fabs(fv.value() - Fv.value());
            double scale = 1.0 + std::fabs(fv.value());
            if (d > 1e-6 * scale) { mismatch = true; break; }
        }
        if (mismatch) {
            wrong.push_back(std::string(e.ref) + " : " + e.input
                            + "  ->  " + fmt.format(F));
        } else if (eval_ok_any) {
            ++n_verified;
        }
    }

    EXPECT_TRUE(wrong.empty())
        << "integrate() returned a non-antiderivative (D(F) != f) for:\n  "
        << [&] { std::string s; for (auto& w : wrong) s += w + "\n  "; return s; }();

    EXPECT_GE(n_verified, 28)
        << "Risch coverage regressed below the established floor (28); "
        << n_verified << "/" << kCorpus.size() << " verified Bronstein ch.6-8.";
}

// Shared contract for the two integrands that previously hung integrate() (S5):
//   * integrate() must TERMINATE — the gtest TIMEOUT is the liveness assertion.
//   * NO SILENT WRONG ANSWER — if a closed form F is returned, D(F) must equal
//     the integrand at every sampled point (Unimplemented is acceptable).
// Split one-per-hang so each stays well under the 60s per-test CTest cap and a
// regression points at the exact integrand.
static void expect_no_hang_no_silent_wrong(
    ExprPtr f, const Symbol& x, symbolic::CASContext& ctx, const char* in) {
    ASSERT_NE(f, nullptr) << in;
    auto r = calculus::integrate(f, x, ctx);  // must terminate (no hang)
    if (r.is_error()) return;                 // Unimplemented: acceptable
    auto dF = calculus::diff(r.value(), x, 1U, ctx);
    ASSERT_TRUE(dF.is_ok()) << in;
    for (double xv : {0.37, 0.81, 1.43, 2.17}) {
        numeric::NumericEnv env{{"x", xv}};
        auto fv = numeric::eval(f, env);
        auto Fv = numeric::eval(dF.value(), env);
        if (fv.is_error() || Fv.is_error()) continue;
        const double d = std::fabs(fv.value() - Fv.value());
        const double scale = 1.0 + std::fabs(fv.value());
        EXPECT_LE(d, 1e-6 * scale) << in << " at x=" << xv;
    }
}

// sin(x)/x — non-elementary (Si).  Weierstrass must NOT apply (integrand is not
// rational in sin/cos), so it returns Unimplemented instead of looping.
TEST_F(BronsteinCorpus, NoHang_SinXOverX_NonElementary) {
    expect_no_hang_no_silent_wrong(parse("sin(x)/x"), x, ctx, "sin(x)/x");
}

// log(x + sqrt(x^2+1)) = asinh(x), elementary.  Now SOLVED in closed form
// (x·log(x+√(x²+1)) − √(x²+1)) via conjugate rationalisation of the parts
// remainder x·g'/g (HC-IBP-RADSUM-RATIONALIZE), enabled by the
// parse_integer_exponent robustness fix.  Asserts the closed form is produced
// and is a genuine antiderivative (D(F) = f).
TEST_F(BronsteinCorpus, SolvesLogXPlusSqrt_Asinh) {
    ExprPtr f = parse("log(x + sqrt(x^2 + 1))");
    ASSERT_NE(f, nullptr);
    auto r = calculus::integrate(f, x, ctx);
    ASSERT_TRUE(r.is_ok()) << "asinh log integral must produce a closed form";
    auto dF = calculus::diff(r.value(), x, 1U, ctx);
    ASSERT_TRUE(dF.is_ok());
    bool any = false;
    for (double xv : {0.37, 0.81, 1.43, 2.17}) {
        numeric::NumericEnv env{{"x", xv}};
        auto fv = numeric::eval(f, env);
        auto Fv = numeric::eval(dF.value(), env);
        if (fv.is_error() || Fv.is_error()) continue;
        any = true;
        EXPECT_LE(std::fabs(fv.value() - Fv.value()), 1e-6 * (1.0 + std::fabs(fv.value())))
            << "D(F) != f at x=" << xv;
    }
    EXPECT_TRUE(any) << "numeric verification could not evaluate any sample point";
}

// General radius: ∫log(x+√(x²+a)) for a≠1.  Needs both the conjugate
// rationalisation AND the T-054 numeric-coefficient lift — without it the parts
// remainder reduces to c·x·(c·√(x²+a))⁻¹ (c=(S·conj)²) whose constant does not
// cancel, leaving NO_STRATEGY.  Asserts the closed form is produced and verified.
TEST_F(BronsteinCorpus, SolvesLogXPlusSqrt_GeneralRadius) {
    for (const char* in : {"log(x + sqrt(x^2 + 4))", "log(x + sqrt(x^2 + 9))"}) {
        ExprPtr f = parse(in);
        ASSERT_NE(f, nullptr) << in;
        auto r = calculus::integrate(f, x, ctx);
        ASSERT_TRUE(r.is_ok()) << in << " must produce a closed form";
        auto dF = calculus::diff(r.value(), x, 1U, ctx);
        ASSERT_TRUE(dF.is_ok()) << in;
        bool any = false;
        for (double xv : {0.4, 0.9, 1.6, 2.3}) {
            numeric::NumericEnv env{{"x", xv}};
            auto fv = numeric::eval(f, env);
            auto Fv = numeric::eval(dF.value(), env);
            if (fv.is_error() || Fv.is_error()) continue;
            any = true;
            EXPECT_LE(std::fabs(fv.value() - Fv.value()), 1e-6 * (1.0 + std::fabs(fv.value())))
                << in << " D(F) != f at x=" << xv;
        }
        EXPECT_TRUE(any) << in;
    }
}

}  // namespace
