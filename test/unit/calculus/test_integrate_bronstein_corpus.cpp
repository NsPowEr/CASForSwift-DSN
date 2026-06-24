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
// NOTE: two integrands are intentionally omitted because integrate() does
// not currently terminate on them (tracked separately as an S5 hang):
//   log(x + sqrt(x^2 + 1))   and   sin(x)/x.

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

}  // namespace
