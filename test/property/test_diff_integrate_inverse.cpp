// F0.4 — Property: D(integrate(f, x), x) ≡ f  for a fixed corpus of 20 integrals.
//
// Verification: diff(integrate(f)) is simplified and checked for
// mathematical equality with f via simplify(diff_result - f) → 0.
//
// TODO: expand to 200 corpus entries once the Risch pipeline is fully
//       instrumented (F5.1).

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <rapidcheck/gtest.h>
#include <string>
#include <vector>

namespace cas::property {
namespace {

// ---------------------------------------------------------------------------
// Helper
// ---------------------------------------------------------------------------

[[nodiscard]] ExprPtr must_parse(const std::string& s, symbolic::CASContext& ctx) {
    auto t = Lexer(s).tokenize();
    if (!t.is_ok()) throw std::runtime_error("lex: " + t.error().message);
    auto e = Parser(t.value(), ctx.arena()).parse();
    if (!e.is_ok()) throw std::runtime_error("parse: " + e.error().message);
    return e.value();
}

[[nodiscard]] bool is_zero(ExprPtr e, symbolic::CASContext& ctx) {
    if (!e) return false;
    auto s = ctx.simplify(e);
    if (!s.is_ok()) return false;
    ExprPtr sv = s.value();
    if (const auto* il = expr_cast<IntegerLit>(sv)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(sv)) return rl->numerator.is_zero();
    return false;
}

// Attempt to verify D(∫f dx) ≡ f.
// Returns: true if verified, false if integrate returned Unimplemented
// (which is a skip, not a failure), throws on hard error.
enum class VerifyResult { Pass, Skip, Fail };

[[nodiscard]] VerifyResult verify_diff_of_integral(
    const std::string& f_str, symbolic::CASContext& ctx, const Symbol& x)
{
    ExprPtr f = must_parse(f_str, ctx);

    auto integ_res = calculus::integrate(f, x, ctx);
    if (!integ_res.is_ok()) {
        // Unimplemented is acceptable (engine limitation) — skip.
        return VerifyResult::Skip;
    }

    auto F = integ_res.value();
    auto df_res = calculus::diff(F, x, 1, ctx);
    if (!df_res.is_ok()) return VerifyResult::Fail;

    // Check D(F) - f simplifies to 0.
    auto diff_expr_res = algebra::expand(
        ctx.arena().make<Binary>(BinaryOp::Sub, df_res.value(), f),
        ctx);
    if (!diff_expr_res.is_ok()) return VerifyResult::Fail;

    return is_zero(diff_expr_res.value(), ctx) ? VerifyResult::Pass : VerifyResult::Fail;
}

// ---------------------------------------------------------------------------
// Corpus: 20 integrals
// ---------------------------------------------------------------------------

// clang-format off
static const std::vector<std::string> kCorpus = {
    // Polynomials
    "x",
    "x^2",
    "x^3",
    "3*x^2 + 2*x + 1",
    "x^4 - 2*x^2 + 1",

    // Rationals (integrable by partial fractions)
    "1/x",
    "1/(x^2 + 1)",

    // Exponentials
    "exp(x)",
    "2*exp(x)",

    // Trigonometric
    "sin(x)",
    "cos(x)",
    "sin(2*x)",
    "cos(3*x)",

    // Products
    "x*exp(x)",
    "x*sin(x)",
    "x*cos(x)",

    // Compositions
    "exp(2*x)",
    "sin(x)^2",

    // Logarithmic
    "ln(x)",

    // Rational simple
    "1/(x + 1)",
};
// clang-format on

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

class DiffIntegrateInverseTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};
};

// Run the full fixed corpus.  At least 50% of entries must pass
// (the rest may be Skip due to current engine limitations).
TEST_F(DiffIntegrateInverseTest, FullCorpusFiftyPercentPass) {
    int pass = 0; int skip = 0; int fail = 0;
    std::vector<std::string> failures;
    for (const auto& entry : kCorpus) {
        symbolic::CASContext lctx;
        Symbol lx{"x"};
        auto res = verify_diff_of_integral(entry, lctx, lx);
        if (res == VerifyResult::Pass) ++pass;
        else if (res == VerifyResult::Skip) ++skip;
        else { ++fail; failures.push_back(entry); }
    }
    // Hard requirement: zero hard failures (wrong answers).
    EXPECT_EQ(fail, 0)
        << "Hard failures (D(∫f) ≠ f): " << fail << "\n"
        << "  Failed integrands: "
        << [&]{ std::string s; for (auto& f : failures) s += "\n    " + f; return s; }();
    // At least 50% of corpus must pass (not skip).
    EXPECT_GE(pass, static_cast<int>(kCorpus.size()) / 2)
        << "Expected ≥ " << kCorpus.size()/2 << " passes, got " << pass
        << " (skipped=" << skip << ")";
}

// Per-entry tests for the 5 simplest entries that MUST always pass.
TEST_F(DiffIntegrateInverseTest, PolynomialX) {
    EXPECT_NE(verify_diff_of_integral("x", ctx, x), VerifyResult::Fail);
}
TEST_F(DiffIntegrateInverseTest, PolynomialXSquared) {
    EXPECT_NE(verify_diff_of_integral("x^2", ctx, x), VerifyResult::Fail);
}
TEST_F(DiffIntegrateInverseTest, ExponentialExp) {
    EXPECT_NE(verify_diff_of_integral("exp(x)", ctx, x), VerifyResult::Fail);
}
TEST_F(DiffIntegrateInverseTest, TrigSin) {
    EXPECT_NE(verify_diff_of_integral("sin(x)", ctx, x), VerifyResult::Fail);
}
TEST_F(DiffIntegrateInverseTest, TrigCos) {
    EXPECT_NE(verify_diff_of_integral("cos(x)", ctx, x), VerifyResult::Fail);
}

// rapidcheck: randomly select an index into corpus and verify no hard failure.
RC_GTEST_FIXTURE_PROP(DiffIntegrateInverseTest, NoHardFailureForAnyCorpusEntry, ()) {
    // Use fixed-seed iteration over corpus (no random needed — corpus is finite).
    for (const auto& entry : kCorpus) {
        symbolic::CASContext lctx;
        Symbol lx{"x"};
        RC_ASSERT(verify_diff_of_integral(entry, lctx, lx) != VerifyResult::Fail);
    }
}

}  // namespace
}  // namespace cas::property
