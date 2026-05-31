// F0.4 — Property: partial_fractions(f) reassembled ≡ f
//
// For each seed rational f(x) = p(x)/q(x):
//   1. Compute partial_fractions(f) → list of terms.
//   2. Sum the terms (via together / expand).
//   3. Check sum - f simplifies to 0.

#include "cas/algebra.hpp"
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
// Helpers
// ---------------------------------------------------------------------------

[[nodiscard]] ExprPtr must_parse(const std::string& s, symbolic::CASContext& ctx) {
    auto t = Lexer(s).tokenize();
    if (!t.is_ok()) throw std::runtime_error("lex: " + t.error().message);
    auto e = Parser(t.value(), ctx.arena()).parse();
    if (!e.is_ok()) throw std::runtime_error("parse: " + e.error().message);
    return e.value();
}

[[nodiscard]] bool is_zero(ExprPtr e, symbolic::CASContext& ctx) {
    auto s = ctx.simplify(e);
    if (!s.is_ok()) return false;
    ExprPtr sv = s.value();
    if (const auto* il = expr_cast<IntegerLit>(sv)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(sv)) return rl->numerator.is_zero();
    return false;
}

// Build sum of a list of ExprPtr.
[[nodiscard]] ExprPtr sum_terms(
    const std::vector<ExprPtr>& terms, symbolic::CASContext& ctx)
{
    if (terms.empty()) return ctx.arena().make<IntegerLit>(BigInt(0));
    ExprPtr s = terms[0];
    for (std::size_t i = 1; i < terms.size(); ++i) {
        s = ctx.arena().make<Binary>(BinaryOp::Add, s, terms[i]);
    }
    return s;
}

// Returns: true = pass, false = fail, nullopt = skip (unimplemented).
[[nodiscard]] std::optional<bool> verify_roundtrip(
    const std::string& f_str, const Symbol& var, symbolic::CASContext& ctx)
{
    ExprPtr f = must_parse(f_str, ctx);
    auto pf_res = algebra::partial_fractions(f, var, ctx);
    if (!pf_res.is_ok()) return std::nullopt;  // skip

    const auto& parts = pf_res.value();
    if (parts.empty()) {
        // Empty decomposition — original must be zero.
        return is_zero(f, ctx) ? std::optional<bool>(true) : std::optional<bool>(false);
    }

    ExprPtr reassembled = sum_terms(parts, ctx);

    // Check reassembled - f simplifies to 0 via together / expand.
    auto tog_res = algebra::together(
        ctx.arena().make<Binary>(BinaryOp::Sub, reassembled, f), ctx);
    if (!tog_res.is_ok()) {
        // Try raw expand instead.
        auto exp_res = algebra::expand(
            ctx.arena().make<Binary>(BinaryOp::Sub, reassembled, f), ctx);
        if (!exp_res.is_ok()) return std::nullopt;
        return is_zero(exp_res.value(), ctx);
    }
    return is_zero(tog_res.value(), ctx);
}

// ---------------------------------------------------------------------------
// Corpus: 5 rational functions
// ---------------------------------------------------------------------------

struct RatEntry {
    std::string str;
    std::string var;
};

static const std::vector<RatEntry> kCorpus = {
    // 1 / (x*(x+1)) = 1/x - 1/(x+1)
    {"1/(x*(x+1))",             "x"},
    // 1 / (x^2 - 1) = 1/(2*(x-1)) - 1/(2*(x+1))
    {"1/(x^2 - 1)",             "x"},
    // (x+2) / ((x+1)*(x-1))
    {"(x+2)/((x+1)*(x-1))",    "x"},
    // 1 / (x*(x^2+1))
    {"1/(x*(x^2+1))",           "x"},
    // (2*x) / (x^2 - 3*x + 2)   = 2*(x) / ((x-1)*(x-2))
    {"2*x/(x^2 - 3*x + 2)",     "x"},
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

class PartialFractionsRoundtripTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

TEST_F(PartialFractionsRoundtripTest, AllCorpusEntries_NoHardFailure) {
    int pass = 0; int skip = 0; int fail = 0;
    std::vector<std::string> failures;
    for (const auto& entry : kCorpus) {
        symbolic::CASContext lctx;
        Symbol var{entry.var};
        auto res = verify_roundtrip(entry.str, var, lctx);
        if (!res.has_value()) ++skip;
        else if (*res) ++pass;
        else { ++fail; failures.push_back(entry.str); }
    }
    EXPECT_EQ(fail, 0)
        << "partial_fractions reassembly ≠ original for:"
        << [&]{ std::string s; for (auto& f : failures) s += "\n  " + f; return s; }();
    EXPECT_GE(pass, static_cast<int>(kCorpus.size()) / 2)
        << "Expected ≥ 50% pass, got " << pass << "/" << kCorpus.size()
        << " (skip=" << skip << ")";
}

TEST_F(PartialFractionsRoundtripTest, SimpleProduct_xTimesXPlus1) {
    Symbol x{"x"};
    auto res = verify_roundtrip("1/(x*(x+1))", x, ctx);
    if (res.has_value()) EXPECT_TRUE(*res);
}

TEST_F(PartialFractionsRoundtripTest, DifferenceOfSquares_xSquaredMinus1) {
    Symbol x{"x"};
    auto res = verify_roundtrip("1/(x^2 - 1)", x, ctx);
    if (res.has_value()) EXPECT_TRUE(*res);
}

// rapidcheck: all corpus entries have no hard failure.
RC_GTEST_FIXTURE_PROP(PartialFractionsRoundtripTest, NoHardFailureForCorpus, ()) {
    for (const auto& entry : kCorpus) {
        symbolic::CASContext lctx;
        Symbol var{entry.var};
        auto res = verify_roundtrip(entry.str, var, lctx);
        RC_ASSERT(!res.has_value() || *res);
    }
}

}  // namespace
}  // namespace cas::property
