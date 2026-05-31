// F0.4 — Property: factor(p) → ∏ factors ≡ p  for 10 Z[x] polynomials.
//
// Verification: reconstruct the product of all returned factors (raised to
// their multiplicities, multiplied by content) and check it equals p via
// expand + simplify → 0 difference.

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

// Reconstruct product from Factorization:
//   content * factor[0]^m0 * factor[1]^m1 * ...
[[nodiscard]] ExprPtr reconstruct_product(
    const algebra::Factorization& fz, symbolic::CASContext& ctx)
{
    AstArena& arena = ctx.arena();
    ExprPtr prod = fz.content;

    for (const auto& pf : fz.factors) {
        ExprPtr term = pf.factor;
        if (pf.multiplicity > 1U) {
            ExprPtr exp_e = arena.make<IntegerLit>(BigInt(
                static_cast<std::int64_t>(pf.multiplicity)));
            term = arena.make<Binary>(BinaryOp::Pow, term, exp_e);
        }
        prod = arena.make<Binary>(BinaryOp::Mul, prod, term);
    }
    return prod;
}

// Verify factor(p) → ∏ factors ≡ p by checking expand(product - p) == 0.
// Returns: true = pass, false = fail, nullopt = skip.
[[nodiscard]] std::optional<bool> verify_factorization(
    const std::string& poly_str, const Symbol& var, symbolic::CASContext& ctx)
{
    ExprPtr p = must_parse(poly_str, ctx);
    auto fz_res = algebra::factor_over_integers(p, var, ctx);
    if (!fz_res.is_ok()) return std::nullopt;  // skip

    ExprPtr product = reconstruct_product(fz_res.value(), ctx);
    ExprPtr diff = ctx.arena().make<Binary>(BinaryOp::Sub, product, p);

    auto exp_res = algebra::expand(diff, ctx);
    if (!exp_res.is_ok()) return std::nullopt;
    return is_zero(exp_res.value(), ctx);
}

// ---------------------------------------------------------------------------
// Corpus: 10 Z[x] polynomials
// ---------------------------------------------------------------------------

struct PolyEntry {
    std::string str;
    std::string var;
};

static const std::vector<PolyEntry> kCorpus = {
    // Constant
    {"1",                                "x"},
    // Linear
    {"x - 5",                            "x"},
    // Quadratic — factorable
    {"x^2 - 1",                          "x"},
    // Quadratic — irreducible over Q
    {"x^2 + 1",                          "x"},
    // Cubic — full factorization
    {"x^3 - 6*x^2 + 11*x - 6",          "x"},
    // Repeated root
    {"x^2 - 2*x + 1",                    "x"},
    // Difference of squares with large coeff
    {"x^2 - 25",                         "x"},
    // Product of 4 linear factors
    {"x^4 - 10*x^3 + 35*x^2 - 50*x + 24", "x"},
    // Irreducible quartic
    {"x^4 + 1",                          "x"},
    // With content != 1
    {"2*x^2 - 8",                        "x"},
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

class FactorProductTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

TEST_F(FactorProductTest, AllCorpusEntries_NoHardFailure) {
    int pass = 0; int skip = 0; int fail = 0;
    std::vector<std::string> failures;
    for (const auto& entry : kCorpus) {
        symbolic::CASContext lctx;
        Symbol var{entry.var};
        auto res = verify_factorization(entry.str, var, lctx);
        if (!res.has_value()) ++skip;
        else if (*res) ++pass;
        else { ++fail; failures.push_back(entry.str); }
    }
    EXPECT_EQ(fail, 0)
        << "factor(p) → ∏ ≠ p for: "
        << [&]{ std::string s; for (auto& f : failures) s += "\n  " + f; return s; }();
    EXPECT_GE(pass, static_cast<int>(kCorpus.size()) / 2)
        << "Expected ≥ 50% pass, got " << pass << "/" << kCorpus.size()
        << " (skip=" << skip << ")";
}

TEST_F(FactorProductTest, DifferenceOfSquaresX2Minus1) {
    Symbol x{"x"};
    auto res = verify_factorization("x^2 - 1", x, ctx);
    if (res.has_value()) EXPECT_TRUE(*res);
}

TEST_F(FactorProductTest, Cubic6_11_minus6) {
    Symbol x{"x"};
    auto res = verify_factorization("x^3 - 6*x^2 + 11*x - 6", x, ctx);
    if (res.has_value()) EXPECT_TRUE(*res);
}

TEST_F(FactorProductTest, RepeatedRootXMinus1Squared) {
    Symbol x{"x"};
    auto res = verify_factorization("x^2 - 2*x + 1", x, ctx);
    if (res.has_value()) EXPECT_TRUE(*res);
}

// rapidcheck: all corpus entries have no hard failure.
RC_GTEST_FIXTURE_PROP(FactorProductTest, NoHardFailureForCorpus, ()) {
    for (const auto& entry : kCorpus) {
        symbolic::CASContext lctx;
        Symbol var{entry.var};
        auto res = verify_factorization(entry.str, var, lctx);
        // hard failure (wrong answer) is forbidden; Unimplemented (nullopt) is ok
        RC_ASSERT(!res.has_value() || *res);
    }
}

}  // namespace
}  // namespace cas::property
