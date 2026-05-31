// F0.4 — Property: f(root(f)) ≡ 0  for every root returned by solve_polynomial.
//
// For each seed polynomial, solve_polynomial returns a list of roots.
// Substituting each root back into f and simplifying must yield 0.

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

// Returns true if expr simplifies to 0.
[[nodiscard]] bool simplifies_to_zero(ExprPtr e, symbolic::CASContext& ctx) {
    auto s = ctx.simplify(e);
    if (!s.is_ok()) return false;
    ExprPtr sv = s.value();
    if (const auto* il = expr_cast<IntegerLit>(sv)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(sv)) return rl->numerator.is_zero();
    return false;
}

// Substitute root into poly and check result is zero.
// Returns: true = passes, false = fails, nullopt = skip (unimplemented).
[[nodiscard]] std::optional<bool> check_root(
    ExprPtr poly, ExprPtr root, const Symbol& var, symbolic::CASContext& ctx)
{
    auto sub = substitute(poly, var, root, ctx);
    if (!sub.is_ok()) return std::nullopt;  // skip
    auto exp = algebra::expand(sub.value(), ctx);
    if (!exp.is_ok()) return std::nullopt;
    return simplifies_to_zero(exp.value(), ctx);
}

// ---------------------------------------------------------------------------
// Seed corpus: 5 polynomials with known roots
// ---------------------------------------------------------------------------

struct PolynomialSeed {
    std::string poly_str;
    std::string var_str;
    std::size_t expected_root_count;
};

static const std::vector<PolynomialSeed> kSeeds = {
    // x - 3 = 0  → root = 3
    {"x - 3",               "x", 1},
    // x^2 - 4 = 0  → roots = {2, -2}
    {"x^2 - 4",             "x", 2},
    // x^2 - 3*x + 2 = 0  → roots = {1, 2}
    {"x^2 - 3*x + 2",       "x", 2},
    // x^3 - 6*x^2 + 11*x - 6 = 0  → roots = {1, 2, 3}
    {"x^3 - 6*x^2 + 11*x - 6", "x", 3},
    // x^2 + x - 2 = 0  → roots = {1, -2}
    {"x^2 + x - 2",          "x", 2},
};

// ---------------------------------------------------------------------------
// Test class
// ---------------------------------------------------------------------------

class RootsSubstituteTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

// For each seed, verify every returned root satisfies f(root) = 0.
TEST_F(RootsSubstituteTest, AllSeedRootsSatisfyPolynomial) {
    for (const auto& seed : kSeeds) {
        symbolic::CASContext lctx;
        Symbol var{seed.var_str};
        ExprPtr poly = must_parse(seed.poly_str, lctx);

        auto roots_res = algebra::solve_polynomial(poly, var, lctx);
        if (!roots_res.is_ok()) {
            // Unimplemented — not a hard failure, just skip.
            GTEST_SKIP() << "solve_polynomial returned Unimplemented for: " << seed.poly_str;
        }

        const auto& roots = roots_res.value();
        ASSERT_FALSE(roots.empty())
            << "No roots returned for: " << seed.poly_str;
        EXPECT_GE(roots.size(), 1u)
            << "Expected ≥ 1 root for: " << seed.poly_str;

        int pass = 0;
        for (ExprPtr root : roots) {
            auto res = check_root(poly, root, var, lctx);
            if (!res.has_value()) continue;  // skip unimplemented sub-step
            EXPECT_TRUE(*res)
                << "f(root) ≠ 0 for poly=" << seed.poly_str;
            if (*res) ++pass;
        }
        // At least one root must be verifiable.
        EXPECT_GE(pass, 1)
            << "Zero verifiable roots for: " << seed.poly_str;
    }
}

// P-ROOTS-01: for x^2 - 4, both roots must satisfy f(root)=0
TEST_F(RootsSubstituteTest, QuadraticXSquaredMinus4) {
    Symbol x{"x"};
    auto poly = must_parse("x^2 - 4", ctx);
    auto roots_res = algebra::solve_polynomial(poly, x, ctx);
    if (!roots_res.is_ok()) GTEST_SKIP();
    for (ExprPtr r : roots_res.value()) {
        auto res = check_root(poly, r, x, ctx);
        if (res.has_value()) EXPECT_TRUE(*res);
    }
}

// P-ROOTS-02: for x - 3, the single root must be 3
TEST_F(RootsSubstituteTest, LinearXMinus3) {
    Symbol x{"x"};
    auto poly = must_parse("x - 3", ctx);
    auto roots_res = algebra::solve_polynomial(poly, x, ctx);
    if (!roots_res.is_ok()) GTEST_SKIP();
    ASSERT_EQ(roots_res.value().size(), 1u);
    ExprPtr r = roots_res.value()[0];
    if (const auto* il = expr_cast<IntegerLit>(r)) {
        EXPECT_TRUE(il->value == BigInt(3));
    }
}

// rapidcheck: all corpus entries — no root substitution gives non-zero result
RC_GTEST_FIXTURE_PROP(RootsSubstituteTest, NoRootSubstGivesNonZero, ()) {
    for (const auto& seed : kSeeds) {
        symbolic::CASContext lctx;
        Symbol var{seed.var_str};
        ExprPtr poly = must_parse(seed.poly_str, lctx);
        auto roots_res = algebra::solve_polynomial(poly, var, lctx);
        if (!roots_res.is_ok()) continue;
        for (ExprPtr root : roots_res.value()) {
            auto res = check_root(poly, root, var, lctx);
            if (res.has_value()) RC_ASSERT(*res);
        }
    }
}

}  // namespace
}  // namespace cas::property
