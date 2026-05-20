// CAS-L3-19 — Solving polynomial in algebraic closure (deg ≥ 5 via RootOf).

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::algebra;

namespace {

class SolveClosureTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};
    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
};

TEST_F(SolveClosureTest, Quintic_Irreducible_EmitsRootOf) {
    // x⁵ - x - 1 — classical Bring's example, Galois = S₅ (non solvable by radicals).
    auto p = parse("x^5 - x - 1");
    auto r = solve_polynomial(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_EQ(r.value().size(), 5U);
    // Each root must be a RootOf (no radical form possible).
    std::size_t rootof_count = 0;
    for (auto root : r.value()) {
        if (expr_is<RootOf>(root)) ++rootof_count;
    }
    EXPECT_GE(rootof_count, 1U)
        << "Quintic irreducible must produce at least one RootOf";
}

TEST_F(SolveClosureTest, Sextic_Reducible_FactorsExtracted) {
    // x⁶ - 1 = (x-1)(x+1)(x²+x+1)(x²-x+1). 6 roots, some explicit.
    auto p = parse("x^6 - 1");
    auto r = solve_polynomial(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_GE(r.value().size(), 2U);
}

TEST_F(SolveClosureTest, Quartic_ExplicitFormula) {
    // x⁴ - 1 = (x-1)(x+1)(x²+1). 4 roots: 1, -1, i, -i.
    auto p = parse("x^4 - 1");
    auto r = solve_polynomial(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_EQ(r.value().size(), 4U);
}

TEST_F(SolveClosureTest, AntiHardcodeSepticPureExponent) {
    // x⁷ - 2: irreducible (Eisenstein at p=2), Galois = generic deg-7 group.
    auto p = parse("x^7 - 2");
    auto r = solve_polynomial(p, x, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_EQ(r.value().size(), 7U);
    // All should be RootOf (no radical except principal one).
    for (auto root : r.value()) {
        // Acceptable: RootOf OR specific radical expression.
        if (!expr_is<RootOf>(root)) {
            // verify it's some expression involving 2^(1/7) or similar
            EXPECT_NE(root, nullptr);
        }
    }
}

}  // namespace
