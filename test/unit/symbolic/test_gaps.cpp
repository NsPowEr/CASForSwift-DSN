// T-053: reactivated Gap-demonstration tests, rewritten to use STRUCTURAL
// assertions (CLAUDE.md: "Mai toString() per validare la logica"). The old
// version compared formatter output strings; this version inspects the AST.

#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas::test {

class GapDemonstrationTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    Result<ExprPtr> parse_str(const std::string& input) {
        auto tokens = Lexer(input).tokenize();
        if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
        Parser parser(tokens.value(), ctx.arena());
        return parser.parse();
    }
};

// GAP #4: Branch cut — sqrt(-1) * sqrt(-1) must collapse to the integer -1
// (via sqrt(-1)^2), NOT erroneously combine into sqrt(1) = 1.
TEST_F(GapDemonstrationTest, GAP4_BranchCuts_SqrtProduct) {
    auto expr = parse_str("sqrt(-1) * sqrt(-1)");
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto simplified = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    const auto* lit = expr_cast<IntegerLit>(simplified.value());
    ASSERT_NE(lit, nullptr) << "sqrt(-1)*sqrt(-1) must simplify to an integer literal";
    EXPECT_EQ(lit->value, BigInt(-1));
}

// GAP #5: Taylor-cancellation limit — lim_{x->0} (cos x - e^{-x^2/2}) / x^4 = -1/12.
TEST_F(GapDemonstrationTest, GAP5_TaylorLimits_Cancellation) {
    auto expr = parse_str("(cos(x) - exp(-x^2/2)) / x^4");
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto lim_res = calculus::limit(expr.value(), Symbol("x"),
        ctx.arena().make<IntegerLit>(BigInt(0)), LimitDirection::Both, ctx);
    ASSERT_TRUE(lim_res.is_ok())
        << (lim_res.is_error() ? lim_res.error().message : "Unknown error");

    // Expected exact value: -1/12.
    const auto* rat = expr_cast<RationalLit>(lim_res.value());
    ASSERT_NE(rat, nullptr) << "limit must be the exact rational -1/12";
    EXPECT_EQ(rat->numerator, BigInt(-1));
    EXPECT_EQ(rat->denominator, BigInt(12));
}

}  // namespace cas::test
