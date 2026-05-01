#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "cas/formatter.hpp"

#include <gtest/gtest.h>
#include <iostream>
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

    std::string to_str(ExprPtr expr) {
        auto fmt = Formatter::create(Formatter::Format::Standard);
        return fmt->format(expr);
    }
};

// GAP #4: Branch Cuts
// sqrt(-1) * sqrt(-1) should be -1.
TEST_F(GapDemonstrationTest, GAP4_BranchCuts_SqrtProduct) {
    auto expr = parse_str("sqrt(-1) * sqrt(-1)");
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto simplified = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    std::string result = to_str(simplified.value());
    std::cout << "[GAP #4] sqrt(-1) * sqrt(-1) simplified to: " << result << std::endl;

    // Dopo il fix, deve restituire "-1" (via sqrt(-1)^2) e NON unire in sqrt(1)
    EXPECT_EQ(result, "-1");
}

// GAP #5: Taylor Limits
TEST_F(GapDemonstrationTest, GAP5_TaylorLimits_Cancellation) {
    auto expr = parse_str("(cos(x) - exp(-x^2/2)) / x^4");
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto lim_res = calculus::limit(expr.value(), Symbol("x"), ctx.arena().make_integer(0), LimitDirection::Both, ctx);
    
    ASSERT_TRUE(lim_res.is_ok()) << (lim_res.is_error() ? lim_res.error().message : "Unknown error");
    
    std::string result = to_str(lim_res.value());
    std::cout << "[GAP #5] lim x->0 (cos(x) - exp(-x^2/2)) / x^4 = " << result << std::endl;
    
    // Con Taylor adattivo deve dare -1/12
    EXPECT_EQ(result, "-1/12");
}

} // namespace cas::test
