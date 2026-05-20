#include <gtest/gtest.h>
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "cas/formatter.hpp"

using namespace cas;
using namespace cas::calculus;

class LimitMrvTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};

    ExprPtr parse(const std::string& input) {
        auto tokens = Lexer(input).tokenize();
        if (tokens.is_error()) {
            throw std::runtime_error("Lex error: " + tokens.error().message);
        }

        Parser parser(tokens.value(), ctx.arena());
        auto res = parser.parse();
        if (res.is_error()) {
            throw std::runtime_error("Parse error: " + res.error().message);
        }
        return res.value();
    }

    void verify_limit(const std::string& expr_str, const std::string& expected_str) {
        auto expr = parse(expr_str);
        auto expected = parse(expected_str);
        auto inf = ctx.arena().make<Constant>(MathConstant::Infinity);
        
        auto res = limit(expr, x, inf, LimitDirection::Both, ctx);
        ASSERT_TRUE(res.is_ok()) << "Limit failed for " << expr_str << ": " << res.error().message;
        
        auto simplified_res = ctx.simplify(res.value());
        ASSERT_TRUE(simplified_res.is_ok());
        
        auto simplified_expected = ctx.simplify(expected);
        ASSERT_TRUE(simplified_expected.is_ok());

        EXPECT_TRUE(structural_equal(simplified_res.value(), simplified_expected.value()))
            << "For " << expr_str << " limit(x->inf) expected " << expected_str 
            << " but got " << formatter::TextFormatter{}.format(simplified_res.value());
    }
};

TEST_F(LimitMrvTest, CancellationSimple) {
    // (x + 1) - x -> 1
    verify_limit("(x + 1) - x", "1");
}

TEST_F(LimitMrvTest, CancellationExponential) {
    // (exp(x) + 1) - exp(x) -> 1
    verify_limit("(exp(x) + 1) - exp(x)", "1");
}

TEST_F(LimitMrvTest, CancellationTower) {
    // exp(x + exp(-x)) - exp(x) -> 1
    verify_limit("exp(x + exp(-x)) - exp(x)", "1");
}

TEST_F(LimitMrvTest, ComplexCancellation) {
    // (exp(x) + x) / (exp(x) + 1) -> 1
    verify_limit("(exp(x) + x) / (exp(x) + 1)", "1");
}
