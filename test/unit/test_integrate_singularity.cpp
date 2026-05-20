#include <gtest/gtest.h>
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "cas/formatter.hpp"

using namespace cas;
using namespace cas::calculus;

class IntegrateSingularityTest : public ::testing::Test {
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

    void verify_singularity_rejected(const std::string& expr_str, const std::string& a_str, const std::string& b_str) {
        auto expr = parse(expr_str);
        auto a = parse(a_str);
        auto b = parse(b_str);
        
        auto res = definite_integral(expr, x, a, b, ctx);
        EXPECT_TRUE(res.is_error()) << "Should have rejected integral of " << expr_str << " from " << a_str << " to " << b_str;
        if (res.is_error()) {
            EXPECT_EQ(res.error().kind, CASErrorKind::Undefined);
        }
    }
    
    void verify_definite(const std::string& expr_str, const std::string& a_str, const std::string& b_str, const std::string& expected_str) {
        auto expr = parse(expr_str);
        auto a = parse(a_str);
        auto b = parse(b_str);
        auto expected = parse(expected_str);
        
        auto res = definite_integral(expr, x, a, b, ctx);
        ASSERT_TRUE(res.is_ok()) << "Integral failed for " << expr_str << ": " << res.error().message;
        
        auto simplified_res = ctx.simplify(res.value());
        ASSERT_TRUE(simplified_res.is_ok());
        
        auto simplified_expected = ctx.simplify(expected);
        ASSERT_TRUE(simplified_expected.is_ok());

        EXPECT_TRUE(structural_equal(simplified_res.value(), simplified_expected.value()))
            << "For " << expr_str << " from " << a_str << " to " << b_str 
            << " expected " << expected_str << " but got " << formatter::TextFormatter{}.format(simplified_res.value());
    }
};

TEST_F(IntegrateSingularityTest, RationalPole) {
    // 1/(x-1) from 0 to 2 cross pole at 1
    verify_singularity_rejected("1/(x-1)", "0", "2");
}

TEST_F(IntegrateSingularityTest, AlgebraicSingularity) {
    // 1/sqrt(x) at 0
    verify_singularity_rejected("1/sqrt(x)", "0", "1");
}

TEST_F(IntegrateSingularityTest, TranscendentalSingularity) {
    // 1/sin(x) from -1 to 1 cross 0
    verify_singularity_rejected("1/sin(x)", "-1", "1");
}

TEST_F(IntegrateSingularityTest, TanSingularity) {
    // tan(x) from 0 to pi
    verify_singularity_rejected("tan(x)", "0", "pi");
}
