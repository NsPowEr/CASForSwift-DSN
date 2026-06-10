// F7.5.C1 — focused test: x^5/120 must NOT collapse to x^(1/24).
//
// Regression for HC-F75-C1-POW-DIV: simplifier or formatter on
// Pow(x, n) / m (n, m integers) produces wrong AST or wrong display.

#include <gtest/gtest.h>

#include "cas/ast.hpp"
#include "cas/formatter.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

TEST(PowDivSimplify, ParserKeepsPowAndDivSeparate) {
    symbolic::CASContext ctx;
    auto toks = Lexer("x^5/120").tokenize();
    ASSERT_TRUE(toks.is_ok());
    Parser p(toks.value(), ctx.arena());
    auto e = p.parse();
    ASSERT_TRUE(e.is_ok());
    std::string rendered = formatter::TextFormatter().format(e.value());
    int kind = static_cast<int>(e.value()->kind);
    const auto* bin = expr_cast<Binary>(e.value());
    EXPECT_EQ(bin->op, BinaryOp::Div) << "parser produced kind=" << kind
        << ", formatted: '" << rendered << "'";
}

TEST(PowDivSimplify, SimplifyCollapsesToCoefficientTimesPow) {
    symbolic::CASContext ctx;
    auto toks = Lexer("x^5/120").tokenize();
    Parser p(toks.value(), ctx.arena());
    auto e = p.parse();
    ASSERT_TRUE(e.is_ok());
    auto s = ctx.simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    std::string rendered = formatter::TextFormatter().format(s.value());
    // Must contain x^5 (not x^(1/24) or sqrt-like form).
    EXPECT_NE(rendered.find("x^5"), std::string::npos)
        << "simplify produced: " << rendered;
    EXPECT_EQ(rendered.find("x^1/24"), std::string::npos)
        << "simplify produced: " << rendered;
}

TEST(PowDivSimplify, MultiTermTaylorSeries) {
    symbolic::CASContext ctx;
    auto toks = Lexer("1+x+x^2/2+x^3/6+x^4/24+x^5/120").tokenize();
    ASSERT_TRUE(toks.is_ok());
    Parser p(toks.value(), ctx.arena());
    auto e = p.parse();
    ASSERT_TRUE(e.is_ok());
    auto s = ctx.simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    std::string rendered = formatter::TextFormatter().format(s.value());
    EXPECT_NE(rendered.find("x^5"), std::string::npos) << rendered;
    EXPECT_NE(rendered.find("x^4"), std::string::npos) << rendered;
    EXPECT_NE(rendered.find("x^3"), std::string::npos) << rendered;
    EXPECT_EQ(rendered.find("sqrt"), std::string::npos)
        << "should not contain sqrt: " << rendered;
}

}  // namespace
