#include "cas/formatter.hpp"
#include "cas/symbolic.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/numeric/sampler.hpp"
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>

namespace cas::formatter {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

} // namespace

TEST(FormatterTest, TextFormatterPrecedence) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("(x + y) * z", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    
    TextFormatter fmt;
    EXPECT_EQ(fmt.format(expr.value()), "(x + y) * z");
    
    auto expr2 = parse_expr("x + y * z", ctx.arena());
    EXPECT_EQ(fmt.format(expr2.value()), "x + y * z");
}

TEST(FormatterTest, LaTeXFormatterBasic) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("1/2 + sqrt(x)", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    
    LaTeXFormatter fmt;
    EXPECT_EQ(fmt.format(expr.value()), "\\frac{1}{2} + \\sqrt{x}");
}

TEST(FormatterTest, Ascii2DNestedFractions) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("1/2 + 3/4", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    
    Ascii2DFormatter fmt;
    std::string res = fmt.format(expr.value());
    
    std::vector<std::string> lines;
    std::stringstream ss(res);
    std::string line;
    while (std::getline(ss, line)) {
        if (!line.empty()) lines.push_back(line);
    }
    
    ASSERT_EQ(lines.size(), 3U);
    EXPECT_TRUE(lines[1].find("+") != std::string::npos);
}

TEST(FormatterTest, Ascii2DDeepNesting) {
    symbolic::CASContext ctx;
    // 1 / (2 / 3)
    auto tokens = Lexer("1 / (2 / 3)").tokenize();
    Parser parser(tokens.value(), ctx.arena());
    auto expr = parser.parse();
    ASSERT_TRUE(expr.is_ok());
    
    Ascii2DFormatter fmt;
    std::string res = fmt.format(expr.value());
    std::cout << "DEEP NESTING:\n" << res << std::endl;
    
    std::vector<std::string> lines;
    std::stringstream ss(res);
    std::string line;
    while (std::getline(ss, line)) {
        if (!line.empty()) lines.push_back(line);
    }
    
    //  1 
    // ---
    //  2 
    // ---
    //  3 
    ASSERT_EQ(lines.size(), 5U);
}

TEST(SamplerTest, AdaptiveSamplerFindsPeaks) {
    symbolic::CASContext ctx;
    auto tokens = Lexer("exp(-(x^2 * 100))").tokenize();
    Parser parser(tokens.value(), ctx.arena());
    auto expr = parser.parse();
    ASSERT_TRUE(expr.is_ok());
    
    numeric::AdaptiveSampler sampler(expr.value(), "x");
    auto res = sampler.sample(-1.0, 1.0);
    ASSERT_TRUE(res.is_ok());
    
    bool found_peak = false;
    for (const auto& p : res.value()) {
        if (std::abs(p.x) < 0.05 && p.y > 0.9) found_peak = true;
    }
    EXPECT_TRUE(found_peak);
    EXPECT_GT(res.value().size(), 20U);
}

} // namespace cas::formatter
