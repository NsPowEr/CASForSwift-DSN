#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

namespace cas::symbolic {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    return Parser(tokens.value(), arena).parse();
}

bool math_eq(const std::string& a, const std::string& b) {
    CASContext ctx;
    auto ea = parse_expr(a, ctx.arena());
    auto eb = parse_expr(b, ctx.arena());
    if (ea.is_error() || eb.is_error()) return false;
    auto r = mathematically_equal(ea.value(), eb.value(), ctx);
    return r.is_ok() && r.value();
}

// --- P0-001 acceptance criteria ---

TEST(MathematicalEquality, ExpandedSquare) {
    EXPECT_TRUE(math_eq("(x+1)^2", "x^2+2*x+1"));
}

TEST(MathematicalEquality, DifferenceOfSquares) {
    EXPECT_TRUE(math_eq("x^2-1", "(x-1)*(x+1)"));
}

TEST(MathematicalEquality, CommutativeAddition) {
    EXPECT_TRUE(math_eq("x+y", "y+x"));
}

TEST(MathematicalEquality, NegativeCaseDistinctPolynomials) {
    EXPECT_FALSE(math_eq("x^2+1", "x^2+2"));
}

// --- P0-001 anti-hardcode tests ---

TEST(MathematicalEquality, VariableIndependence_z) {
    EXPECT_TRUE(math_eq("(z+1)^2", "z^2+2*z+1"));
}

TEST(MathematicalEquality, LargeCoefficients) {
    EXPECT_TRUE(math_eq("(x+100)^2", "x^2+200*x+10000"));
}

TEST(MathematicalEquality, NegativeCaseOffByOne) {
    EXPECT_FALSE(math_eq("(x+1)^2", "x^2+x+1"));
}

TEST(MathematicalEquality, MultivariateExpansion) {
    EXPECT_TRUE(math_eq("(a+b)^2", "a^2+2*a*b+b^2"));
}

TEST(MathematicalEquality, HigherDegreePolynomial) {
    EXPECT_TRUE(math_eq("(x+1)^3", "x^3+3*x^2+3*x+1"));
}

TEST(MathematicalEquality, ProductOfLinears) {
    EXPECT_TRUE(math_eq("(x+2)*(x-3)", "x^2-x-6"));
}

TEST(MathematicalEquality, SymmetryMultivariate) {
    EXPECT_TRUE(math_eq("(x+y)*(x-y)", "x^2-y^2"));
}

TEST(MathematicalEquality, ScaledEquality) {
    EXPECT_TRUE(math_eq("2*(x+1)", "2*x+2"));
}

TEST(MathematicalEquality, NegativeCase_MultivariateDistinct) {
    EXPECT_FALSE(math_eq("x^2+y^2", "(x+y)^2"));
}

}  // namespace
}  // namespace cas::symbolic
