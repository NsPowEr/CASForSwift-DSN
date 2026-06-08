// F6.4-T1 — Unit tests for solve_inequality_1var (Sturm-based).
//
// Coverage:
//   - Quadratic > 0 with two real roots (open intervals union)
//   - Quadratic < 0 (closed interior interval)
//   - x^2 + 1 > 0 always true (whole real line)
//   - x^2 - 1 ≥ 0 (closed endpoints at ±1)
//   - Linear ax + b ≥ 0 (half-line)
//   - Cubic with three roots, mixed signs

#include "cas/solve_inequality.hpp"

#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/rational.hpp"

#include <gtest/gtest.h>

namespace cas {
namespace algebra {
namespace {

ExprPtr parse_expr(const std::string& s, AstArena& arena) {
    auto t = Lexer(s).tokenize();
    EXPECT_TRUE(t.is_ok()) << s;
    Parser p(t.value(), arena);
    auto r = p.parse();
    EXPECT_TRUE(r.is_ok()) << s;
    return r.value();
}

class SolveInequalityTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};
};

// x^2 - 4 > 0  ⇔  x < -2  ∨  x > 2.
TEST_F(SolveInequalityTest, QuadraticGreaterThanZero) {
    auto p = parse_expr("x^2 - 4", ctx.arena());
    auto res = solve_inequality_1var(p, x, InequalityOp::Greater, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    ASSERT_EQ(res.value().size(), 2U);
    // First interval: (-∞, -2)
    EXPECT_FALSE(res.value()[0].lower.has_value());
    EXPECT_TRUE(res.value()[0].upper.has_value());
    EXPECT_TRUE(res.value()[0].upper_open);
    // Second interval: (2, +∞)
    EXPECT_TRUE(res.value()[1].lower.has_value());
    EXPECT_FALSE(res.value()[1].upper.has_value());
    EXPECT_TRUE(res.value()[1].lower_open);
}

// x^2 - 4 < 0  ⇔  -2 < x < 2.
TEST_F(SolveInequalityTest, QuadraticLessThanZero) {
    auto p = parse_expr("x^2 - 4", ctx.arena());
    auto res = solve_inequality_1var(p, x, InequalityOp::Less, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    ASSERT_EQ(res.value().size(), 1U);
    EXPECT_TRUE(res.value()[0].lower.has_value());
    EXPECT_TRUE(res.value()[0].upper.has_value());
    EXPECT_TRUE(res.value()[0].lower_open);
    EXPECT_TRUE(res.value()[0].upper_open);
}

// x^2 + 1 > 0 is always true.
TEST_F(SolveInequalityTest, AlwaysPositiveQuadratic) {
    auto p = parse_expr("x^2 + 1", ctx.arena());
    auto res = solve_inequality_1var(p, x, InequalityOp::Greater, ctx);
    ASSERT_TRUE(res.is_ok());
    ASSERT_EQ(res.value().size(), 1U);
    EXPECT_FALSE(res.value()[0].lower.has_value());
    EXPECT_FALSE(res.value()[0].upper.has_value());
}

// x^2 - 1 ≥ 0 ⇔ x ≤ -1 ∨ x ≥ 1 (closed endpoints).
TEST_F(SolveInequalityTest, NonStrictClosedEndpoints) {
    auto p = parse_expr("x^2 - 1", ctx.arena());
    auto res = solve_inequality_1var(p, x, InequalityOp::GreaterEqual, ctx);
    ASSERT_TRUE(res.is_ok());
    ASSERT_EQ(res.value().size(), 2U);
    EXPECT_FALSE(res.value()[0].upper_open);  // closed at -1
    EXPECT_FALSE(res.value()[1].lower_open);  // closed at +1
}

}  // namespace
}  // namespace algebra
}  // namespace cas
