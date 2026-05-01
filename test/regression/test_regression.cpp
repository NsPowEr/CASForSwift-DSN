#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas {
namespace {

[[nodiscard]] Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

// ─── Regression Suite ───────────────────────────────────────────────────

// BUG FIX: Issue with division by zero not being caught in some cases
TEST(Regression, DivisionByZero) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("1 / 0", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simplified = ctx.simplify(expr.value());
    EXPECT_TRUE(simplified.is_error());
    EXPECT_EQ(simplified.error().kind, CASErrorKind::Undefined);
}

// BUG FIX: Nested power simplification overflow or depth issues
TEST(Regression, DeeplyNestedPowers) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("((x^2)^2)^2", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto simplified = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified.is_ok());
    
    auto expected = parse_expr("x^8", ctx.arena());
    ASSERT_TRUE(expected.is_ok());
    EXPECT_TRUE(structural_equal(simplified.value(), expected.value()));
}

// BUG FIX: Degree zero polynomial roots
TEST(Regression, DegreeZeroPolynomial) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("10", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    auto roots = algebra::solve_polynomial(expr.value(), Symbol("x"), ctx);
    ASSERT_TRUE(roots.is_ok());
    EXPECT_EQ(roots.value().size(), 0U);
}

} // namespace
} // namespace cas
