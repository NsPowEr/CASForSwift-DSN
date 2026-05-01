#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/linalg/Matrix.hpp"
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

static ExprPtr integer(symbolic::CASContext& ctx, long long value) {
    return ctx.arena().make<IntegerLit>(BigInt(value));
}

// ─── Regression Suite: Algebra ──────────────────────────────────────────

TEST(AlgebraRegression, DegreeZeroPolynomial) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("5", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");

    auto roots = algebra::solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok());
    EXPECT_EQ(roots.value().size(), 0U);

    auto expr2 = parse_expr("x + 5", ctx.arena());
    auto gcd = algebra::polynomial_gcd(expr.value(), expr2.value(), x, ctx);
    ASSERT_TRUE(gcd.is_ok());
    EXPECT_TRUE(expr_is<IntegerLit>(gcd.value()));
}

TEST(AlgebraRegression, FactorZeroPolynomial) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("0", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");

    (void)algebra::factor_over_integers(expr.value(), x, ctx);
}

TEST(AlgebraRegression, DivisionByZeroInAlgebra) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("x / 0", ctx.arena());
    ASSERT_TRUE(expr.is_ok());

    auto simplified = ctx.simplify(expr.value());
    EXPECT_TRUE(simplified.is_error());
    EXPECT_EQ(simplified.error().kind, CASErrorKind::Undefined);
}

TEST(AlgebraRegression, MatrixNonSquareOperations) {
    symbolic::CASContext context;
    linalg::MatrixExpr m(2U, 3U, {
        integer(context, 1), integer(context, 2), integer(context, 3),
        integer(context, 4), integer(context, 5), integer(context, 6),
    });

    auto det = linalg::determinant(m, context);
    EXPECT_TRUE(det.is_error());
    EXPECT_EQ(det.error().kind, CASErrorKind::InvalidArgument);

    auto inv = linalg::inverse(m, context);
    EXPECT_TRUE(inv.is_error());
    EXPECT_EQ(inv.error().kind, CASErrorKind::InvalidArgument);
}

TEST(AlgebraRegression, DeeplyNestedExpressions) {
    symbolic::CASContext ctx;
    std::string deep = "x";
    for(int i = 0; i < 100; ++i) {
        deep = "sin(" + deep + ")";
    }
    auto expr = parse_expr(deep, ctx.arena());
    ASSERT_TRUE(expr.is_ok());

    auto simplified = ctx.simplify(expr.value());
    ASSERT_TRUE(simplified.is_ok());
}

} // namespace
} // namespace cas
