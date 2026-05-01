#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace cas::linalg {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

ExprPtr integer(symbolic::CASContext& ctx, long long value) {
    return ctx.arena().make<IntegerLit>(BigInt(value));
}

ExprPtr symbol(symbolic::CASContext& ctx, std::string name) {
    return ctx.arena().make<Symbol>(std::move(name));
}

void expect_equivalent(ExprPtr actual, const std::string& expected_text, symbolic::CASContext& context) {
    AstArena expected_arena;
    auto expected = parse_expr(expected_text, expected_arena);
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    auto equal = symbolic::mathematically_equal(actual, expected.value(), context);
    ASSERT_TRUE(equal.is_ok()) << equal.error().message;
    EXPECT_TRUE(equal.value()) << "Expected " << expected_text;
}

}  // namespace

TEST(MatrixAdvancedTest, JordanNormalForm3x3NonDiagonalizable) {
    symbolic::CASContext context;
    // Matrix with eigenvalue 2 of multiplicity 3, but only one Jordan block of size 3
    // A = [[2, 1, 0], [0, 2, 1], [0, 0, 2]]
    MatrixExpr A(3U, 3U, {
        integer(context, 2), integer(context, 1), integer(context, 0),
        integer(context, 0), integer(context, 2), integer(context, 1),
        integer(context, 0), integer(context, 0), integer(context, 2)
    });

    auto res = jordan_normal_form(A, context);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    
    const auto& J = res.value().J;
    const auto& P = res.value().P;

    // J should be the same as A (already in Jordan form)
    expect_equivalent(J(0, 0), "2", context);
    expect_equivalent(J(0, 1), "1", context);
    expect_equivalent(J(1, 1), "2", context);
    expect_equivalent(J(1, 2), "1", context);
    expect_equivalent(J(2, 2), "2", context);

    // Verify A * P = P * J
    auto AP = multiply(A, P, context);
    auto PJ = multiply(P, J, context);
    ASSERT_TRUE(AP.is_ok());
    ASSERT_TRUE(PJ.is_ok());

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            auto equal = symbolic::mathematically_equal(AP.value()(i, j), PJ.value()(i, j), context);
            ASSERT_TRUE(equal.is_ok());
            EXPECT_TRUE(equal.value());
        }
    }
}

TEST(MatrixAdvancedTest, JordanNormalFormDefective) {
    symbolic::CASContext context;
    // A = [[3, 1], [0, 3]]
    MatrixExpr A(2U, 2U, {
        integer(context, 3), integer(context, 1),
        integer(context, 0), integer(context, 3)
    });
    auto res = jordan_normal_form(A, context);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_equivalent(res.value().J(0, 1), "1", context);
}

TEST(MatrixAdvancedTest, ModularDeterminant10x10Symbolic) {
    symbolic::CASContext context;
    // Create a 10x10 matrix with symbols and some large integers
    const size_t n = 10;
    MatrixExpr A(n, n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i == j) {
                A(i, j) = symbol(context, "x" + std::to_string(i));
            } else if (i == j + 1 || j == i + 1) {
                A(i, j) = integer(context, 1);
            } else {
                A(i, j) = integer(context, 0);
            }
        }
    }

    // Determinant of tridiagonal matrix with x_i on diagonal and 1 on off-diagonals.
    // This will be symbolic. Our determinant_modular falls back to Bareiss for symbolic.
    auto det_mod = determinant_modular(A, context);
    ASSERT_TRUE(det_mod.is_ok()) << det_mod.error().message;

    auto det_std = determinant(A, context);
    ASSERT_TRUE(det_std.is_ok()) << det_std.error().message;

    auto equal = symbolic::mathematically_equal(det_mod.value(), det_std.value(), context);
    ASSERT_TRUE(equal.is_ok());
    EXPECT_TRUE(equal.value());
}

TEST(MatrixAdvancedTest, ModularDeterminantGiantCoefficients) {
    symbolic::CASContext context;
    // 2x2 matrix with large integers to test CRT
    // Max entry ~ 10^20
    BigInt giant = BigInt::parse("100000000000000000000").value(); // 10^20
    MatrixExpr A(2U, 2U, {
        context.arena().make<IntegerLit>(giant), integer(context, 1),
        integer(context, 1), context.arena().make<IntegerLit>(giant)
    });

    auto det_mod = determinant_modular(A, context);
    ASSERT_TRUE(det_mod.is_ok()) << det_mod.error().message;
    
    // det = giant^2 - 1
    BigInt expected = giant * giant - BigInt(1);
    expect_equivalent(det_mod.value(), expected.decimal(), context);
}

}  // namespace cas::linalg
