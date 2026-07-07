#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "cas/symbolic.hpp"
#include <gtest/gtest.h>
#include <string>

namespace cas::linalg {
namespace {

ExprPtr symbol(symbolic::CASContext& ctx, std::string name) {
    return ctx.arena().make<Symbol>(std::move(name));
}

ExprPtr pow_expr(symbolic::CASContext& ctx, ExprPtr base, long long exp) {
    return ctx.arena().make<Binary>(BinaryOp::Pow, base, integer(ctx, exp));
}

ExprPtr add_expr_raw(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    return ctx.arena().make<Binary>(BinaryOp::Add, lhs, rhs);
}

}  // namespace

TEST(MatrixBareissTest, Vandermonde4x4) {
    symbolic::CASContext context;
    MatrixExpr m(4U, 4U);
    for (int i = 0; i < 4; ++i) {
        ExprPtr xi = symbol(context, "x" + std::to_string(i));
        for (int j = 0; j < 4; ++j) {
            if (j == 0) m(i, j) = integer(context, 1);
            else if (j == 1) m(i, j) = xi;
            else m(i, j) = pow_expr(context, xi, j);
        }
    }

    auto det_res = determinant(m, context);
    ASSERT_TRUE(det_res.is_ok());

    ExprPtr det = det_res.value();
    EXPECT_NE(det, nullptr);
}

TEST(MatrixBareissTest, Symbolic2x2) {
    symbolic::CASContext context;
    MatrixExpr m(2U, 2U);
    ExprPtr a = symbol(context, "a");
    ExprPtr b = symbol(context, "b");
    ExprPtr c = symbol(context, "c");
    ExprPtr d = symbol(context, "d");

    m(0, 0) = a; m(0, 1) = b;
    m(1, 0) = c; m(1, 1) = d;

    auto det_res = determinant(m, context);
    ASSERT_TRUE(det_res.is_ok());
    EXPECT_NE(det_res.value(), nullptr);
}

TEST(MatrixBareissTest, ZeroPivotHandling) {
    symbolic::CASContext ctx;
    MatrixExpr m(4U, 4U);
    ExprPtr x = symbol(ctx, "x");
    ExprPtr y = symbol(ctx, "y");
    ExprPtr z = symbol(ctx, "z");
    ExprPtr w = symbol(ctx, "w");

    // Row 0 has 0 in (0,0) position
    m(0, 0) = integer(ctx, 0); m(0, 1) = x; m(0, 2) = integer(ctx, 2); m(0, 3) = integer(ctx, 1);
    m(1, 0) = integer(ctx, 0); m(1, 1) = y; m(1, 2) = integer(ctx, 3); m(1, 3) = integer(ctx, 1);
    m(2, 0) = integer(ctx, 2); m(2, 1) = integer(ctx, 0); m(2, 2) = z; m(2, 3) = integer(ctx, 4);
    m(3, 0) = integer(ctx, 1); m(3, 1) = integer(ctx, 2); m(3, 2) = integer(ctx, 0); m(3, 3) = w;

    auto res = determinant(m, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NE(res.value(), nullptr);
}

TEST(MatrixBareissTest, SingleElement) {
    symbolic::CASContext ctx;
    MatrixExpr m(1U, 1U);
    m(0, 0) = symbol(ctx, "alpha");

    auto res = determinant(m, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_EQ(res.value(), m(0, 0));
}

TEST(MatrixBareissTest, IdentityMatrix) {
    symbolic::CASContext ctx;
    MatrixExpr m(3U, 3U);
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            if (i == j) {
                m(i, j) = integer(ctx, 1);
            } else {
                m(i, j) = integer(ctx, 0);
            }
        }
    }

    auto res = determinant(m, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NE(res.value(), nullptr);
}

TEST(MatrixBareissTest, ComplexNumericEvaluation) {
    symbolic::CASContext ctx;
    ExprPtr complex_sum = add_expr_raw(ctx, integer(ctx, 3), integer(ctx, 4));
    MatrixExpr m(2U, 2U);
    m(0, 0) = complex_sum; m(0, 1) = integer(ctx, 1);
    m(1, 0) = integer(ctx, 1); m(1, 1) = integer(ctx, 1);

    auto res = determinant(m, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NE(res.value(), nullptr);
}

TEST(MatrixBareissTest, SymbolicElimination4x4CorrectnessAndSize) {
    symbolic::CASContext ctx;
    auto x = symbol(ctx, "x");
    auto y = symbol(ctx, "y");
    auto z = symbol(ctx, "z");
    auto w = symbol(ctx, "w");

    MatrixExpr mat(4U, 4U, {
        x, integer(ctx, 1), integer(ctx, 0), integer(ctx, 2),
        integer(ctx, 0), y, integer(ctx, 3), integer(ctx, 1),
        integer(ctx, 2), integer(ctx, 0), z, integer(ctx, 4),
        integer(ctx, 1), integer(ctx, 2), integer(ctx, 0), w
    });

    auto res = bareiss(mat, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    const auto& res_mat = res.value();

    for (std::size_t i = 1; i < 4; ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            EXPECT_TRUE(is_zero_expr(res_mat(i, j)));
        }
    }

    std::size_t cpx_corner = estimate_complexity(res_mat(3, 3));
    EXPECT_GT(cpx_corner, 0U);
    EXPECT_LT(cpx_corner, 1000U);
}

TEST(MatrixBareissTest, SymbolicElimination6x6CorrectnessAndSize) {
    symbolic::CASContext ctx;
    MatrixExpr mat(6U, 6U);
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            if (i == j) {
                mat(i, j) = symbol(ctx, "s" + std::to_string(i));
            } else if (i + 1 == j || j + 1 == i) {
                mat(i, j) = integer(ctx, 1);
            } else {
                mat(i, j) = integer(ctx, 0);
            }
        }
    }

    auto res = bareiss(mat, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    const auto& res_mat = res.value();

    for (std::size_t i = 1; i < 6; ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            EXPECT_TRUE(is_zero_expr(res_mat(i, j)));
        }
    }

    std::size_t cpx_corner = estimate_complexity(res_mat(5, 5));
    EXPECT_GT(cpx_corner, 0U);
    EXPECT_LT(cpx_corner, 5000U);
}

TEST(MatrixBareissTest, ContextualPivotSelectionPrefersMonomialsOverSums) {
    symbolic::CASContext ctx;
    auto a = symbol(ctx, "a");
    auto b = symbol(ctx, "b");
    auto c = symbol(ctx, "c");

    ExprPtr complex_sum = add_expr_raw(ctx, add_expr_raw(ctx, a, b), c);
    ExprPtr simple_sym = a;

    MatrixExpr mat(2U, 2U, {
        complex_sum, integer(ctx, 2),
        simple_sym, integer(ctx, 3)
    });

    auto res = bareiss(mat, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
}

}  // namespace cas::linalg
