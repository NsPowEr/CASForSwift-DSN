// F7.5.A2 / HC-F75-A2-MATRIX-SCALAR-OP — unit coverage for the new
// top-level expression evaluator in test/golden/matrix_adapter.hpp.
//
// The adapter handles inputs that mix scalar subexpressions with matrix
// literals `[[…]]` via standard +/-/* (and / by scalar) trees, recursing
// into operands and dispatching to cas::linalg::add/subtract/multiply
// plus an element-wise scalar multiplier. The tests below validate the
// dispatch surface — they intentionally bypass the Maxima oracle so the
// coverage is independent of any external golden run.

#include <gtest/gtest.h>

#include "matrix_adapter.hpp"

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class MatrixAdapterD2Test : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] cas::golden::MatrixOrScalar eval(const std::string& s) {
        auto r = cas::golden::evaluate_matrix_expression(s, ctx);
        EXPECT_TRUE(r.is_ok()) << s << " — " << (r.is_ok() ? "" : r.error().message);
        return r.is_ok() ? r.value() : cas::golden::MatrixOrScalar{};
    }

    [[nodiscard]] bool integer_value(ExprPtr e, long long n) {
        if (!e) return false;
        const auto* lit = expr_cast<IntegerLit>(e);
        return lit != nullptr && lit->value == BigInt(n);
    }
};

TEST_F(MatrixAdapterD2Test, LeafMatrixLiteral) {
    auto v = eval("[[1,2],[3,4]]");
    ASSERT_TRUE(v.is_matrix);
    ASSERT_EQ(v.matrix.rows(), 2U);
    ASSERT_EQ(v.matrix.cols(), 2U);
    EXPECT_TRUE(integer_value(v.matrix(0, 0), 1));
    EXPECT_TRUE(integer_value(v.matrix(1, 1), 4));
}

TEST_F(MatrixAdapterD2Test, LeafScalar) {
    auto v = eval("42");
    ASSERT_FALSE(v.is_matrix);
    EXPECT_TRUE(integer_value(v.scalar, 42));
}

TEST_F(MatrixAdapterD2Test, ScalarTimesMatrix) {
    // 2 * [[1,2],[3,4]] = [[2,4],[6,8]]
    auto v = eval("2 * [[1,2],[3,4]]");
    ASSERT_TRUE(v.is_matrix);
    EXPECT_TRUE(integer_value(v.matrix(0, 0), 2));
    EXPECT_TRUE(integer_value(v.matrix(0, 1), 4));
    EXPECT_TRUE(integer_value(v.matrix(1, 0), 6));
    EXPECT_TRUE(integer_value(v.matrix(1, 1), 8));
}

TEST_F(MatrixAdapterD2Test, MatrixTimesScalar) {
    auto v = eval("[[1,2],[3,4]] * 3");
    ASSERT_TRUE(v.is_matrix);
    EXPECT_TRUE(integer_value(v.matrix(0, 0), 3));
    EXPECT_TRUE(integer_value(v.matrix(1, 1), 12));
}

TEST_F(MatrixAdapterD2Test, MatrixPlusMatrix) {
    auto v = eval("[[1,2],[3,4]] + [[5,6],[7,8]]");
    ASSERT_TRUE(v.is_matrix);
    EXPECT_TRUE(integer_value(v.matrix(0, 0), 6));
    EXPECT_TRUE(integer_value(v.matrix(0, 1), 8));
    EXPECT_TRUE(integer_value(v.matrix(1, 0), 10));
    EXPECT_TRUE(integer_value(v.matrix(1, 1), 12));
}

TEST_F(MatrixAdapterD2Test, MatrixMinusMatrix) {
    auto v = eval("[[5,6],[7,8]] - [[1,2],[3,4]]");
    ASSERT_TRUE(v.is_matrix);
    EXPECT_TRUE(integer_value(v.matrix(0, 0), 4));
    EXPECT_TRUE(integer_value(v.matrix(1, 1), 4));
}

TEST_F(MatrixAdapterD2Test, MatrixTimesMatrixIdentity) {
    auto v = eval("[[1,0],[0,1]] * [[7,8],[9,10]]");
    ASSERT_TRUE(v.is_matrix);
    EXPECT_TRUE(integer_value(v.matrix(0, 0), 7));
    EXPECT_TRUE(integer_value(v.matrix(0, 1), 8));
    EXPECT_TRUE(integer_value(v.matrix(1, 0), 9));
    EXPECT_TRUE(integer_value(v.matrix(1, 1), 10));
}

TEST_F(MatrixAdapterD2Test, MatrixTimesMatrixGeneral) {
    // [[1,2],[3,4]] * [[5,6],[7,8]] = [[19,22],[43,50]]
    auto v = eval("[[1,2],[3,4]] * [[5,6],[7,8]]");
    ASSERT_TRUE(v.is_matrix);
    EXPECT_TRUE(integer_value(v.matrix(0, 0), 19));
    EXPECT_TRUE(integer_value(v.matrix(0, 1), 22));
    EXPECT_TRUE(integer_value(v.matrix(1, 0), 43));
    EXPECT_TRUE(integer_value(v.matrix(1, 1), 50));
}

TEST_F(MatrixAdapterD2Test, ScalarMatrixPlusMatrix) {
    // 2*[[1,2],[3,4]] + [[1,0],[0,1]] = [[3,4],[6,9]]
    auto v = eval("2*[[1,2],[3,4]] + [[1,0],[0,1]]");
    ASSERT_TRUE(v.is_matrix);
    EXPECT_TRUE(integer_value(v.matrix(0, 0), 3));
    EXPECT_TRUE(integer_value(v.matrix(0, 1), 4));
    EXPECT_TRUE(integer_value(v.matrix(1, 0), 6));
    EXPECT_TRUE(integer_value(v.matrix(1, 1), 9));
}

TEST_F(MatrixAdapterD2Test, UnaryMinusMatrix) {
    auto v = eval("-[[1,2],[3,4]]");
    ASSERT_TRUE(v.is_matrix);
    EXPECT_TRUE(integer_value(v.matrix(0, 0), -1));
    EXPECT_TRUE(integer_value(v.matrix(1, 1), -4));
}

TEST_F(MatrixAdapterD2Test, ParenthesizedExpression) {
    // (2 * [[1,1],[1,1]]) + [[0,1],[1,0]] = [[2,3],[3,2]]
    auto v = eval("(2 * [[1,1],[1,1]]) + [[0,1],[1,0]]");
    ASSERT_TRUE(v.is_matrix);
    EXPECT_TRUE(integer_value(v.matrix(0, 0), 2));
    EXPECT_TRUE(integer_value(v.matrix(0, 1), 3));
}

TEST_F(MatrixAdapterD2Test, ScalarPlusMatrixIsRejected) {
    auto r = cas::golden::evaluate_matrix_expression("1 + [[1,2],[3,4]]", ctx);
    ASSERT_FALSE(r.is_ok());
    EXPECT_EQ(r.error().kind, CASErrorKind::Unimplemented);
}

TEST_F(MatrixAdapterD2Test, DivisionByMatrixIsRejected) {
    auto r = cas::golden::evaluate_matrix_expression("[[1,2],[3,4]] / [[1,0],[0,1]]", ctx);
    ASSERT_FALSE(r.is_ok());
    EXPECT_EQ(r.error().kind, CASErrorKind::Unimplemented);
}

TEST_F(MatrixAdapterD2Test, MatrixDividedByScalar) {
    auto v = eval("[[2,4],[6,8]] / 2");
    ASSERT_TRUE(v.is_matrix);
    EXPECT_TRUE(integer_value(v.matrix(0, 0), 1));
    EXPECT_TRUE(integer_value(v.matrix(1, 1), 4));
}

TEST_F(MatrixAdapterD2Test, MattraceWrapperEvaluatesTrace) {
    // HC-F75-A2-MAXIMA-MATTRACE: Maxima leaves mattrace(matrix(...))
    // unevaluated; helper must compute trace on the CAS side.
    auto r = cas::golden::try_evaluate_mattrace_wrapper(
        "mattrace(matrix([1,2],[3,4]))", ctx);
    ASSERT_TRUE(r.has_value());
    ASSERT_TRUE(r->is_ok()) << r->error().message;
    EXPECT_TRUE(integer_value(r->value(), 5));
}

TEST_F(MatrixAdapterD2Test, MattraceWrapperWithWhitespaceAndTerminator) {
    auto r = cas::golden::try_evaluate_mattrace_wrapper(
        "  mattrace(matrix([7,0],[0,3])) ;", ctx);
    ASSERT_TRUE(r.has_value());
    ASSERT_TRUE(r->is_ok());
    EXPECT_TRUE(integer_value(r->value(), 10));
}

TEST_F(MatrixAdapterD2Test, MattraceWrapperSkipsNonMattrace) {
    auto r = cas::golden::try_evaluate_mattrace_wrapper(
        "matrix([1,2],[3,4])", ctx);
    EXPECT_FALSE(r.has_value());
}

TEST_F(MatrixAdapterD2Test, MattraceWrapperRejectsMalformed) {
    auto r = cas::golden::try_evaluate_mattrace_wrapper(
        "mattrace(not_a_matrix(1,2))", ctx);
    ASSERT_TRUE(r.has_value());
    EXPECT_FALSE(r->is_ok());
}

TEST_F(MatrixAdapterD2Test, NestedScalarExpression) {
    // (2+1)*[[1,1],[1,1]] = [[3,3],[3,3]]
    auto v = eval("(2+1)*[[1,1],[1,1]]");
    ASSERT_TRUE(v.is_matrix);
    EXPECT_TRUE(integer_value(v.matrix(0, 0), 3));
    EXPECT_TRUE(integer_value(v.matrix(1, 1), 3));
}

}  // namespace
