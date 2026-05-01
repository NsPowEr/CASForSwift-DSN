#pragma once
#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include <cstddef>
#include <vector>

namespace cas::linalg {
template <typename T>
class Matrix {
public:
    Matrix(size_t rows, size_t cols) : rows_(rows), cols_(cols), data_(rows * cols) {}
    Matrix(size_t rows, size_t cols, const std::vector<T>& data) : rows_(rows), cols_(cols), data_(data) {
        data_.resize(rows * cols);
    }
    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }
    const std::vector<T>& elements() const { return data_; }
    T& operator()(size_t r, size_t c) { return data_[r * cols_ + c]; }
    const T& operator()(size_t r, size_t c) const { return data_[r * cols_ + c]; }
private:
    size_t rows_, cols_;
    std::vector<T> data_;
};
using MatrixExpr = Matrix<ExprPtr>;
Result<MatrixExpr> add(const MatrixExpr& a, const MatrixExpr& b, symbolic::CASContext& ctx);
Result<MatrixExpr> subtract(const MatrixExpr& a, const MatrixExpr& b, symbolic::CASContext& ctx);
Result<MatrixExpr> multiply(const MatrixExpr& a, const MatrixExpr& b, symbolic::CASContext& ctx);
Result<MatrixExpr> transpose(const MatrixExpr& matrix);
Result<MatrixExpr> identity(std::size_t n, symbolic::CASContext& ctx);
Result<ExprPtr> determinant(const MatrixExpr& matrix, symbolic::CASContext& ctx);
Result<MatrixExpr> inverse(const MatrixExpr& matrix, symbolic::CASContext& ctx);
Result<MatrixExpr> rref(const MatrixExpr& matrix, symbolic::CASContext& ctx);
Result<std::vector<ExprPtr>> linsolve(const MatrixExpr& a, const std::vector<ExprPtr>& b, symbolic::CASContext& ctx);
Result<std::size_t> rank(const MatrixExpr& matrix, symbolic::CASContext& ctx);
Result<ExprPtr> trace(const MatrixExpr& matrix, symbolic::CASContext& ctx);

struct Eigenpair {
    ExprPtr eigenvalue;
    std::vector<ExprPtr> eigenvector;
};

Result<ExprPtr> characteristic_polynomial(const MatrixExpr& matrix, const Symbol& lambda_var, symbolic::CASContext& ctx);
Result<std::vector<ExprPtr>> eigenvalues(const MatrixExpr& matrix, symbolic::CASContext& ctx);
Result<std::vector<Eigenpair>> eigenvectors(const MatrixExpr& matrix, symbolic::CASContext& ctx);

struct SmithNormalForm {
    MatrixExpr S; // Smith Normal Form matrix
    MatrixExpr U; // Left transformation matrix (U*A*V = S)
    MatrixExpr V; // Right transformation matrix
};

Result<SmithNormalForm> smith_normal_form(const MatrixExpr& matrix, symbolic::CASContext& ctx);
}
