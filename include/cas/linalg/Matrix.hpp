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
    std::vector<T>& elements() { return data_; }
    T& operator()(size_t r, size_t c) { return data_[r * cols_ + c]; }
    const T& operator()(size_t r, size_t c) const { return data_[r * cols_ + c]; }
    
    void fill(const T& value) {
        std::fill(data_.begin(), data_.end(), value);
    }
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
Result<MatrixExpr> bareiss(const MatrixExpr& matrix, symbolic::CASContext& ctx);

// CAS-L3-17 — Symbolic LU decomposition (Doolittle form).
// Returns (L, U) such that A = L · U (no pivoting). Fails if a zero
// pivot is encountered. For numerical stability/robustness use the
// permuted variant LU_with_pivoting (TODO follow-up).
//   L: lower triangular with unit diagonal (size n×n)
//   U: upper triangular (size n×n)
struct LUDecomposition {
    MatrixExpr L;
    MatrixExpr U;
};
[[nodiscard]] Result<LUDecomposition> lu_decompose(const MatrixExpr& matrix,
                                                    symbolic::CASContext& ctx);

// Solve A·x = b via existing LU factorization: forward-substitute L·y=b
// then back-substitute U·x=y. Returns x as vector of n ExprPtr.
[[nodiscard]] Result<std::vector<ExprPtr>> lu_solve(
    const LUDecomposition& lu, const std::vector<ExprPtr>& b,
    symbolic::CASContext& ctx);

// L3-17 — Permuted LU with partial pivoting: P·A = L·U.
// Risolve matrici con zero (o near-zero) pivots scegliendo riga via
// PivotScore (preferisce IntegerLit/RationalLit, poi nonzero certo,
// poi simbolico). P = permutation as vector of row indices.
struct PLUDecomposition {
    std::vector<std::size_t> P;  // row permutation: P[i] = original row at position i
    MatrixExpr L;
    MatrixExpr U;
};
[[nodiscard]] Result<PLUDecomposition> lu_decompose_pivoted(
    const MatrixExpr& matrix, symbolic::CASContext& ctx);

// L3-17 — QR decomposition via classical Gram-Schmidt (symbolic).
// A = Q · R where Q has orthonormal columns and R upper triangular.
// MVP: supporta solo matrici a coefficienti razionali; le entry di Q
// possono coinvolgere sqrt() simbolico per le norme.
struct QRDecomposition {
    MatrixExpr Q;
    MatrixExpr R;
};
[[nodiscard]] Result<QRDecomposition> qr_decompose(
    const MatrixExpr& matrix, symbolic::CASContext& ctx);

// F4.1c — Cholesky LDL^T per matrici simmetriche.
// A = L · D · L^T, con L unit-lower triangular e D diagonale (storata come
// vettore di n entries). Non richiede sqrt (rispetto a Cholesky classico L·L^T)
// e funziona anche per matrici indefinite finché tutti i D_j != 0 (no
// Bunch-Kaufman pivoting in questa versione: fallisce esplicitamente con
// reason code se D_j = 0).  Requisito: A simmetrica (verificato).
struct LDLTDecomposition {
    MatrixExpr L;
    std::vector<ExprPtr> D;
};
[[nodiscard]] Result<LDLTDecomposition> cholesky_ldlt(
    const MatrixExpr& matrix, symbolic::CASContext& ctx);

// F4.2d — Companion matrix di p(x) = c_0 + c_1 x + ... + c_{n-1} x^{n-1} + x^n.
// Restituisce la matrice n×n con 1 sulla sub-diagonale e [-c_0, -c_1, ...,
// -c_{n-1}] nell'ultima colonna (forma di Frobenius standard).  Richiede
// che `polynomial` sia monico in `var` (lc = 1 dopo normalizzazione).
[[nodiscard]] Result<MatrixExpr> companion_matrix(
    ExprPtr polynomial, const Symbol& var, symbolic::CASContext& ctx);

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
Result<std::vector<std::vector<ExprPtr>>> null_space(const MatrixExpr& matrix, symbolic::CASContext& ctx);

// Compute null space of `matrix` treating `alpha_expr` (a RootOf with rational
// minimal polynomial) as an algebraic generator.  Each matrix entry is reduced
// in Q(alpha) before RREF, allowing exact kernel computation when the entries
// are polynomial expressions in alpha that ordinary symbolic RREF cannot pivot
// on (because zero-equivalence cannot be decided structurally).
// Falls back to ordinary null_space() when entries are not in Q(alpha) or the
// supplied minimal polynomial is reducible.
Result<std::vector<std::vector<ExprPtr>>> null_space_over_extension(
    const MatrixExpr& matrix,
    ExprPtr alpha_expr,
    symbolic::CASContext& ctx);

struct JordanDecomposition {
    MatrixExpr J; // Jordan Normal Form matrix
    MatrixExpr P; // Transformation matrix such that A = P*J*P^-1
};

Result<JordanDecomposition> jordan_normal_form(const MatrixExpr& matrix, symbolic::CASContext& ctx);
Result<ExprPtr> determinant_modular(const MatrixExpr& matrix, symbolic::CASContext& ctx);

struct SmithNormalForm {
    MatrixExpr S; // Smith Normal Form matrix
    MatrixExpr U; // Left transformation matrix (U*A*V = S)
    MatrixExpr V; // Right transformation matrix
};

Result<SmithNormalForm> smith_normal_form(const MatrixExpr& matrix, symbolic::CASContext& ctx);

// F4.2c — Hermite Normal Form (HNF). Per A m×n su PID R esiste U unimodulare
// m×m tale che U·A = H è upper-triangular con entries "ridotte":
//   - su Z:    0 ≤ H[i][j] < H[i][i] per j > i (modulo pivot)
//   - su Q[x]: deg(H[i][j]) < deg(H[i][i]) per j > i
// Routing: matrice tutta-integer → Z path; entries polinomiali in 1 var → Q[x].
struct HermiteNormalForm {
    MatrixExpr H;
    MatrixExpr U;
};
[[nodiscard]] Result<HermiteNormalForm> hermite_normal_form(
    const MatrixExpr& matrix, symbolic::CASContext& ctx);
}
