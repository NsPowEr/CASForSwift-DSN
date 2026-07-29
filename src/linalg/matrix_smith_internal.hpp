#ifndef CAS_LINALG_MATRIX_SMITH_INTERNAL_HPP
#define CAS_LINALG_MATRIX_SMITH_INTERNAL_HPP

#include "cas/linalg/Matrix.hpp"

#include <cstddef>
#include <utility>

// Shared internals of the Smith normal form implementation, split (T-047
// anti-monolith) across two translation units:
//   - matrix_smith.cpp     : integer (Z) path + the public dispatcher
//   - matrix_smith_qx.cpp   : single-variable polynomial Q[x] path
namespace cas::linalg::smith_detail {

// Row/column swaps applied in lock-step to the working matrix S and its
// accumulated unimodular transform (U for rows, V for columns). Shared by
// both the Z and Q[x] reductions, hence inline in this header.
inline void swap_rows(MatrixExpr& S, MatrixExpr& U, std::size_t i, std::size_t j) {
    if (i == j) return;
    for (std::size_t k = 0; k < S.cols(); ++k) std::swap(S(i, k), S(j, k));
    for (std::size_t k = 0; k < U.cols(); ++k) std::swap(U(i, k), U(j, k));
}

inline void swap_cols(MatrixExpr& S, MatrixExpr& V, std::size_t i, std::size_t j) {
    if (i == j) return;
    for (std::size_t k = 0; k < S.rows(); ++k) std::swap(S(k, i), S(k, j));
    for (std::size_t k = 0; k < V.rows(); ++k) std::swap(V(k, i), V(k, j));
}

// F4.2b — Smith Normal Form over Q[x] (univariate polynomial PID). Defined in
// matrix_smith_qx.cpp; called by the dispatcher in matrix_smith.cpp.
Result<SmithNormalForm> smith_normal_form_qx(const MatrixExpr& matrix,
                                             const Symbol& var,
                                             symbolic::CASContext& ctx);

}  // namespace cas::linalg::smith_detail

#endif  // CAS_LINALG_MATRIX_SMITH_INTERNAL_HPP
