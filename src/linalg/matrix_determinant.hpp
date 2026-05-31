#pragma once

#include "cas/linalg/Matrix.hpp"

namespace cas::linalg::detail {

// F4.1d — Bareiss fraction-free Gaussian elimination per il calcolo del
// determinante di una matrice quadrata simbolica. PivotScore selection
// (preferisce letterali numerici, poi simbolici strutturalmente nonzero
// penalizzati per complessità). Restituisce 0 se la matrice è singolare
// (pivot zero in qualche colonna senza alternativa valida).
[[nodiscard]] Result<ExprPtr> bareiss_determinant(
    const MatrixExpr& matrix, symbolic::CASContext& ctx);

}  // namespace cas::linalg::detail
