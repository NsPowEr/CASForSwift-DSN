#pragma once

#include "cas/linalg/Matrix.hpp"

#include <cstddef>
#include <optional>

namespace cas::linalg::detail {

[[nodiscard]] Result<std::optional<ExprPtr>> determinant_tridiagonal_if_applicable(
    const MatrixExpr& matrix,
    symbolic::CASContext& ctx);

// F4.3 — Vandermonde: M[i][j] = x_i^j (colonna 0 = 1). Riconosce la
// struttura e applica formula chiusa det = ∏_{i<j} (x_j − x_i).
[[nodiscard]] Result<std::optional<ExprPtr>> determinant_vandermonde_if_applicable(
    const MatrixExpr& matrix,
    symbolic::CASContext& ctx);

// F4.3 — Banda generale di bandwidth bw: M[i][j] = 0 per |i−j| > bw.
// Esegue Bareiss specializzato che salta entry esterne alla banda.
[[nodiscard]] Result<std::optional<ExprPtr>> determinant_banded_if_applicable(
    const MatrixExpr& matrix,
    symbolic::CASContext& ctx);

// F4.3 — Circulant: M[i][j] = c_{(j−i) mod n}. det = ∏_{k=0}^{n-1} P(ω^k),
// con ω = e^{2πi/n} e P(x) = Σ c_j x^j.  Per n ≤ 4 produce forma chiusa
// senza nested DFT, per n > 4 fallisce esplicitamente (debito documentato).
[[nodiscard]] Result<std::optional<ExprPtr>> determinant_circulant_if_applicable(
    const MatrixExpr& matrix,
    symbolic::CASContext& ctx);

// F4.3 — Toeplitz: M[i][j] = t_{i−j} (costante lungo diagonali).
// Algoritmo Levinson-Durbin O(n²) tramite predizione lineare; per ora
// implementato direttamente come det(LU) via pivoting Toeplitz-aware.
[[nodiscard]] Result<std::optional<ExprPtr>> determinant_toeplitz_if_applicable(
    const MatrixExpr& matrix,
    symbolic::CASContext& ctx);

}  // namespace cas::linalg::detail
