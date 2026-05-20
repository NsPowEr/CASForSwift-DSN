#pragma once

#include "cas/linalg/Matrix.hpp"

#include <optional>

namespace cas::linalg::detail {

[[nodiscard]] Result<std::optional<ExprPtr>> determinant_tridiagonal_if_applicable(
    const MatrixExpr& matrix,
    symbolic::CASContext& ctx);

}  // namespace cas::linalg::detail
