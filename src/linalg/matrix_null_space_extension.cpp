// Null space of a symbolic matrix whose entries lie in an algebraic extension
// Q(alpha) of Q, where alpha is the abstract generator of a RootOf expression.
//
// Strategy:
//   1. Convert every matrix entry to an AlgebraicNumber over the supplied
//      generator's minimal polynomial.
//   2. Perform Gauss-Jordan elimination using exact field arithmetic in Q(alpha).
//   3. Read off the kernel basis and render each component back to an ExprPtr.
//
// When any entry is not expressible in Q(alpha), or the minimal polynomial
// turns out to be reducible (inversion fails during RREF), we fall back to the
// generic null_space() routine that operates structurally on ExprPtr.

#include "cas/algebraic_number_bridge.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace cas::linalg {

namespace {

using cas::algebra::AlgebraicNumber;
using MinPoly = AlgebraicNumber::CoeffVec;

struct AlgMatrix {
    std::size_t rows{0U};
    std::size_t cols{0U};
    std::vector<AlgebraicNumber> data;  // row-major

    AlgMatrix(std::size_t r, std::size_t c, const AlgebraicNumber& zero)
        : rows(r), cols(c), data(r * c, zero) {}

    [[nodiscard]] AlgebraicNumber& at(std::size_t r, std::size_t c) noexcept {
        return data[r * cols + c];
    }
    [[nodiscard]] const AlgebraicNumber& at(std::size_t r, std::size_t c) const noexcept {
        return data[r * cols + c];
    }
};

[[nodiscard]] AlgebraicNumber make_alg_zero(const MinPoly& mp) {
    return AlgebraicNumber(AlgebraicNumber::CoeffVec{Rational(BigInt(0))}, mp);
}

[[nodiscard]] Result<std::optional<AlgMatrix>> build_alg_matrix(
    const MatrixExpr& matrix,
    ExprPtr generator,
    const MinPoly& min_poly,
    symbolic::CASContext& ctx) {
    const std::size_t r = matrix.rows();
    const std::size_t c = matrix.cols();
    AlgMatrix result(r, c, make_alg_zero(min_poly));

    for (std::size_t i = 0; i < r; ++i) {
        for (std::size_t j = 0; j < c; ++j) {
            auto conv = cas::algebra::try_express_in_q_alpha(matrix(i, j), generator, min_poly, ctx);
            if (conv.is_error() || !conv.value().has_value()) return ok(std::optional<AlgMatrix>{}); // Fallback
            result.at(i, j) = std::move(conv.value().value());
        }
    }
    return ok(std::optional<AlgMatrix>{std::move(result)});
}

}  // namespace

Result<std::vector<std::vector<ExprPtr>>> null_space_over_extension(
    const MatrixExpr& matrix,
    ExprPtr alpha_expr,
    symbolic::CASContext& ctx) {
    // 1. Estrai polinomio minimo dal generator (deve essere un RootOf).
    const auto* root_of = expr_cast<RootOf>(alpha_expr);
    if (!root_of) {
        return null_space(matrix, ctx); // Fallback
    }

    auto mp_res = cas::algebra::rootof_min_poly(*root_of, ctx);
    if (mp_res.is_error()) return null_space(matrix, ctx);
    const MinPoly min_poly = std::move(mp_res.value());

    // 2. Tenta conversione a matrice algebrica.
    auto alg_res = build_alg_matrix(matrix, alpha_expr, min_poly, ctx);
    if (alg_res.is_error() || !alg_res.value().has_value()) {
        return null_space(matrix, ctx);
    }
    AlgMatrix mat = std::move(*alg_res.value());

    // 3. RREF in Q(alpha).
    const std::size_t r = mat.rows;
    const std::size_t c = mat.cols;
    std::vector<std::size_t> pivot_cols;
    std::size_t curr_row = 0U;

    for (std::size_t j = 0; j < c && curr_row < r; ++j) {
        std::size_t sel = r;
        for (std::size_t i = curr_row; i < r; ++i) {
            if (!mat.at(i, j).is_zero()) { sel = i; break; }
        }
        if (sel == r) continue;

        if (sel != curr_row) {
            for (std::size_t k = j; k < c; ++k) std::swap(mat.at(curr_row, k), mat.at(sel, k));
        }

        auto pivot_inv = mat.at(curr_row, j).inverse();
        if (pivot_inv.is_error()) return null_space(matrix, ctx); // Reducible min-poly fallback

        for (std::size_t k = j; k < c; ++k) mat.at(curr_row, k) = mat.at(curr_row, k) * pivot_inv.value();
        
        for (std::size_t i = 0; i < r; ++i) {
            if (i == curr_row || mat.at(i, j).is_zero()) continue;
            auto factor = mat.at(i, j);
            for (std::size_t k = j; k < c; ++k) {
                mat.at(i, k) = mat.at(i, k) - (factor * mat.at(curr_row, k));
            }
        }
        pivot_cols.push_back(j);
        curr_row++;
    }

    // 4. Kernel basis.
    std::vector<bool> is_pivot(c, false);
    for (auto pc : pivot_cols) is_pivot[pc] = true;

    std::vector<std::vector<ExprPtr>> basis;
    for (std::size_t j = 0; j < c; ++j) {
        if (is_pivot[j]) continue;

        std::vector<ExprPtr> vec(c);
        for (std::size_t k = 0; k < c; ++k) vec[k] = integer(ctx, 0);
        vec[j] = integer(ctx, 1);

        for (std::size_t i = 0; i < pivot_cols.size(); ++i) {
            auto comp = -mat.at(i, j);
            auto expr = cas::algebra::algebraic_number_to_expr(comp, alpha_expr, ctx);
            if (expr.is_error()) return fail<std::vector<std::vector<ExprPtr>>>(expr.error());
            vec[pivot_cols[i]] = expr.value();
        }
        basis.push_back(std::move(vec));
    }

    return ok(std::move(basis));
}

}  // namespace cas::linalg
