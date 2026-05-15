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
    ExprPtr alpha_expr,
    const MinPoly& mp,
    symbolic::CASContext& ctx) {
    AlgMatrix out(matrix.rows(), matrix.cols(), make_alg_zero(mp));
    for (std::size_t r = 0; r < matrix.rows(); ++r) {
        for (std::size_t c = 0; c < matrix.cols(); ++c) {
            auto opt = cas::algebra::try_express_in_q_alpha(matrix(r, c), alpha_expr, mp, ctx);
            if (opt.is_error()) {
                return fail<std::optional<AlgMatrix>>(opt.error());
            }
            if (!opt.value().has_value()) {
                return ok(std::optional<AlgMatrix>{});
            }
            out.at(r, c) = std::move(opt.value().value());
        }
    }
    return ok(std::optional<AlgMatrix>(std::move(out)));
}

// Reduced row echelon form, in place, over Q(alpha).
// Returns an error if a pivot inversion fails (which indicates the supplied
// minimal polynomial is reducible: alpha generates a ring with zero divisors,
// not a field).
[[nodiscard]] Result<void> rref_alg(AlgMatrix& m) {
    std::size_t pivot_row = 0U;
    for (std::size_t lead = 0U; pivot_row < m.rows && lead < m.cols; ++lead) {
        std::size_t candidate = pivot_row;
        while (candidate < m.rows && m.at(candidate, lead).is_zero()) {
            ++candidate;
        }
        if (candidate == m.rows) continue;
        if (candidate != pivot_row) {
            for (std::size_t k = 0; k < m.cols; ++k) {
                std::swap(m.at(pivot_row, k), m.at(candidate, k));
            }
        }
        auto inv_pivot = m.at(pivot_row, lead).inverse();
        if (inv_pivot.is_error()) {
            return Result<void>(inv_pivot.error());
        }
        // Scale row.
        for (std::size_t k = 0; k < m.cols; ++k) {
            m.at(pivot_row, k) = m.at(pivot_row, k) * inv_pivot.value();
        }
        // Eliminate other rows.
        for (std::size_t r2 = 0; r2 < m.rows; ++r2) {
            if (r2 == pivot_row) continue;
            AlgebraicNumber factor = m.at(r2, lead);
            if (factor.is_zero()) continue;
            for (std::size_t k = 0; k < m.cols; ++k) {
                AlgebraicNumber product = factor * m.at(pivot_row, k);
                m.at(r2, k) = m.at(r2, k) - product;
            }
        }
        ++pivot_row;
    }
    return ok();
}

}  // namespace

Result<std::vector<std::vector<ExprPtr>>> null_space_over_extension(
    const MatrixExpr& matrix,
    ExprPtr alpha_expr,
    symbolic::CASContext& ctx) {
    if (!alpha_expr) {
        return null_space(matrix, ctx);
    }
    // Canonicalize alpha_expr through the simplifier so its structural form
    // matches whatever appears inside matrix entries that themselves passed
    // through ctx.simplify().  Without this normalization a fresh RootOf node
    // built ad-hoc by the caller can fail structural_equal against the
    // simplifier-canonical RootOf that ends up inside matrix entries.
    {
        auto canon = ctx.simplify(alpha_expr);
        if (canon.is_ok()) alpha_expr = canon.value();
    }
    const auto* root = expr_cast<RootOf>(alpha_expr);
    if (!root) {
        return null_space(matrix, ctx);
    }

    auto mp_res = cas::algebra::rootof_min_poly(*root, ctx);
    if (mp_res.is_error()) {
        return null_space(matrix, ctx);
    }
    const auto& mp = mp_res.value();

    auto am_res = build_alg_matrix(matrix, alpha_expr, mp, ctx);
    if (am_res.is_error()) {
        return fail<std::vector<std::vector<ExprPtr>>>(am_res.error());
    }
    if (!am_res.value().has_value()) {
        return null_space(matrix, ctx);
    }
    AlgMatrix m = std::move(am_res.value().value());

    auto rref_res = rref_alg(m);
    if (rref_res.is_error()) {
        // Minimal polynomial likely reducible.  Fall back to structural path.
        return null_space(matrix, ctx);
    }

    // Identify pivot columns.
    std::vector<std::optional<std::size_t>> pivot_row_for_col(m.cols, std::nullopt);
    std::vector<bool> is_pivot(m.cols, false);
    std::size_t current_row = 0U;
    for (std::size_t c = 0; c < m.cols && current_row < m.rows; ++c) {
        if (!m.at(current_row, c).is_zero()) {
            is_pivot[c] = true;
            pivot_row_for_col[c] = current_row;
            ++current_row;
        }
    }

    AstArena& arena = ctx.arena();
    std::vector<std::vector<ExprPtr>> basis;
    for (std::size_t free_col = 0; free_col < m.cols; ++free_col) {
        if (is_pivot[free_col]) continue;

        std::vector<ExprPtr> vec(m.cols);
        for (auto& slot : vec) slot = arena.make<IntegerLit>(BigInt(0));
        vec[free_col] = arena.make<IntegerLit>(BigInt(1));

        for (std::size_t pivot_c = 0; pivot_c < m.cols; ++pivot_c) {
            if (!is_pivot[pivot_c] || !pivot_row_for_col[pivot_c].has_value()) continue;
            const std::size_t pr = *pivot_row_for_col[pivot_c];
            const AlgebraicNumber& coeff = m.at(pr, free_col);
            if (coeff.is_zero()) continue;
            AlgebraicNumber neg = -coeff;
            ExprPtr raw = cas::algebra::algebraic_number_to_expr_raw(neg, alpha_expr, arena);
            auto simp = ctx.simplify(raw);
            if (simp.is_error()) {
                return fail<std::vector<std::vector<ExprPtr>>>(simp.error());
            }
            vec[pivot_c] = simp.value();
        }
        basis.push_back(std::move(vec));
    }
    return ok(std::move(basis));
}

}  // namespace cas::linalg
