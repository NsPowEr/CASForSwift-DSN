#include "polynomial_groebner_fglm_internal.hpp"
#include <algorithm>
#include <utility>

namespace cas::algebra {

[[nodiscard]] std::vector<Rational> mat_vec_mul(
    const std::vector<std::vector<Rational>>& M,
    const std::vector<Rational>& v)
{
    const std::size_t D = v.size();
    std::vector<Rational> result(D, Rational(0));
    for (std::size_t i = 0; i < D; ++i) {
        for (std::size_t j = 0; j < D; ++j) {
            result[i] = result[i] + M[i][j] * v[j];
        }
    }
    return result;
}

[[nodiscard]] std::optional<std::vector<Rational>> linear_dependency(
    const std::vector<std::vector<Rational>>& basis_vecs,
    const std::vector<Rational>& vec)
{
    const std::size_t k = basis_vecs.size();
    const std::size_t D = vec.size();
    if (k == 0) return std::nullopt;

    std::vector<std::vector<Rational>> aug(D, std::vector<Rational>(k + 1, Rational(0)));
    for (std::size_t i = 0; i < D; ++i) {
        for (std::size_t j = 0; j < k; ++j) aug[i][j] = basis_vecs[j][i];
        aug[i][k] = vec[i];
    }

    std::vector<std::size_t> pivot_col(k, D);
    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < k && pivot_row < D; ++col) {
        std::size_t sel = D;
        for (std::size_t row = pivot_row; row < D; ++row) {
            if (!aug[row][col].numerator().is_zero()) { sel = row; break; }
        }
        if (sel == D) continue;
        std::swap(aug[pivot_row], aug[sel]);
        Rational piv = aug[pivot_row][col];
        for (std::size_t c = col; c <= k; ++c) aug[pivot_row][c] = aug[pivot_row][c] / piv;
        for (std::size_t row = 0; row < D; ++row) {
            if (row == pivot_row) continue;
            if (aug[row][col].numerator().is_zero()) continue;
            Rational factor = aug[row][col];
            for (std::size_t c = col; c <= k; ++c) {
                aug[row][c] = aug[row][c] - factor * aug[pivot_row][c];
            }
        }
        pivot_col[col] = pivot_row;
        pivot_row++;
    }

    for (std::size_t row = pivot_row; row < D; ++row) {
        if (!aug[row][k].numerator().is_zero()) return std::nullopt;
    }

    std::vector<Rational> coeffs(k, Rational(0));
    std::size_t pr = 0;
    for (std::size_t col = 0; col < k && pr < D; ++col) {
        if (pivot_col[col] != D) {
            coeffs[col] = aug[pivot_col[col]][k];
            ++pr;
        }
    }
    return coeffs;
}

} // namespace cas::algebra
