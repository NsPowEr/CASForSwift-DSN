#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "cas/error_helpers.hpp"
#include "cas/algebra.hpp"
#include "matrix_structured_determinant.hpp"
#include <cmath>

namespace cas::linalg::detail {

[[nodiscard]] Result<std::optional<ExprPtr>> determinant_tridiagonal_if_applicable(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n == 0) return ok(std::optional<ExprPtr>{integer(ctx, 1)});
    
    // Check if tridiagonal: A(i,j) != 0 only if |i-j| <= 1
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (std::abs(static_cast<int>(i) - static_cast<int>(j)) > 1) {
                if (!is_zero_expr(matrix(i, j))) return ok(std::optional<ExprPtr>{});
            }
        }
    }

    // f_i = a_i * f_{i-1} - b_i * c_{i-1} * f_{i-2}
    std::vector<ExprPtr> f(n + 1);
    f[0] = integer(ctx, 1);
    
    auto a0_res = simplify(ctx, matrix(0, 0));
    if (a0_res.is_error()) return fail<std::optional<ExprPtr>>(a0_res.error());
    f[1] = a0_res.value();

    for (std::size_t i = 2; i <= n; ++i) {
        auto term1_res = mul_expr(ctx, matrix(i - 1, i - 1), f[i - 1]);
        if (term1_res.is_error()) return fail<std::optional<ExprPtr>>(term1_res.error());
        
        auto bc_res = mul_expr(ctx, matrix(i - 1, i - 2), matrix(i - 2, i - 1));
        if (bc_res.is_error()) return fail<std::optional<ExprPtr>>(bc_res.error());
        
        auto term2_res = mul_expr(ctx, bc_res.value(), f[i - 2]);
        if (term2_res.is_error()) return fail<std::optional<ExprPtr>>(term2_res.error());
        
        auto fi_res = sub_expr(ctx, term1_res.value(), term2_res.value());
        if (fi_res.is_error()) return fail<std::optional<ExprPtr>>(fi_res.error());
        f[i] = fi_res.value();
    }

    return ok(std::optional<ExprPtr>{f[n]});
}

[[nodiscard]] Result<std::optional<ExprPtr>> determinant_vandermonde_if_applicable(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n < 2) return ok(std::optional<ExprPtr>{});

    // Check if Vandermonde: A(i, j) = x_i^j (or x_i^{j-1})
    std::vector<ExprPtr> x(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (!is_one_expr(matrix(i, 0))) return ok(std::optional<ExprPtr>{});
        x[i] = matrix(i, 1);
    }

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 2; j < n; ++j) {
            auto prod_res = mul_expr(ctx, matrix(i, j - 1), x[i]);
            if (prod_res.is_error() || prod_res.value() != matrix(i, j)) return ok(std::optional<ExprPtr>{});
        }
    }

    // Det = ∏_{0<=i<j<n} (x_j - x_i)
    ExprPtr det = integer(ctx, 1);
    for (std::size_t j = 1; j < n; ++j) {
        for (std::size_t i = 0; i < j; ++i) {
            auto diff_res = sub_expr(ctx, x[j], x[i]);
            if (diff_res.is_error()) return fail<std::optional<ExprPtr>>(diff_res.error());
            auto next_det = mul_expr(ctx, det, diff_res.value());
            if (next_det.is_error()) return fail<std::optional<ExprPtr>>(next_det.error());
            det = next_det.value();
        }
    }

    return ok(std::optional<ExprPtr>{det});
}

[[nodiscard]] Result<std::optional<ExprPtr>> determinant_circulant_if_applicable(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n < 2) return ok(std::optional<ExprPtr>{});
    
    std::vector<ExprPtr> c(n);
    for (std::size_t j = 0; j < n; ++j) c[j] = matrix(0, j);

    for (std::size_t i = 1; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (matrix(i, j) != c[(j - i + n) % n]) return ok(std::optional<ExprPtr>{});
        }
    }

    if (n == 2) {
        auto diff = sub_expr(ctx, c[0], c[1]);
        auto sum = add_expr(ctx, c[0], c[1]);
        if (diff.is_error() || sum.is_error()) return fail<std::optional<ExprPtr>>(CASError{CASErrorKind::InternalError, "Circulant internal error", std::nullopt});
        auto prod = mul_expr(ctx, diff.value(), sum.value());
        if (prod.is_error()) return fail<std::optional<ExprPtr>>(prod.error());
        return ok(std::optional<ExprPtr>{prod.value()});
    }
    
    if (n == 3) {
        auto s1 = add_expr(ctx, c[0], c[1]);
        if (s1.is_error()) return fail<std::optional<ExprPtr>>(s1.error());
        auto s = add_expr(ctx, s1.value(), c[2]);
        if (s.is_error()) return fail<std::optional<ExprPtr>>(s.error());
        
        auto m01 = mul_expr(ctx, c[0], c[1]);
        auto m12 = mul_expr(ctx, c[1], c[2]);
        auto m02 = mul_expr(ctx, c[0], c[2]);
        
        auto sq0 = mul_expr(ctx, c[0], c[0]);
        auto sq1 = mul_expr(ctx, c[1], c[1]);
        auto sq2 = mul_expr(ctx, c[2], c[2]);
        
        if (m01.is_error() || m12.is_error() || m02.is_error() || sq0.is_error() || sq1.is_error() || sq2.is_error())
             return fail<std::optional<ExprPtr>>(CASError{CASErrorKind::InternalError, "Circulant internal error", std::nullopt});
             
        auto q1 = add_expr(ctx, sq0.value(), sq1.value());
        if (q1.is_error()) return fail<std::optional<ExprPtr>>(q1.error());
        auto q2 = add_expr(ctx, q1.value(), sq2.value());
        if (q2.is_error()) return fail<std::optional<ExprPtr>>(q2.error());
        auto q3 = sub_expr(ctx, q2.value(), m01.value());
        if (q3.is_error()) return fail<std::optional<ExprPtr>>(q3.error());
        auto q4 = sub_expr(ctx, q3.value(), m12.value());
        if (q4.is_error()) return fail<std::optional<ExprPtr>>(q4.error());
        auto q5 = sub_expr(ctx, q4.value(), m02.value());
        if (q5.is_error()) return fail<std::optional<ExprPtr>>(q5.error());
        
        auto prod = mul_expr(ctx, s.value(), q5.value());
        if (prod.is_error()) return fail<std::optional<ExprPtr>>(prod.error());
        return ok(std::optional<ExprPtr>{prod.value()});
    }

    if (n == 4) {
        auto d02 = sub_expr(ctx, c[0], c[2]);
        auto d13 = sub_expr(ctx, c[1], c[3]);
        if (d02.is_error() || d13.is_error()) return fail<std::optional<ExprPtr>>(CASError{CASErrorKind::InternalError, "Circulant internal error", std::nullopt});

        auto s_pos_a = add_expr(ctx, c[0], c[1]);
        auto s_pos_b = add_expr(ctx, c[2], c[3]);
        if (s_pos_a.is_error() || s_pos_b.is_error()) return fail<std::optional<ExprPtr>>(CASError{CASErrorKind::InternalError, "Circulant internal error", std::nullopt});
        auto s_pos = add_expr(ctx, s_pos_a.value(), s_pos_b.value());
        
        auto s_alt_a = sub_expr(ctx, c[0], c[1]);
        auto s_alt_b = sub_expr(ctx, c[2], c[3]);
        if (s_alt_a.is_error() || s_alt_b.is_error()) return fail<std::optional<ExprPtr>>(CASError{CASErrorKind::InternalError, "Circulant internal error", std::nullopt});
        auto s_alt = add_expr(ctx, s_alt_a.value(), s_alt_b.value()); 
        
        if (s_pos.is_error() || s_alt.is_error()) return fail<std::optional<ExprPtr>>(CASError{CASErrorKind::InternalError, "Circulant internal error", std::nullopt});
        
        auto d02_sq = mul_expr(ctx, d02.value(), d02.value());
        auto d13_sq = mul_expr(ctx, d13.value(), d13.value());
        if (d02_sq.is_error() || d13_sq.is_error()) return fail<std::optional<ExprPtr>>(CASError{CASErrorKind::InternalError, "Circulant internal error", std::nullopt});
        
        auto q = add_expr(ctx, d02_sq.value(), d13_sq.value());
        if (q.is_error()) return fail<std::optional<ExprPtr>>(q.error());
        
        auto h1 = mul_expr(ctx, s_pos.value(), s_alt.value());
        if (h1.is_error()) return fail<std::optional<ExprPtr>>(h1.error());
        
        auto res = mul_expr(ctx, h1.value(), q.value());
        if (res.is_error()) return fail<std::optional<ExprPtr>>(res.error());
        return ok(std::optional<ExprPtr>{res.value()});
    }

    return ok(std::optional<ExprPtr>{});
}

[[nodiscard]] Result<std::optional<ExprPtr>> determinant_toeplitz_if_applicable(const MatrixExpr&, symbolic::CASContext&) {
    // HARDCODE-OF-PASSAGE: HC-F43-TOEPLITZ
    // Fast-path Toeplitz non implementato (richiede algoritmo Trench/Levinson
    // simbolico, O(n²) vs O(n³) Bareiss). Caller cade su general path.
    // Vedi HARDCODE_LEDGER.md HC-F43-TOEPLITZ.
    return ok(std::optional<ExprPtr>{});
}

[[nodiscard]] Result<std::optional<ExprPtr>> determinant_banded_if_applicable(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    // Bandwidth detection: cerca il più piccolo k tale che A(i,j)=0 per |i-j|>k.
    // Caso k=0 (diagonale) e k=1 (tridiagonale) sono gestiti dai detector dedicati;
    // qui copriamo k>=2.
    const std::size_t n = matrix.rows();
    if (n != matrix.cols() || n < 3) return ok(std::optional<ExprPtr>{});

    std::size_t bandwidth = 0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const std::size_t d = (i > j) ? (i - j) : (j - i);
            if (d > bandwidth && !is_zero_expr(matrix(i, j))) {
                bandwidth = d;
            }
        }
    }
    // k<=1 → cade ai detector specializzati (gestiti prima nel dispatcher).
    // k>=n-1 → matrice "piena", nessuna specializzazione utile.
    if (bandwidth <= 1U || bandwidth >= n - 1U) return ok(std::optional<ExprPtr>{});

    // HARDCODE-OF-PASSAGE: HC-F43-BANDED
    // Fast-path banded (k>=2) non implementato. Approccio corretto:
    // LU specializzato per matrici banded O(n·k²) anziché O(n³) Bareiss generale.
    // Per ora ritorniamo nullopt → caller cade su general path (Bareiss/cofactor)
    // che produce risultato corretto ma in O(n³).
    // Vedi HARDCODE_LEDGER.md HC-F43-BANDED.
    (void)ctx;
    return ok(std::optional<ExprPtr>{});
}

} // namespace cas::linalg::detail
