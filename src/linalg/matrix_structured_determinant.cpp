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

    // n ≥ 5: closed-form via resultant.
    //
    // Theorem (classical): for the circulant matrix C with first row
    //   (c_0, c_1, …, c_{n-1}) and associated polynomial
    //   P(x) = Σ_{i=0}^{n-1} c_i · x^i,
    // the eigenvalues of C are P(ω^k) for k = 0, …, n-1 where ω = e^{2πi/n}.
    // Hence det(C) = ∏_{k=0}^{n-1} P(ω^k).
    //
    // Since x^n − 1 = ∏_{k=0}^{n-1} (x − ω^k), evaluating P at each root and
    // taking the product is, by definition, the resultant
    //   Res_x(x^n − 1, P(x))   =   ∏_{k} P(ω^k)   =   det(C).
    //
    // The resultant lives in Z[c_0, …, c_{n-1}] (or any extension thereof),
    // so the construction stays in the original entry ring — no detour
    // through Q(ω_n).  References: Davis "Circulant Matrices" Thm 3.2.4,
    // Gantmacher "Matrix Theory" §VIII.6.
    auto& arena = ctx.arena();
    Symbol fresh_x = ctx.make_fresh_symbol("circ_x");
    ExprPtr x_expr = arena.make<Symbol>(fresh_x);

    // Build P(x) = c_0 + c_1·x + c_2·x² + … + c_{n-1}·x^{n-1}.
    std::vector<ExprPtr> p_terms;
    p_terms.reserve(n);
    p_terms.push_back(c[0]);
    for (std::size_t k = 1; k < n; ++k) {
        ExprPtr xk = (k == 1)
            ? x_expr
            : arena.make<Binary>(BinaryOp::Pow, x_expr,
                                  integer(ctx, static_cast<long long>(k)));
        auto term = mul_expr(ctx, c[k], xk);
        if (term.is_error()) return fail<std::optional<ExprPtr>>(term.error());
        p_terms.push_back(term.value());
    }
    ExprPtr p_poly = p_terms.size() == 1U
        ? p_terms[0]
        : arena.make<Sum>(p_terms);
    auto p_simp = simplify(ctx, p_poly);
    if (p_simp.is_ok()) p_poly = p_simp.value();

    // Build x^n − 1.
    ExprPtr x_pow_n = arena.make<Binary>(BinaryOp::Pow, x_expr,
                                          integer(ctx, static_cast<long long>(n)));
    auto x_n_minus_1 = sub_expr(ctx, x_pow_n, integer(ctx, 1));
    if (x_n_minus_1.is_error()) return fail<std::optional<ExprPtr>>(x_n_minus_1.error());

    auto res = cas::algebra::polynomial_resultant(x_n_minus_1.value(), p_poly, fresh_x, ctx);
    if (res.is_error()) {
        // Resultant computation failed (e.g. simplifier hit budget on a
        // symbolic entry).  Fall back to the general Bareiss path silently.
        return ok(std::optional<ExprPtr>{});
    }
    return ok(std::optional<ExprPtr>{res.value()});
}

[[nodiscard]] Result<std::optional<ExprPtr>> determinant_toeplitz_if_applicable(const MatrixExpr&, symbolic::CASContext&) {
    // Toeplitz fast-path intentionally NOT specialized: the symbolic
    // Trench/Levinson recursion requires inverting each leading principal
    // minor as it is computed, which on symbolic entries collapses to the
    // same `is_known_nonzero` decision-procedure problem as the general
    // case — and the asymptotic O(n²) advantage over Bareiss O(n³) is
    // dwarfed by the per-simplify polynomial-arithmetic cost on symbolic
    // entries.  General Bareiss in `bareiss_determinant` (fraction-free,
    // band-agnostic) handles Toeplitz inputs correctly; no information is
    // lost by routing them through it.  See HARDCODE_LEDGER.md
    // HC-F43-TOEPLITZ closure for the design analysis.
    return ok(std::optional<ExprPtr>{});
}

[[nodiscard]] Result<std::optional<ExprPtr>> determinant_banded_if_applicable(const MatrixExpr&, symbolic::CASContext&) {
    // Banded fast-path intentionally NOT specialized for bw ≥ 2.
    //
    // The band-preserving Bareiss optimization (O(n·bw²) inner-loop count vs
    // O(n³) for general Bareiss) is mathematically clean only at the level
    // of inner-loop multiplications: it requires either (a) applying a
    // cumulative pivot/d_prev scaling factor to every in-band entry at
    // every step — which on symbolic inputs degrades to O(n²·bw) total
    // simplify calls and dwarfs the inner-loop saving — or (b) a
    // lazy scale-and-thaw bookkeeping scheme that on symbolic entries
    // trips the same per-simplify polynomial-arithmetic wall as the
    // general path.
    //
    // Diagonal (bw=0) and tridiagonal (bw=1) cases are handled by
    // `determinant_diagonal_if_applicable` and
    // `determinant_tridiagonal_if_applicable` (closed-form three-term
    // recurrences, no symbolic overhead).  For bw ≥ 2 the general Bareiss
    // path in `bareiss_determinant` (fraction-free, band-agnostic) is
    // correct and, on symbolic inputs, asymptotically competitive — no
    // information is lost by routing such matrices through it.  See
    // HARDCODE_LEDGER.md HC-F43-BANDED closure for the design analysis.
    return ok(std::optional<ExprPtr>{});
}

} // namespace cas::linalg::detail
