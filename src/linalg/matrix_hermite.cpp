// CAS-F4.2c — Hermite Normal Form per Z e Q[x].
//
// Per A m×n su un PID R, calcola U unimodulare m×m e H upper-triangular m×n
// con: U·A = H e entries "ridotte" sopra ogni pivot.
//
// Algoritmo (column-wise sweep):
//   for c = 0..n-1:
//     pivot_row = trova prima riga r ≥ c con H[r][c] ≠ 0
//     se non esiste, salta colonna (defective rank)
//     swap_rows(pivot_row, c)
//     per ogni i in (c+1)..m-1 con H[i][c] ≠ 0:
//       (g, s, t) = bezout(H[c][c], H[i][c])
//       (a, b) = (H[c][c]/g, H[i][c]/g)
//       row_op (c,i) ← ((s, t; -b, a)) · (row c, row i)   [unimodulare det=1]
//     riduci H[i][c] per i < c modulo H[c][c] (modulo polinomiale su Q[x]):
//       q = H[i][c] div_exact H[c][c] (Z: floor div centrato; Q[x]: divisione esatta del quoziente)
//       row_op: row i := row i − q · row c

#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/error.hpp"
#include "cas/error_helpers.hpp"
#include "cas/numtheory.hpp"
#include "cas/symbolic.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace cas::linalg {

namespace {

void collect_symbols(ExprPtr e, std::set<std::string>& out) {
    if (!e) return;
    switch (e->kind) {
    case ExprKind::Symbol: out.insert(expr_ref<Symbol>(e).name); return;
    case ExprKind::Unary: collect_symbols(expr_cast<Unary>(e)->operand, out); return;
    case ExprKind::Binary: {
        const auto* b = expr_cast<Binary>(e);
        collect_symbols(b->left, out);
        collect_symbols(b->right, out);
        return;
    }
    case ExprKind::Sum:
        for (auto t : expr_cast<Sum>(e)->terms) collect_symbols(t, out);
        return;
    case ExprKind::Product:
        for (auto f : expr_cast<Product>(e)->factors) collect_symbols(f, out);
        return;
    case ExprKind::FuncCall:
        for (auto a : expr_cast<FuncCall>(e)->args) collect_symbols(a, out);
        return;
    default: return;
    }
}

[[nodiscard]] Result<void> row_op_symbolic(
    MatrixExpr& H, MatrixExpr& U, std::size_t i, std::size_t j,
    ExprPtr a, ExprPtr b, ExprPtr c, ExprPtr d, symbolic::CASContext& ctx) {
    auto apply = [&](MatrixExpr& M) -> Result<void> {
        for (std::size_t k = 0; k < M.cols(); ++k) {
            ExprPtr mi = M(i, k), mj = M(j, k);
            auto prod1 = mul_expr(ctx, a, mi);
            if (prod1.is_error()) return fail<void>(prod1.error());
            auto prod2 = mul_expr(ctx, b, mj);
            if (prod2.is_error()) return fail<void>(prod2.error());
            auto sum1 = add_expr(ctx, prod1.value(), prod2.value());
            if (sum1.is_error()) return fail<void>(sum1.error());
            
            auto prod3 = mul_expr(ctx, c, mi);
            if (prod3.is_error()) return fail<void>(prod3.error());
            auto prod4 = mul_expr(ctx, d, mj);
            if (prod4.is_error()) return fail<void>(prod4.error());
            auto sum2 = add_expr(ctx, prod3.value(), prod4.value());
            if (sum2.is_error()) return fail<void>(sum2.error());

            auto ni = algebra::expand(sum1.value(), ctx);
            auto nj = algebra::expand(sum2.value(), ctx);
            if (ni.is_error()) return fail<void>(ni.error());
            if (nj.is_error()) return fail<void>(nj.error());
            M(i, k) = ni.value();
            M(j, k) = nj.value();
        }
        return ok();
    };
    auto r1 = apply(H);
    if (r1.is_error()) return r1;
    return apply(U);
}

void swap_rows_HU(MatrixExpr& H, MatrixExpr& U, std::size_t i, std::size_t j) {
    if (i == j) return;
    for (std::size_t k = 0; k < H.cols(); ++k) std::swap(H(i, k), H(j, k));
    for (std::size_t k = 0; k < U.cols(); ++k) std::swap(U(i, k), U(j, k));
}

[[nodiscard]] Result<HermiteNormalForm> hermite_normal_form_z(
    const MatrixExpr& A, symbolic::CASContext& ctx) {
    const std::size_t m = A.rows();
    const std::size_t n = A.cols();
    MatrixExpr H(m, n, A.elements());
    auto Ui = identity(m, ctx);
    if (Ui.is_error()) return fail<HermiteNormalForm>(Ui.error());
    MatrixExpr U = std::move(Ui.value());

    std::size_t r = 0;
    for (std::size_t c = 0; c < n && r < m; ++c) {
        std::size_t pivot = m;
        for (std::size_t i = r; i < m; ++i) {
            auto v = try_get_bigint(H(i, c));
            if (v && !v->is_zero()) { pivot = i; break; }
        }
        if (pivot == m) continue;
        swap_rows_HU(H, U, r, pivot);

        for (std::size_t i = r + 1; i < m; ++i) {
            auto bi_opt = try_get_bigint(H(i, c));
            if (!bi_opt || bi_opt->is_zero()) continue;
            auto a_opt = try_get_bigint(H(r, c));
            if (!a_opt) continue;

            auto xg = numtheory::extended_gcd(*a_opt, *bi_opt);
            if (xg.is_error()) return fail<HermiteNormalForm>(xg.error());
            auto [g, s, t] = xg.value();
            BigInt a_red = (*a_opt) / g;
            BigInt b_red = (*bi_opt) / g;
            auto rop = row_op_symbolic(H, U, r, i,
                integer(ctx, s), integer(ctx, t),
                integer(ctx, -b_red), integer(ctx, a_red), ctx);
            if (rop.is_error()) return fail<HermiteNormalForm>(rop.error());
        }
        
        auto piv_opt = try_get_bigint(H(r, c));
        if (piv_opt && piv_opt->is_negative()) {
            for (std::size_t k = 0; k < n; ++k) {
                auto v = try_get_bigint(H(r, k));
                if (v) H(r, k) = integer(ctx, -*v);
            }
            for (std::size_t k = 0; k < m; ++k) {
                auto v = try_get_bigint(U(r, k));
                if (v) U(r, k) = integer(ctx, -*v);
            }
        }
        
        auto piv_v = try_get_bigint(H(r, c));
        if (piv_v && !piv_v->is_zero()) {
            BigInt piv = *piv_v;
            for (std::size_t i = 0; i < r; ++i) {
                auto vi_opt = try_get_bigint(H(i, c));
                if (!vi_opt) continue;
                BigInt vi = *vi_opt;
                BigInt q = vi / piv;
                if ((vi % piv).is_negative()) q = q - BigInt(1);
                if (q.is_zero()) continue;
                
                auto rop = row_op_symbolic(H, U, i, r, integer(ctx, 1), integer(ctx, -q), integer(ctx, 0), integer(ctx, 1), ctx);
                if (rop.is_error()) return fail<HermiteNormalForm>(rop.error());
            }
        }
        ++r;
    }
    return ok(HermiteNormalForm{std::move(H), std::move(U)});
}

[[nodiscard]] Result<HermiteNormalForm> hermite_normal_form_qx(
    const MatrixExpr& A, const Symbol& var, symbolic::CASContext& ctx) {
    const std::size_t m = A.rows();
    const std::size_t n = A.cols();
    MatrixExpr H(m, n, A.elements());
    auto Ui = identity(m, ctx);
    if (Ui.is_error()) return fail<HermiteNormalForm>(Ui.error());
    MatrixExpr U = std::move(Ui.value());

    std::size_t r = 0;
    for (std::size_t c = 0; c < n && r < m; ++c) {
        std::size_t pivot = m;
        std::size_t pivot_deg = 0;
        for (std::size_t i = r; i < m; ++i) {
            if (is_zero_expr(H(i, c))) continue;
            auto d = algebra::polynomial_degree(H(i, c), var, ctx);
            if (d.is_error()) return fail<HermiteNormalForm>(d.error());
            if (pivot == m || d.value() < pivot_deg) {
                pivot = i;
                pivot_deg = d.value();
            }
        }
        if (pivot == m) continue;
        swap_rows_HU(H, U, r, pivot);

        for (std::size_t i = r + 1; i < m; ++i) {
            if (is_zero_expr(H(i, c))) continue;
            auto bz = algebra::polynomial_bezout(H(r, c), H(i, c), var, ctx);
            if (bz.is_error()) return fail<HermiteNormalForm>(bz.error());
            auto a_red = algebra::polynomial_exact_divide(H(r, c), bz.value().gcd, var, ctx);
            auto b_red = algebra::polynomial_exact_divide(H(i, c), bz.value().gcd, var, ctx);
            if (a_red.is_error()) return fail<HermiteNormalForm>(a_red.error());
            if (b_red.is_error()) return fail<HermiteNormalForm>(b_red.error());
            
            auto neg_b = negate_expr(ctx, b_red.value());
            if (neg_b.is_error()) return fail<HermiteNormalForm>(neg_b.error());

            auto rop = row_op_symbolic(H, U, r, i,
                bz.value().s, bz.value().t, neg_b.value(), a_red.value(), ctx);
            if (rop.is_error()) return fail<HermiteNormalForm>(rop.error());
        }

        // C3-FIX: use polynomial_divmod for reduction
        for (std::size_t i = 0; i < r; ++i) {
            if (is_zero_expr(H(i, c))) continue;
            auto dm = algebra::polynomial_divmod(H(i, c), H(r, c), var, ctx);
            if (dm.is_error()) continue; // skip if not polynomial
            
            ExprPtr q = dm.value().quotient;
            if (is_zero_expr(q)) continue;
            
            auto neg_q = negate_expr(ctx, q);
            if (neg_q.is_error()) return fail<HermiteNormalForm>(neg_q.error());

            auto rop = row_op_symbolic(H, U, i, r, integer(ctx, 1), neg_q.value(), integer(ctx, 0), integer(ctx, 1), ctx);
            if (rop.is_error()) return fail<HermiteNormalForm>(rop.error());
        }
        ++r;
    }
    return ok(HermiteNormalForm{std::move(H), std::move(U)});
}

}  // namespace

Result<HermiteNormalForm> hermite_normal_form(const MatrixExpr& matrix,
                                                symbolic::CASContext& ctx) {
    if (matrix.rows() == 0U || matrix.cols() == 0U) {
        return fail<HermiteNormalForm>(CASError{CASErrorKind::InvalidArgument,
            "hermite_normal_form: empty matrix", std::nullopt});
    }
    std::set<std::string> syms;
    bool all_int = true;
    for (std::size_t i = 0; i < matrix.rows(); ++i) {
        for (std::size_t j = 0; j < matrix.cols(); ++j) {
            collect_symbols(matrix(i, j), syms);
            if (!try_get_bigint(matrix(i, j))) all_int = false;
        }
    }
    if (all_int) return hermite_normal_form_z(matrix, ctx);
    if (syms.size() == 1U) {
        Symbol var(*syms.begin());
        return hermite_normal_form_qx(matrix, var, ctx);
    }
    return make_unimplemented<HermiteNormalForm>(
        "linalg", "hermite_normal_form",
        "entries must be integer or polynomial in a single variable",
        error::reason_codes::LINALG_SMITH_NON_INTEGER,
        "Convert entries to integer or single-var polynomial form", "F4.2c");
}

}  // namespace cas::linalg
