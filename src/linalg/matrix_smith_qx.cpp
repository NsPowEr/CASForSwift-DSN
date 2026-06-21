#include "matrix_smith_internal.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/error_helpers.hpp"

#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

namespace cas::linalg {
namespace smith_detail {

// F4.2b — Smith Normal Form su Q[x] (PID polinomiale univariato).
[[nodiscard]] Result<SmithNormalForm> smith_normal_form_qx(
    const MatrixExpr& matrix, const Symbol& var, symbolic::CASContext& ctx) {
    const std::size_t rows = matrix.rows();
    const std::size_t cols = matrix.cols();
    MatrixExpr S(rows, cols, matrix.elements());
    auto U_res = identity(rows, ctx);
    auto V_res = identity(cols, ctx);
    if (U_res.is_error()) return fail<SmithNormalForm>(U_res.error());
    if (V_res.is_error()) return fail<SmithNormalForm>(V_res.error());
    MatrixExpr U = std::move(U_res.value());
    MatrixExpr V = std::move(V_res.value());

    auto row_op = [&](std::size_t i, std::size_t j, ExprPtr a, ExprPtr b, ExprPtr c, ExprPtr d) -> Result<void> {
        for (std::size_t k = 0; k < S.cols(); ++k) {
            ExprPtr si = S(i, k), sj = S(j, k);
            auto ai_res = mul_expr(ctx, a, si);
            if (ai_res.is_error()) return fail<void>(ai_res.error());
            auto bj_res = mul_expr(ctx, b, sj);
            if (bj_res.is_error()) return fail<void>(bj_res.error());
            auto ci_res = mul_expr(ctx, c, si);
            if (ci_res.is_error()) return fail<void>(ci_res.error());
            auto dj_res = mul_expr(ctx, d, sj);
            if (dj_res.is_error()) return fail<void>(dj_res.error());

            auto sum_i = add_expr(ctx, ai_res.value(), bj_res.value());
            if (sum_i.is_error()) return fail<void>(sum_i.error());
            auto ni = algebra::expand(sum_i.value(), ctx);
            if (ni.is_error()) return fail<void>(ni.error());

            auto sum_j = add_expr(ctx, ci_res.value(), dj_res.value());
            if (sum_j.is_error()) return fail<void>(sum_j.error());
            auto nj = algebra::expand(sum_j.value(), ctx);
            if (nj.is_error()) return fail<void>(nj.error());

            S(i, k) = ni.value();
            S(j, k) = nj.value();
        }
        for (std::size_t k = 0; k < U.cols(); ++k) {
            ExprPtr ui = U(i, k), uj = U(j, k);
            auto ai_res = mul_expr(ctx, a, ui);
            if (ai_res.is_error()) return fail<void>(ai_res.error());
            auto bj_res = mul_expr(ctx, b, uj);
            if (bj_res.is_error()) return fail<void>(bj_res.error());
            auto ci_res = mul_expr(ctx, c, ui);
            if (ci_res.is_error()) return fail<void>(ci_res.error());
            auto dj_res = mul_expr(ctx, d, uj);
            if (dj_res.is_error()) return fail<void>(dj_res.error());

            auto sum_i = add_expr(ctx, ai_res.value(), bj_res.value());
            if (sum_i.is_error()) return fail<void>(sum_i.error());
            auto ni = algebra::expand(sum_i.value(), ctx);
            if (ni.is_error()) return fail<void>(ni.error());

            auto sum_j = add_expr(ctx, ci_res.value(), dj_res.value());
            if (sum_j.is_error()) return fail<void>(sum_j.error());
            auto nj = algebra::expand(sum_j.value(), ctx);
            if (nj.is_error()) return fail<void>(nj.error());

            U(i, k) = ni.value();
            U(j, k) = nj.value();
        }
        return ok();
    };

    auto col_op = [&](std::size_t i, std::size_t j, ExprPtr a, ExprPtr b, ExprPtr c, ExprPtr d) -> Result<void> {
        for (std::size_t k = 0; k < S.rows(); ++k) {
            ExprPtr si = S(k, i), sj = S(k, j);
            auto ai_res = mul_expr(ctx, a, si);
            if (ai_res.is_error()) return fail<void>(ai_res.error());
            auto bj_res = mul_expr(ctx, b, sj);
            if (bj_res.is_error()) return fail<void>(bj_res.error());
            auto ci_res = mul_expr(ctx, c, si);
            if (ci_res.is_error()) return fail<void>(ci_res.error());
            auto dj_res = mul_expr(ctx, d, sj);
            if (dj_res.is_error()) return fail<void>(dj_res.error());

            auto sum_i = add_expr(ctx, ai_res.value(), bj_res.value());
            if (sum_i.is_error()) return fail<void>(sum_i.error());
            auto ni = algebra::expand(sum_i.value(), ctx);
            if (ni.is_error()) return fail<void>(ni.error());

            auto sum_j = add_expr(ctx, ci_res.value(), dj_res.value());
            if (sum_j.is_error()) return fail<void>(sum_j.error());
            auto nj = algebra::expand(sum_j.value(), ctx);
            if (nj.is_error()) return fail<void>(nj.error());

            S(k, i) = ni.value();
            S(k, j) = nj.value();
        }
        for (std::size_t k = 0; k < V.rows(); ++k) {
            ExprPtr vi = V(k, i), vj = V(k, j);
            auto ai_res = mul_expr(ctx, a, vi);
            if (ai_res.is_error()) return fail<void>(ai_res.error());
            auto bj_res = mul_expr(ctx, b, vj);
            if (bj_res.is_error()) return fail<void>(bj_res.error());
            auto ci_res = mul_expr(ctx, c, vi);
            if (ci_res.is_error()) return fail<void>(ci_res.error());
            auto dj_res = mul_expr(ctx, d, vj);
            if (dj_res.is_error()) return fail<void>(dj_res.error());

            auto sum_i = add_expr(ctx, ai_res.value(), bj_res.value());
            if (sum_i.is_error()) return fail<void>(sum_i.error());
            auto ni = algebra::expand(sum_i.value(), ctx);
            if (ni.is_error()) return fail<void>(ni.error());

            auto sum_j = add_expr(ctx, ci_res.value(), dj_res.value());
            if (sum_j.is_error()) return fail<void>(sum_j.error());
            auto nj = algebra::expand(sum_j.value(), ctx);
            if (nj.is_error()) return fail<void>(nj.error());

            V(k, i) = ni.value();
            V(k, j) = nj.value();
        }
        return ok();
    };

    auto poly_degree = [&](ExprPtr e) -> Result<std::size_t> {
        if (is_zero_expr(e)) return ok(static_cast<std::size_t>(0U));
        auto deg_res = algebra::polynomial_degree(e, var, ctx);
        if (deg_res.is_error()) return fail<std::size_t>(deg_res.error());
        return ok(static_cast<std::size_t>(deg_res.value()));
    };

    std::size_t p = 0;
    while (p < std::min(rows, cols)) {
        if (auto chk = ctx.check_interrupt(); chk.is_error()) return fail<SmithNormalForm>(chk.error());
        std::optional<std::pair<std::size_t, std::size_t>> pivot;
        std::size_t min_deg = 0;
        for (std::size_t i = p; i < rows; ++i) {
            for (std::size_t j = p; j < cols; ++j) {
                if (is_zero_expr(S(i, j))) continue;
                auto d = poly_degree(S(i, j));
                if (d.is_error()) return fail<SmithNormalForm>(d.error());
                if (!pivot || d.value() < min_deg) {
                    pivot = {i, j};
                    min_deg = d.value();
                }
            }
        }
        if (!pivot) break;

        swap_rows(S, U, p, pivot->first);
        swap_cols(S, V, p, pivot->second);

        bool changed = true;
        std::size_t guard = 0;
        const std::size_t guard_max = (rows + cols) * ctx.smith_stabilization_multiplier();
        while (changed) {
            if (auto chk = ctx.check_interrupt(); chk.is_error()) return fail<SmithNormalForm>(chk.error());
            changed = false;
            if (++guard > guard_max) {
                return fail<SmithNormalForm>(CASError{CASErrorKind::Unimplemented,
                    "smith_normal_form_qx: pivot stabilization loop exceeded guard",
                    std::nullopt});
            }
            // Row elimination via Bezout(S[p][p], S[i][p])
            for (std::size_t i = p + 1; i < rows; ++i) {
                if (is_zero_expr(S(i, p))) continue;
                auto bz = algebra::polynomial_bezout(S(p, p), S(i, p), var, ctx);
                if (bz.is_error()) return fail<SmithNormalForm>(bz.error());
                auto a_red = algebra::polynomial_exact_divide(S(p, p), bz.value().gcd, var, ctx);
                auto b_red = algebra::polynomial_exact_divide(S(i, p), bz.value().gcd, var, ctx);
                if (a_red.is_error()) return fail<SmithNormalForm>(a_red.error());
                if (b_red.is_error()) return fail<SmithNormalForm>(b_red.error());
                
                auto neg_b_res = negate_expr(ctx, b_red.value());
                if (neg_b_res.is_error()) return fail<SmithNormalForm>(neg_b_res.error());
                
                auto r = row_op(p, i, bz.value().s, bz.value().t, neg_b_res.value(), a_red.value());
                if (r.is_error()) return fail<SmithNormalForm>(r.error());
                changed = true;
            }
            for (std::size_t j = p + 1; j < cols; ++j) {
                if (is_zero_expr(S(p, j))) continue;
                auto bz = algebra::polynomial_bezout(S(p, p), S(p, j), var, ctx);
                if (bz.is_error()) return fail<SmithNormalForm>(bz.error());
                auto a_red = algebra::polynomial_exact_divide(S(p, p), bz.value().gcd, var, ctx);
                auto b_red = algebra::polynomial_exact_divide(S(p, j), bz.value().gcd, var, ctx);
                if (a_red.is_error()) return fail<SmithNormalForm>(a_red.error());
                if (b_red.is_error()) return fail<SmithNormalForm>(b_red.error());
                
                auto neg_b_res = negate_expr(ctx, b_red.value());
                if (neg_b_res.is_error()) return fail<SmithNormalForm>(neg_b_res.error());

                auto r = col_op(p, j, bz.value().s, bz.value().t, neg_b_res.value(), a_red.value());
                if (r.is_error()) return fail<SmithNormalForm>(r.error());
                changed = true;
            }
        }

        // Divisibility: S[p][p] divides all S[i,j] for i,j > p
        for (std::size_t i = p + 1; i < rows; ++i) {
            for (std::size_t j = p + 1; j < cols; ++j) {
                if (is_zero_expr(S(i, j))) continue;
                // Bezout again to make S(p,p) = gcd(S(p,p), S(i,j))
                auto bz = algebra::polynomial_bezout(S(p, p), S(i, j), var, ctx);
                if (bz.is_error()) return fail<SmithNormalForm>(bz.error());
                
                // If gcd is same degree as S(p,p), it divides (assuming normalized)
                auto deg_p = poly_degree(S(p, p));
                auto deg_g = poly_degree(bz.value().gcd);
                if (deg_p.is_ok() && deg_g.is_ok() && deg_g.value() == deg_p.value()) continue;

                // Add row i to row p and restart
                for (std::size_t k = 0; k < cols; ++k) {
                    auto res = add_expr(ctx, S(p, k), S(i, k));
                    if (res.is_error()) return fail<SmithNormalForm>(res.error());
                    S(p, k) = res.value();
                }
                for (std::size_t k = 0; k < rows; ++k) {
                    auto res = add_expr(ctx, U(p, k), U(i, k));
                    if (res.is_error()) return fail<SmithNormalForm>(res.error());
                    U(p, k) = res.value();
                }
                changed = true;
                break;
            }
            if (changed) break;
        }
        if (changed) continue;

        // Normalize monic
        if (!is_zero_expr(S(p, p))) {
            auto coeffs = algebra::univariate_coefficients(S(p, p), var, ctx);
            if (coeffs.is_ok() && !coeffs.value().empty()) {
                auto lc = coeffs.value().back();
                if (!is_one_expr(lc)) {
                    auto inv_lc = div_expr(ctx, integer(ctx, 1), lc);
                    if (inv_lc.is_ok()) {
                        for (std::size_t k = 0; k < S.cols(); ++k) {
                            auto res = mul_expr(ctx, inv_lc.value(), S(p, k));
                            if (res.is_ok()) {
                                auto expanded = algebra::expand(res.value(), ctx);
                                S(p, k) = expanded.is_ok() ? expanded.value() : res.value();
                            }
                        }
                        for (std::size_t k = 0; k < U.cols(); ++k) {
                            auto res = mul_expr(ctx, inv_lc.value(), U(p, k));
                            if (res.is_ok()) {
                                auto expanded = algebra::expand(res.value(), ctx);
                                U(p, k) = expanded.is_ok() ? expanded.value() : res.value();
                            }
                        }
                    }
                }
            }
        }
        ++p;
    }

    return ok(SmithNormalForm{std::move(S), std::move(U), std::move(V)});
}

}  // namespace smith_detail
}  // namespace cas::linalg
