#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/error_helpers.hpp"
#include "cas/numtheory.hpp"

#include "matrix_smith_internal.hpp"

#include <algorithm>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace cas::linalg {
namespace {

Result<void> combine_rows(MatrixExpr& S, MatrixExpr& U, std::size_t i, std::size_t j, 
                         const BigInt& a, const BigInt& b, const BigInt& c, const BigInt& d,
                         symbolic::CASContext& ctx) {
    auto expr_a = integer(ctx, a);
    auto expr_b = integer(ctx, b);
    auto expr_c = integer(ctx, c);
    auto expr_d = integer(ctx, d);

    for (std::size_t k = 0; k < S.cols(); ++k) {
        auto si = S(i, k);
        auto sj = S(j, k);
        auto prod1 = mul_expr(ctx, expr_a, si);
        if (prod1.is_error()) return fail<void>(prod1.error());
        auto prod2 = mul_expr(ctx, expr_b, sj);
        if (prod2.is_error()) return fail<void>(prod2.error());
        auto val_i = add_expr(ctx, prod1.value(), prod2.value());
        if (val_i.is_error()) return fail<void>(val_i.error());

        auto prod3 = mul_expr(ctx, expr_c, si);
        if (prod3.is_error()) return fail<void>(prod3.error());
        auto prod4 = mul_expr(ctx, expr_d, sj);
        if (prod4.is_error()) return fail<void>(prod4.error());
        auto val_j = add_expr(ctx, prod3.value(), prod4.value());
        if (val_j.is_error()) return fail<void>(val_j.error());

        S(i, k) = val_i.value();
        S(j, k) = val_j.value();
    }
    for (std::size_t k = 0; k < U.cols(); ++k) {
        auto ui = U(i, k);
        auto uj = U(j, k);
        auto prod1 = mul_expr(ctx, expr_a, ui);
        if (prod1.is_error()) return fail<void>(prod1.error());
        auto prod2 = mul_expr(ctx, expr_b, uj);
        if (prod2.is_error()) return fail<void>(prod2.error());
        auto val_i = add_expr(ctx, prod1.value(), prod2.value());
        if (val_i.is_error()) return fail<void>(val_i.error());

        auto prod3 = mul_expr(ctx, expr_c, ui);
        if (prod3.is_error()) return fail<void>(prod3.error());
        auto prod4 = mul_expr(ctx, expr_d, uj);
        if (prod4.is_error()) return fail<void>(prod4.error());
        auto val_j = add_expr(ctx, prod3.value(), prod4.value());
        if (val_j.is_error()) return fail<void>(val_j.error());

        U(i, k) = val_i.value();
        U(j, k) = val_j.value();
    }
    return ok();
}

Result<void> combine_cols(MatrixExpr& S, MatrixExpr& V, std::size_t i, std::size_t j, 
                         const BigInt& a, const BigInt& b, const BigInt& c, const BigInt& d,
                         symbolic::CASContext& ctx) {
    auto expr_a = integer(ctx, a);
    auto expr_b = integer(ctx, b);
    auto expr_c = integer(ctx, c);
    auto expr_d = integer(ctx, d);

    for (std::size_t k = 0; k < S.rows(); ++k) {
        auto si = S(k, i);
        auto sj = S(k, j);
        auto prod1 = mul_expr(ctx, expr_a, si);
        if (prod1.is_error()) return fail<void>(prod1.error());
        auto prod2 = mul_expr(ctx, expr_b, sj);
        if (prod2.is_error()) return fail<void>(prod2.error());
        auto val_i = add_expr(ctx, prod1.value(), prod2.value());
        if (val_i.is_error()) return fail<void>(val_i.error());

        auto prod3 = mul_expr(ctx, expr_c, si);
        if (prod3.is_error()) return fail<void>(prod3.error());
        auto prod4 = mul_expr(ctx, expr_d, sj);
        if (prod4.is_error()) return fail<void>(prod4.error());
        auto val_j = add_expr(ctx, prod3.value(), prod4.value());
        if (val_j.is_error()) return fail<void>(val_j.error());

        S(k, i) = val_i.value();
        S(k, j) = val_j.value();
    }
    for (std::size_t k = 0; k < V.rows(); ++k) {
        auto vi = V(k, i);
        auto vj = V(k, j);
        auto prod1 = mul_expr(ctx, expr_a, vi);
        if (prod1.is_error()) return fail<void>(prod1.error());
        auto prod2 = mul_expr(ctx, expr_b, vj);
        if (prod2.is_error()) return fail<void>(prod2.error());
        auto val_i = add_expr(ctx, prod1.value(), prod2.value());
        if (val_i.is_error()) return fail<void>(val_i.error());

        auto prod3 = mul_expr(ctx, expr_c, vi);
        if (prod3.is_error()) return fail<void>(prod3.error());
        auto prod4 = mul_expr(ctx, expr_d, vj);
        if (prod4.is_error()) return fail<void>(prod4.error());
        auto val_j = add_expr(ctx, prod3.value(), prod4.value());
        if (val_j.is_error()) return fail<void>(val_j.error());

        V(k, i) = val_i.value();
        V(k, j) = val_j.value();
    }
    return ok();
}

void collect_symbols(ExprPtr e, std::set<std::string>& out) {
    if (!e) return;
    switch (e->kind) {
    case ExprKind::Symbol:
        out.insert(expr_ref<Symbol>(e).name);
        return;
    case ExprKind::Unary:
        collect_symbols(expr_cast<Unary>(e)->operand, out);
        return;
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

} // namespace

Result<SmithNormalForm> smith_normal_form(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    std::set<std::string> all_syms;
    bool all_integer = true;
    for (std::size_t i = 0; i < matrix.rows(); ++i) {
        for (std::size_t j = 0; j < matrix.cols(); ++j) {
            collect_symbols(matrix(i, j), all_syms);
            if (!try_get_bigint(matrix(i, j))) all_integer = false;
        }
    }
    if (!all_integer && all_syms.size() == 1U) {
        Symbol var(*all_syms.begin());
        return smith_detail::smith_normal_form_qx(matrix, var, ctx);
    }
    if (!all_integer) {
        return make_unimplemented<SmithNormalForm>(
            "linalg", "smith_normal_form",
            "matrix entries must be either integer or polynomial in a single variable",
            error::reason_codes::LINALG_SMITH_NON_INTEGER,
            "Convert symbolic entries to a single-variable polynomial form, or expand to numerics",
            "F4.2b");
    }

    const std::size_t rows = matrix.rows();
    const std::size_t cols = matrix.cols();
    MatrixExpr S(rows, cols, matrix.elements());
    auto U_res = identity(rows, ctx);
    auto V_res = identity(cols, ctx);
    if (U_res.is_error()) return fail<SmithNormalForm>(U_res.error());
    if (V_res.is_error()) return fail<SmithNormalForm>(V_res.error());
    MatrixExpr U = std::move(U_res.value());
    MatrixExpr V = std::move(V_res.value());

    std::size_t p = 0;
    while (p < std::min(rows, cols)) {
        if (auto chk = ctx.check_interrupt(); chk.is_error()) return fail<SmithNormalForm>(chk.error());
        std::optional<std::pair<std::size_t, std::size_t>> pivot;
        BigInt min_val;

        for (std::size_t i = p; i < rows; ++i) {
            for (std::size_t j = p; j < cols; ++j) {
                auto val = try_get_bigint(S(i, j));
                if (!val) continue; // Should not happen in Z-path
                if (val->is_zero()) continue;
                
                BigInt abs_val = val->is_negative() ? -*val : *val;
                if (!pivot || abs_val < min_val) {
                    pivot = {i, j};
                    min_val = abs_val;
                }
            }
        }

        if (!pivot) break;

        smith_detail::swap_rows(S, U, p, pivot->first);
        smith_detail::swap_cols(S, V, p, pivot->second);

        bool changed = true;
        std::size_t guard = 0;
        const std::size_t guard_max = (rows + cols) * ctx.smith_stabilization_multiplier();
        while (changed) {
            if (auto chk = ctx.check_interrupt(); chk.is_error()) return fail<SmithNormalForm>(chk.error());
            changed = false;
            if (++guard > guard_max) {
                 return fail<SmithNormalForm>(CASError{CASErrorKind::Unimplemented,
                    "smith_normal_form: pivot stabilization loop exceeded guard",
                    std::nullopt});
            }
            for (std::size_t i = p + 1; i < rows; ++i) {
                auto a = try_get_bigint(S(p, p)).value();
                auto b = try_get_bigint(S(i, p)).value();
                if (b.is_zero()) continue;

                auto gcd_res = numtheory::extended_gcd(a, b);
                if (gcd_res.is_error()) return fail<SmithNormalForm>(gcd_res.error());
                auto [g, s, t] = gcd_res.value();
                
                auto a_red = a / g;
                auto b_red = b / g;
                auto res = combine_rows(S, U, p, i, s, t, -b_red, a_red, ctx);
                if (res.is_error()) return fail<SmithNormalForm>(res.error());
                changed = true;
            }
            for (std::size_t j = p + 1; j < cols; ++j) {
                auto a = try_get_bigint(S(p, p)).value();
                auto b = try_get_bigint(S(p, j)).value();
                if (b.is_zero()) continue;

                auto gcd_res = numtheory::extended_gcd(a, b);
                if (gcd_res.is_error()) return fail<SmithNormalForm>(gcd_res.error());
                auto [g, s, t] = gcd_res.value();
                
                auto a_red = a / g;
                auto b_red = b / g;
                auto res = combine_cols(S, V, p, j, s, t, -b_red, a_red, ctx);
                if (res.is_error()) return fail<SmithNormalForm>(res.error());
                changed = true;
            }
        }

        // Ensure S[p,p] divides all S[i,j] for i,j > p
        for (std::size_t i = p + 1; i < rows; ++i) {
            for (std::size_t j = p + 1; j < cols; ++j) {
                auto a = try_get_bigint(S(p, p)).value();
                auto b = try_get_bigint(S(i, j)).value();
                if ((b % a).is_zero()) continue;

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

        // Normalize S[p,p] positive
        auto val = try_get_bigint(S(p, p)).value();
        if (val.is_negative()) {
            for (std::size_t k = 0; k < cols; ++k) S(p, k) = integer(ctx, -try_get_bigint(S(p, k)).value());
            for (std::size_t k = 0; k < rows; ++k) U(p, k) = integer(ctx, -try_get_bigint(U(p, k)).value());
        }
        ++p;
    }

    return ok(SmithNormalForm{std::move(S), std::move(U), std::move(V)});
}

} // namespace cas::linalg
