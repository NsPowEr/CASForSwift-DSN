#include "cas/linalg/Matrix.hpp"
#include "cas/numtheory.hpp"
#include "cas/ast.hpp"
#include "cas/error_helpers.hpp"

#include <algorithm>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::linalg {
namespace {

[[nodiscard]] ExprPtr integer(symbolic::CASContext& ctx, const BigInt& value) {
    return ctx.arena().make<IntegerLit>(value);
}

[[nodiscard]] std::optional<BigInt> try_get_bigint(ExprPtr expr) {
    if (const auto* il = expr_cast<IntegerLit>(expr)) return il->value;
    return std::nullopt;
}

[[nodiscard]] Result<ExprPtr> add_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    return ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Add, lhs, rhs));
}

[[nodiscard]] Result<ExprPtr> mul_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    return ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Mul, lhs, rhs));
}

// Row operations on S and U
void swap_rows(MatrixExpr& S, MatrixExpr& U, std::size_t i, std::size_t j) {
    if (i == j) return;
    for (std::size_t k = 0; k < S.cols(); ++k) std::swap(S(i, k), S(j, k));
    for (std::size_t k = 0; k < U.cols(); ++k) std::swap(U(i, k), U(j, k));
}

// Column operations on S and V
void swap_cols(MatrixExpr& S, MatrixExpr& V, std::size_t i, std::size_t j) {
    if (i == j) return;
    for (std::size_t k = 0; k < S.rows(); ++k) std::swap(S(k, i), S(k, j));
    for (std::size_t k = 0; k < V.rows(); ++k) std::swap(V(k, i), V(k, j));
}

Result<void> combine_rows(MatrixExpr& S, MatrixExpr& U, std::size_t i, std::size_t j, 
                         const BigInt& a, const BigInt& b, const BigInt& c, const BigInt& d,
                         symbolic::CASContext& ctx) {
    // S[i,:] = a*S[i,:] + b*S[j,:]
    // S[j,:] = c*S[i,:] + d*S[j,:]
    for (std::size_t k = 0; k < S.cols(); ++k) {
        auto si = S(i, k);
        auto sj = S(j, k);
        auto val_i = add_expr(ctx, mul_expr(ctx, integer(ctx, a), si).value(), mul_expr(ctx, integer(ctx, b), sj).value());
        auto val_j = add_expr(ctx, mul_expr(ctx, integer(ctx, c), si).value(), mul_expr(ctx, integer(ctx, d), sj).value());
        if (val_i.is_error()) return fail<void>(val_i.error());
        if (val_j.is_error()) return fail<void>(val_j.error());
        S(i, k) = val_i.value();
        S(j, k) = val_j.value();
    }
    for (std::size_t k = 0; k < U.cols(); ++k) {
        auto ui = U(i, k);
        auto uj = U(j, k);
        auto val_i = add_expr(ctx, mul_expr(ctx, integer(ctx, a), ui).value(), mul_expr(ctx, integer(ctx, b), uj).value());
        auto val_j = add_expr(ctx, mul_expr(ctx, integer(ctx, c), ui).value(), mul_expr(ctx, integer(ctx, d), uj).value());
        if (val_i.is_error()) return fail<void>(val_i.error());
        if (val_j.is_error()) return fail<void>(val_j.error());
        U(i, k) = val_i.value();
        U(j, k) = val_j.value();
    }
    return ok();
}

Result<void> combine_cols(MatrixExpr& S, MatrixExpr& V, std::size_t i, std::size_t j, 
                         const BigInt& a, const BigInt& b, const BigInt& c, const BigInt& d,
                         symbolic::CASContext& ctx) {
    for (std::size_t k = 0; k < S.rows(); ++k) {
        auto si = S(k, i);
        auto sj = S(k, j);
        auto val_i = add_expr(ctx, mul_expr(ctx, integer(ctx, a), si).value(), mul_expr(ctx, integer(ctx, b), sj).value());
        auto val_j = add_expr(ctx, mul_expr(ctx, integer(ctx, c), si).value(), mul_expr(ctx, integer(ctx, d), sj).value());
        if (val_i.is_error()) return fail<void>(val_i.error());
        if (val_j.is_error()) return fail<void>(val_j.error());
        S(k, i) = val_i.value();
        S(k, j) = val_j.value();
    }
    for (std::size_t k = 0; k < V.rows(); ++k) {
        auto vi = V(k, i);
        auto vj = V(k, j);
        auto val_i = add_expr(ctx, mul_expr(ctx, integer(ctx, a), vi).value(), mul_expr(ctx, integer(ctx, b), vj).value());
        auto val_j = add_expr(ctx, mul_expr(ctx, integer(ctx, c), vi).value(), mul_expr(ctx, integer(ctx, d), vj).value());
        if (val_i.is_error()) return fail<void>(val_i.error());
        if (val_j.is_error()) return fail<void>(val_j.error());
        V(k, i) = val_i.value();
        V(k, j) = val_j.value();
    }
    return ok();
}

} // namespace

Result<SmithNormalForm> smith_normal_form(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
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
        // Find smallest non-zero element in S[p:, p:]
        std::optional<std::pair<std::size_t, std::size_t>> pivot;
        BigInt min_val;

        for (std::size_t i = p; i < rows; ++i) {
            for (std::size_t j = p; j < cols; ++j) {
                auto val = try_get_bigint(S(i, j));
                if (!val) {
                    // F0.8-MIGRATED
                    return make_unimplemented<SmithNormalForm>(
                        "linalg", "smith_normal_form",
                        "matrix entry at (" + std::to_string(i) + "," + std::to_string(j)
                            + ") is not an integer literal",
                        error::reason_codes::LINALG_SMITH_NON_INTEGER,
                        "Smith Normal Form requires integer matrices; "
                        "convert symbolic entries to integers first (e.g., evaluate constants)",
                        "L3-13");
                }
                if (val->is_zero()) continue;
                
                BigInt abs_val = val->is_negative() ? -*val : *val;
                if (!pivot || abs_val < min_val) {
                    pivot = {i, j};
                    min_val = abs_val;
                }
            }
        }

        if (!pivot) break; // All remaining are zero

        swap_rows(S, U, p, pivot->first);
        swap_cols(S, V, p, pivot->second);

        bool changed = true;
        while (changed) {
            changed = false;
            // Row eliminations
            for (std::size_t i = p + 1; i < rows; ++i) {
                auto a = try_get_bigint(S(p, p)).value();
                auto b = try_get_bigint(S(i, p)).value();
                if (b.is_zero()) continue;

                auto gcd_res = numtheory::extended_gcd(a, b);
                if (gcd_res.is_error()) return fail<SmithNormalForm>(gcd_res.error());
                auto [g, s, t] = gcd_res.value();
                
                // New S[p,p] will be g = s*a + t*b
                // New S[i,p] will be 0 = (-b/g)*a + (a/g)*b
                auto a_red = a / g;
                auto b_red = b / g;
                auto res = combine_rows(S, U, p, i, s, t, -b_red, a_red, ctx);
                if (res.is_error()) return fail<SmithNormalForm>(res.error());
                changed = true;
            }

            // Column eliminations
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

                // S[p,:] = S[p,:] + S[i,:]
                for (std::size_t k = 0; k < cols; ++k) {
                    auto val = add_expr(ctx, S(p, k), S(i, k));
                    if (val.is_error()) return fail<SmithNormalForm>(val.error());
                    S(p, k) = val.value();
                }
                for (std::size_t k = 0; k < rows; ++k) {
                    auto val = add_expr(ctx, U(p, k), U(i, k));
                    if (val.is_error()) return fail<SmithNormalForm>(val.error());
                    U(p, k) = val.value();
                }
                
                // Restart eliminations for this p
                changed = true;
                break;
            }
            if (changed) break;
        }
        if (!changed) {
            // Normalize S[p,p] to be positive
            auto val = try_get_bigint(S(p, p)).value();
            if (val.is_negative()) {
                for (std::size_t k = 0; k < cols; ++k) S(p, k) = integer(ctx, -try_get_bigint(S(p, k)).value());
                for (std::size_t k = 0; k < rows; ++k) U(p, k) = integer(ctx, -try_get_bigint(U(p, k)).value());
            }
            ++p;
        }
    }

    return ok(SmithNormalForm{std::move(S), std::move(U), std::move(V)});
}

} // namespace cas::linalg
