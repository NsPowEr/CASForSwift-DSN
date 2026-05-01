#include "cas/linalg/Matrix.hpp"
#include "cas/algebra.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::linalg {
namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] ExprPtr integer(symbolic::CASContext& ctx, long long value) {
    return ctx.arena().make<IntegerLit>(BigInt(value));
}

[[nodiscard]] bool is_zero_expr(ExprPtr expr) {
    if (const auto* il = expr_cast<IntegerLit>(expr)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(expr)) return rl->numerator.is_zero();
    return false;
}

[[nodiscard]] Result<ExprPtr> add_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    return ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Add, lhs, rhs));
}

[[nodiscard]] Result<ExprPtr> sub_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    return ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Sub, lhs, rhs));
}

[[nodiscard]] Result<ExprPtr> mul_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    return ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Mul, lhs, rhs));
}

// Division-free cofactor expansion — correct for symbolic entries where Bareiss pivots
// may not be decidably nonzero (e.g., entries of the form a_ii - lambda).
// DEPRECATED for large matrices; kept as fallback for non-square or small sizes if needed.
[[nodiscard]] Result<ExprPtr> cofactor_det(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n == 0U) return ok(integer(ctx, 1));
    if (n == 1U) return ok(matrix(0U, 0U));

    if (n == 2U) {
        auto ad = mul_expr(ctx, matrix(0U, 0U), matrix(1U, 1U));
        if (ad.is_error()) return ad;
        auto bc = mul_expr(ctx, matrix(0U, 1U), matrix(1U, 0U));
        if (bc.is_error()) return bc;
        return sub_expr(ctx, ad.value(), bc.value());
    }

    // Expand along first row.
    ExprPtr result = integer(ctx, 0);
    for (std::size_t col = 0; col < n; ++col) {
        ExprPtr entry = matrix(0U, col);
        if (is_zero_expr(entry)) continue;

        MatrixExpr minor(n - 1U, n - 1U);
        std::size_t out_col = 0U;
        for (std::size_t c = 0; c < n; ++c) {
            if (c == col) continue;
            for (std::size_t r = 1U; r < n; ++r) {
                minor(r - 1U, out_col) = matrix(r, c);
            }
            ++out_col;
        }

        auto minor_val = cofactor_det(minor, ctx);
        if (minor_val.is_error()) return minor_val;

        auto term = mul_expr(ctx, entry, minor_val.value());
        if (term.is_error()) return term;

        ExprPtr signed_term = term.value();
        if ((col % 2U) != 0U) {
            auto negated = ctx.simplify(ctx.arena().make<Unary>(UnaryOp::Neg, signed_term));
            if (negated.is_error()) return negated;
            signed_term = negated.value();
        }

        auto next = add_expr(ctx, result, signed_term);
        if (next.is_error()) return next;
        result = next.value();
    }
    return ok(result);
}

[[nodiscard]] Result<MatrixExpr> mat_scale_id(ExprPtr s, std::size_t n, symbolic::CASContext& ctx) {
    MatrixExpr res(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            res(i, j) = (i == j) ? s : integer(ctx, 0);
        }
    }
    return ok(res);
}

// Faddeev-Leverrier algorithm for characteristic polynomial.
// p(lambda) = det(lambda*I - A) = lambda^n + c_{n-1}*lambda^{n-1} + ... + c_0
// Note: our signature returns det(A - lambda*I) to match standard CAS conventions.
[[nodiscard]] Result<ExprPtr> faddeev_leverrier(const MatrixExpr& matrix, const Symbol& lambda_var, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n == 0U) return ok(integer(ctx, 1));

    std::vector<ExprPtr> coeffs(n + 1, ExprPtr{});
    coeffs[0] = integer(ctx, 1); // lambda^n coefficient

    MatrixExpr M = identity(n, ctx).value();
    for (std::size_t k = 1; k <= n; ++k) {
        auto AM = multiply(matrix, M, ctx);
        if (AM.is_error()) return fail<ExprPtr>(AM.error());
        
        auto tr = trace(AM.value(), ctx);
        if (tr.is_error()) return tr;
        
        auto neg_tr = ctx.simplify(ctx.arena().make<Unary>(UnaryOp::Neg, tr.value()));
        if (neg_tr.is_error()) return neg_tr;
        
        auto ck = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Div, neg_tr.value(), integer(ctx, k)));
        if (ck.is_error()) return ck;
        
        coeffs[k] = ck.value();
        
        if (k < n) {
            auto ckI = mat_scale_id(ck.value(), n, ctx);
            if (ckI.is_error()) return fail<ExprPtr>(ckI.error());
            auto next_M = add(AM.value(), ckI.value(), ctx);
            if (next_M.is_error()) return fail<ExprPtr>(next_M.error());
            M = std::move(next_M.value());
        }
    }
    // This is det(lambda*I - A).
    // To get det(A - lambda*I), we multiply by (-1)^n.
    ExprPtr lambda_expr = ctx.arena().make<Symbol>(lambda_var);
    std::vector<ExprPtr> terms;
    for (std::size_t k = 0; k <= n; ++k) {
        ExprPtr coeff = coeffs[k];
        if (is_zero_expr(coeff)) continue;
        
        ExprPtr pwr;
        if (n - k == 0) {
            pwr = coeff;
        } else if (n - k == 1) {
            if (const auto* il = expr_cast<IntegerLit>(coeff); il && il->value == BigInt(1)) {
                pwr = lambda_expr;
            } else {
                auto res = mul_expr(ctx, coeff, lambda_expr);
                if (res.is_error()) return res;
                pwr = res.value();
            }
        } else {
            auto exp = ctx.arena().make<Binary>(BinaryOp::Pow, lambda_expr, integer(ctx, n - k));
            if (const auto* il = expr_cast<IntegerLit>(coeff); il && il->value == BigInt(1)) {
                pwr = exp;
            } else {
                auto res = mul_expr(ctx, coeff, exp);
                if (res.is_error()) return res;
                pwr = res.value();
            }
        }
        terms.push_back(pwr);
    }

    auto poly = ctx.simplify(ctx.arena().make<Sum>(std::move(terms)));
    if (poly.is_error()) return poly;

    if (n % 2 != 0) {
        return ctx.simplify(ctx.arena().make<Unary>(UnaryOp::Neg, poly.value()));
    }
    return poly;
}

// Computes the null space basis of a matrix via RREF.
// Each basis vector is a column that spans ker(matrix).
[[nodiscard]] Result<std::vector<std::vector<ExprPtr>>> null_space(
    const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t n_cols = matrix.cols();

    auto reduced = rref(matrix, ctx);
    if (reduced.is_error()) {
        return fail<std::vector<std::vector<ExprPtr>>>(reduced.error());
    }

    // Identify pivot columns and their rows.
    std::vector<std::optional<std::size_t>> pivot_row_for_col(n_cols, std::nullopt);
    std::vector<bool> is_pivot_col(n_cols, false);
    
    std::size_t current_row = 0;
    for (std::size_t c = 0; c < n_cols && current_row < reduced.value().rows(); ++c) {
        if (!is_zero_expr(reduced.value()(current_row, c))) {
            is_pivot_col[c] = true;
            pivot_row_for_col[c] = current_row;
            ++current_row;
        }
    }

    // For each free column, build a null space basis vector.
    std::vector<std::vector<ExprPtr>> basis;
    for (std::size_t free_col = 0; free_col < n_cols; ++free_col) {
        if (is_pivot_col[free_col]) continue;

        std::vector<ExprPtr> vec(n_cols, integer(ctx, 0));
        vec[free_col] = integer(ctx, 1);

        for (std::size_t pivot_c = 0; pivot_c < n_cols; ++pivot_c) {
            if (!is_pivot_col[pivot_c] || !pivot_row_for_col[pivot_c].has_value()) {
                continue;
            }
            const std::size_t r = *pivot_row_for_col[pivot_c];
            ExprPtr coeff = reduced.value()(r, free_col);
            if (is_zero_expr(coeff)) continue;
            auto neg_coeff = ctx.simplify(ctx.arena().make<Unary>(UnaryOp::Neg, coeff));
            if (neg_coeff.is_error()) {
                return fail<std::vector<std::vector<ExprPtr>>>(neg_coeff.error());
            }
            vec[pivot_c] = neg_coeff.value();
        }

        basis.push_back(std::move(vec));
    }

    return ok(std::move(basis));
}

}  // namespace

Result<ExprPtr> characteristic_polynomial(
    const MatrixExpr& matrix, const Symbol& lambda_var, symbolic::CASContext& ctx) {
    if (matrix.rows() != matrix.cols()) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument,
            "Characteristic polynomial requires a square matrix"));
    }

    const std::size_t n = matrix.rows();
    if (n <= 3) {
        // For very small matrices, cofactor is fine and often produces cleaner terms
        ExprPtr lambda_expr = ctx.arena().make<Symbol>(lambda_var);
        MatrixExpr char_matrix(n, n);
        for (std::size_t r = 0; r < n; ++r) {
            for (std::size_t c = 0; c < n; ++c) {
                if (r == c) {
                    auto entry = sub_expr(ctx, matrix(r, c), lambda_expr);
                    if (entry.is_error()) return entry;
                    char_matrix(r, c) = entry.value();
                } else {
                    char_matrix(r, c) = matrix(r, c);
                }
            }
        }
        return cofactor_det(char_matrix, ctx);
    }

    return faddeev_leverrier(matrix, lambda_var, ctx);
}

Result<std::vector<ExprPtr>> eigenvalues(
    const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    if (matrix.rows() != matrix.cols()) {
        return fail<std::vector<ExprPtr>>(make_error(CASErrorKind::InvalidArgument,
            "Eigenvalues require a square matrix"));
    }

    const Symbol lambda_var("_lambda_");
    auto char_poly = characteristic_polynomial(matrix, lambda_var, ctx);
    if (char_poly.is_error()) {
        return fail<std::vector<ExprPtr>>(char_poly.error());
    }

    return algebra::solve_polynomial(char_poly.value(), lambda_var, ctx);
}

Result<std::vector<Eigenpair>> eigenvectors(
    const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    if (matrix.rows() != matrix.cols()) {
        return fail<std::vector<Eigenpair>>(make_error(CASErrorKind::InvalidArgument,
            "Eigenvectors require a square matrix"));
    }

    const std::size_t n = matrix.rows();
    auto eigenvals = eigenvalues(matrix, ctx);
    if (eigenvals.is_error()) {
        return fail<std::vector<Eigenpair>>(eigenvals.error());
    }

    // Group eigenvalues by structural equality to avoid redundant work
    std::vector<ExprPtr> unique_eigenvals;
    for (ExprPtr val : eigenvals.value()) {
        bool found = false;
        for (ExprPtr uval : unique_eigenvals) {
            if (structural_equal(val, uval)) {
                found = true;
                break;
            }
        }
        if (!found) {
            unique_eigenvals.push_back(val);
        }
    }

    std::vector<Eigenpair> result;
    for (ExprPtr eigenval : unique_eigenvals) {
        MatrixExpr shifted(n, n);
        bool build_ok = true;
        for (std::size_t r = 0; r < n && build_ok; ++r) {
            for (std::size_t c = 0; c < n && build_ok; ++c) {
                if (r == c) {
                    auto entry = sub_expr(ctx, matrix(r, c), eigenval);
                    if (entry.is_error()) { build_ok = false; break; }
                    shifted(r, c) = entry.value();
                } else {
                    shifted(r, c) = matrix(r, c);
                }
            }
        }

        if (build_ok) {
            auto kernel = null_space(shifted, ctx);
            if (kernel.is_ok()) {
                for (auto& vec : kernel.value()) {
                    result.push_back(Eigenpair{.eigenvalue = eigenval, .eigenvector = std::move(vec)});
                }
            }
        }
    }

    return ok(std::move(result));
}

}  // namespace cas::linalg
