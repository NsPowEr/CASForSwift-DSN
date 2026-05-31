#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"

#include <optional>
#include <string>
#include <utility>

namespace cas::linalg {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

void swap_rows(MatrixExpr& matrix, MatrixExpr& identity, std::size_t r1, std::size_t r2) {
    if (r1 == r2) return;
    for (std::size_t c = 0; c < matrix.cols(); ++c) std::swap(matrix(r1, c), matrix(r2, c));
    for (std::size_t c = 0; c < identity.cols(); ++c) std::swap(identity(r1, c), identity(r2, c));
}

[[nodiscard]] Result<void> scale_row(MatrixExpr& matrix, MatrixExpr& identity, std::size_t row, ExprPtr divisor, symbolic::CASContext& ctx) {
    if (is_one_expr(divisor)) return ok();
    if (!is_known_nonzero(divisor, ctx)) {
        return make_unimplemented<void>(
            "linalg", "inverse",
            "pivot expression cannot be decided nonzero",
            error::reason_codes::LINALG_RREF_UNDECIDABLE_PIVOT,
            "Add assumptions about pivot variables", "L3-13");
    }
    for (std::size_t c = 0; c < matrix.cols(); ++c) {
        auto res = div_expr(ctx, matrix(row, c), divisor);
        if (res.is_error()) return fail<void>(res.error());
        matrix(row, c) = res.value();
    }
    for (std::size_t c = 0; c < identity.cols(); ++c) {
        auto res = div_expr(ctx, identity(row, c), divisor);
        if (res.is_error()) return fail<void>(res.error());
        identity(row, c) = res.value();
    }
    return ok();
}

[[nodiscard]] Result<void> eliminate_row(MatrixExpr& matrix, MatrixExpr& identity, std::size_t target, std::size_t pivot_row, std::size_t col, symbolic::CASContext& ctx) {
    ExprPtr factor = matrix(target, col);
    if (is_zero_expr(factor)) return ok();

    for (std::size_t c = 0; c < matrix.cols(); ++c) {
        auto prod = mul_expr(ctx, factor, matrix(pivot_row, c));
        if (prod.is_error()) return fail<void>(prod.error());
        auto diff = sub_expr(ctx, matrix(target, c), prod.value());
        if (diff.is_error()) return fail<void>(diff.error());
        auto simp = simplify(ctx, diff.value());
        matrix(target, c) = simp.is_ok() ? simp.value() : diff.value();
    }
    for (std::size_t c = 0; c < identity.cols(); ++c) {
        auto prod = mul_expr(ctx, factor, identity(pivot_row, c));
        if (prod.is_error()) return fail<void>(prod.error());
        auto diff = sub_expr(ctx, identity(target, c), prod.value());
        if (diff.is_error()) return fail<void>(diff.error());
        auto simp = simplify(ctx, diff.value());
        identity(target, c) = simp.is_ok() ? simp.value() : diff.value();
    }
    return ok();
}

} // namespace

Result<MatrixExpr> inverse(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n != matrix.cols()) {
        return fail<MatrixExpr>(make_error(CASErrorKind::InvalidArgument, "Inverse requires a square matrix"));
    }
    if (n == 0) {
        return fail<MatrixExpr>(make_error(CASErrorKind::InvalidArgument, "Inverse of empty matrix"));
    }

    if (n == 2) {
        auto det_res = determinant(matrix, ctx);
        if (det_res.is_error()) return fail<MatrixExpr>(det_res.error());
        ExprPtr det = det_res.value();
        if (is_zero_expr(det)) return fail<MatrixExpr>(make_error(CASErrorKind::Undefined, "inverse: matrix is singular"));

        MatrixExpr res(2, 2);
        auto d_res = div_expr(ctx, matrix(1, 1), det);
        auto mb_r = negate_expr(ctx, matrix(0, 1));
        auto mc_r = negate_expr(ctx, matrix(1, 0));
        if (mb_r.is_error() || mc_r.is_error()) return fail<MatrixExpr>(make_error(CASErrorKind::InternalError, "negate failed"));
        
        auto mb_res = div_expr(ctx, mb_r.value(), det);
        auto mc_res = div_expr(ctx, mc_r.value(), det);
        auto a_res = div_expr(ctx, matrix(0, 0), det);
        
        if (d_res.is_error() || mb_res.is_error() || mc_res.is_error() || a_res.is_error())
            return fail<MatrixExpr>(make_error(CASErrorKind::InternalError, "2x2 inverse failed"));
            
        res(0, 0) = d_res.value();
        res(0, 1) = mb_res.value();
        res(1, 0) = mc_res.value();
        res(1, 1) = a_res.value();
        
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                auto s = simplify(ctx, res(i, j));
                if (s.is_ok()) res(i, j) = s.value();
            }
        }
        return ok(std::move(res));
    }

    MatrixExpr a(n, n, matrix.elements());
    auto id_res = identity(n, ctx);
    if (id_res.is_error()) return id_res;
    MatrixExpr inv = std::move(id_res.value());

    for (std::size_t i = 0; i < n; ++i) {
        // Find pivot
        std::size_t sel = n;
        PivotScore best{-1, 0, 0};
        for (std::size_t r = i; r < n; ++r) {
            if (is_zero_expr(a(r, i))) continue;
            PivotScore s = make_pivot_score(a(r, i), ctx);
            if (sel == n || s > best) {
                best = s;
                sel = r;
                if (s.certainty == 3) break;
            }
        }

        if (sel == n || best.certainty < 0) {
            return fail<MatrixExpr>(make_error(CASErrorKind::Undefined, "Matrix is singular (or pivot cannot be decided nonzero)"));
        }

        swap_rows(a, inv, i, sel);

        auto scale_res = scale_row(a, inv, i, a(i, i), ctx);
        if (scale_res.is_error()) return fail<MatrixExpr>(scale_res.error());

        for (std::size_t r = 0; r < n; ++r) {
            if (r != i) {
                auto elim_res = eliminate_row(a, inv, r, i, i, ctx);
                if (elim_res.is_error()) return fail<MatrixExpr>(elim_res.error());
            }
        }
    }

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            auto s = simplify(ctx, inv(i, j));
            if (s.is_ok()) {
                auto exp = algebra::expand(s.value(), ctx);
                auto val = exp.is_ok() ? exp.value() : s.value();
                auto tg = algebra::together(val, ctx);
                inv(i, j) = tg.is_ok() ? tg.value() : val;
            }
        }
    }
    
    // DEBUG
    // if (n == 2) std::cout << "INV(0,0): " << formatter::format(inv(0,0)) << std::endl;

    return ok(std::move(inv));
}

} // namespace cas::linalg
