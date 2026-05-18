#include "cas/linalg/Matrix.hpp"
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::linalg {
namespace {

[[nodiscard]] ExprPtr integer(symbolic::CASContext& ctx, long long value) {
    return ctx.arena().make<IntegerLit>(BigInt(value));
}

[[nodiscard]] Result<ExprPtr> simplify(symbolic::CASContext& ctx, ExprPtr expr) {
    return ctx.simplify(expr);
}

[[nodiscard]] Result<ExprPtr> sub_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Sub, lhs, rhs));
}

[[nodiscard]] Result<ExprPtr> mul_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Mul, lhs, rhs));
}

[[nodiscard]] Result<ExprPtr> div_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Div, lhs, rhs));
}

[[nodiscard]] bool is_zero_expr(ExprPtr expr) {
    if (const auto* integer_lit = expr_cast<IntegerLit>(expr)) {
        return integer_lit->value.is_zero();
    }
    if (const auto* rational_lit = expr_cast<RationalLit>(expr)) {
        return rational_lit->numerator.is_zero();
    }
    return false;
}

[[nodiscard]] std::size_t estimate_complexity(ExprPtr expr) {
    if (!expr) return 0;
    switch (expr->kind) {
    case ExprKind::IntegerLit:
    case ExprKind::RationalLit:
    case ExprKind::DecimalLit:
    case ExprKind::Constant:
        return 1;
    case ExprKind::Symbol:
        return 2;
    case ExprKind::Unary:
        return 1 + estimate_complexity(expr_cast<Unary>(expr)->operand);
    case ExprKind::Binary: {
        const auto* b = expr_cast<Binary>(expr);
        return 1 + estimate_complexity(b->left) + estimate_complexity(b->right);
    }
    case ExprKind::Sum: {
        std::size_t c = 1;
        for (auto t : expr_cast<Sum>(expr)->terms) c += estimate_complexity(t);
        return c;
    }
    case ExprKind::Product: {
        std::size_t c = 1;
        for (auto f : expr_cast<Product>(expr)->factors) c += estimate_complexity(f);
        return c;
    }
    case ExprKind::FuncCall: {
        std::size_t c = 1;
        for (auto a : expr_cast<FuncCall>(expr)->args) c += estimate_complexity(a);
        return c;
    }
    default:
        return 5;
    }
}

[[nodiscard]] bool is_structurally_nonzero(ExprPtr expr) {
    if (!expr) return false;
    if (const auto* i = expr_cast<IntegerLit>(expr)) return !i->value.is_zero();
    if (const auto* r = expr_cast<RationalLit>(expr)) return !r->numerator.is_zero();
    if (expr_is<Constant>(expr)) return true;
    if (expr_is<Symbol>(expr)) return true;
    return !is_zero_expr(expr);
}

// Worst-case total degree of `expr` viewed as a polynomial in all symbolic
// indeterminates. Used as a tiebreaker in pivot selection so the Bareiss
// formula  (a*d - b*c)/prev_pivot  has a numerator of minimal symbolic degree.
// Conservative (upper bound) for non-polynomial nodes.
[[nodiscard]] std::size_t total_degree(ExprPtr expr) {
    if (!expr) return 0;
    switch (expr->kind) {
    case ExprKind::IntegerLit:
    case ExprKind::RationalLit:
    case ExprKind::DecimalLit:
    case ExprKind::Constant:
        return 0;
    case ExprKind::Symbol:
        return 1;
    case ExprKind::Unary:
        return total_degree(expr_cast<Unary>(expr)->operand);
    case ExprKind::Binary: {
        const auto* b = expr_cast<Binary>(expr);
        switch (b->op) {
        case BinaryOp::Add:
        case BinaryOp::Sub:
            return std::max(total_degree(b->left), total_degree(b->right));
        case BinaryOp::Mul:
            return total_degree(b->left) + total_degree(b->right);
        case BinaryOp::Div:
            return total_degree(b->left);
        case BinaryOp::Pow: {
            if (const auto* exp_lit = expr_cast<IntegerLit>(b->right)) {
                if (exp_lit->value.is_negative()) return total_degree(b->left);
                auto d = exp_lit->value.to_double();
                if (d > 0 && d < 1e6) {
                    return total_degree(b->left) * static_cast<std::size_t>(d);
                }
            }
            return total_degree(b->left) + 1;
        }
        default:
            return std::max(total_degree(b->left), total_degree(b->right));
        }
    }
    case ExprKind::Sum: {
        std::size_t d = 0;
        for (auto t : expr_cast<Sum>(expr)->terms) d = std::max(d, total_degree(t));
        return d;
    }
    case ExprKind::Product: {
        std::size_t d = 0;
        for (auto f : expr_cast<Product>(expr)->factors) d += total_degree(f);
        return d;
    }
    case ExprKind::FuncCall: {
        // Transcendental wrappers contribute their argument's degree as an
        // upper-bound proxy: prefer simpler-argument pivots when possible.
        std::size_t d = 0;
        for (auto a : expr_cast<FuncCall>(expr)->args) d = std::max(d, total_degree(a));
        return d + 1;
    }
    default:
        return 1;
    }
}

// Lexicographic pivot quality, strictly orderable.
//
// Components (compared in order, larger is better):
//
//   1. `certainty`:    how confidently we can prove the candidate is nonzero,
//                      which is the *correctness* dimension of pivot choice:
//                        3 = literal numeric (Int/Rat) — provably nonzero
//                        2 = ctx.assumptions().is_positive — provably nonzero
//                        1 = ctx.assumptions().is_nonzero — provably nonzero
//                        0 = only structural nonzero (heuristic)
//                      RootOf elements demote one level since their nonzero
//                      status depends on the defining minimal polynomial
//                      (handled via `rootof_penalty`).
//
//   2. `neg_total_degree`:
//                      -total_degree(val).  Lower symbolic degree → smaller
//                      coefficient blow-up in the Bareiss recurrence.
//
//   3. `neg_complexity`:
//                      -estimate_complexity(val). Tiebreaker on AST size.
//
// All components are signed integers so a plain operator< on the tuple yields
// the correct "higher is better" ordering when compared with `<`. Pivot
// selection picks the maximum.
struct PivotScore {
    int certainty;
    int neg_total_degree;
    int neg_complexity;

    [[nodiscard]] auto as_tuple() const noexcept {
        return std::tie(certainty, neg_total_degree, neg_complexity);
    }
    [[nodiscard]] bool operator<(const PivotScore& other) const noexcept {
        return as_tuple() < other.as_tuple();
    }
    [[nodiscard]] bool operator==(const PivotScore& other) const noexcept {
        return as_tuple() == other.as_tuple();
    }
};

[[nodiscard]] PivotScore make_pivot_score(ExprPtr val, symbolic::CASContext& ctx) {
    int certainty;
    if (expr_is<IntegerLit>(val) || expr_is<RationalLit>(val)) {
        certainty = 3;
    } else if (ctx.assumptions().is_positive(val)) {
        certainty = 2;
    } else if (ctx.assumptions().is_nonzero(val)) {
        certainty = 1;
    } else {
        certainty = 0;  // structurally nonzero but not certified
    }
    // RootOf: nonzero status depends on minimal polynomial — demote one level.
    if (expr_is<RootOf>(val) && certainty > 0) {
        certainty -= 1;
    }
    const std::size_t deg = total_degree(val);
    const std::size_t cpx = estimate_complexity(val);
    // Clamp to int range to avoid overflow in negation.
    const int deg_i = static_cast<int>(std::min<std::size_t>(deg, 1'000'000));
    const int cpx_i = static_cast<int>(std::min<std::size_t>(cpx, 1'000'000));
    return PivotScore{certainty, -deg_i, -cpx_i};
}

void swap_rows(MatrixExpr& matrix, std::size_t lhs, std::size_t rhs) {
    if (lhs == rhs) return;
    for (std::size_t col = 0; col < matrix.cols(); ++col) {
        std::swap(matrix(lhs, col), matrix(rhs, col));
    }
}

}  // namespace

Result<MatrixExpr> bareiss(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    MatrixExpr result(matrix.rows(), matrix.cols(), matrix.elements());
    std::size_t r = 0;
    ExprPtr prev_pivot = integer(ctx, 1);

    for (std::size_t c = 0; c < result.cols() && r < result.rows(); ++c) {
        std::size_t pivot_row = result.rows();
        std::optional<PivotScore> best_score;

        for (std::size_t i = r; i < result.rows(); ++i) {
            ExprPtr val = result(i, c);
            if (is_zero_expr(val)) continue;
            if (!is_structurally_nonzero(val)) continue;

            PivotScore score = make_pivot_score(val, ctx);
            if (!best_score.has_value() || *best_score < score) {
                best_score = score;
                pivot_row = i;
            }
            // Early exit: a literal-nonzero, degree-0 pivot dominates any
            // other candidate this column can offer. (Higher certainty
            // would be impossible; lower degree would be impossible too.)
            if (score.certainty == 3 && score.neg_total_degree == 0) break;
        }

        if (pivot_row == result.rows()) {
            continue;
        }

        if (pivot_row != r) {
            swap_rows(result, r, pivot_row);
        }

        ExprPtr current_pivot = result(r, c);

        for (std::size_t i = r + 1; i < result.rows(); ++i) {
            for (std::size_t j = c + 1; j < result.cols(); ++j) {
                auto left_res = mul_expr(ctx, current_pivot, result(i, j));
                if (left_res.is_error()) return fail<MatrixExpr>(left_res.error());
                
                auto right_res = mul_expr(ctx, result(i, c), result(r, j));
                if (right_res.is_error()) return fail<MatrixExpr>(right_res.error());
                
                auto num_res = sub_expr(ctx, left_res.value(), right_res.value());
                if (num_res.is_error()) return fail<MatrixExpr>(num_res.error());
                
                auto val_res = div_expr(ctx, num_res.value(), prev_pivot);
                if (val_res.is_error()) return fail<MatrixExpr>(val_res.error());
                
                result(i, j) = val_res.value();
            }
            result(i, c) = integer(ctx, 0);
        }

        prev_pivot = current_pivot;
        ++r;
    }

    return ok(std::move(result));
}

}  // namespace cas::linalg
