#include "cas/linalg/matrix_expr_helpers.hpp"
#include <algorithm>

namespace cas::linalg {

ExprPtr integer(symbolic::CASContext& ctx, long long value) {
    return ctx.arena().make<IntegerLit>(BigInt(value));
}

ExprPtr integer(symbolic::CASContext& ctx, const BigInt& value) {
    return ctx.arena().make<IntegerLit>(value);
}

Result<ExprPtr> simplify(symbolic::CASContext& ctx, ExprPtr expr) {
    return ctx.simplify(expr);
}

Result<ExprPtr> add_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    if (is_zero_expr(lhs)) return ok(rhs);
    if (is_zero_expr(rhs)) return ok(lhs);
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Add, lhs, rhs));
}

Result<ExprPtr> negate_expr(symbolic::CASContext& ctx, ExprPtr expr) {
    if (is_zero_expr(expr)) return ok(expr);
    return simplify(ctx, ctx.arena().make<Unary>(UnaryOp::Neg, expr));
}

Result<ExprPtr> sub_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    if (is_zero_expr(rhs)) return ok(lhs);
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Sub, lhs, rhs));
}

Result<ExprPtr> mul_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    if (is_zero_expr(lhs) || is_zero_expr(rhs)) return ok(integer(ctx, 0));
    if (is_one_expr(lhs)) return ok(rhs);
    if (is_one_expr(rhs)) return ok(lhs);
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Mul, lhs, rhs));
}

Result<ExprPtr> div_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    if (is_zero_expr(lhs)) return ok(integer(ctx, 0));
    if (is_one_expr(rhs)) return ok(lhs);
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Div, lhs, rhs));
}

bool is_zero_expr(ExprPtr expr) {
    if (!expr) return false;
    if (const auto* il = expr_cast<IntegerLit>(expr)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(expr)) return rl->numerator.is_zero();
    return false;
}

bool is_one_expr(ExprPtr expr) {
    if (!expr) return false;
    if (const auto* il = expr_cast<IntegerLit>(expr)) return il->value == BigInt(1);
    if (const auto* rl = expr_cast<RationalLit>(expr)) return rl->numerator == BigInt(1) && rl->denominator == BigInt(1);
    return false;
}

std::optional<BigInt> try_get_bigint(ExprPtr expr) {
    if (const auto* il = expr_cast<IntegerLit>(expr)) return il->value;
    return std::nullopt;
}

std::size_t estimate_complexity(ExprPtr expr) {
    if (!expr) return 0;
    switch (expr->kind) {
    case ExprKind::IntegerLit: {
        const auto* il = expr_cast<IntegerLit>(expr);
        if (il->value == BigInt(1) || il->value == BigInt(-1)) return 1;
        return 2;
    }
    case ExprKind::RationalLit: {
        const auto* rl = expr_cast<RationalLit>(expr);
        if (rl->numerator == BigInt(1) && rl->denominator == BigInt(1)) return 1;
        if (rl->numerator == BigInt(-1) && rl->denominator == BigInt(1)) return 1;
        return 2;
    }
    case ExprKind::DecimalLit:
    case ExprKind::Constant:
        return 2;
    case ExprKind::Symbol:
        return 4;
    case ExprKind::Unary:
        return 1 + estimate_complexity(expr_cast<Unary>(expr)->operand);
    case ExprKind::Binary: {
        const auto* b = expr_cast<Binary>(expr);
        switch (b->op) {
        case BinaryOp::Mul:
            return 2 + estimate_complexity(b->left) + estimate_complexity(b->right);
        case BinaryOp::Div:
            return 3 + estimate_complexity(b->left) + estimate_complexity(b->right);
        case BinaryOp::Add:
        case BinaryOp::Sub:
            return 10 + 2 * (estimate_complexity(b->left) + estimate_complexity(b->right));
        case BinaryOp::Pow:
            return 4 + estimate_complexity(b->left) + estimate_complexity(b->right);
        default:
            return 5 + estimate_complexity(b->left) + estimate_complexity(b->right);
        }
    }
    case ExprKind::Sum: {
        std::size_t c = 10;
        for (auto t : expr_cast<Sum>(expr)->terms) c += 2 * estimate_complexity(t);
        return c;
    }
    case ExprKind::Product: {
        std::size_t c = 2;
        for (auto f : expr_cast<Product>(expr)->factors) c += estimate_complexity(f);
        return c;
    }
    case ExprKind::FuncCall: {
        std::size_t c = 5;
        for (auto a : expr_cast<FuncCall>(expr)->args) c += estimate_complexity(a);
        return c;
    }
    default:
        return 10;
    }
}

bool is_structurally_nonzero(ExprPtr expr) {
    if (!expr) return false;
    if (const auto* i = expr_cast<IntegerLit>(expr)) return !i->value.is_zero();
    if (const auto* r = expr_cast<RationalLit>(expr)) return !r->numerator.is_zero();
    if (expr_is<Constant>(expr)) return true;
    if (expr_is<Symbol>(expr)) return true;
    return !is_zero_expr(expr);
}

bool is_known_nonzero(ExprPtr expr, symbolic::CASContext& ctx) {
    if (is_structurally_nonzero(expr)) return true;
    return ctx.assumptions().is_nonzero(expr) || ctx.assumptions().is_positive(expr) || ctx.assumptions().is_negative(expr);
}

std::size_t total_degree(ExprPtr expr) {
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
        std::size_t d = 0;
        for (auto a : expr_cast<FuncCall>(expr)->args) d = std::max(d, total_degree(a));
        return d + 1;
    }
    default:
        return 1;
    }
}

Result<ExprPtr> sym_norm_sq(const std::vector<ExprPtr>& v, symbolic::CASContext& ctx) {
    ExprPtr sq = integer(ctx, 0);
    for (auto& vi : v) {
        auto prod = mul_expr(ctx, vi, vi);
        if (prod.is_error()) return fail<ExprPtr>(prod.error());
        auto sum = add_expr(ctx, sq, prod.value());
        if (sum.is_error()) return fail<ExprPtr>(sum.error());
        sq = sum.value();
    }
    return ok(sq);
}

Result<ExprPtr> sym_norm(const std::vector<ExprPtr>& v, symbolic::CASContext& ctx) {
    auto sq_res = sym_norm_sq(v, ctx);
    if (sq_res.is_error()) return sq_res;
    // Do NOT simplify(sqrt(...)) here: on a rational sum-of-squares with large
    // numerator/denominator (typical for random rational inputs of size n≥6),
    // simplify attempts integer factorization to detect perfect squares and
    // hits exponential cost. We leave sqrt as an unevaluated algebraic block;
    // downstream mul/sub/expand handle it lazily via structural sharing.
    return ok(ctx.arena().make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{sq_res.value()}));
}

PivotScore make_pivot_score(ExprPtr val, symbolic::CASContext& ctx) {
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
    const int deg_i = static_cast<int>(std::min<std::size_t>(deg, 1'000'000));
    const int cpx_i = static_cast<int>(std::min<std::size_t>(cpx, 1'000'000));
    return PivotScore{certainty, -cpx_i, -deg_i};
}

} // namespace cas::linalg
