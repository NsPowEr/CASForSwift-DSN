#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/symbolic.hpp"
#include "cas/normal_form.hpp"
#include "algebra_internal.hpp"

#include <type_traits>
#include <utility>
#include <vector>

namespace cas::symbolic {

Result<bool> mathematically_equal(ExprPtr lhs, ExprPtr rhs, CASContext& context) {
    const bool owns_operation = !context.operation_active_;
    if (owns_operation) {
        context.operation_active_ = true;
        context.trace_capture_active_ = false;
        context.trace_.clear();
    }

    auto finalize = [&]() {
        if (owns_operation) context.operation_active_ = false;
    };

    auto lhs_s = context.simplify(lhs);
    if (lhs_s.is_error()) { finalize(); return fail<bool>(lhs_s.error()); }

    auto rhs_s = context.simplify(rhs);
    if (rhs_s.is_error()) { finalize(); return fail<bool>(rhs_s.error()); }

    // Fast path: structural equality after simplify
    if (structural_equal(lhs_s.value(), rhs_s.value())) {
        finalize();
        return ok(true);
    }

    // Mathematical equality using polynomial normal form: normal_form(expand(lhs - rhs)) == 0
    auto diff_expr = context.arena().make<Binary>(BinaryOp::Sub, lhs_s.value(), rhs_s.value());
    auto normal_diff = polynomial_normal_form(diff_expr, context);
    if (normal_diff.is_ok()) {
        if (expr_is<IntegerLit>(normal_diff.value()) && expr_cast<IntegerLit>(normal_diff.value())->value.is_zero()) {
            finalize();
            return ok(true);
        }
    }

    // Rational equality fallback: check num_L * den_R - num_R * den_L == 0 using normal form
    auto lhs_parts = algebra::split_num_den(lhs_s.value(), context);
    auto rhs_parts = algebra::split_num_den(rhs_s.value(), context);
    if (lhs_parts.is_ok() && rhs_parts.is_ok()) {
        auto cross_l = algebra::multiply_exprs(lhs_parts.value().numerator, rhs_parts.value().denominator, context);
        auto cross_r = algebra::multiply_exprs(rhs_parts.value().numerator, lhs_parts.value().denominator, context);
        if (cross_l.is_ok() && cross_r.is_ok()) {
            auto cross_diff_expr = context.arena().make<Binary>(BinaryOp::Sub, cross_l.value(), cross_r.value());
            auto cross_diff = polynomial_normal_form(cross_diff_expr, context);
            if (cross_diff.is_ok() && expr_is<IntegerLit>(cross_diff.value()) && expr_cast<IntegerLit>(cross_diff.value())->value.is_zero()) {
                finalize();
                return ok(true);
            }
        }
    }

    finalize();
    return ok(false);
}

namespace {

// L2-19 helpers.  All transformations preserve mathematical equality and
// fire only when the relevant positivity / integrality assumptions are
// derivable from CASContext::assumptions().

[[nodiscard]] bool is_known_positive(ExprPtr expr, const Assumptions& a) {
    if (!expr) return false;
    if (const auto* lit = expr_cast<IntegerLit>(expr)) {
        return !lit->value.is_zero() && !lit->value.is_negative();
    }
    if (const auto* lit = expr_cast<RationalLit>(expr)) {
        const bool num_neg = lit->numerator.is_negative();
        const bool den_neg = lit->denominator.is_negative();
        if (lit->numerator.is_zero()) return false;
        return num_neg == den_neg;
    }
    return a.is_positive(expr);
}

[[nodiscard]] bool is_known_integer(ExprPtr expr, const Assumptions& a) {
    if (!expr) return false;
    if (expr_is<IntegerLit>(expr)) return true;
    return a.is_integer(expr);
}

ExprPtr expand_log_walker(ExprPtr expr, CASContext& ctx);
ExprPtr expand_exp_walker(ExprPtr expr, CASContext& ctx);

ExprPtr rebuild_with_children(
    ExprPtr expr,
    CASContext& ctx,
    ExprPtr (*walker)(ExprPtr, CASContext&)) {
    if (!expr) return expr;
    return visit_expr(expr, [&](const auto& node) -> ExprPtr {
        using Node = std::decay_t<decltype(node)>;
        AstArena& arena = ctx.arena();
        if constexpr (std::is_same_v<Node, Unary>) {
            return arena.make<Unary>(node.op, walker(node.operand, ctx));
        } else if constexpr (std::is_same_v<Node, Binary>) {
            return arena.make<Binary>(node.op, walker(node.left, ctx), walker(node.right, ctx));
        } else if constexpr (std::is_same_v<Node, FuncCall>) {
            std::vector<ExprPtr> args;
            args.reserve(node.args.size());
            for (ExprPtr a : node.args) args.push_back(walker(a, ctx));
            return arena.make<FuncCall>(node.func_id, std::move(args));
        } else if constexpr (std::is_same_v<Node, Sum>) {
            std::vector<ExprPtr> terms;
            terms.reserve(node.terms.size());
            for (ExprPtr t : node.terms) terms.push_back(walker(t, ctx));
            return arena.make<Sum>(std::move(terms));
        } else if constexpr (std::is_same_v<Node, Product>) {
            std::vector<ExprPtr> factors;
            factors.reserve(node.factors.size());
            for (ExprPtr f : node.factors) factors.push_back(walker(f, ctx));
            return arena.make<Product>(std::move(factors));
        } else {
            return expr;
        }
    });
}

[[nodiscard]] bool is_log_funccall(const FuncCall& fc) {
    return (fc.func_id == BuiltinOp::Ln || fc.func_id == BuiltinOp::Log) && fc.args.size() == 1U;
}

[[nodiscard]] bool is_exp_funccall(const FuncCall& fc) {
    return fc.func_id == BuiltinOp::Exp && fc.args.size() == 1U;
}

ExprPtr expand_log_walker(ExprPtr expr, CASContext& ctx) {
    if (!expr) return expr;
    ExprPtr transformed = rebuild_with_children(expr, ctx, expand_log_walker);

    const auto* fc = expr_cast<FuncCall>(transformed);
    if (!fc || !is_log_funccall(*fc)) return transformed;

    ExprPtr arg = fc->args[0];
    const Assumptions& a = ctx.assumptions();
    AstArena& arena = ctx.arena();

    if (const auto* prod = expr_cast<Product>(arg)) {
        bool all_positive = true;
        for (ExprPtr f : prod->factors) {
            if (!is_known_positive(f, a)) { all_positive = false; break; }
        }
        if (all_positive && prod->factors.size() >= 2U) {
            std::vector<ExprPtr> terms;
            terms.reserve(prod->factors.size());
            for (ExprPtr f : prod->factors) {
                terms.push_back(arena.make<FuncCall>(fc->func_id, std::vector<ExprPtr>{f}));
            }
            return arena.make<Sum>(std::move(terms));
        }
    }

    if (const auto* bin = expr_cast<Binary>(arg)) {
        if (bin->op == BinaryOp::Pow && is_known_positive(bin->left, a)) {
            ExprPtr inner_log = arena.make<FuncCall>(fc->func_id, std::vector<ExprPtr>{bin->left});
            return arena.make<Binary>(BinaryOp::Mul, bin->right, inner_log);
        }
        if (bin->op == BinaryOp::Div
            && is_known_positive(bin->left, a)
            && is_known_positive(bin->right, a)) {
            ExprPtr log_num = arena.make<FuncCall>(fc->func_id, std::vector<ExprPtr>{bin->left});
            ExprPtr log_den = arena.make<FuncCall>(fc->func_id, std::vector<ExprPtr>{bin->right});
            return arena.make<Binary>(BinaryOp::Sub, log_num, log_den);
        }
        if (bin->op == BinaryOp::Mul
            && is_known_positive(bin->left, a)
            && is_known_positive(bin->right, a)) {
            ExprPtr log_l = arena.make<FuncCall>(fc->func_id, std::vector<ExprPtr>{bin->left});
            ExprPtr log_r = arena.make<FuncCall>(fc->func_id, std::vector<ExprPtr>{bin->right});
            return arena.make<Binary>(BinaryOp::Add, log_l, log_r);
        }
    }
    return transformed;
}

// Detect a single factor of the form n * ln(x) inside an exponent and
// return the (n, ln_arg) pair on success.
[[nodiscard]] bool try_match_integer_times_log(
    ExprPtr expr,
    const Assumptions& a,
    ExprPtr& out_n,
    ExprPtr& out_log_arg) {
    if (const auto* bin = expr_cast<Binary>(expr); bin && bin->op == BinaryOp::Mul) {
        ExprPtr lhs = bin->left;
        ExprPtr rhs = bin->right;
        for (int trial = 0; trial < 2; ++trial) {
            if (const auto* fc = expr_cast<FuncCall>(rhs); fc && is_log_funccall(*fc)) {
                if (is_known_positive(fc->args[0], a) && is_known_integer(lhs, a)) {
                    out_n = lhs;
                    out_log_arg = fc->args[0];
                    return true;
                }
            }
            std::swap(lhs, rhs);
        }
    }
    if (const auto* prod = expr_cast<Product>(expr); prod && prod->factors.size() == 2U) {
        ExprPtr lhs = prod->factors[0];
        ExprPtr rhs = prod->factors[1];
        for (int trial = 0; trial < 2; ++trial) {
            if (const auto* fc = expr_cast<FuncCall>(rhs); fc && is_log_funccall(*fc)) {
                if (is_known_positive(fc->args[0], a) && is_known_integer(lhs, a)) {
                    out_n = lhs;
                    out_log_arg = fc->args[0];
                    return true;
                }
            }
            std::swap(lhs, rhs);
        }
    }
    return false;
}

ExprPtr expand_exp_walker(ExprPtr expr, CASContext& ctx) {
    if (!expr) return expr;
    ExprPtr transformed = rebuild_with_children(expr, ctx, expand_exp_walker);

    const auto* fc = expr_cast<FuncCall>(transformed);
    if (!fc || !is_exp_funccall(*fc)) return transformed;

    ExprPtr arg = fc->args[0];
    const Assumptions& a = ctx.assumptions();
    AstArena& arena = ctx.arena();

    // exp(sum t_i) -> prod exp(t_i)  (always valid; no positivity needed).
    if (const auto* sum = expr_cast<Sum>(arg); sum && sum->terms.size() >= 2U) {
        std::vector<ExprPtr> factors;
        factors.reserve(sum->terms.size());
        for (ExprPtr t : sum->terms) {
            factors.push_back(arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{t}));
        }
        return arena.make<Product>(std::move(factors));
    }

    // exp(ln(x)) -> x  when x > 0.
    if (const auto* inner = expr_cast<FuncCall>(arg); inner && is_log_funccall(*inner)) {
        if (is_known_positive(inner->args[0], a)) {
            return inner->args[0];
        }
    }

    // exp(n * ln(x)) -> x^n  when x > 0 and n integer (handles the common
    // case where the original expression was 2*ln(x) etc).
    ExprPtr n_factor{};
    ExprPtr log_arg{};
    if (try_match_integer_times_log(arg, a, n_factor, log_arg)) {
        return arena.make<Binary>(BinaryOp::Pow, log_arg, n_factor);
    }

    return transformed;
}

}  // namespace

Result<bool> mathematically_equal_subset_risch(ExprPtr lhs, ExprPtr rhs, CASContext& context) {
    if (!lhs || !rhs) {
        return fail<bool>(CASError{
            CASErrorKind::InvalidArgument,
            "mathematically_equal_subset_risch: null operand",
            std::nullopt});
    }

    // Apply log expansion first, then exp expansion.  The expansions
    // commute on the decidable subset (positive arguments / integer
    // exponents) so a single pass each suffices.
    ExprPtr lhs_log = expand_log_walker(lhs, context);
    ExprPtr rhs_log = expand_log_walker(rhs, context);
    ExprPtr lhs_norm = expand_exp_walker(lhs_log, context);
    ExprPtr rhs_norm = expand_exp_walker(rhs_log, context);

    return mathematically_equal(lhs_norm, rhs_norm, context);
}

} // namespace cas::symbolic
