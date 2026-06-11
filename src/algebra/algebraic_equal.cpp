#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/symbolic.hpp"
#include "cas/normal_form.hpp"
#include "algebra_internal.hpp"

#include <cstddef>
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

    // F7.5.A4: rewrite sech/csch/coth/tanh to canonical cosh/sinh quotients
    // BEFORE structural / algebraic comparison so notational mismatches
    // (CAS cosh(x)^-2 vs Maxima sech(x)^2 etc.) collapse.
    {
        auto lhs_h = algebra::hyperbolic_normalize(lhs_s.value(), context.arena());
        auto rhs_h = algebra::hyperbolic_normalize(rhs_s.value(), context.arena());
        if (lhs_h.get() != lhs_s.value().get() || rhs_h.get() != rhs_s.value().get()) {
            auto lhs_h_s = context.simplify(lhs_h);
            auto rhs_h_s = context.simplify(rhs_h);
            if (lhs_h_s.is_ok()) lhs_s = lhs_h_s;
            if (rhs_h_s.is_ok()) rhs_s = rhs_h_s;
        }
    }

    if (structural_equal(lhs_s.value(), rhs_s.value())) {
        finalize();
        return ok(true);
    }

    // F7.5.A1 / HC-F75-CYCLOTOMIC-ROOTOF: RootOf-specific decisions.
    // See src/algebra/algebraic_equal_cyclotomic.cpp::try_rootof_decision.
    if (auto rootof_dec = algebra::try_rootof_decision(
            lhs_s.value(), rhs_s.value(), context);
        rootof_dec.has_value()) {
        finalize();
        return ok(*rootof_dec);
    }

    // F1.6 bridge: Constant::I ↔ ComplexLit(0,1) (and other purely-imaginary
    // ComplexLit forms) live on different AST kinds but represent the same
    // algebraic value.  A direct simplify of (lhs - rhs) drives both forms
    // through the same Sum coefficient pool (ComplexRational accumulator) and
    // collapses to IntegerLit(0) when equal.
    {
        auto diff_simple = context.simplify(
            context.arena().make<Binary>(BinaryOp::Sub, lhs_s.value(), rhs_s.value()));
        if (diff_simple.is_ok()) {
            if (const auto* il = expr_cast<IntegerLit>(diff_simple.value());
                il != nullptr && il->value.is_zero()) {
                finalize();
                return ok(true);
            }
            if (const auto* rl = expr_cast<RationalLit>(diff_simple.value());
                rl != nullptr && rl->numerator.is_zero()) {
                finalize();
                return ok(true);
            }
            if (const auto* cl = expr_cast<ComplexLit>(diff_simple.value());
                cl != nullptr && cl->re_num.is_zero() && cl->im_num.is_zero()) {
                finalize();
                return ok(true);
            }
        }
    }

    auto diff_expr = context.arena().make<Binary>(BinaryOp::Sub, lhs_s.value(), rhs_s.value());
    auto normal_diff = polynomial_normal_form(diff_expr, context);
    if (normal_diff.is_ok()) {
        if (expr_is<IntegerLit>(normal_diff.value()) && expr_cast<IntegerLit>(normal_diff.value())->value.is_zero()) {
            finalize();
            return ok(true);
        }
    }

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

// L2-19 subset Risch normalisation helpers.
//
// All transformations preserve mathematical equality and fire only when
// the relevant positivity / numerical assumptions are derivable from
// CASContext::assumptions() (extended with structural inference for
// `exp`, `cosh`, `sqrt`, even-power sums and constants).

// ---- Positivity oracle -------------------------------------------------

[[nodiscard]] bool is_positive_constant_literal(ExprPtr expr) {
    if (const auto* lit = expr_cast<IntegerLit>(expr)) {
        return !lit->value.is_zero() && !lit->value.is_negative();
    }
    if (const auto* lit = expr_cast<RationalLit>(expr)) {
        if (lit->numerator.is_zero()) return false;
        return lit->numerator.is_negative() == lit->denominator.is_negative();
    }
    return false;
}

[[nodiscard]] bool infer_positive(ExprPtr expr, const Assumptions& a);
[[nodiscard]] bool infer_zero(ExprPtr expr, const Assumptions& a);

[[nodiscard]] bool infer_nonnegative(ExprPtr expr, const Assumptions& a) {
    if (!expr) return false;
    if (infer_positive(expr, a)) return true;
    if (const auto* lit = expr_cast<IntegerLit>(expr)) return !lit->value.is_negative();
    if (const auto* lit = expr_cast<RationalLit>(expr)) {
        return lit->numerator.is_zero()
            || lit->numerator.is_negative() == lit->denominator.is_negative();
    }
    if (a.is_nonnegative(expr)) return true;
    if (const auto* bin = expr_cast<Binary>(expr); bin && bin->op == BinaryOp::Pow) {
        if (const auto* exp = expr_cast<IntegerLit>(bin->right)) {
            if (!exp->value.is_negative() && exp->value.to_u64() % 2U == 0U) return true;
        }
    }
    return false;
}

[[nodiscard]] bool infer_positive(ExprPtr expr, const Assumptions& a) {
    if (!expr) return false;
    if (is_positive_constant_literal(expr)) return true;
    if (a.is_positive(expr)) return true;

    if (const auto* fc = expr_cast<FuncCall>(expr)) {
        switch (fc->func_id) {
            case BuiltinOp::Exp:
                if (fc->args.size() == 1U && a.is_real(fc->args[0])) return true;
                return false;
            case BuiltinOp::Cosh:
                if (fc->args.size() == 1U && a.is_real(fc->args[0])) return true;
                return false;
            case BuiltinOp::Sqrt:
                if (fc->args.size() == 1U && infer_positive(fc->args[0], a)) return true;
                return false;
            default:
                break;
        }
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr f : prod->factors) if (!infer_positive(f, a)) return false;
        return !prod->factors.empty();
    }
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Mul) {
            return infer_positive(bin->left, a) && infer_positive(bin->right, a);
        }
        if (bin->op == BinaryOp::Div) {
            return infer_positive(bin->left, a) && infer_positive(bin->right, a);
        }
        if (bin->op == BinaryOp::Pow) {
            if (infer_positive(bin->left, a)) return true;
            if (const auto* exp = expr_cast<IntegerLit>(bin->right)) {
                if (!exp->value.is_negative() && exp->value.to_u64() % 2U == 0U) {
                    return !infer_zero(bin->left, a);
                }
            }
        }
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        // Sum of strictly positives is strictly positive.  Sum of
        // non-negatives with at least one strictly positive is positive.
        bool any_strict = false;
        for (ExprPtr t : sum->terms) {
            if (infer_positive(t, a)) { any_strict = true; continue; }
            if (!infer_nonnegative(t, a)) return false;
        }
        return any_strict;
    }
    return false;
}

[[nodiscard]] bool infer_zero(ExprPtr expr, const Assumptions& a) {
    if (!expr) return false;
    if (const auto* lit = expr_cast<IntegerLit>(expr)) return lit->value.is_zero();
    if (const auto* lit = expr_cast<RationalLit>(expr)) return lit->numerator.is_zero();
    (void)a;
    return false;
}

[[nodiscard]] bool is_known_positive(ExprPtr expr, const Assumptions& a) {
    return infer_positive(expr, a);
}

// ---- Identity-preserving recursion (REGOLA 2: Structural Sharing) ------

using WalkerFn = ExprPtr (*)(ExprPtr, CASContext&);

[[nodiscard]] ExprPtr rebuild_with_children(ExprPtr expr, CASContext& ctx, WalkerFn walker) {
    if (!expr) return expr;
    return visit_expr(expr, [&](const auto& node) -> ExprPtr {
        using Node = std::decay_t<decltype(node)>;
        AstArena& arena = ctx.arena();
        if constexpr (std::is_same_v<Node, Unary>) {
            ExprPtr child = walker(node.operand, ctx);
            if (child == node.operand) return expr;
            return arena.make<Unary>(node.op, child);
        } else if constexpr (std::is_same_v<Node, Binary>) {
            ExprPtr l = walker(node.left, ctx);
            ExprPtr r = walker(node.right, ctx);
            if (l == node.left && r == node.right) return expr;
            return arena.make<Binary>(node.op, l, r);
        } else if constexpr (std::is_same_v<Node, FuncCall>) {
            std::vector<ExprPtr> args;
            args.reserve(node.args.size());
            bool changed = false;
            for (ExprPtr a : node.args) {
                ExprPtr a2 = walker(a, ctx);
                if (a2 != a) changed = true;
                args.push_back(a2);
            }
            if (!changed) return expr;
            return arena.make<FuncCall>(node.func_id, std::move(args));
        } else if constexpr (std::is_same_v<Node, Sum>) {
            std::vector<ExprPtr> terms;
            terms.reserve(node.terms.size());
            bool changed = false;
            for (ExprPtr t : node.terms) {
                ExprPtr t2 = walker(t, ctx);
                if (t2 != t) changed = true;
                terms.push_back(t2);
            }
            if (!changed) return expr;
            return arena.make<Sum>(std::move(terms));
        } else if constexpr (std::is_same_v<Node, Product>) {
            std::vector<ExprPtr> factors;
            factors.reserve(node.factors.size());
            bool changed = false;
            for (ExprPtr f : node.factors) {
                ExprPtr f2 = walker(f, ctx);
                if (f2 != f) changed = true;
                factors.push_back(f2);
            }
            if (!changed) return expr;
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

ExprPtr expand_log_walker(ExprPtr expr, CASContext& ctx);
ExprPtr expand_exp_walker(ExprPtr expr, CASContext& ctx);

// One-shot local rewrite for a log(arg).  Returns nullptr if no rewrite
// applies.  The expand_log_walker calls this AND then re-applies itself
// to the rewritten subtree to honour the recursive contract (B2 fix).
[[nodiscard]] ExprPtr try_local_log_rewrite(
    const FuncCall& fc,
    ExprPtr arg,
    CASContext& ctx) {
    const Assumptions& a = ctx.assumptions();
    AstArena& arena = ctx.arena();

    if (const auto* prod = expr_cast<Product>(arg)) {
        bool all_positive = !prod->factors.empty();
        for (ExprPtr f : prod->factors) {
            if (!is_known_positive(f, a)) { all_positive = false; break; }
        }
        if (all_positive && prod->factors.size() >= 2U) {
            std::vector<ExprPtr> terms;
            terms.reserve(prod->factors.size());
            for (ExprPtr f : prod->factors) {
                terms.push_back(arena.make<FuncCall>(fc.func_id, std::vector<ExprPtr>{f}));
            }
            return arena.make<Sum>(std::move(terms));
        }
    }

    if (const auto* bin = expr_cast<Binary>(arg)) {
        if (bin->op == BinaryOp::Pow && is_known_positive(bin->left, a)) {
            ExprPtr inner_log = arena.make<FuncCall>(fc.func_id, std::vector<ExprPtr>{bin->left});
            return arena.make<Binary>(BinaryOp::Mul, bin->right, inner_log);
        }
        if (bin->op == BinaryOp::Div
            && is_known_positive(bin->left, a)
            && is_known_positive(bin->right, a)) {
            ExprPtr log_num = arena.make<FuncCall>(fc.func_id, std::vector<ExprPtr>{bin->left});
            ExprPtr log_den = arena.make<FuncCall>(fc.func_id, std::vector<ExprPtr>{bin->right});
            return arena.make<Binary>(BinaryOp::Sub, log_num, log_den);
        }
        if (bin->op == BinaryOp::Mul
            && is_known_positive(bin->left, a)
            && is_known_positive(bin->right, a)) {
            ExprPtr log_l = arena.make<FuncCall>(fc.func_id, std::vector<ExprPtr>{bin->left});
            ExprPtr log_r = arena.make<FuncCall>(fc.func_id, std::vector<ExprPtr>{bin->right});
            return arena.make<Binary>(BinaryOp::Add, log_l, log_r);
        }
    }
    return ExprPtr{};
}

ExprPtr expand_log_walker(ExprPtr expr, CASContext& ctx) {
    if (!expr) return expr;
    ExprPtr transformed = rebuild_with_children(expr, ctx, expand_log_walker);

    const auto* fc = expr_cast<FuncCall>(transformed);
    if (!fc || !is_log_funccall(*fc)) return transformed;

    ExprPtr local = try_local_log_rewrite(*fc, fc->args[0], ctx);
    if (!local) return transformed;
    // Re-apply walker to the rewritten subtree so newly-produced log(...)
    // subterms (e.g. log(x^a) from log(x^a · y^b)) collapse further.
    return expand_log_walker(local, ctx);
}

// Walk a Mul/Product node trying to factor it as  scalar * ln(arg)  where
// scalar is anything (no integer constraint — for x>0 the identity
// exp(c * ln(x)) = x^c holds for any real c on the principal branch).
[[nodiscard]] bool try_match_scalar_times_log(
    ExprPtr expr,
    const Assumptions& a,
    AstArena& arena,
    ExprPtr& out_scalar,
    ExprPtr& out_log_arg) {
    std::vector<ExprPtr> factors;
    if (const auto* prod = expr_cast<Product>(expr)) {
        factors = prod->factors;
    } else if (const auto* bin = expr_cast<Binary>(expr); bin && bin->op == BinaryOp::Mul) {
        factors = {bin->left, bin->right};
    } else {
        return false;
    }
    if (factors.size() < 2U) return false;

    int log_index = -1;
    for (std::size_t i = 0; i < factors.size(); ++i) {
        if (const auto* fc = expr_cast<FuncCall>(factors[i]); fc && is_log_funccall(*fc)) {
            if (is_known_positive(fc->args[0], a)) {
                if (log_index != -1) return false;  // ambiguous: > 1 log factor
                log_index = static_cast<int>(i);
            }
        }
    }
    if (log_index == -1) return false;

    const auto& log_fc = expr_ref<FuncCall>(factors[static_cast<std::size_t>(log_index)]);
    out_log_arg = log_fc.args[0];

    std::vector<ExprPtr> rest;
    rest.reserve(factors.size() - 1U);
    for (std::size_t i = 0; i < factors.size(); ++i) {
        if (static_cast<int>(i) == log_index) continue;
        rest.push_back(factors[i]);
    }
    if (rest.size() == 1U) {
        out_scalar = rest.front();
    } else {
        out_scalar = arena.make<Product>(std::move(rest));
    }
    return true;
}

[[nodiscard]] ExprPtr try_local_exp_rewrite(
    const FuncCall& fc,
    ExprPtr arg,
    CASContext& ctx) {
    const Assumptions& a = ctx.assumptions();
    AstArena& arena = ctx.arena();

    if (const auto* sum = expr_cast<Sum>(arg); sum && sum->terms.size() >= 2U) {
        std::vector<ExprPtr> factors;
        factors.reserve(sum->terms.size());
        for (ExprPtr t : sum->terms) {
            factors.push_back(arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{t}));
        }
        return arena.make<Product>(std::move(factors));
    }

    if (const auto* inner = expr_cast<FuncCall>(arg); inner && is_log_funccall(*inner)) {
        if (is_known_positive(inner->args[0], a)) {
            return inner->args[0];
        }
    }

    ExprPtr scalar{};
    ExprPtr log_arg{};
    if (try_match_scalar_times_log(arg, a, arena, scalar, log_arg)) {
        return arena.make<Binary>(BinaryOp::Pow, log_arg, scalar);
    }

    (void)fc;
    return ExprPtr{};
}

ExprPtr expand_exp_walker(ExprPtr expr, CASContext& ctx) {
    if (!expr) return expr;
    ExprPtr transformed = rebuild_with_children(expr, ctx, expand_exp_walker);

    const auto* fc = expr_cast<FuncCall>(transformed);
    if (!fc || !is_exp_funccall(*fc)) return transformed;

    ExprPtr local = try_local_exp_rewrite(*fc, fc->args[0], ctx);
    if (!local) return transformed;
    return expand_exp_walker(local, ctx);
}

}  // namespace

Result<bool> mathematically_equal_subset_risch(ExprPtr lhs, ExprPtr rhs, CASContext& context) {
    if (!lhs || !rhs) {
        return fail<bool>(CASError{
            CASErrorKind::InvalidArgument,
            "mathematically_equal_subset_risch: null operand",
            std::nullopt});
    }

    ExprPtr lhs_log = expand_log_walker(lhs, context);
    ExprPtr rhs_log = expand_log_walker(rhs, context);
    ExprPtr lhs_norm = expand_exp_walker(lhs_log, context);
    ExprPtr rhs_norm = expand_exp_walker(rhs_log, context);

    return mathematically_equal(lhs_norm, rhs_norm, context);
}

} // namespace cas::symbolic
