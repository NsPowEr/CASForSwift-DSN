// limit_mrv_set.cpp — MRV set computation and expression rewriting for the
// Gruntz algorithm (Gruntz 1996 §3.5).
//
// Provides:
//   - collect_mrv_candidates  (file-local)
//   - mrv_set                 (external linkage, declared in calculus_internal.hpp)
//   - rewrite_mrv             (external linkage, declared in calculus_internal.hpp)

#include "limit_mrv_internal.hpp"

#include <vector>

namespace cas::calculus {

namespace {

void collect_mrv_candidates(
    ExprPtr e,
    const Symbol& var,
    MRVSet& candidates,
    symbolic::CASContext& ctx) {
    if (!depends_on(e, var)) return;

    if (const auto* sym = expr_cast<Symbol>(e)) {
        if (sym->name == var.name) candidates.insert(e);
        return;
    }

    if (const auto* call = expr_cast<FuncCall>(e)) {
        if (call->func_id == BuiltinOp::Exp) {
            auto lim_u = limit(
                call->args.front(),
                var,
                ctx.arena().make<Constant>(MathConstant::Infinity),
                LimitDirection::Both,
                ctx);
            if (lim_u.is_ok() && limit_is_infinity(lim_u.value())) {
                if (expr_is<Unary>(lim_u.value())) {
                    // exp(-inf) → 0, but its reciprocal exp(inf) is in the MRV set.
                    candidates.insert(ctx.arena().make<FuncCall>(
                        BuiltinOp::Exp,
                        std::vector<ExprPtr>{ctx.arena().make<Unary>(
                            UnaryOp::Neg, call->args.front())}));
                } else {
                    candidates.insert(e);
                }
            }
        } else if (call->func_id == BuiltinOp::Ln) {
            candidates.insert(e);
        }
        for (auto arg : call->args)
            collect_mrv_candidates(arg, var, candidates, ctx);
    } else if (const auto* unary = expr_cast<Unary>(e)) {
        collect_mrv_candidates(unary->operand, var, candidates, ctx);
    } else if (const auto* binary = expr_cast<Binary>(e)) {
        collect_mrv_candidates(binary->left, var, candidates, ctx);
        collect_mrv_candidates(binary->right, var, candidates, ctx);
    } else if (const auto* sum = expr_cast<Sum>(e)) {
        for (auto term : sum->terms)
            collect_mrv_candidates(term, var, candidates, ctx);
    } else if (const auto* product = expr_cast<Product>(e)) {
        for (auto factor : product->factors)
            collect_mrv_candidates(factor, var, candidates, ctx);
    }
}

} // namespace

MRVSet mrv_set(ExprPtr e, const Symbol& var, symbolic::CASContext& ctx) {
    MRVSet candidates;
    collect_mrv_candidates(e, var, candidates, ctx);
    if (candidates.empty()) return {};

    MRVSet mrv;
    for (auto c : candidates) {
        if (mrv.empty()) {
            mrv.insert(c);
        } else {
            int cmp = compare_growth(c, *mrv.begin(), var, ctx);
            if (cmp > 0) {
                mrv.clear();
                mrv.insert(c);
            } else if (cmp == 0) {
                mrv.insert(c);
            }
        }
    }
    return mrv;
}

Result<ExprPtr> rewrite_mrv(
    ExprPtr e,
    const MRVSet& mrv,
    ExprPtr w,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    for (ExprPtr m : mrv) {
        if (structural_equal(e, m)) return ok(w);
    }

    if (expr_is<Symbol>(e) || expr_is<IntegerLit>(e) || expr_is<RationalLit>(e) ||
        expr_is<Constant>(e) || expr_is<DecimalLit>(e)) {
        return ok(e);
    }

    if (const auto* call = expr_cast<FuncCall>(e)) {
        if (call->func_id == BuiltinOp::Exp) {
            ExprPtr neg_arg = ctx.arena().make<Unary>(UnaryOp::Neg, call->args.front());
            auto s_neg = ctx.simplify(neg_arg);
            if (s_neg.is_ok()) {
                ExprPtr exp_neg = ctx.arena().make<FuncCall>(
                    BuiltinOp::Exp, std::vector<ExprPtr>{s_neg.value()});
                for (ExprPtr m : mrv) {
                    if (structural_equal(exp_neg, m)) {
                        // Reciprocal MRV element exp(-u) maps to 1/w,
                        // so exp(u) maps to the denominator of w.
                        if (const auto* bin = expr_cast<Binary>(w)) {
                            if (bin->op == BinaryOp::Div) return ok(bin->right);
                        }
                    }
                }
            }
        }
        std::vector<ExprPtr> new_args;
        for (auto arg : call->args) {
            auto res = rewrite_mrv(arg, mrv, w, var, ctx);
            if (res.is_error()) return res;
            new_args.push_back(res.value());
        }
        return ok(ctx.arena().make<FuncCall>(call->name, std::move(new_args)));
    }
    if (const auto* unary = expr_cast<Unary>(e)) {
        auto res = rewrite_mrv(unary->operand, mrv, w, var, ctx);
        if (res.is_error()) return res;
        return ok(ctx.arena().make<Unary>(unary->op, res.value()));
    }
    if (const auto* binary = expr_cast<Binary>(e)) {
        auto l = rewrite_mrv(binary->left, mrv, w, var, ctx);
        if (l.is_error()) return l;
        auto r = rewrite_mrv(binary->right, mrv, w, var, ctx);
        if (r.is_error()) return r;
        return ok(ctx.arena().make<Binary>(binary->op, l.value(), r.value()));
    }
    if (const auto* sum = expr_cast<Sum>(e)) {
        std::vector<ExprPtr> new_terms;
        for (auto t : sum->terms) {
            auto res = rewrite_mrv(t, mrv, w, var, ctx);
            if (res.is_error()) return res;
            new_terms.push_back(res.value());
        }
        return ok(ctx.arena().make<Sum>(std::move(new_terms)));
    }
    if (const auto* product = expr_cast<Product>(e)) {
        std::vector<ExprPtr> new_factors;
        for (auto f : product->factors) {
            auto res = rewrite_mrv(f, mrv, w, var, ctx);
            if (res.is_error()) return res;
            new_factors.push_back(res.value());
        }
        return ok(ctx.arena().make<Product>(std::move(new_factors)));
    }

    return ok(e);
}

} // namespace cas::calculus
