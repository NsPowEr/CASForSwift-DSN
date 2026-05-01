#include "calculus_internal.hpp"
#include "cas/ast_debug.hpp"
#include "cas/error.hpp"
#include "cas/formatter.hpp"
#include "cas/symbolic.hpp"

#include <algorithm>
#include <set>
#include <vector>

namespace cas::calculus {

bool MRVCompare::operator()(ExprPtr lhs, ExprPtr rhs) const noexcept {
    return symbolic::canonical_compare(lhs, rhs) < 0;
}

namespace {

void collect_mrv_candidates(ExprPtr e, const Symbol& var, MRVSet& candidates, symbolic::CASContext& ctx) {
    if (!depends_on(e, var)) return;

    if (const auto* sym = expr_cast<Symbol>(e)) {
        if (sym->name == var.name) candidates.insert(e);
        return;
    }

    if (const auto* call = expr_cast<FuncCall>(e)) {
        if (call->func_id == BuiltinOp::Exp) {
            candidates.insert(e);
        }
        for (auto arg : call->args) collect_mrv_candidates(arg, var, candidates, ctx);
    } else if (const auto* unary = expr_cast<Unary>(e)) {
        collect_mrv_candidates(unary->operand, var, candidates, ctx);
    } else if (const auto* binary = expr_cast<Binary>(e)) {
        collect_mrv_candidates(binary->left, var, candidates, ctx);
        collect_mrv_candidates(binary->right, var, candidates, ctx);
    } else if (const auto* sum = expr_cast<Sum>(e)) {
        for (auto term : sum->terms) collect_mrv_candidates(term, var, candidates, ctx);
    } else if (const auto* product = expr_cast<Product>(e)) {
        for (auto factor : product->factors) collect_mrv_candidates(factor, var, candidates, ctx);
    }
}

int compare_growth(ExprPtr a, ExprPtr b, const Symbol&, symbolic::CASContext&) {
    if (structural_equal(a, b)) return 0;

    auto get_exp_depth = [](ExprPtr e) {
        int depth = 0;
        while (const auto* call = expr_cast<FuncCall>(e)) {
            if (call->func_id == BuiltinOp::Exp) {
                depth++;
                e = call->args[0];
            } else break;
        }
        return depth;
    };

    int da = get_exp_depth(a);
    int db = get_exp_depth(b);
    if (da != db) return da > db ? 1 : -1;

    return 0; 
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

Result<ExprPtr> rewrite_mrv(ExprPtr e, const MRVSet& mrv, ExprPtr w, const Symbol& var, symbolic::CASContext& ctx) {
    if (mrv.count(e)) return ok(w);

    if (expr_is<Symbol>(e) || expr_is<IntegerLit>(e) || expr_is<RationalLit>(e) || 
        expr_is<Constant>(e) || expr_is<DecimalLit>(e)) {
        return ok(e);
    }

    if (const auto* call = expr_cast<FuncCall>(e)) {
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

Result<ExprPtr> compute_limit_mrv(ExprPtr expr, const Symbol& var, ExprPtr point, symbolic::CASContext& ctx) {
    if (!limit_is_infinity(point)) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "MRV solo per x -> inf",
            .hint = std::nullopt
        });
    }

    MRVSet mrv = mrv_set(expr, var, ctx);
    if (mrv.empty()) return ok(expr);

    AstArena& arena = ctx.arena();
    Symbol w_var("w");
    ExprPtr w_sym = arena.make<Symbol>("w");

    ExprPtr replacement = arena.make<Binary>(BinaryOp::Div, limit_make_integer(arena, 1), w_sym);
    auto rewritten = rewrite_mrv(expr, mrv, replacement, var, ctx);
    if (rewritten.is_error()) return rewritten;

    auto simplified = ctx.simplify(rewritten.value());
    if (simplified.is_error()) return simplified;

    // Use Taylor series on w -> 0
    auto series = taylor_series(simplified.value(), w_var, limit_make_integer(arena, 0), 4U, ctx);
    if (series.is_error()) return fail<ExprPtr>(series.error());

    auto res = ctx.substitute(series.value().polynomial, w_var, limit_make_integer(arena, 0));
    if (res.is_ok()) {
        auto final_res = ctx.simplify(res.value());
        if (final_res.is_ok() && !depends_on(final_res.value(), w_var)) return final_res;
    }

    return ok(series.value().polynomial);
}

} // namespace cas::calculus
