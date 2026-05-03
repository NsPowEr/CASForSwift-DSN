#include "calculus_internal.hpp"
#include "cas/error.hpp"
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
        if (call->func_id == BuiltinOp::Exp || call->func_id == BuiltinOp::Ln) {
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

// Returns polynomial degree of e w.r.t. var, or nullopt if not a polynomial.
std::optional<int> poly_degree_wrt(ExprPtr e, const Symbol& var) {
    if (!depends_on(e, var)) return 0;
    if (const auto* sym = expr_cast<Symbol>(e)) {
        return (sym->name == var.name) ? 1 : 0;
    }
    if (const auto* sum = expr_cast<Sum>(e)) {
        int deg = 0;
        for (auto t : sum->terms) {
            auto d = poly_degree_wrt(t, var);
            if (!d) return std::nullopt;
            deg = std::max(deg, *d);
        }
        return deg;
    }
    if (const auto* product = expr_cast<Product>(e)) {
        int deg = 0;
        for (auto f : product->factors) {
            auto d = poly_degree_wrt(f, var);
            if (!d) return std::nullopt;
            deg += *d;
        }
        return deg;
    }
    if (const auto* binary = expr_cast<Binary>(e)) {
        if (binary->op == BinaryOp::Pow) {
            auto base_deg = poly_degree_wrt(binary->left, var);
            if (!base_deg) return std::nullopt;
            if (*base_deg == 0) return 0;
            // base depends on var — need non-negative integer exponent
            if (const auto* exp_lit = expr_cast<IntegerLit>(binary->right)) {
                if (!exp_lit->value.is_negative()) {
                    auto eu = exp_lit->value.to_u64();
                    if (eu > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) return std::nullopt;
                    return *base_deg * static_cast<int>(eu);
                }
            }
            return std::nullopt;
        }
        if (binary->op == BinaryOp::Add) {
            auto l = poly_degree_wrt(binary->left, var);
            auto r = poly_degree_wrt(binary->right, var);
            if (!l || !r) return std::nullopt;
            return std::max(*l, *r);
        }
        if (binary->op == BinaryOp::Mul) {
            auto l = poly_degree_wrt(binary->left, var);
            auto r = poly_degree_wrt(binary->right, var);
            if (!l || !r) return std::nullopt;
            return *l + *r;
        }
    }
    return std::nullopt;  // FuncCall or unrecognized — not polynomial
}

// Growth comparison for x → ∞.
// Returns: +1 if a grows faster, -1 if b grows faster, 0 if same rate, 0 for INCOMPARABLE.
// Ranks: constant(0) < ln(...)(1) < polynomial(2, degree as tiebreak) < exp(...)(3)
// For exp vs exp: compare arguments recursively (one level only to avoid infinite recursion).
int compare_growth(ExprPtr a, ExprPtr b, const Symbol& var, symbolic::CASContext& ctx) {
    (void)ctx;
    if (structural_equal(a, b)) return 0;

    // Assigns a coarse rank: 0=constant, 1=log, 2=polynomial, 3=exp
    auto get_growth_rank = [&](ExprPtr e) -> int {
        if (!depends_on(e, var)) return 0;
        if (const auto* call = expr_cast<FuncCall>(e)) {
            if (call->func_id == BuiltinOp::Ln) return 1;
            if (call->func_id == BuiltinOp::Exp) return 3;
        }
        return 2;  // polynomial or mixed expression depending on var
    };

    int ra = get_growth_rank(a);
    int rb = get_growth_rank(b);

    if (ra != rb) return ra > rb ? 1 : -1;

    // Same coarse rank — use polynomial degree as tiebreak
    if (ra == 2) {
        auto da = poly_degree_wrt(a, var);
        auto db = poly_degree_wrt(b, var);
        if (da && db) {
            if (*da != *db) return *da > *db ? 1 : -1;
            return 0;
        }
    }

    // Both exp: compare arguments (non-recursive — one level only)
    if (ra == 3) {
        const auto* ca = expr_cast<FuncCall>(a);
        const auto* cb = expr_cast<FuncCall>(b);
        if (ca && cb && !ca->args.empty() && !cb->args.empty()) {
            auto darg_a = poly_degree_wrt(ca->args[0], var);
            auto darg_b = poly_degree_wrt(cb->args[0], var);
            if (darg_a && darg_b && *darg_a != *darg_b)
                return *darg_a > *darg_b ? 1 : -1;
        }
    }

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

    auto series = taylor_series(simplified.value(), w_var, limit_make_integer(arena, 0), 4U, ctx);
    if (series.is_error()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "MRV taylor series failed: " + series.error().message,
            .hint = std::nullopt,
        });
    }

    auto res = ctx.substitute(series.value().polynomial, w_var, limit_make_integer(arena, 0));
    if (res.is_ok()) {
        auto final_res = ctx.simplify(res.value());
        if (final_res.is_ok() && !depends_on(final_res.value(), w_var)) return final_res;
    }

    return ok(series.value().polynomial);
}

} // namespace cas::calculus
