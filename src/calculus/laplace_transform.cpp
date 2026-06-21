// CAS-L3-07 — Laplace transform L{f(t)}(s) MVP.
//
// Pattern-based implementation covering elementary functions:
//   L{1}      = 1/s
//   L{t^n}    = n! / s^(n+1)        (n ∈ N)
//   L{exp(a·t)} = 1/(s-a)
//   L{sin(a·t)} = a / (s²+a²)
//   L{cos(a·t)} = s / (s²+a²)
//   L{a·f + b·g} = a·L{f} + b·L{g}  (linearity)
//   L{exp(a·t)·f(t)} = F(s-a)        (frequency shift) — for elementary f
//
// Returns Unimplemented for non-elementary integrands. Inverse transform
// deferred follow-up (requires partial-fraction decomp on F(s)).

#include "calculus_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <vector>

#include "laplace_transform_internal.hpp"

namespace cas::calculus {

namespace {

using laplace_detail::make_int;
using laplace_detail::is_sym;
using laplace_detail::factorial_expr;
using laplace_detail::depends_on_t;

// L{c} = c/s when c independent of t.
[[nodiscard]] ExprPtr laplace_constant(ExprPtr c, const Symbol& s, symbolic::CASContext& ctx) {
    ExprPtr s_expr = ctx.arena().make<Symbol>(s);
    return ctx.arena().make<Binary>(BinaryOp::Div, c, s_expr);
}

}  // namespace

[[nodiscard]] Result<ExprPtr> laplace_transform(
    ExprPtr expr, const Symbol& t, const Symbol& s,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr s_expr = arena.make<Symbol>(s);

    // Constant (t-independent): L{c} = c/s
    if (!depends_on_t(expr, t)) {
        return ctx.simplify(laplace_constant(expr, s, ctx));
    }

    // L{t} = 1/s²
    if (is_sym(expr, t)) {
        ExprPtr two = make_int(ctx, 2);
        ExprPtr s_sq = arena.make<Binary>(BinaryOp::Pow, s_expr, two);
        return ctx.simplify(arena.make<Binary>(BinaryOp::Div, make_int(ctx, 1), s_sq));
    }

    // L{t^n} = n!/s^(n+1), n ∈ N
    if (const auto* bin = expr_cast<Binary>(expr); bin && bin->op == BinaryOp::Pow) {
        if (is_sym(bin->left, t)) {
            if (const auto* n_lit = expr_cast<IntegerLit>(bin->right);
                n_lit && !n_lit->value.is_negative() && n_lit->value.bit_length() < 32) {
                const long long n = static_cast<long long>(n_lit->value.to_u64());
                ExprPtr fact = factorial_expr(ctx, n);
                ExprPtr s_pow = arena.make<Binary>(BinaryOp::Pow, s_expr, make_int(ctx, n + 1));
                return ctx.simplify(arena.make<Binary>(BinaryOp::Div, fact, s_pow));
            }
        }
    }

    // L{exp(a·t)} = 1/(s-a)
    if (const auto* call = expr_cast<FuncCall>(expr);
        call && call->func_id == BuiltinOp::Exp && call->args.size() == 1U) {
        // Extract a·t pattern: arg = a · t (Product) or t alone.
        ExprPtr arg = call->args[0];
        ExprPtr a = nullptr;
        if (is_sym(arg, t)) {
            a = make_int(ctx, 1);
        } else if (auto* prod = expr_cast<Product>(arg)) {
            std::vector<ExprPtr> non_t;
            bool found_t = false;
            for (ExprPtr f : prod->factors) {
                if (is_sym(f, t)) {
                    if (found_t) { a = nullptr; break; }
                    found_t = true;
                } else if (!depends_on_t(f, t)) {
                    non_t.push_back(f);
                } else {
                    a = nullptr; break;
                }
            }
            if (found_t) {
                a = non_t.empty()
                    ? make_int(ctx, 1)
                    : (non_t.size() == 1 ? non_t[0]
                       : static_cast<ExprPtr>(arena.make<Product>(std::move(non_t))));
            }
        }
        if (a) {
            ExprPtr s_minus_a = arena.make<Binary>(BinaryOp::Sub, s_expr, a);
            return ctx.simplify(arena.make<Binary>(BinaryOp::Div, make_int(ctx, 1), s_minus_a));
        }
    }

    // L{sin(a·t)} = a / (s² + a²);  L{cos(a·t)} = s / (s² + a²)
    if (const auto* call = expr_cast<FuncCall>(expr);
        call && (call->func_id == BuiltinOp::Sin || call->func_id == BuiltinOp::Cos)
        && call->args.size() == 1U) {
        ExprPtr arg = call->args[0];
        ExprPtr a = nullptr;
        if (is_sym(arg, t)) {
            a = make_int(ctx, 1);
        } else if (auto* prod = expr_cast<Product>(arg)) {
            std::vector<ExprPtr> non_t;
            bool found_t = false;
            for (ExprPtr f : prod->factors) {
                if (is_sym(f, t)) { found_t = true; }
                else if (!depends_on_t(f, t)) non_t.push_back(f);
                else { a = nullptr; break; }
            }
            if (found_t) {
                a = non_t.empty() ? make_int(ctx, 1)
                    : (non_t.size() == 1 ? non_t[0]
                       : static_cast<ExprPtr>(arena.make<Product>(std::move(non_t))));
            }
        }
        if (a) {
            ExprPtr two = make_int(ctx, 2);
            ExprPtr s_sq = arena.make<Binary>(BinaryOp::Pow, s_expr, two);
            ExprPtr a_sq = arena.make<Binary>(BinaryOp::Pow, a, two);
            ExprPtr den = arena.make<Binary>(BinaryOp::Add, s_sq, a_sq);
            ExprPtr num = (call->func_id == BuiltinOp::Sin) ? a : s_expr;
            return ctx.simplify(arena.make<Binary>(BinaryOp::Div, num, den));
        }
    }

    // Linearity: L{Σ c_i · f_i(t)} = Σ c_i · L{f_i(t)}
    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> terms;
        for (ExprPtr term : sum->terms) {
            auto l = laplace_transform(term, t, s, ctx);
            if (l.is_error()) return l;
            terms.push_back(l.value());
        }
        return ctx.simplify(arena.make<Sum>(std::move(terms)));
    }

    // Constant scalar product: c · f(t) → c · L{f(t)} when c indep of t.
    if (const auto* prod = expr_cast<Product>(expr)) {
        std::vector<ExprPtr> consts, vars;
        ExprPtr exp_factor = nullptr;
        ExprPtr exp_arg = nullptr;
        std::vector<ExprPtr> non_exp_vars;
        for (ExprPtr f : prod->factors) {
            if (depends_on_t(f, t)) {
                vars.push_back(f);
                // Detect exp(a·t) factor for frequency-shift rule.
                if (auto* call = expr_cast<FuncCall>(f);
                    call && call->func_id == BuiltinOp::Exp && call->args.size() == 1U) {
                    if (!exp_factor) {
                        exp_factor = f;
                        exp_arg = call->args[0];
                        continue;
                    }
                }
                non_exp_vars.push_back(f);
            }
            else consts.push_back(f);
        }
        // Frequency shift: L{exp(a·t)·g(t)} = G(s - a)
        // where a must be constant in t and a = coefficient of t in exp_arg.
        if (exp_factor && !non_exp_vars.empty()) {
            // Extract a from a·t pattern.
            ExprPtr a = nullptr;
            if (is_sym(exp_arg, t)) {
                a = make_int(ctx, 1);
            } else if (auto* p_arg = expr_cast<Product>(exp_arg)) {
                std::vector<ExprPtr> a_factors;
                bool found_t = false;
                bool ok = true;
                for (auto pf : p_arg->factors) {
                    if (is_sym(pf, t)) {
                        if (found_t) { ok = false; break; }
                        found_t = true;
                    } else if (!depends_on_t(pf, t)) {
                        a_factors.push_back(pf);
                    } else { ok = false; break; }
                }
                if (ok && found_t) {
                    a = a_factors.empty() ? make_int(ctx, 1)
                        : (a_factors.size() == 1 ? a_factors[0]
                           : static_cast<ExprPtr>(arena.make<Product>(std::move(a_factors))));
                }
            }
            if (a) {
                // Compute L{g(t)} = G(s), then substitute s → s - a.
                ExprPtr g_only = non_exp_vars.size() == 1
                    ? non_exp_vars[0]
                    : static_cast<ExprPtr>(arena.make<Product>(std::move(non_exp_vars)));
                auto Gs = laplace_transform(g_only, t, s, ctx);
                if (Gs.is_ok()) {
                    ExprPtr s_minus_a = arena.make<Binary>(BinaryOp::Sub,
                        arena.make<Symbol>(s), a);
                    auto shifted = symbolic::substitute(Gs.value(), s, s_minus_a, ctx);
                    if (shifted.is_ok()) {
                        consts.push_back(shifted.value());
                        return ctx.simplify(consts.size() == 1
                            ? consts[0]
                            : static_cast<ExprPtr>(arena.make<Product>(std::move(consts))));
                    }
                }
            }
        }
        if (!consts.empty() && !vars.empty()) {
            ExprPtr f_only = vars.size() == 1 ? vars[0]
                : static_cast<ExprPtr>(arena.make<Product>(std::move(vars)));
            auto l = laplace_transform(f_only, t, s, ctx);
            if (l.is_error()) return l;
            consts.push_back(l.value());
            return ctx.simplify(arena.make<Product>(std::move(consts)));
        }
    }

    // F0.8-MIGRATED
    return make_unimplemented<ExprPtr>(
        "calculus", "laplace_transform",
        "expression not in elementary Laplace table",
        cas::error::reason_codes::LAPLACE_UNKNOWN_FORM,
        "Extend the Laplace table or implement algorithmic transform (e.g. Meijer G-functions)",
        "F0.8");
}

}  // namespace cas::calculus
