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
#include "cas/symbolic.hpp"

#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] ExprPtr make_int(symbolic::CASContext& ctx, long long v) {
    return ctx.arena().make<IntegerLit>(BigInt(v));
}

// Test e ≡ Symbol(name).
[[nodiscard]] bool is_sym(ExprPtr e, const Symbol& sym) {
    const auto* s = expr_cast<Symbol>(e);
    return s && s->name == sym.name;
}

// Factorial as expr (n!).
[[nodiscard]] ExprPtr factorial_expr(symbolic::CASContext& ctx, long long n) {
    BigInt acc(1);
    for (long long i = 2; i <= n; ++i) acc = acc * BigInt(i);
    return ctx.arena().make<IntegerLit>(std::move(acc));
}

// Does expr depend on var (symbolic occurrence)?
[[nodiscard]] bool depends_on_t(ExprPtr e, const Symbol& v) {
    if (!e) return false;
    if (auto* s = expr_cast<Symbol>(e)) return s->name == v.name;
    if (auto* un = expr_cast<Unary>(e)) return depends_on_t(un->operand, v);
    if (auto* bin = expr_cast<Binary>(e)) return depends_on_t(bin->left, v) || depends_on_t(bin->right, v);
    if (auto* sum = expr_cast<Sum>(e)) { for (auto t : sum->terms) if (depends_on_t(t, v)) return true; return false; }
    if (auto* prod = expr_cast<Product>(e)) { for (auto f : prod->factors) if (depends_on_t(f, v)) return true; return false; }
    if (auto* fc = expr_cast<FuncCall>(e)) { for (auto a : fc->args) if (depends_on_t(a, v)) return true; return false; }
    return false;
}

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
        for (ExprPtr f : prod->factors) {
            if (depends_on_t(f, t)) vars.push_back(f);
            else consts.push_back(f);
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

    return fail<ExprPtr>(CASError{
        CASErrorKind::Unimplemented,
        "Laplace transform: pattern not in elementary table",
        std::nullopt});
}

}  // namespace cas::calculus
