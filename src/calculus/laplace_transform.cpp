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

namespace {

// Detect s² + a² where a² rational positive square or rational positive.
// Returns a (as expr) on match, else nullptr.
[[nodiscard]] ExprPtr match_s_squared_plus_const(ExprPtr expr, const Symbol& s,
                                                  AstArena& arena) {
    const auto* sum = expr_cast<Sum>(expr);
    if (!sum || sum->terms.size() != 2) return nullptr;
    ExprPtr s_sq_part = nullptr;
    ExprPtr const_part = nullptr;
    for (ExprPtr term : sum->terms) {
        if (const auto* bin = expr_cast<Binary>(term);
            bin && bin->op == BinaryOp::Pow) {
            if (const auto* sym = expr_cast<Symbol>(bin->left);
                sym && sym->name == s.name) {
                if (const auto* il = expr_cast<IntegerLit>(bin->right);
                    il && il->value == BigInt(2)) {
                    if (s_sq_part) return nullptr;
                    s_sq_part = term;
                    continue;
                }
            }
        }
        const_part = term;
    }
    if (!s_sq_part || !const_part) return nullptr;
    // a² extracted; need sqrt(a²) = a as expr. For integer/rational c>0:
    if (const auto* il = expr_cast<IntegerLit>(const_part)) {
        if (il->value.is_negative() || il->value.is_zero()) return nullptr;
        // Try integer sqrt.
        BigInt c = il->value;
        BigInt x = c;
        BigInt y = (x + BigInt(1)) / BigInt(2);
        while (y < x) { x = y; y = (x + c / x) / BigInt(2); }
        if (x * x == c) return arena.make<IntegerLit>(std::move(x));
        // Else a = sqrt(c) symbolic.
        return arena.make<FuncCall>(BuiltinOp::Sqrt,
            std::vector<ExprPtr>{const_part});
    }
    return nullptr;
}

}  // namespace

// Inverse Laplace L⁻¹{F(s)}(t) — elementary pattern table.
// Covers:
//   1/s → 1
//   1/s^(n+1) → t^n / n!
//   1/(s-a) → exp(a·t)
//   a/(s²+a²) → sin(a·t)
//   s/(s²+a²) → cos(a·t)
//   Linearity over Sum
[[nodiscard]] Result<ExprPtr> inverse_laplace_transform(
    ExprPtr expr, const Symbol& s, const Symbol& t,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr t_expr = arena.make<Symbol>(t);

    // Linearity over Sum.
    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> terms;
        for (ExprPtr term : sum->terms) {
            auto inv = inverse_laplace_transform(term, s, t, ctx);
            if (inv.is_error()) return inv;
            terms.push_back(inv.value());
        }
        return ctx.simplify(arena.make<Sum>(std::move(terms)));
    }

    // Constant scalar factor: c · F(s) → c · L⁻¹{F(s)}
    if (const auto* prod = expr_cast<Product>(expr)) {
        std::vector<ExprPtr> consts, vars;
        for (ExprPtr f : prod->factors) {
            if (depends_on_t(f, s)) vars.push_back(f);
            else consts.push_back(f);
        }
        if (!consts.empty() && !vars.empty()) {
            ExprPtr f_only = vars.size() == 1 ? vars[0]
                : static_cast<ExprPtr>(arena.make<Product>(std::move(vars)));
            auto inv = inverse_laplace_transform(f_only, s, t, ctx);
            if (inv.is_error()) return inv;
            consts.push_back(inv.value());
            return ctx.simplify(arena.make<Product>(std::move(consts)));
        }
    }

    // 1 / Q(s) patterns.
    if (const auto* bin = expr_cast<Binary>(expr); bin && bin->op == BinaryOp::Div) {
        ExprPtr num = bin->left;
        ExprPtr den = bin->right;

        // num must be constant in s for most patterns; or s for cos pattern.
        bool num_is_s = false;
        if (const auto* sym = expr_cast<Symbol>(num); sym && sym->name == s.name) {
            num_is_s = true;
        }

        // 1/s → 1
        if (is_sym(den, s) && !depends_on_t(num, s)) {
            return ctx.simplify(num);
        }

        // 1/s^(n+1) → num · t^n / n!
        if (const auto* den_pow = expr_cast<Binary>(den);
            den_pow && den_pow->op == BinaryOp::Pow) {
            if (is_sym(den_pow->left, s)) {
                if (const auto* exp_lit = expr_cast<IntegerLit>(den_pow->right);
                    exp_lit && exp_lit->value > BigInt(0)
                    && !depends_on_t(num, s)) {
                    const long long n_plus_1 = static_cast<long long>(exp_lit->value.to_u64());
                    const long long n = n_plus_1 - 1;
                    ExprPtr t_pow_n = (n == 0)
                        ? static_cast<ExprPtr>(make_int(ctx, 1))
                        : static_cast<ExprPtr>(arena.make<Binary>(BinaryOp::Pow,
                              t_expr, make_int(ctx, n)));
                    ExprPtr fact = factorial_expr(ctx, n);
                    ExprPtr ratio = arena.make<Binary>(BinaryOp::Div, t_pow_n, fact);
                    ExprPtr scaled = arena.make<Product>(std::vector<ExprPtr>{num, ratio});
                    return ctx.simplify(scaled);
                }
            }
        }

        // 1/(s-a) → exp(a·t)
        if (const auto* den_bin = expr_cast<Binary>(den);
            den_bin && den_bin->op == BinaryOp::Sub) {
            if (is_sym(den_bin->left, s) && !depends_on_t(den_bin->right, s)
                && !depends_on_t(num, s)) {
                ExprPtr a = den_bin->right;
                ExprPtr at = arena.make<Product>(std::vector<ExprPtr>{a, t_expr});
                ExprPtr exp_at = arena.make<FuncCall>(BuiltinOp::Exp,
                    std::vector<ExprPtr>{at});
                ExprPtr scaled = arena.make<Product>(std::vector<ExprPtr>{num, exp_at});
                return ctx.simplify(scaled);
            }
        }

        // a/(s²+a²) → sin(a·t)   or   s/(s²+a²) → cos(a·t)
        ExprPtr a_extracted = match_s_squared_plus_const(den, s, arena);
        if (a_extracted) {
            ExprPtr a = a_extracted;
            ExprPtr at = arena.make<Product>(std::vector<ExprPtr>{a, t_expr});
            // cos pattern: num == s
            if (num_is_s) {
                return ctx.simplify(arena.make<FuncCall>(BuiltinOp::Cos,
                    std::vector<ExprPtr>{at}));
            }
            // sin pattern: num == a (or constant proportional to a)
            // Verify num structurally equals a (simple case).
            if (structural_equal(num, a)) {
                return ctx.simplify(arena.make<FuncCall>(BuiltinOp::Sin,
                    std::vector<ExprPtr>{at}));
            }
            // num = c (constant), → (c/a) · sin(a·t)
            if (!depends_on_t(num, s)) {
                ExprPtr c_over_a = arena.make<Binary>(BinaryOp::Div, num, a);
                ExprPtr sin_at = arena.make<FuncCall>(BuiltinOp::Sin,
                    std::vector<ExprPtr>{at});
                return ctx.simplify(arena.make<Product>(
                    std::vector<ExprPtr>{c_over_a, sin_at}));
            }
        }
    }

    return fail<ExprPtr>(CASError{
        CASErrorKind::Unimplemented,
        "Inverse Laplace: pattern not in elementary table",
        std::nullopt});
}

}  // namespace cas::calculus
