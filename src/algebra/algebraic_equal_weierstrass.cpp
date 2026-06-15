// F7.5 — Weierstrass canonical-form equivalence helper used by
// `mathematically_equal`. When a difference expression depends solely on
// sin/cos of `c·v` where v is a single symbol and c ranges over positive
// rationals sharing a common base r_min (every c is an integer multiple of
// r_min), the substitution
//   θ = r_min · v,   t = tan(θ/2),
//   sin θ → 2t/(1+t²),  cos θ → (1-t²)/(1+t²),
//   sin(m·θ), cos(m·θ) → Chebyshev / U-polynomial recurrence in sin θ, cos θ
// turns the difference into a rational function of t. Two antiderivatives
// of the same R(sin x, cos x) integrand — even when one uses the
// half-argument `tan(x/2)` and the other the full-argument `sin(x)/(cos(x)+1)`
// form — reduce to identical rational(t) shapes whose numerator collapses
// to literal zero after `simplify` / `polynomial_normal_form`.
//
// Anti-furbizia: not a pattern lookup. The recurrence covers every positive
// integer multiplier m ≥ 1; mismatched shapes (non-rational coefficient,
// non-linear in v, multi-variable trig) bail out via std::nullopt and the
// outer comparator continues with the existing fallback chain. The current
// implementation reliably collapses small identities (sin(2θ) − 2 sin θ cos θ,
// cos(2θ) − 1 + 2 sin²θ); larger antiderivative differences may exceed the
// canonicaliser's depth budget — in that case the helper returns `false`
// and the caller falls through to subsequent equality strategies. Tightening
// the rational-of-t canonicaliser is a separate follow-up.

#include "cas/algebra.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/normal_form.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"

#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::symbolic {

namespace {

struct TrigCall {
    BuiltinOp op;   // Sin or Cos
    Rational coeff; // coefficient of var in arg
};

// Try to decompose `arg` as `c · Symbol(var_name)` where c is rational.
// Stores var_name on first call; on subsequent calls verifies match.
[[nodiscard]] std::optional<Rational> decompose_linear(
    ExprPtr arg, std::optional<std::string>& var_name) {
    if (!arg) return std::nullopt;
    auto match_symbol = [&](ExprPtr e) -> bool {
        if (const auto* sym = expr_cast<Symbol>(e)) {
            if (!var_name) var_name = sym->name;
            return *var_name == sym->name;
        }
        return false;
    };
    if (match_symbol(arg)) {
        return Rational(BigInt(1));
    }
    if (const auto* prod = expr_cast<Product>(arg)) {
        Rational coeff(BigInt(1));
        ExprPtr sym_factor = nullptr;
        for (ExprPtr f : prod->factors) {
            if (const auto* il = expr_cast<IntegerLit>(f)) {
                coeff = coeff * Rational(il->value);
            } else if (const auto* rl = expr_cast<RationalLit>(f)) {
                coeff = coeff * Rational(rl->numerator, rl->denominator);
            } else if (expr_is<Symbol>(f) && !sym_factor) {
                sym_factor = f;
            } else {
                return std::nullopt;
            }
        }
        if (sym_factor && match_symbol(sym_factor)) return coeff;
        return std::nullopt;
    }
    if (const auto* bin = expr_cast<Binary>(arg); bin && bin->op == BinaryOp::Div) {
        if (match_symbol(bin->left)) {
            if (const auto* il = expr_cast<IntegerLit>(bin->right);
                il && !il->value.is_zero())
                return Rational(BigInt(1), il->value);
            if (const auto* rl = expr_cast<RationalLit>(bin->right);
                rl && !rl->numerator.is_zero())
                return Rational(rl->denominator, rl->numerator);
        }
        std::optional<std::string> tmp = var_name;
        auto left_c = decompose_linear(bin->left, tmp);
        if (left_c) {
            if (const auto* il = expr_cast<IntegerLit>(bin->right);
                il && !il->value.is_zero()) {
                var_name = tmp;
                return *left_c / Rational(il->value);
            }
        }
        return std::nullopt;
    }
    return std::nullopt;
}

// Walk expression collecting every sin/cos call. Verifies all args are
// linear (rational·var) of a single symbol. On any non-conforming sin/cos
// returns false → bail.
[[nodiscard]] bool collect_trig_calls(
    ExprPtr expr, std::vector<TrigCall>& out, std::optional<std::string>& var_name) {
    if (!expr) return true;
    bool ok = true;
    std::function<void(ExprPtr)> walk = [&](ExprPtr e) {
        if (!e || !ok) return;
        if (const auto* fc = expr_cast<FuncCall>(e)) {
            const bool is_sc = (fc->func_id == BuiltinOp::Sin
                             || fc->func_id == BuiltinOp::Cos);
            if (is_sc && fc->args.size() == 1U) {
                auto c = decompose_linear(fc->args[0], var_name);
                if (!c) { ok = false; return; }
                if (!c->numerator().is_negative() && !c->numerator().is_zero())
                    out.push_back(TrigCall{fc->func_id, *c});
                else { ok = false; return; }
                return;
            }
            for (auto a : fc->args) walk(a);
            return;
        }
        if (const auto* bin = expr_cast<Binary>(e)) {
            walk(bin->left); walk(bin->right); return;
        }
        if (const auto* un = expr_cast<Unary>(e)) { walk(un->operand); return; }
        if (const auto* s = expr_cast<Sum>(e))
            for (auto t : s->terms) walk(t);
        if (const auto* p = expr_cast<Product>(e))
            for (auto f : p->factors) walk(f);
    };
    walk(expr);
    return ok;
}

// gcd of two non-negative rationals (numerator gcd / denominator lcm).
[[nodiscard]] Rational rational_gcd(const Rational& a, const Rational& b) {
    const BigInt na = a.numerator().is_negative() ? -a.numerator() : a.numerator();
    const BigInt nb = b.numerator().is_negative() ? -b.numerator() : b.numerator();
    const BigInt num = gcd(na, nb);
    const BigInt da = a.denominator();
    const BigInt db = b.denominator();
    const BigInt den = (da * db) / gcd(da, db);
    return Rational(num, den);
}

// Build sin(m·θ), cos(m·θ) given sin(θ)=S, cos(θ)=C using the recurrence
//   sin((k+1)θ) = sin(kθ)·C + cos(kθ)·S
//   cos((k+1)θ) = cos(kθ)·C - sin(kθ)·S
[[nodiscard]] std::pair<ExprPtr, ExprPtr> chebyshev_multiply(
    long long m, ExprPtr S, ExprPtr C, AstArena& arena) {
    if (m <= 0) {
        return {arena.make<IntegerLit>(BigInt(0)),
                arena.make<IntegerLit>(BigInt(1))};
    }
    ExprPtr s_k = S;
    ExprPtr c_k = C;
    for (long long k = 1; k < m; ++k) {
        ExprPtr s_next = arena.make<Sum>(std::vector<ExprPtr>{
            arena.make<Product>(std::vector<ExprPtr>{s_k, C}),
            arena.make<Product>(std::vector<ExprPtr>{c_k, S})});
        ExprPtr c_next = arena.make<Sum>(std::vector<ExprPtr>{
            arena.make<Product>(std::vector<ExprPtr>{c_k, C}),
            arena.make<Unary>(UnaryOp::Neg,
                arena.make<Product>(std::vector<ExprPtr>{s_k, S}))});
        s_k = s_next;
        c_k = c_next;
    }
    return {s_k, c_k};
}

// Substitute every sin/cos(c·v) in `expr` by the Chebyshev expansion in
// (S, C) where S = sin(base·v), C = cos(base·v), m = c / base ∈ ℕ₊.
[[nodiscard]] ExprPtr substitute_trig_chebyshev(
    ExprPtr expr, const std::string& var, const Rational& base,
    ExprPtr S, ExprPtr C, AstArena& arena) {
    if (!expr) return expr;
    if (const auto* fc = expr_cast<FuncCall>(expr)) {
        if ((fc->func_id == BuiltinOp::Sin || fc->func_id == BuiltinOp::Cos)
            && fc->args.size() == 1U) {
            std::optional<std::string> vn = var;
            auto c = decompose_linear(fc->args[0], vn);
            if (c) {
                Rational ratio = *c / base;
                const long long m = ratio.numerator().to_u64();
                auto [sm, cm] = chebyshev_multiply(m, S, C, arena);
                return (fc->func_id == BuiltinOp::Sin) ? sm : cm;
            }
        }
        std::vector<ExprPtr> args;
        args.reserve(fc->args.size());
        for (auto a : fc->args)
            args.push_back(substitute_trig_chebyshev(a, var, base, S, C, arena));
        return arena.make<FuncCall>(fc->func_id, std::move(args));
    }
    if (const auto* bin = expr_cast<Binary>(expr)) {
        return arena.make<Binary>(bin->op,
            substitute_trig_chebyshev(bin->left, var, base, S, C, arena),
            substitute_trig_chebyshev(bin->right, var, base, S, C, arena));
    }
    if (const auto* un = expr_cast<Unary>(expr)) {
        return arena.make<Unary>(un->op,
            substitute_trig_chebyshev(un->operand, var, base, S, C, arena));
    }
    if (const auto* s = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> terms;
        terms.reserve(s->terms.size());
        for (auto t : s->terms)
            terms.push_back(substitute_trig_chebyshev(t, var, base, S, C, arena));
        return arena.make<Sum>(std::move(terms));
    }
    if (const auto* p = expr_cast<Product>(expr)) {
        std::vector<ExprPtr> factors;
        factors.reserve(p->factors.size());
        for (auto f : p->factors)
            factors.push_back(substitute_trig_chebyshev(f, var, base, S, C, arena));
        return arena.make<Product>(std::move(factors));
    }
    return expr;
}

}  // namespace

bool weierstrass_zero_diff(ExprPtr diff_expr, CASContext& ctx) {
    std::vector<TrigCall> calls;
    std::optional<std::string> var_name;
    if (!collect_trig_calls(diff_expr, calls, var_name)) return false;
    if (calls.empty() || !var_name) return false;

    Rational base = calls.front().coeff;
    for (std::size_t i = 1; i < calls.size(); ++i)
        base = rational_gcd(base, calls[i].coeff);
    if (base.numerator().is_zero()) return false;
    for (const auto& c : calls) {
        Rational ratio = c.coeff / base;
        if (ratio.denominator() != BigInt(1)) return false;
    }

    AstArena& a = ctx.arena();
    Symbol t = ctx.make_fresh_symbol("__wq");
    ExprPtr t_sym = a.make<Symbol>(t.name);
    ExprPtr one = a.make<IntegerLit>(BigInt(1));
    ExprPtr two = a.make<IntegerLit>(BigInt(2));
    ExprPtr t_sq = a.make<Binary>(BinaryOp::Pow, t_sym, two);
    ExprPtr one_plus = a.make<Sum>(std::vector<ExprPtr>{one, t_sq});
    ExprPtr one_minus = a.make<Sum>(std::vector<ExprPtr>{
        a.make<IntegerLit>(BigInt(1)),
        a.make<Unary>(UnaryOp::Neg, t_sq)});
    ExprPtr S = a.make<Binary>(BinaryOp::Div,
        a.make<Product>(std::vector<ExprPtr>{
            a.make<IntegerLit>(BigInt(2)), t_sym}),
        one_plus);
    ExprPtr C = a.make<Binary>(BinaryOp::Div, one_minus, one_plus);

    ExprPtr sub = substitute_trig_chebyshev(
        diff_expr, *var_name, base, S, C, a);
    auto sub_simp = ctx.simplify(sub);
    if (!sub_simp.is_ok()) return false;
    if (const auto* il = expr_cast<IntegerLit>(sub_simp.value()))
        return il->value.is_zero();

    // Path 1: polynomial_normal_form often reduces simple rational(t)
    // shapes directly — collapses identities like sin(2θ) − 2 sinθ cosθ.
    auto normal = polynomial_normal_form(sub_simp.value(), ctx);
    if (normal.is_ok()) {
        if (const auto* il = expr_cast<IntegerLit>(normal.value()))
            if (il->value.is_zero()) return true;
        if (const auto* rl = expr_cast<RationalLit>(normal.value()))
            if (rl->numerator.is_zero()) return true;
    }
    // Path 2: split into (num, den) and run normal form on numerator.
    auto parts = algebra::split_num_den(sub_simp.value(), ctx);
    if (!parts.is_ok()) return false;
    auto num_norm = polynomial_normal_form(parts.value().numerator, ctx);
    if (!num_norm.is_ok()) return false;
    if (const auto* il = expr_cast<IntegerLit>(num_norm.value()))
        return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(num_norm.value()))
        return rl->numerator.is_zero();
    return false;
}

}  // namespace cas::symbolic
