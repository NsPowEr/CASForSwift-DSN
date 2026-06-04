// F5.1 / B9-Task#22 — Risch Bronstein cap.8 estensioni trascendentali.
//
// Implementa la Risch differential equation
//
//     y' + f·y = g
//
// per due classi di estensioni trascendentali del campo k = Q(var):
//
//   1. solve_risch_de_logarithmic_q — θ = log(u), θ' = u'/u ∈ k.
//      Bronstein §6.4.1 "PolyRischDE_logarithmic".  Algoritmo top-down
//      sui gradi di θ con accoppiamento via shift y_{i+1}:
//         y_i' + f₀·y_i = g_i − (i+1)·(u'/u)·y_{i+1}
//
//   2. solve_risch_de_exponential_q — θ = exp(u), θ' = u'·θ.
//      Bronstein §6.4.2 "PolyRischDE_exponential".  Equazioni disaccoppiate
//      grazie a θ^(i-1)·θ' = u'·θ^i:
//         y_i' + (i·u' + f₀)·y_i = g_i
//
// In entrambi i casi:
//   - sub-case f costante in θ;  caso f polinomiale in θ marcato
//     Unimplemented diagnostico esplicito (non debt: scope estensione
//     futuro è cap.8 generale §6.4.1/§6.4.2 con shifts non-locali).
//   - coefficienti razionali in var per parsing polynomial.
//   - dispatch interno a solve_risch_de_q per ciascuna equazione su Q(var).

#include "calculus_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

// Walker che sostituisce un sotto-espressione `target` con `replacement`
// ovunque in `e`.  Identica per spirito al `deep_replace_expr` di
// integrate_risch.cpp (HPP-007), qui replicata per evitare cross-TU
// coupling.  TODO future refactor: estrarre in shared ast_helpers.
[[nodiscard]] bool deep_struct_equal(ExprPtr a, ExprPtr b);
[[nodiscard]] bool deep_struct_equal(ExprPtr a, ExprPtr b) {
    if (a == b) return true;
    if (!a || !b) return false;
    if (const auto* sa = expr_cast<Symbol>(a))
        if (const auto* sb = expr_cast<Symbol>(b)) return sa->name == sb->name;
    if (const auto* ia = expr_cast<IntegerLit>(a))
        if (const auto* ib = expr_cast<IntegerLit>(b)) return ia->value == ib->value;
    if (const auto* ra = expr_cast<RationalLit>(a))
        if (const auto* rb = expr_cast<RationalLit>(b))
            return ra->numerator == rb->numerator && ra->denominator == rb->denominator;
    if (const auto* ca = expr_cast<Constant>(a))
        if (const auto* cb = expr_cast<Constant>(b)) return ca->value == cb->value;
    if (const auto* fa = expr_cast<FuncCall>(a))
        if (const auto* fb = expr_cast<FuncCall>(b)) {
            if (fa->func_id != fb->func_id || fa->args.size() != fb->args.size()) return false;
            for (std::size_t i = 0; i < fa->args.size(); ++i)
                if (!deep_struct_equal(fa->args[i], fb->args[i])) return false;
            return true;
        }
    if (const auto* ba = expr_cast<Binary>(a))
        if (const auto* bb = expr_cast<Binary>(b))
            return ba->op == bb->op &&
                   deep_struct_equal(ba->left, bb->left) &&
                   deep_struct_equal(ba->right, bb->right);
    if (const auto* ua = expr_cast<Unary>(a))
        if (const auto* ub = expr_cast<Unary>(b))
            return ua->op == ub->op && deep_struct_equal(ua->operand, ub->operand);
    if (const auto* su_a = expr_cast<Sum>(a))
        if (const auto* su_b = expr_cast<Sum>(b)) {
            if (su_a->terms.size() != su_b->terms.size()) return false;
            for (std::size_t i = 0; i < su_a->terms.size(); ++i)
                if (!deep_struct_equal(su_a->terms[i], su_b->terms[i])) return false;
            return true;
        }
    if (const auto* pa = expr_cast<Product>(a))
        if (const auto* pb = expr_cast<Product>(b)) {
            if (pa->factors.size() != pb->factors.size()) return false;
            for (std::size_t i = 0; i < pa->factors.size(); ++i)
                if (!deep_struct_equal(pa->factors[i], pb->factors[i])) return false;
            return true;
        }
    return false;
}

[[nodiscard]] ExprPtr deep_replace(
    ExprPtr e, ExprPtr target, ExprPtr replacement, AstArena& arena) {
    if (!e) return e;
    if (deep_struct_equal(e, target)) return replacement;
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        std::vector<ExprPtr> na;
        na.reserve(fc->args.size());
        for (ExprPtr a : fc->args) na.push_back(deep_replace(a, target, replacement, arena));
        return arena.make<FuncCall>(fc->func_id, std::move(na));
    }
    if (const auto* bin = expr_cast<Binary>(e))
        return arena.make<Binary>(bin->op,
            deep_replace(bin->left, target, replacement, arena),
            deep_replace(bin->right, target, replacement, arena));
    if (const auto* un = expr_cast<Unary>(e))
        return arena.make<Unary>(un->op,
            deep_replace(un->operand, target, replacement, arena));
    if (const auto* sum = expr_cast<Sum>(e)) {
        std::vector<ExprPtr> nt;
        nt.reserve(sum->terms.size());
        for (ExprPtr t : sum->terms) nt.push_back(deep_replace(t, target, replacement, arena));
        return arena.make<Sum>(std::move(nt));
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        std::vector<ExprPtr> nf;
        nf.reserve(prod->factors.size());
        for (ExprPtr f : prod->factors) nf.push_back(deep_replace(f, target, replacement, arena));
        return arena.make<Product>(std::move(nf));
    }
    return e;
}

}  // namespace

Result<ExprPtr> solve_risch_de_logarithmic_q(
    ExprPtr f_expr,
    ExprPtr g_expr,
    ExprPtr u_arg,
    const Symbol& var,
    symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();

    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<ExprPtr>(
            "calculus", "solve_risch_de_logarithmic_q",
            msg,
            cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
            "Logarithmic Risch DE: estendere a f polinomiale in θ (Bronstein "
            "§6.4.1 caso non-costante) o cap.8 generale",
            "F0.8");
    };

    // 1. Costruisci θ = log(u) e θ' = u'/u in k = Q(var).
    ExprPtr theta = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{u_arg});
    auto u_prime_res = diff(u_arg, var, 1U, ctx);
    if (u_prime_res.is_error()) return fail_unimpl("cannot differentiate u(x)");
    ExprPtr theta_prime = arena.make<Binary>(BinaryOp::Div, u_prime_res.value(), u_arg);
    auto theta_prime_simp = ctx.simplify(theta_prime);
    if (theta_prime_simp.is_error()) return fail_unimpl("cannot simplify θ'");
    theta_prime = theta_prime_simp.value();

    // 2. Sostituisci θ con simbolo fresco s in f, g per il parsing polinomiale.
    Symbol theta_sym = ctx.make_fresh_symbol("riscθ");
    ExprPtr theta_sym_e = arena.make<Symbol>(theta_sym);
    ExprPtr f_sub = deep_replace(f_expr, theta, theta_sym_e, arena);
    ExprPtr g_sub = deep_replace(g_expr, theta, theta_sym_e, arena);

    // 3. f deve essere costante in θ (sub-case Bronstein §6.4.1 base).  Se f
    //    contiene θ_sym dopo la sostituzione, è polinomiale in θ e cade fuori
    //    scope di questa routine (cap.8 caso generale, non implementato).
    if (depends_on(f_sub, theta_sym))
        return fail_unimpl("f is polynomial in θ; Bronstein §6.4.1 general case "
                           "not yet implemented");

    // 4. Parse g come polinomio in θ_sym; estrai coefficienti (ExprPtr in k).
    auto g_poly_res = algebra::parse_polynomial(g_sub, theta_sym, ctx);
    if (g_poly_res.is_error())
        return fail_unimpl("g is not polynomial in θ");
    const std::vector<ExprPtr>& g_coeffs_raw = g_poly_res.value().coefficients();
    std::vector<ExprPtr> g_coeffs(g_coeffs_raw.begin(), g_coeffs_raw.end());
    // Strip trailing zeros: in PolyExpr i trailing zeros sono possibili.
    while (!g_coeffs.empty()) {
        auto simp = ctx.simplify(g_coeffs.back());
        if (simp.is_ok()) {
            if (const auto* il = expr_cast<IntegerLit>(simp.value());
                il && il->value.is_zero()) {
                g_coeffs.pop_back();
                continue;
            }
        }
        break;
    }
    if (g_coeffs.empty()) {
        // g ≡ 0: y = 0 è soluzione banale.
        return ok(arena.make<IntegerLit>(BigInt(0)));
    }

    const std::size_t n = g_coeffs.size() - 1U;  // grado massimo di g (e quindi di y).

    // 5. Risolvi top-down per y_n, y_{n-1}, ..., y_0.
    //    Bronstein §6.4.1: per i da n a 0,
    //       y_i' + f · y_i  =  g_i − (i+1) · θ' · y_{i+1}
    //    dove y_{n+1} = 0.
    std::vector<ExprPtr> y_coeffs(n + 1U, arena.make<IntegerLit>(BigInt(0)));

    for (std::size_t step = 0; step <= n; ++step) {
        const std::size_t i = n - step;
        // rhs = g_i − (i+1) · θ' · y_{i+1}
        ExprPtr rhs = g_coeffs[i];
        if (i + 1U <= n) {
            ExprPtr coef = arena.make<IntegerLit>(
                BigInt(static_cast<std::int64_t>(i + 1U)));
            ExprPtr term = arena.make<Product>(
                std::vector<ExprPtr>{coef, theta_prime, y_coeffs[i + 1U]});
            rhs = arena.make<Binary>(BinaryOp::Sub, rhs, term);
        }
        // together + simplify: il simplifier alone non sempre cancella
        // (1/x)·x → 1 quando i fattori sono distribuiti in Product opaco.
        auto rhs_tog  = algebra::together(rhs, ctx);
        if (rhs_tog.is_ok()) rhs = rhs_tog.value();
        auto rhs_simp = ctx.simplify(rhs);
        if (rhs_simp.is_ok()) rhs = rhs_simp.value();

        auto y_i = solve_risch_de_q(f_sub, rhs, var, ctx);
        if (y_i.is_error()) return fail_unimpl(
            "Risch DE su Q(var) non risolto per qualche coefficiente y_i");
        y_coeffs[i] = y_i.value();
    }

    // 6. Ricostruisci y = Σ y_i · θ^i.
    std::vector<ExprPtr> terms;
    terms.reserve(n + 1U);
    for (std::size_t i = 0; i <= n; ++i) {
        auto y_simp = ctx.simplify(y_coeffs[i]);
        ExprPtr y_i = y_simp.is_ok() ? y_simp.value() : y_coeffs[i];
        if (const auto* il = expr_cast<IntegerLit>(y_i); il && il->value.is_zero()) continue;
        ExprPtr term;
        if (i == 0U) {
            term = y_i;
        } else if (i == 1U) {
            term = arena.make<Binary>(BinaryOp::Mul, y_i, theta);
        } else {
            ExprPtr theta_pow = arena.make<Binary>(BinaryOp::Pow, theta,
                arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(i))));
            term = arena.make<Binary>(BinaryOp::Mul, y_i, theta_pow);
        }
        terms.push_back(term);
    }
    if (terms.empty()) return ok(arena.make<IntegerLit>(BigInt(0)));
    if (terms.size() == 1U) {
        auto s = ctx.simplify(terms.front());
        return s.is_ok() ? s : ok(terms.front());
    }
    ExprPtr sum = arena.make<Sum>(std::move(terms));
    auto s = ctx.simplify(sum);
    return s.is_ok() ? s : ok(sum);
}

// ============================================================================
// Bronstein §6.4.2 — caso esponenziale θ = exp(u).
// ============================================================================
Result<ExprPtr> solve_risch_de_exponential_q(
    ExprPtr f_expr,
    ExprPtr g_expr,
    ExprPtr u_arg,
    const Symbol& var,
    symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();

    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<ExprPtr>(
            "calculus", "solve_risch_de_exponential_q",
            msg,
            cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
            "Exponential Risch DE: estendere a f polinomiale in θ o cap.8 generale",
            "F0.8");
    };

    // 1. θ = exp(u);  θ' = u' · θ ∈ k[θ].
    //    Tuttavia per i ≥ 0,  i·θ^(i-1)·θ' = i·u'·θ^i, quindi le equazioni
    //    per i diversi disaccoppiano (coefficiente di θ^i dipende solo da y_i).
    ExprPtr theta = arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{u_arg});
    auto u_prime_res = diff(u_arg, var, 1U, ctx);
    if (u_prime_res.is_error()) return fail_unimpl("cannot differentiate u(x)");
    ExprPtr u_prime = u_prime_res.value();
    auto u_prime_simp = ctx.simplify(u_prime);
    if (u_prime_simp.is_ok()) u_prime = u_prime_simp.value();

    // 2. Sostituisci θ con simbolo fresco per parsing polinomiale.
    Symbol theta_sym = ctx.make_fresh_symbol("riscτ");
    ExprPtr theta_sym_e = arena.make<Symbol>(theta_sym);
    ExprPtr f_sub = deep_replace(f_expr, theta, theta_sym_e, arena);
    ExprPtr g_sub = deep_replace(g_expr, theta, theta_sym_e, arena);

    // 3. f costante in θ richiesto (sub-case base §6.4.2).
    if (depends_on(f_sub, theta_sym))
        return fail_unimpl("f is polynomial in θ; Bronstein §6.4.2 general case "
                           "not yet implemented");

    // 4. Parse g come polinomio in θ.
    auto g_poly_res = algebra::parse_polynomial(g_sub, theta_sym, ctx);
    if (g_poly_res.is_error())
        return fail_unimpl("g is not polynomial in θ");
    const std::vector<ExprPtr>& g_coeffs_raw = g_poly_res.value().coefficients();
    std::vector<ExprPtr> g_coeffs(g_coeffs_raw.begin(), g_coeffs_raw.end());
    while (!g_coeffs.empty()) {
        auto simp = ctx.simplify(g_coeffs.back());
        if (simp.is_ok()) {
            if (const auto* il = expr_cast<IntegerLit>(simp.value());
                il && il->value.is_zero()) {
                g_coeffs.pop_back();
                continue;
            }
        }
        break;
    }
    if (g_coeffs.empty())
        return ok(arena.make<IntegerLit>(BigInt(0)));

    const std::size_t n = g_coeffs.size() - 1U;

    // 5. Per ciascun grado i (indipendentemente), risolvi
    //       y_i' + (i · u' + f) · y_i = g_i
    std::vector<ExprPtr> y_coeffs(n + 1U, arena.make<IntegerLit>(BigInt(0)));
    for (std::size_t i = 0; i <= n; ++i) {
        // f_eff_i = i · u' + f_sub
        ExprPtr f_eff;
        if (i == 0U) {
            f_eff = f_sub;
        } else {
            ExprPtr i_coef = arena.make<IntegerLit>(
                BigInt(static_cast<std::int64_t>(i)));
            ExprPtr term = arena.make<Binary>(BinaryOp::Mul, i_coef, u_prime);
            f_eff = arena.make<Binary>(BinaryOp::Add, term, f_sub);
        }
        auto f_eff_tog = algebra::together(f_eff, ctx);
        if (f_eff_tog.is_ok()) f_eff = f_eff_tog.value();
        auto f_eff_simp = ctx.simplify(f_eff);
        if (f_eff_simp.is_ok()) f_eff = f_eff_simp.value();

        auto y_i = solve_risch_de_q(f_eff, g_coeffs[i], var, ctx);
        if (y_i.is_error()) return fail_unimpl(
            "Risch DE non risolto per qualche coefficiente y_i (caso exp)");
        y_coeffs[i] = y_i.value();
    }

    // 6. Ricostruisci y = Σ y_i · θ^i.
    std::vector<ExprPtr> terms;
    terms.reserve(n + 1U);
    for (std::size_t i = 0; i <= n; ++i) {
        auto y_simp = ctx.simplify(y_coeffs[i]);
        ExprPtr y_i = y_simp.is_ok() ? y_simp.value() : y_coeffs[i];
        if (const auto* il = expr_cast<IntegerLit>(y_i); il && il->value.is_zero()) continue;
        ExprPtr term;
        if (i == 0U) {
            term = y_i;
        } else if (i == 1U) {
            term = arena.make<Binary>(BinaryOp::Mul, y_i, theta);
        } else {
            ExprPtr theta_pow = arena.make<Binary>(BinaryOp::Pow, theta,
                arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(i))));
            term = arena.make<Binary>(BinaryOp::Mul, y_i, theta_pow);
        }
        terms.push_back(term);
    }
    if (terms.empty()) return ok(arena.make<IntegerLit>(BigInt(0)));
    if (terms.size() == 1U) {
        auto s = ctx.simplify(terms.front());
        return s.is_ok() ? s : ok(terms.front());
    }
    ExprPtr sum = arena.make<Sum>(std::move(terms));
    auto s = ctx.simplify(sum);
    return s.is_ok() ? s : ok(sum);
}

}  // namespace cas::calculus
