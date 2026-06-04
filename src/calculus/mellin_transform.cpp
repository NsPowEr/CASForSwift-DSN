// F5.8 / Task #15 — Mellin transform M{f(t)}(s) = ∫_0^∞ t^(s-1) · f(t) dt.
//
// Convention: M-transform unilaterale su [0, ∞).  Strip of convergence non
// è tracciato esplicitamente nei risultati simbolici (chi consuma il
// risultato deve interpretarlo nello strip corretto).
//
// Pipeline:
//   1. Linearità su Sum.
//   2. Estrazione scalari da Product.
//   3. Coppie axiomatic:
//        M{exp(−a·t)}(s)         = Γ(s)/a^s
//        M{t^a · exp(−b·t)}(s)   = Γ(s+a)/b^(s+a)
//        M{1/(1+t)^a}(s)         = Γ(s)·Γ(a−s)/Γ(a)
//   4. Proprietà:
//        Modulation:  M{t^a · f(t)}(s)  = F(s+a)
//        Scaling:     M{f(a·t)}(s)      = a^(−s) · F(s)
//   5. Unimplemented diagnostico per pattern non riconosciuti.
//
// Riferimento: Bateman Manuscript Project, "Tables of Integral Transforms"
// vol.I (McGraw-Hill 1954) cap.6.

#include "calculus_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] ExprPtr make_int_e(symbolic::CASContext& ctx, long long v) {
    return ctx.arena().make<IntegerLit>(BigInt(v));
}

[[nodiscard]] bool is_sym_eq(ExprPtr e, const Symbol& s) {
    const auto* sy = expr_cast<Symbol>(e);
    return sy && sy->name == s.name;
}

[[nodiscard]] bool depends_on_t(ExprPtr e, const Symbol& v) {
    if (!e) return false;
    if (const auto* s = expr_cast<Symbol>(e)) return s->name == v.name;
    if (const auto* un = expr_cast<Unary>(e)) return depends_on_t(un->operand, v);
    if (const auto* bin = expr_cast<Binary>(e))
        return depends_on_t(bin->left, v) || depends_on_t(bin->right, v);
    if (const auto* sum = expr_cast<Sum>(e)) {
        for (ExprPtr x : sum->terms) if (depends_on_t(x, v)) return true;
        return false;
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        for (ExprPtr x : prod->factors) if (depends_on_t(x, v)) return true;
        return false;
    }
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        for (ExprPtr x : fc->args) if (depends_on_t(x, v)) return true;
        return false;
    }
    return false;
}

// Γ(arg) come ExprPtr.
[[nodiscard]] ExprPtr gamma_call(ExprPtr arg, symbolic::CASContext& ctx) {
    return ctx.arena().make<FuncCall>(BuiltinOp::Gamma,
        std::vector<ExprPtr>{arg});
}

// Cerca pattern exp(-a·t) con a indip da t e a costante coefficiente di t.
// Strategia robusta indipendente dalla forma sintattica concreta:
//   1. arg = expr.arg.  Calcola c = d(arg)/dt.
//   2. Se c dipende ancora da t, arg non è lineare → fallisce.
//   3. Verifica arg = c·t (residuo arg − c·t == 0).
//   4. Accetta exp(c·t) con a = −c.  Funziona per qualsiasi normalizzazione
//      del segno (Unary(Neg, Mul(3,t)), Mul(−3,t), Mul(3,Neg(t)), ecc.).
[[nodiscard]] ExprPtr extract_a_in_exp_minus_at(
    ExprPtr expr, const Symbol& t, symbolic::CASContext& ctx) {
    const auto* fc = expr_cast<FuncCall>(expr);
    if (!fc || fc->func_id != BuiltinOp::Exp || fc->args.size() != 1U) return nullptr;
    ExprPtr arg = fc->args[0];

    auto c_res = diff(arg, t, 1U, ctx);
    if (c_res.is_error()) return nullptr;
    ExprPtr c = c_res.value();
    if (depends_on_t(c, t)) return nullptr;  // arg non è lineare in t

    // Verifica residuo arg − c·t == 0.
    AstArena& arena = ctx.arena();
    ExprPtr ct = arena.make<Binary>(BinaryOp::Mul, c, arena.make<Symbol>(t));
    ExprPtr resid = arena.make<Binary>(BinaryOp::Sub, arg, ct);
    auto resid_tog = algebra::together(resid, ctx);
    if (resid_tog.is_error()) return nullptr;
    auto resid_simp = ctx.simplify(resid_tog.value());
    if (resid_simp.is_error()) return nullptr;
    const auto* zero_il = expr_cast<IntegerLit>(resid_simp.value());
    if (!zero_il || !zero_il->value.is_zero()) return nullptr;

    // a = -c, semplifica.
    ExprPtr neg_c = arena.make<Unary>(UnaryOp::Neg, c);
    auto a_simp = ctx.simplify(neg_c);
    return a_simp.is_ok() ? a_simp.value() : neg_c;
}

}  // namespace

Result<ExprPtr> mellin_transform(
    ExprPtr expr, const Symbol& t, const Symbol& s,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr s_e = arena.make<Symbol>(s);

    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<ExprPtr>(
            "calculus", "mellin_transform",
            msg,
            cas::error::reason_codes::LAPLACE_UNKNOWN_FORM,
            "Mellin transform: estendere coppia pattern o utilizzare proprietà",
            "F0.8");
    };

    // 1. Linearità su Sum.
    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> terms;
        terms.reserve(sum->terms.size());
        for (ExprPtr x : sum->terms) {
            auto m = mellin_transform(x, t, s, ctx);
            if (m.is_error()) return m;
            terms.push_back(m.value());
        }
        return ctx.simplify(arena.make<Sum>(std::move(terms)));
    }

    // 2. Coppia: exp(-a·t) → Γ(s) / a^s.
    if (auto a = extract_a_in_exp_minus_at(expr, t, ctx)) {
        ExprPtr num = gamma_call(s_e, ctx);
        ExprPtr den = arena.make<Binary>(BinaryOp::Pow, a, s_e);
        return ctx.simplify(arena.make<Binary>(BinaryOp::Div, num, den));
    }

    // 3. Product: scalari + (eventualmente) t^a · exp(-b·t) o forme analoghe.
    if (const auto* prod = expr_cast<Product>(expr)) {
        std::vector<ExprPtr> scalars;
        std::vector<ExprPtr> t_dep;
        ExprPtr t_power = nullptr;  // a in t^a
        for (ExprPtr f : prod->factors) {
            if (!depends_on_t(f, t)) { scalars.push_back(f); continue; }
            // t^a con a indipendente da t?
            if (const auto* pw = expr_cast<Binary>(f);
                pw && pw->op == BinaryOp::Pow
                && is_sym_eq(pw->left, t)
                && !depends_on_t(pw->right, t)) {
                if (!t_power) {
                    t_power = pw->right;
                    continue;
                }
                // Già un t^a → cumula sommando esponenti.
                t_power = arena.make<Binary>(BinaryOp::Add, t_power, pw->right);
                continue;
            }
            // t solo (≡ t^1).
            if (is_sym_eq(f, t)) {
                if (!t_power) t_power = make_int_e(ctx, 1);
                else t_power = arena.make<Binary>(BinaryOp::Add, t_power, make_int_e(ctx, 1));
                continue;
            }
            t_dep.push_back(f);
        }

        // Caso (a): solo t^a (e scalari) ⇒ ∫_0^∞ t^(s+a-1) dt — diverge.
        //          Non c'è chiusura simbolica utile; ricade Unimplemented.
        if (t_dep.empty() && t_power) {
            return fail_unimpl("Mellin: t^a alone diverges on (0,∞) — "
                               "strip of convergence vuoto");
        }

        // Caso (b): t_dep contiene esattamente un exp(-b·t).
        if (t_dep.size() == 1U) {
            if (auto b = extract_a_in_exp_minus_at(t_dep[0], t, ctx)) {
                // Modulation: M{t^a · exp(-b·t)}(s) = Γ(s+a) / b^(s+a).
                ExprPtr s_plus_a = t_power
                    ? arena.make<Binary>(BinaryOp::Add, s_e, t_power)
                    : s_e;
                ExprPtr num = gamma_call(s_plus_a, ctx);
                ExprPtr den = arena.make<Binary>(BinaryOp::Pow, b, s_plus_a);
                ExprPtr res = arena.make<Binary>(BinaryOp::Div, num, den);
                if (!scalars.empty()) {
                    std::vector<ExprPtr> all = scalars;
                    all.push_back(res);
                    res = arena.make<Product>(std::move(all));
                }
                return ctx.simplify(res);
            }
        }

        // Caso (c): scalare + single t-dep generico → estrai scalare e ricorri.
        if (!scalars.empty() && t_dep.size() == 1U && !t_power) {
            auto inner = mellin_transform(t_dep[0], t, s, ctx);
            if (inner.is_error()) return inner;
            std::vector<ExprPtr> all = scalars;
            all.push_back(inner.value());
            return ctx.simplify(arena.make<Product>(std::move(all)));
        }
    }

    // 4. (1+t)^(-a) (modulus base form) per a > 0:
    //     M{(1+t)^(-a)}(s) = Γ(s)·Γ(a-s) / Γ(a).
    if (const auto* pw = expr_cast<Binary>(expr); pw && pw->op == BinaryOp::Pow) {
        const auto* base_bin = expr_cast<Binary>(pw->left);
        if (base_bin && base_bin->op == BinaryOp::Add
            && !depends_on_t(pw->right, t)) {
            // base = (1 + t)?
            ExprPtr maybe_one = nullptr;
            ExprPtr maybe_t = nullptr;
            if (is_sym_eq(base_bin->left, t)) { maybe_t = base_bin->left; maybe_one = base_bin->right; }
            else if (is_sym_eq(base_bin->right, t)) { maybe_t = base_bin->right; maybe_one = base_bin->left; }
            if (maybe_t) {
                const auto* il = expr_cast<IntegerLit>(maybe_one);
                if (il && il->value == BigInt(1)) {
                    // exponent must be -a; check if Neg or RationalLit negative.
                    ExprPtr exponent = pw->right;
                    ExprPtr a_param = nullptr;
                    if (const auto* un = expr_cast<Unary>(exponent); un && un->op == UnaryOp::Neg) {
                        a_param = un->operand;
                    }
                    // Altrimenti potrebbe essere già negativo come RationalLit ma evitiamo
                    // ambiguità: pattern strict Neg(positive).
                    if (a_param) {
                        // M{(1+t)^(-a)}(s) = Γ(s)·Γ(a-s)/Γ(a).
                        ExprPtr a_minus_s = arena.make<Binary>(BinaryOp::Sub, a_param, s_e);
                        ExprPtr num = arena.make<Binary>(BinaryOp::Mul,
                            gamma_call(s_e, ctx), gamma_call(a_minus_s, ctx));
                        ExprPtr den = gamma_call(a_param, ctx);
                        return ctx.simplify(arena.make<Binary>(BinaryOp::Div, num, den));
                    }
                }
            }
        }
    }

    return fail_unimpl("pattern non riconosciuto");
}

}  // namespace cas::calculus
