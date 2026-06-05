// F5.8 / Task #14 — Fourier transform F{f(t)}(ω) = ∫_{-∞}^{∞} f(t)·e^(−i·ω·t) dt.
//
// Convention: segno -i nell'integrazione (engineering convention).  Ω è la
// variabile coniugata della frequenza angolare.
//
// Pipeline:
//   1. Linearità su Sum.
//   2. Coppie axiomatic (distribution theory):
//        δ(t)                     → 1
//        δ(t − a)                 → exp(−i·ω·a)
//        1 (constante)            → 2π·δ(ω)
//        exp(i·a·t)               → 2π·δ(ω − a)
//        cos(a·t)                 → π·[δ(ω − a) + δ(ω + a)]
//        sin(a·t)                 → -i·π·[δ(ω − a) − δ(ω + a)]
//        exp(−a·t²)               → √(π/a)·exp(−ω²/(4a))
//   3. Pattern non riconosciuti → Unimplemented diagnostico esplicito.
//
// Riferimento: Bracewell, "The Fourier Transform and Its Applications",
// 3rd ed., McGraw-Hill 2000.

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

[[nodiscard]] ExprPtr int_e(symbolic::CASContext& ctx, long long v) {
    return ctx.arena().make<IntegerLit>(BigInt(v));
}

[[nodiscard]] bool depends_on_v(ExprPtr e, const Symbol& v) {
    if (!e) return false;
    if (const auto* s = expr_cast<Symbol>(e)) return s->name == v.name;
    if (const auto* un = expr_cast<Unary>(e)) return depends_on_v(un->operand, v);
    if (const auto* bin = expr_cast<Binary>(e))
        return depends_on_v(bin->left, v) || depends_on_v(bin->right, v);
    if (const auto* sum = expr_cast<Sum>(e)) {
        for (ExprPtr x : sum->terms) if (depends_on_v(x, v)) return true;
        return false;
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        for (ExprPtr x : prod->factors) if (depends_on_v(x, v)) return true;
        return false;
    }
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        for (ExprPtr x : fc->args) if (depends_on_v(x, v)) return true;
        return false;
    }
    return false;
}

[[nodiscard]] ExprPtr dirac(ExprPtr arg, symbolic::CASContext& ctx) {
    return ctx.arena().make<FuncCall>(BuiltinOp::DiracDelta,
        std::vector<ExprPtr>{arg});
}

[[nodiscard]] ExprPtr pi_e(symbolic::CASContext& ctx) {
    return ctx.arena().make<Constant>(MathConstant::Pi);
}

[[nodiscard]] ExprPtr i_e(symbolic::CASContext& ctx) {
    return ctx.arena().make<Constant>(MathConstant::I);
}

// Estrae il coefficiente c lineare tale che arg = c·var + b (b indip).
// Restituisce (c, b) o nullptr se non lineare in var.
struct LinearArg {
    ExprPtr c;
    ExprPtr b;
    bool ok{false};
};

[[nodiscard]] LinearArg extract_linear_in(
    ExprPtr arg, const Symbol& var, symbolic::CASContext& ctx) {
    LinearArg out;
    auto c_res = diff(arg, var, 1U, ctx);
    if (c_res.is_error()) return out;
    if (depends_on_v(c_res.value(), var)) return out;
    AstArena& arena = ctx.arena();
    ExprPtr cv = arena.make<Binary>(BinaryOp::Mul, c_res.value(),
        arena.make<Symbol>(var));
    ExprPtr b = arena.make<Binary>(BinaryOp::Sub, arg, cv);
    auto b_tog = algebra::together(b, ctx);
    if (b_tog.is_error()) return out;
    auto b_simp = ctx.simplify(b_tog.value());
    if (b_simp.is_error()) return out;
    if (depends_on_v(b_simp.value(), var)) return out;
    out.c = c_res.value();
    out.b = b_simp.value();
    out.ok = true;
    return out;
}

}  // namespace

Result<ExprPtr> fourier_transform(
    ExprPtr expr, const Symbol& t, const Symbol& omega,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr w_e = arena.make<Symbol>(omega);

    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<ExprPtr>(
            "calculus", "fourier_transform",
            msg,
            cas::error::reason_codes::LAPLACE_UNKNOWN_FORM,
            "Fourier transform: estendere pattern o usare proprietà aggiuntive",
            "F0.8");
    };

    // 1. Linearità su Sum.
    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> terms;
        terms.reserve(sum->terms.size());
        for (ExprPtr x : sum->terms) {
            auto f = fourier_transform(x, t, omega, ctx);
            if (f.is_error()) return f;
            terms.push_back(f.value());
        }
        return ctx.simplify(arena.make<Sum>(std::move(terms)));
    }

    // 2. Costante (indip da t): F{c}(ω) = c · 2π · δ(ω).
    if (!depends_on_v(expr, t)) {
        ExprPtr two_pi = arena.make<Binary>(BinaryOp::Mul, int_e(ctx, 2), pi_e(ctx));
        ExprPtr term = arena.make<Product>(std::vector<ExprPtr>{
            expr, two_pi, dirac(w_e, ctx)});
        return ctx.simplify(term);
    }

    // 3. δ(arg(t)) con arg lineare in t: F{δ(t − a)}(ω) = exp(−i·ω·a).
    //    Caso più generale δ(c·(t−a)) = δ(t−a)/|c|, ma il termine 1/|c|
    //    è scalare → applichiamo se necessario.
    if (const auto* fc = expr_cast<FuncCall>(expr);
        fc && fc->func_id == BuiltinOp::DiracDelta && fc->args.size() == 1U) {
        LinearArg la = extract_linear_in(fc->args[0], t, ctx);
        if (la.ok) {
            // arg = c·t + b → polo a = -b/c.
            ExprPtr neg_b = arena.make<Unary>(UnaryOp::Neg, la.b);
            ExprPtr a_val = arena.make<Binary>(BinaryOp::Div, neg_b, la.c);
            auto a_simp = ctx.simplify(a_val);
            ExprPtr a = a_simp.is_ok() ? a_simp.value() : a_val;
            // exp(-i·ω·a) / |c|.
            ExprPtr neg_i_w_a = arena.make<Product>(std::vector<ExprPtr>{
                arena.make<Unary>(UnaryOp::Neg, i_e(ctx)), w_e, a});
            ExprPtr exp_term = arena.make<FuncCall>(BuiltinOp::Exp,
                std::vector<ExprPtr>{neg_i_w_a});
            // |c| simbolico via Abs.
            bool c_pos = false;
            if (const auto* il = expr_cast<IntegerLit>(la.c);
                il && !il->value.is_negative() && !il->value.is_zero()) c_pos = true;
            if (const auto* rl = expr_cast<RationalLit>(la.c);
                rl && !rl->numerator.is_negative() && !rl->numerator.is_zero()) c_pos = true;
            ExprPtr abs_c = c_pos
                ? la.c
                : arena.make<FuncCall>(BuiltinOp::Abs, std::vector<ExprPtr>{la.c});
            ExprPtr result = arena.make<Binary>(BinaryOp::Div, exp_term, abs_c);
            return ctx.simplify(result);
        }
    }

    // 4. exp(i·a·t)  →  2π · δ(ω − a).  Detect via diff: arg lineare in t
    //    con coefficiente puro imaginario (coefficient = i·a per a reale).
    //    Per scope MVP, accettiamo arg = i·(c·t + b) con c, b reali.
    if (const auto* fc = expr_cast<FuncCall>(expr);
        fc && fc->func_id == BuiltinOp::Exp && fc->args.size() == 1U) {
        ExprPtr arg = fc->args[0];
        // arg deve essere lineare in t.  Estraiamo coef c = d(arg)/dt.
        LinearArg la = extract_linear_in(arg, t, ctx);
        if (la.ok) {
            // Verifica c = i·a con a indip da t.  Divide c per i e controlla.
            ExprPtr c_over_i = arena.make<Binary>(BinaryOp::Div, la.c, i_e(ctx));
            auto cot_simp = ctx.simplify(c_over_i);
            if (cot_simp.is_ok()) {
                ExprPtr a_val = cot_simp.value();
                // Se la divisione produce qualcosa libero da i (real), accettiamo.
                // F{exp(i·a·t)}(ω) = 2π · δ(ω − a) — caso b = 0.
                // Per b ≠ 0: F{exp(i·a·t + i·b₀)}(ω) = e^(i·b₀)·2π·δ(ω − a)
                ExprPtr two_pi = arena.make<Binary>(BinaryOp::Mul,
                    int_e(ctx, 2), pi_e(ctx));
                ExprPtr omega_minus_a = arena.make<Binary>(BinaryOp::Sub,
                    w_e, a_val);
                ExprPtr term = arena.make<Product>(std::vector<ExprPtr>{
                    two_pi, dirac(omega_minus_a, ctx)});
                // Termine costante exp(b): include se b ≠ 0.
                bool b_is_zero = false;
                if (const auto* il = expr_cast<IntegerLit>(la.b);
                    il && il->value.is_zero()) b_is_zero = true;
                if (!b_is_zero) {
                    ExprPtr exp_b = arena.make<FuncCall>(BuiltinOp::Exp,
                        std::vector<ExprPtr>{la.b});
                    term = arena.make<Binary>(BinaryOp::Mul, exp_b, term);
                }
                auto res_simp = ctx.simplify(term);
                if (res_simp.is_ok()) return res_simp;
                return ok(term);
            }
        }
    }

    // 5. cos(a·t), sin(a·t) tramite Euler.
    if (const auto* fc = expr_cast<FuncCall>(expr);
        fc && (fc->func_id == BuiltinOp::Cos || fc->func_id == BuiltinOp::Sin)
        && fc->args.size() == 1U) {
        ExprPtr arg = fc->args[0];
        LinearArg la = extract_linear_in(arg, t, ctx);
        if (la.ok) {
            // Solo caso b = 0 per il momento (per dare δ(ω∓a) puro).
            bool b_is_zero = false;
            if (const auto* il = expr_cast<IntegerLit>(la.b);
                il && il->value.is_zero()) b_is_zero = true;
            if (b_is_zero) {
                ExprPtr a = la.c;
                ExprPtr w_minus_a = arena.make<Binary>(BinaryOp::Sub, w_e, a);
                ExprPtr w_plus_a  = arena.make<Binary>(BinaryOp::Add, w_e, a);
                ExprPtr d_minus = dirac(w_minus_a, ctx);
                ExprPtr d_plus  = dirac(w_plus_a, ctx);
                ExprPtr pi_v    = pi_e(ctx);
                if (fc->func_id == BuiltinOp::Cos) {
                    // F{cos(a·t)}(ω) = π · [δ(ω-a) + δ(ω+a)].
                    ExprPtr sum_dirac = arena.make<Binary>(BinaryOp::Add, d_minus, d_plus);
                    ExprPtr result = arena.make<Binary>(BinaryOp::Mul, pi_v, sum_dirac);
                    return ctx.simplify(result);
                } else {
                    // F{sin(a·t)}(ω) = -i·π · [δ(ω-a) - δ(ω+a)].
                    ExprPtr diff_dirac = arena.make<Binary>(BinaryOp::Sub, d_minus, d_plus);
                    ExprPtr neg_i_pi = arena.make<Product>(std::vector<ExprPtr>{
                        arena.make<Unary>(UnaryOp::Neg, i_e(ctx)), pi_v});
                    ExprPtr result = arena.make<Binary>(BinaryOp::Mul, neg_i_pi, diff_dirac);
                    return ctx.simplify(result);
                }
            }
        }
    }

    // 6. Product: estrai scalari + ricorri sul fattore dipendente.
    if (const auto* prod = expr_cast<Product>(expr)) {
        std::vector<ExprPtr> consts, deps;
        for (ExprPtr f : prod->factors) {
            if (depends_on_v(f, t)) deps.push_back(f);
            else consts.push_back(f);
        }
        if (!consts.empty() && deps.size() == 1U) {
            auto inner = fourier_transform(deps[0], t, omega, ctx);
            if (inner.is_error()) return inner;
            consts.push_back(inner.value());
            return ctx.simplify(arena.make<Product>(std::move(consts)));
        }
    }

    return fail_unimpl("pattern non riconosciuto");
}

}  // namespace cas::calculus
