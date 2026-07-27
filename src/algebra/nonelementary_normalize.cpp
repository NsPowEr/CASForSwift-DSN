// A43 §4 — Identità di riduzione per la famiglia delle antiderivate non
// elementari (Ei, Si, Ci, li, Li₂, erfi, Shi, Chi).
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Nonelementary_Antiderivatives.md
// (ogni identità verificata numericamente con mpmath a 30 cifre, §7).
//
// Il problema che risolve: la stessa quantità ha più ortografie legittime, e
// senza un rappresentante unico il confronto strutturale si biforca (esattamente
// come `x!` vs `Γ(x+1)` in A44 — 17 WARN e un golden FAIL). Le riscritture:
//
//   li(u)   → Ei(ln u)                    (spec §4, verificata)
//   Shi(u)  → (Ei(u) − Ei(−u))/2          (spec §4, verificata)
//   Chi(u)  → (Ei(u) + Ei(−u))/2          (spec §4, verificata)
//   erfi(u) → −i·erf(i·u)                 (spec §4, verificata)
//
// Rappresentante scelto: `Ei` per la famiglia degli integrali esponenziali —
// è il generatore (li, Shi, Chi si esprimono tramite Ei, non viceversa) — ed
// `erf` per la coppia erf/erfi, perché è il nodo che il motore già possiede e
// la forma in cui Maxima emette il risultato (`erf(%i*x)`).
//
// `Si` e `Ci` NON vengono riscritti: le identità che li legano a Ei passano per
// `Ci(x) ± i·Si(x) = Ei(±ix)` (spec §4, ultima riga) e introdurrebbero una
// unità immaginaria spuria in un risultato reale.
//
// Collocazione: applicato ai due lati di `mathematically_equal`, accanto a
// `hyperbolic_normalize` e `factorial_gamma_normalize` — NON una
// canonicalizzazione globale in `simplify`, che sarebbe churn strutturale di
// massa (stessa valutazione motivata in A44 e ribadita dalla spec §4).
//
// Structural sharing preservato: se in un sottoalbero non scatta nessuna
// riscrittura viene restituito il puntatore originale.

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/builtin_functions.hpp"
#include "algebra_internal.hpp"

#include <string>
#include <utility>
#include <vector>

namespace cas::algebra {

namespace {

[[nodiscard]] ExprPtr make_ei(ExprPtr u, AstArena& arena) {
    std::vector<ExprPtr> args;
    args.push_back(u);
    return arena.make<FuncCall>(BuiltinOp::ExpIntegralEi, std::move(args));
}

// (Ei(u) + sign·Ei(−u))/2 — la forma comune a Shi (sign = −1) e Chi (sign = +1).
[[nodiscard]] ExprPtr make_ei_half_combination(ExprPtr u, bool negate_reflected, AstArena& arena) {
    ExprPtr reflected = make_ei(arena.make<Unary>(UnaryOp::Neg, u), arena);
    if (negate_reflected) {
        reflected = arena.make<Unary>(UnaryOp::Neg, reflected);
    }
    ExprPtr sum = arena.make<Sum>(std::vector<ExprPtr>{make_ei(u, arena), reflected});
    return arena.make<Binary>(BinaryOp::Div, sum, arena.make<IntegerLit>(BigInt(2)));
}

}  // namespace

ExprPtr nonelementary_normalize(ExprPtr expr, AstArena& arena) {
    if (!expr) return expr;
    switch (expr->kind) {
        case ExprKind::FuncCall: {
            const auto* fc = expr_cast<FuncCall>(expr);
            if (fc == nullptr) return expr;
            std::vector<ExprPtr> new_args;
            new_args.reserve(fc->args.size());
            bool changed = false;
            for (auto& a : fc->args) {
                auto na = nonelementary_normalize(a, arena);
                if (na.get() != a.get()) changed = true;
                new_args.push_back(na);
            }
            if (new_args.size() == 1U) {
                ExprPtr u = new_args[0];
                switch (fc->func_id) {
                    case BuiltinOp::LogIntegral: {
                        std::vector<ExprPtr> ln_args;
                        ln_args.push_back(u);
                        return make_ei(arena.make<FuncCall>(BuiltinOp::Ln, std::move(ln_args)), arena);
                    }
                    case BuiltinOp::SinhIntegral:
                        return make_ei_half_combination(u, /*negate_reflected=*/true, arena);
                    case BuiltinOp::CoshIntegral:
                        return make_ei_half_combination(u, /*negate_reflected=*/false, arena);
                    case BuiltinOp::Erfi: {
                        // erfi(u) = −i·erf(i·u)
                        ExprPtr imaginary_unit = arena.make<Constant>(MathConstant::I);
                        std::vector<ExprPtr> erf_args;
                        erf_args.push_back(arena.make<Product>(std::vector<ExprPtr>{imaginary_unit, u}));
                        ExprPtr erf_call = arena.make<FuncCall>(BuiltinOp::Erf, std::move(erf_args));
                        return arena.make<Unary>(UnaryOp::Neg,
                            arena.make<Product>(std::vector<ExprPtr>{imaginary_unit, erf_call}));
                    }
                    default:
                        break;
                }
            }
            if (!changed) return expr;
            return arena.make<FuncCall>(fc->name, std::move(new_args));
        }
        case ExprKind::Unary: {
            const auto* un = expr_cast<Unary>(expr);
            if (un == nullptr) return expr;
            ExprPtr inner = nonelementary_normalize(un->operand, arena);
            if (inner.get() == un->operand.get()) return expr;
            return arena.make<Unary>(un->op, inner);
        }
        case ExprKind::Binary: {
            const auto* bin = expr_cast<Binary>(expr);
            if (bin == nullptr) return expr;
            ExprPtr l = nonelementary_normalize(bin->left, arena);
            ExprPtr r = nonelementary_normalize(bin->right, arena);
            if (l.get() == bin->left.get() && r.get() == bin->right.get()) return expr;
            return arena.make<Binary>(bin->op, l, r);
        }
        case ExprKind::Sum: {
            const auto* sm = expr_cast<Sum>(expr);
            if (sm == nullptr) return expr;
            std::vector<ExprPtr> nt;
            nt.reserve(sm->terms.size());
            bool changed = false;
            for (auto& t : sm->terms) {
                auto nx = nonelementary_normalize(t, arena);
                if (nx.get() != t.get()) changed = true;
                nt.push_back(nx);
            }
            if (!changed) return expr;
            return arena.make<Sum>(std::move(nt));
        }
        case ExprKind::Product: {
            const auto* pr = expr_cast<Product>(expr);
            if (pr == nullptr) return expr;
            std::vector<ExprPtr> nf;
            nf.reserve(pr->factors.size());
            bool changed = false;
            for (auto& f : pr->factors) {
                auto nx = nonelementary_normalize(f, arena);
                if (nx.get() != f.get()) changed = true;
                nf.push_back(nx);
            }
            if (!changed) return expr;
            return arena.make<Product>(std::move(nf));
        }
        default:
            return expr;
    }
}

}  // namespace cas::algebra
