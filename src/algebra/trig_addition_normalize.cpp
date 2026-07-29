// A50 (identità 2) — formule di addizione trigonometriche/iperboliche.
//
// Bottom-up rewrite: `sin/cos/sinh/cosh(u)` con `u` un `Sum` che mescola
// termini SENZA simboli liberi (fase/traslazione costante, es. `1`, `-2`,
// `sqrt(2)`) e termini CON simboli, si espande via la formula di addizione
// esatta:
//
//   sin(P+C)  = sin(P)·cos(C) + cos(P)·sin(C)
//   cos(P+C)  = cos(P)·cos(C) − sin(P)·sin(C)
//   sinh(P+C) = sinh(P)·cosh(C) + cosh(P)·sinh(C)
//   cosh(P+C) = cosh(P)·cosh(C) + sinh(P)·sinh(C)
//
// dove P = somma dei termini con simboli, C = somma dei termini senza.
// Se `u` è interamente costante o interamente simbolico non c'è nulla da
// separare: nessuna riscrittura (altrimenti sin(x) diventerebbe
// sin(x)·cos(0)+cos(x)·sin(0), lavoro sprecato e churn di forma).
//
// Non è una canonicalizzazione globale (stessa cautela di A43/A44): vive
// solo nel confronto di `mathematically_equal`, applicata a entrambi i lati.
// Chiude `D(F)=f` sulla famiglia trigonometrica traslata (A50): sia
// `sin(x+c)/x` (l'argomento del trig ha fase) sia `sin(x)/(x-c)` (il POLO è
// traslato, e la verifica ricombina `cos(c)·sin(x-c)+sin(c)·cos(x-c)` — Si/Ci
// shiftate nella derivata di F — di nuovo in `sin(x)` via la stessa identità,
// applicata all'argomento `x-c` che appare nella derivata).

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "algebra_internal.hpp"

#include <string>
#include <vector>

namespace cas::algebra {

namespace {

// Nessun simbolo libero al suo interno: una fase/traslazione costante può
// essere un letterale, una costante nominata (pi, e, ...), o una loro
// combinazione — non solo un IntegerLit. Nodi non enumerati esplicitamente
// (Integral/Derivative/Limit/RootOf/Matrix/SeriesExp/Quantity/Null) sono
// trattati come NON costanti: conservativo, al più la riscrittura non parte,
// mai un falso "è costante".
[[nodiscard]] bool has_no_free_symbols(ExprPtr e) {
    if (!e) return true;
    switch (e->kind) {
        case ExprKind::IntegerLit:
        case ExprKind::RationalLit:
        case ExprKind::DecimalLit:
        case ExprKind::Constant:
        case ExprKind::ComplexLit:
            return true;
        case ExprKind::Symbol:
            return false;
        case ExprKind::Unary:
            return has_no_free_symbols(expr_cast<Unary>(e)->operand);
        case ExprKind::Binary: {
            const auto* b = expr_cast<Binary>(e);
            return has_no_free_symbols(b->left) && has_no_free_symbols(b->right);
        }
        case ExprKind::Sum: {
            for (auto t : expr_cast<Sum>(e)->terms) if (!has_no_free_symbols(t)) return false;
            return true;
        }
        case ExprKind::Product: {
            for (auto f : expr_cast<Product>(e)->factors) if (!has_no_free_symbols(f)) return false;
            return true;
        }
        case ExprKind::FuncCall: {
            for (auto a : expr_cast<FuncCall>(e)->args) if (!has_no_free_symbols(a)) return false;
            return true;
        }
        default:
            return false;
    }
}

[[nodiscard]] ExprPtr make_call(AstArena& arena, const char* name, ExprPtr arg) {
    std::vector<ExprPtr> args{arg};
    return arena.make<FuncCall>(std::string(name), std::move(args));
}

// Costruisce `name(c)` canonicalizzando il SEGNO di `c` per parità, così che
// `cos(-2)` e `cos(2)` diventino lo stesso nodo strutturale (e non due atomi
// distinti che nessuna regola confronta) — necessario perché un letterale
// negativo è un `IntegerLit`/`RationalLit` col valore già negativo, non un
// `Unary(Neg, ...)`: le regole di parità esistenti (`builtin_rewrite_trig.cpp`)
// pattern-matchano sulla forma sintattica `Neg`, che un letterale negativo
// non ha mai. Riconosce letterali negativi ed espressioni già in forma
// `Unary(Neg, ...)`; qualunque altra forma passa invariata (nessun rischio di
// scorrettezza, solo un mancato collasso in più).
[[nodiscard]] ExprPtr phase_call(AstArena& arena, const char* name, bool odd_parity, ExprPtr c) {
    ExprPtr magnitude = c;
    bool negate = false;
    if (const auto* il = expr_cast<IntegerLit>(c)) {
        if (il->value.is_negative()) {
            magnitude = arena.make<IntegerLit>(BigInt(0) - il->value);
            negate = true;
        }
    } else if (const auto* rl = expr_cast<RationalLit>(c)) {
        if (rl->numerator.is_negative()) {
            magnitude = arena.make<RationalLit>(BigInt(0) - rl->numerator, rl->denominator);
            negate = true;
        }
    } else if (const auto* un = expr_cast<Unary>(c); un && un->op == UnaryOp::Neg) {
        magnitude = un->operand;
        negate = true;
    }
    ExprPtr call = make_call(arena, name, magnitude);
    return (negate && odd_parity) ? static_cast<ExprPtr>(arena.make<Unary>(UnaryOp::Neg, call)) : call;
}

[[nodiscard]] ExprPtr sum_of(AstArena& arena, std::vector<ExprPtr> terms) {
    return terms.size() == 1U ? terms.front()
                               : static_cast<ExprPtr>(arena.make<Sum>(std::move(terms)));
}

// Tenta la separazione P (con simboli) + C (senza) dell'argomento `u` di una
// FuncCall trig/iperbolica. Restituisce nullopt se `u` non è un Sum o se una
// delle due parti è vuota (nulla da separare).
struct PhaseSplit { ExprPtr p; ExprPtr c; };

[[nodiscard]] std::optional<PhaseSplit> split_constant_phase(ExprPtr u, AstArena& arena) {
    const auto* sum = expr_cast<Sum>(u);
    if (!sum) return std::nullopt;
    std::vector<ExprPtr> var_terms, const_terms;
    for (auto t : sum->terms) {
        (has_no_free_symbols(t) ? const_terms : var_terms).push_back(t);
    }
    if (var_terms.empty() || const_terms.empty()) return std::nullopt;
    return PhaseSplit{sum_of(arena, std::move(var_terms)), sum_of(arena, std::move(const_terms))};
}

}  // namespace

ExprPtr trig_addition_normalize(ExprPtr expr, AstArena& arena) {
    if (!expr) return expr;
    switch (expr->kind) {
        case ExprKind::FuncCall: {
            const auto* fc = expr_cast<FuncCall>(expr);
            if (fc == nullptr) return expr;
            if (fc->args.size() != 1U) {
                std::vector<ExprPtr> new_args;
                new_args.reserve(fc->args.size());
                bool changed = false;
                for (auto& a : fc->args) {
                    auto na = trig_addition_normalize(a, arena);
                    if (na.get() != a.get()) changed = true;
                    new_args.push_back(na);
                }
                if (!changed) return expr;
                return arena.make<FuncCall>(fc->name, std::move(new_args));
            }

            ExprPtr u = trig_addition_normalize(fc->args[0], arena);
            const std::string& nm = fc->name;
            const bool is_sin = (nm == "sin");
            const bool is_cos = (nm == "cos");
            const bool is_sinh = (nm == "sinh");
            const bool is_cosh = (nm == "cosh");
            if (!is_sin && !is_cos && !is_sinh && !is_cosh) {
                if (u.get() == fc->args[0].get()) return expr;
                return make_call(arena, nm.c_str(), u);
            }

            auto split = split_constant_phase(u, arena);
            if (!split.has_value()) {
                if (u.get() == fc->args[0].get()) return expr;
                return make_call(arena, nm.c_str(), u);
            }
            ExprPtr p = split->p;
            ExprPtr c = split->c;

            if (is_sin) {
                // sin(P)cos(C) + cos(P)sin(C)
                return arena.make<Sum>(std::vector<ExprPtr>{
                    arena.make<Product>(std::vector<ExprPtr>{make_call(arena, "sin", p), phase_call(arena, "cos", false, c)}),
                    arena.make<Product>(std::vector<ExprPtr>{make_call(arena, "cos", p), phase_call(arena, "sin", true, c)}),
                });
            }
            if (is_cos) {
                // cos(P)cos(C) - sin(P)sin(C)
                ExprPtr cross = arena.make<Product>(std::vector<ExprPtr>{make_call(arena, "sin", p), phase_call(arena, "sin", true, c)});
                return arena.make<Sum>(std::vector<ExprPtr>{
                    arena.make<Product>(std::vector<ExprPtr>{make_call(arena, "cos", p), phase_call(arena, "cos", false, c)}),
                    arena.make<Unary>(UnaryOp::Neg, cross),
                });
            }
            if (is_sinh) {
                // sinh(P)cosh(C) + cosh(P)sinh(C)
                return arena.make<Sum>(std::vector<ExprPtr>{
                    arena.make<Product>(std::vector<ExprPtr>{make_call(arena, "sinh", p), phase_call(arena, "cosh", false, c)}),
                    arena.make<Product>(std::vector<ExprPtr>{make_call(arena, "cosh", p), phase_call(arena, "sinh", true, c)}),
                });
            }
            // cosh: cosh(P)cosh(C) + sinh(P)sinh(C)
            return arena.make<Sum>(std::vector<ExprPtr>{
                arena.make<Product>(std::vector<ExprPtr>{make_call(arena, "cosh", p), phase_call(arena, "cosh", false, c)}),
                arena.make<Product>(std::vector<ExprPtr>{make_call(arena, "sinh", p), phase_call(arena, "sinh", true, c)}),
            });
        }
        case ExprKind::Unary: {
            const auto* un = expr_cast<Unary>(expr);
            if (un == nullptr) return expr;
            auto inner = trig_addition_normalize(un->operand, arena);
            if (inner.get() == un->operand.get()) return expr;
            return arena.make<Unary>(un->op, inner);
        }
        case ExprKind::Binary: {
            const auto* bi = expr_cast<Binary>(expr);
            if (bi == nullptr) return expr;
            auto l = trig_addition_normalize(bi->left, arena);
            auto r = trig_addition_normalize(bi->right, arena);
            if (l.get() == bi->left.get() && r.get() == bi->right.get()) return expr;
            return arena.make<Binary>(bi->op, l, r);
        }
        case ExprKind::Sum: {
            const auto* su = expr_cast<Sum>(expr);
            if (su == nullptr) return expr;
            std::vector<ExprPtr> nt;
            nt.reserve(su->terms.size());
            bool changed = false;
            for (auto& t : su->terms) {
                auto nx = trig_addition_normalize(t, arena);
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
                auto nx = trig_addition_normalize(f, arena);
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
