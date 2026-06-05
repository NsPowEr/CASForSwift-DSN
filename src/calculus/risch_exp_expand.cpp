// F5.1 / B9-Task#23 — Cap.9 structure theorem, sub-case exp-decomposition.
//
// Bronstein "Symbolic Integration I" cap.9.  Per torre trascendentale
// K(t_1,...,t_n) con t_i = exp(u_i), l'indipendenza algebrica richiede che
// gli {u_i} siano Z-linearmente indipendenti modulo costanti.  In caso
// contrario, exp(u_i) e exp(c·u_j) sono algebricamente dipendenti
// (exp(c·u) = exp(u)^c per c razionale ⇒ algebrico se c ∈ Z, altrimenti
// dipendenza via radicale).
//
// Questa routine implementa due decomposizioni canonical che catturano i
// casi pratici più frequenti:
//
//   (1) exp(Sum(a_1, ..., a_k)) → Product(exp(a_1), ..., exp(a_k))
//       Identità formale del derivation: D(exp(Σa_i)) = (Σa_i)' · exp(Σa_i)
//       = Σa_i' · Π exp(a_i) = D(Π exp(a_i)).  Le due forme producono lo
//       stesso elemento del campo differenziale; preferiamo la forma
//       prodotto perché esibisce esplicitamente i generatori indipendenti.
//
//   (2) exp(n · u) → (exp(u))^n  per n ∈ Z literal
//       Caso speciale di (1) iterato; estrae la molteplicità Z dal generatore
//       primitivo.  Riconosce Product(IntegerLit, ...) come argomento di exp
//       e fattorizza il coefficiente intero al di fuori della trascendentale.
//
// Esempio integrale:
//   ∫ (exp(x) + exp(2x)) dx     — generatori distinti exp(x), exp(2x)
//   → expand: exp(x) + exp(x)²  — generatore primitivo unico exp(x)
//   → integrate via cap.8 esponenziale (post-substitution polynomial in t).

#include "calculus_internal.hpp"

#include "cas/symbolic.hpp"

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

// Estrae da `arg` la struttura  n · rest  dove n ∈ Z letterale (positivo o
// negativo) e `rest` è il residuo.  Ritorna std::nullopt se non è un Product
// con esattamente un fattore IntegerLit non unitario.  Match conservativo:
// non gestiamo razionali (Q-coefficienti) né scomposizione di Sum nidificati
// nell'argomento del Product.
struct IntScalarExtract {
    BigInt n;
    ExprPtr rest;
};

[[nodiscard]] std::optional<IntScalarExtract> extract_integer_scalar(
    ExprPtr arg, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    const auto* prod = expr_cast<Product>(arg);
    if (!prod) return std::nullopt;

    std::vector<ExprPtr> non_int_factors;
    non_int_factors.reserve(prod->factors.size());
    BigInt scalar(1);
    bool found_scalar = false;
    for (ExprPtr f : prod->factors) {
        if (const auto* il = expr_cast<IntegerLit>(f)) {
            scalar = scalar * il->value;
            found_scalar = true;
            continue;
        }
        non_int_factors.push_back(f);
    }
    if (!found_scalar) return std::nullopt;
    if (scalar == BigInt(1) || scalar == BigInt(-1)) {
        // Coefficient ±1: nessuna estrazione utile (lasciamo segno al simplifier).
        return std::nullopt;
    }
    if (non_int_factors.empty()) {
        // exp(n) puro: nessuna decomposizione moltiplicativa utile (è costante).
        return std::nullopt;
    }

    ExprPtr rest;
    if (non_int_factors.size() == 1U) rest = non_int_factors.front();
    else rest = arena.make<Product>(std::move(non_int_factors));
    return IntScalarExtract{ .n = scalar, .rest = rest };
}

[[nodiscard]] ExprPtr decompose_exp_arg(
    ExprPtr arg, const Symbol& /*var*/, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    if (!arg) return arena.make<IntegerLit>(BigInt(1));

    // Caso (1): exp(Sum(a_1,...,a_k)) → Π exp(a_i).
    if (const auto* sum = expr_cast<Sum>(arg)) {
        std::vector<ExprPtr> factors;
        factors.reserve(sum->terms.size());
        for (ExprPtr t : sum->terms) {
            // Recurse su ciascun termine: exp(2x + ln(y)) → exp(x)²·exp(ln(y))
            // (poi il simplifier collassa exp(ln(y)) → y).
            ExprPtr piece = decompose_exp_arg(t, /*var*/Symbol{""}, ctx);
            factors.push_back(piece);
        }
        if (factors.empty()) return arena.make<IntegerLit>(BigInt(1));
        if (factors.size() == 1U) return factors.front();
        return arena.make<Product>(std::move(factors));
    }

    // Caso (2): exp(n·u) → exp(u)^n  per n ∈ Z letterale, |n| ≥ 2.
    if (auto extracted = extract_integer_scalar(arg, ctx); extracted.has_value()) {
        ExprPtr base = arena.make<FuncCall>(BuiltinOp::Exp,
            std::vector<ExprPtr>{extracted->rest});
        ExprPtr exponent = arena.make<IntegerLit>(extracted->n);
        return arena.make<Binary>(BinaryOp::Pow, base, exponent);
    }

    // Default: lascia exp(arg) come tale.
    return arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{arg});
}

}  // namespace

ExprPtr expand_exp_args_via_decomposition(
    ExprPtr e, const Symbol& var, symbolic::CASContext& ctx) {
    if (!e) return e;
    AstArena& arena = ctx.arena();

    if (const auto* fc = expr_cast<FuncCall>(e)) {
        std::vector<ExprPtr> new_args;
        new_args.reserve(fc->args.size());
        for (ExprPtr a : fc->args)
            new_args.push_back(expand_exp_args_via_decomposition(a, var, ctx));
        if (fc->func_id == BuiltinOp::Exp && new_args.size() == 1U) {
            return decompose_exp_arg(new_args[0], var, ctx);
        }
        return arena.make<FuncCall>(fc->func_id, std::move(new_args));
    }
    if (const auto* bin = expr_cast<Binary>(e))
        return arena.make<Binary>(bin->op,
            expand_exp_args_via_decomposition(bin->left,  var, ctx),
            expand_exp_args_via_decomposition(bin->right, var, ctx));
    if (const auto* un = expr_cast<Unary>(e))
        return arena.make<Unary>(un->op,
            expand_exp_args_via_decomposition(un->operand, var, ctx));
    if (const auto* sum = expr_cast<Sum>(e)) {
        std::vector<ExprPtr> nt;
        nt.reserve(sum->terms.size());
        for (ExprPtr t : sum->terms)
            nt.push_back(expand_exp_args_via_decomposition(t, var, ctx));
        return arena.make<Sum>(std::move(nt));
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        std::vector<ExprPtr> nf;
        nf.reserve(prod->factors.size());
        for (ExprPtr f : prod->factors)
            nf.push_back(expand_exp_args_via_decomposition(f, var, ctx));
        return arena.make<Product>(std::move(nf));
    }
    return e;
}

}  // namespace cas::calculus
