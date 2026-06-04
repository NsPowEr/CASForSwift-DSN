// F5.1 / B9-Task#23 (partial) — Cap.9 structure theorem, caso log-factorization.
//
// Bronstein "Symbolic Integration I" cap.9 (structure theorem) determina
// indipendenza algebrica di generatori in torre trascendentale.  Un sub-case
// concreto e ben definito è la decomposizione
//
//     log(u(x)) = Σ_i m_i · log(p_i(x)) + log(c)
//
// dove u(x) = c · Π_i p_i(x)^m_i è la fattorizzazione in irreducibili di
// u ∈ Q[var].  Questa identità è una conseguenza diretta del structure
// theorem applicata ai logaritmi: log(u) appartiene al Q-span dei
// {log(p_i)} + log della costante, quindi è "algebricamente dipendente"
// dai logaritmi degli irreducibili.
//
// La routine `expand_log_args_via_factorization(expr, var, ctx)` walks
// l'espressione, riconosce ogni FuncCall(Ln, u) con u polinomio in var,
// e lo sostituisce con la decomposizione Σ m_i · log(p_i) + log(content).
// I generatori di log irreducibili condivisi si fondono in un singolo
// generatore, abilitando il fast-path cap.8 (wiring section 2b-bis) su
// integranda che altrimenti avrebbe generatori multipli.
//
// Esempio:
//   ∫ (log(x) + log(x²)) dx                  — generatori distinti, no wiring
//   →  expand: log(x) + 2·log(x)             — generatore unico log(x)
//   →  integrate via cap.8 logaritmico:  3·(x·log(x) - x)

#include "calculus_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] bool is_polynomial_in(ExprPtr e, const Symbol& var, symbolic::CASContext& ctx) {
    // Crude check: parse_polynomial succeeds.  Avoids cost of full
    // factorization when the argument is clearly not polynomial.
    auto r = algebra::univariate_coefficients(e, var, ctx);
    return r.is_ok();
}

[[nodiscard]] ExprPtr decompose_log_arg(
    ExprPtr u, const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    if (!u) return arena.make<IntegerLit>(BigInt(0));

    if (!is_polynomial_in(u, var, ctx)) {
        // Non polinomiale: lascia log(u) intatto.
        return arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{u});
    }

    auto fact_res = algebra::factor_over_integers(u, var, ctx);
    if (fact_res.is_error()) {
        return arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{u});
    }
    const auto& fact = fact_res.value();

    std::vector<ExprPtr> terms;
    terms.reserve(fact.factors.size() + 1U);

    // log(content) — solo se ≠ 1.
    bool content_is_one = false;
    if (const auto* il = expr_cast<IntegerLit>(fact.content);
        il && il->value == BigInt(1)) content_is_one = true;
    if (!content_is_one) {
        terms.push_back(arena.make<FuncCall>(BuiltinOp::Ln,
            std::vector<ExprPtr>{fact.content}));
    }

    for (const auto& pf : fact.factors) {
        ExprPtr log_p = arena.make<FuncCall>(BuiltinOp::Ln,
            std::vector<ExprPtr>{pf.factor});
        if (pf.multiplicity == 1U) {
            terms.push_back(log_p);
        } else {
            ExprPtr mult = arena.make<IntegerLit>(
                BigInt(static_cast<std::int64_t>(pf.multiplicity)));
            terms.push_back(arena.make<Binary>(BinaryOp::Mul, mult, log_p));
        }
    }

    if (terms.empty()) return arena.make<IntegerLit>(BigInt(0));
    if (terms.size() == 1U) return terms.front();
    return arena.make<Sum>(std::move(terms));
}

}  // namespace

ExprPtr expand_log_args_via_factorization(
    ExprPtr e, const Symbol& var, symbolic::CASContext& ctx) {
    if (!e) return e;
    AstArena& arena = ctx.arena();

    if (const auto* fc = expr_cast<FuncCall>(e)) {
        // Espandi args ricorsivamente prima di decidere.
        std::vector<ExprPtr> new_args;
        new_args.reserve(fc->args.size());
        for (ExprPtr a : fc->args)
            new_args.push_back(expand_log_args_via_factorization(a, var, ctx));
        if (fc->func_id == BuiltinOp::Ln && new_args.size() == 1U) {
            return decompose_log_arg(new_args[0], var, ctx);
        }
        return arena.make<FuncCall>(fc->func_id, std::move(new_args));
    }
    if (const auto* bin = expr_cast<Binary>(e))
        return arena.make<Binary>(bin->op,
            expand_log_args_via_factorization(bin->left,  var, ctx),
            expand_log_args_via_factorization(bin->right, var, ctx));
    if (const auto* un = expr_cast<Unary>(e))
        return arena.make<Unary>(un->op,
            expand_log_args_via_factorization(un->operand, var, ctx));
    if (const auto* sum = expr_cast<Sum>(e)) {
        std::vector<ExprPtr> nt;
        nt.reserve(sum->terms.size());
        for (ExprPtr t : sum->terms)
            nt.push_back(expand_log_args_via_factorization(t, var, ctx));
        return arena.make<Sum>(std::move(nt));
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        std::vector<ExprPtr> nf;
        nf.reserve(prod->factors.size());
        for (ExprPtr f : prod->factors)
            nf.push_back(expand_log_args_via_factorization(f, var, ctx));
        return arena.make<Product>(std::move(nf));
    }
    return e;
}

}  // namespace cas::calculus
