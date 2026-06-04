// F5.8 / Task #17 — Inverse Laplace via Bronstein residue.
//
// Formula del residuo per la trasformata inversa di Laplace di F(s) razionale:
//
//     L⁻¹{F(s)}(t) = Σ_k Res_{s = p_k} [F(s) · e^(s·t)]
//
// dove {p_k} sono i poli di F(s) nel piano complesso.  Bronstein "Symbolic
// Integration I" §3.4–§3.6 fornisce l'algoritmo per calcolare i residui in
// modo algoritmico (non tabulare); questo modulo riusa la routine `residue`
// già presente nel sub-engine F5.6 (residue theorem driver).
//
// Pipeline:
//   1. Decompone F = N(s)/D(s) via algebra::apart_num_den.
//   2. Trova i poli risolvendo D(s) = 0 via algebra::solve_polynomial.
//   3. Per ciascun polo p calcola Res(F(s)·e^(s·t), s, p).
//   4. Somma i residui; semplifica.
//
// Differenze rispetto al pattern-based inverse_laplace_transform:
//   - Funziona per F(s) razionale con polos qualsiasi (anche multipli,
//     anche oltre il pattern table elementare).
//   - Polos complessi coniugati restano nella forma e^(p·t); il
//     riconoscimento di cos/sin tramite Euler è scope follow-up.
//   - Limiti: polos irrazionali (es. radici di polinomi di grado ≥ 3
//     senza Q-roots) ricadono su Unimplemented esplicito.

#include "calculus_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace cas::calculus {

Result<ExprPtr> inverse_laplace_residue_q(
    ExprPtr F, const Symbol& s, const Symbol& t,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<ExprPtr>(
            "calculus", "inverse_laplace_residue_q",
            msg,
            cas::error::reason_codes::LAPLACE_UNKNOWN_FORM,
            "Bronstein residue: poli irrazionali / RootOf richiedono dispatch dedicato",
            "F0.8");
    };

    // 1. F = N/D.
    auto parts = algebra::apart_num_den(F, ctx);
    if (parts.is_error()) return fail_unimpl("F non decomponibile in N/D");
    ExprPtr D = parts.value().denominator;

    // 2. Risolvi D(s) = 0.
    auto poles_res = algebra::solve_polynomial(D, s, ctx);
    if (poles_res.is_error()) return fail_unimpl(
        "denominatore non risolvibile via solve_polynomial");
    const auto& poles = poles_res.value();
    if (poles.empty()) return fail_unimpl("nessun polo trovato");

    // 3. Per ogni polo, calcola residuo di F(s)·exp(s·t).
    ExprPtr t_e = arena.make<Symbol>(t);
    ExprPtr s_e = arena.make<Symbol>(s);
    ExprPtr st = arena.make<Binary>(BinaryOp::Mul, s_e, t_e);
    ExprPtr exp_st = arena.make<FuncCall>(BuiltinOp::Exp,
        std::vector<ExprPtr>{st});
    ExprPtr integrand = arena.make<Binary>(BinaryOp::Mul, F, exp_st);

    std::vector<ExprPtr> residues;
    residues.reserve(poles.size());
    // Deduplicazione: lo stesso polo non va sommato due volte (residue()
    // gestisce molteplicità internamente via Laurent).
    std::vector<ExprPtr> seen_poles;
    for (ExprPtr p : poles) {
        bool dup = false;
        for (ExprPtr q : seen_poles) {
            ExprPtr diff_pq = arena.make<Binary>(BinaryOp::Sub, p, q);
            auto simp = ctx.simplify(diff_pq);
            if (simp.is_ok()) {
                if (const auto* il = expr_cast<IntegerLit>(simp.value());
                    il && il->value.is_zero()) { dup = true; break; }
                if (const auto* rl = expr_cast<RationalLit>(simp.value());
                    rl && rl->numerator.is_zero()) { dup = true; break; }
            }
        }
        if (dup) continue;
        seen_poles.push_back(p);
        auto r = residue(integrand, s, p, ctx);
        if (r.is_error()) return fail_unimpl(
            "residue computation failed at some pole");
        residues.push_back(r.value());
    }

    if (residues.empty()) return fail_unimpl("no residues computed");
    ExprPtr sum;
    if (residues.size() == 1U) sum = residues.front();
    else sum = arena.make<Sum>(std::move(residues));
    auto s_res = ctx.simplify(sum);
    return s_res.is_ok() ? s_res : ok(sum);
}

}  // namespace cas::calculus
