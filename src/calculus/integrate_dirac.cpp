// F5.9-pre — Sifting property ∫ δ(t - a) · f(t) dt = f(a).
//
// Routine try_integrate_dirac_sifting(expr, var, ctx):
//   - Riconosce expr = Π_i f_i(t) dove esattamente un fattore è
//     FuncCall(DiracDelta, arg) con arg lineare in var.
//   - Estrae il polo a risolvendo arg(a) = 0 (cap.6 caso sub-lineare).
//   - Valuta il resto dei fattori in var = a → f(a).
//   - Per arg(t) = t − a (più frequente) o arg(t) = c·(t − a):
//     sifting estesa  ∫ δ(c·(t−a)) · f(t) dt = f(a) / |c|.
//     Il simbolic /|c| è rappresentato come abs(c)^(-1) per evitare
//     ramificazioni di segno; il consumatore può specializzare via
//     assumptions su c.
//
// Caso non-lineare (arg = g(t) con g'(a) ≠ 0): scope follow-up.

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

[[nodiscard]] bool depends_on_v(ExprPtr e, const Symbol& v) {
    if (!e) return false;
    if (const auto* s = expr_cast<Symbol>(e)) return s->name == v.name;
    if (const auto* un = expr_cast<Unary>(e)) return depends_on_v(un->operand, v);
    if (const auto* bin = expr_cast<Binary>(e))
        return depends_on_v(bin->left, v) || depends_on_v(bin->right, v);
    if (const auto* sum = expr_cast<Sum>(e)) {
        for (ExprPtr t : sum->terms) if (depends_on_v(t, v)) return true;
        return false;
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        for (ExprPtr f : prod->factors) if (depends_on_v(f, v)) return true;
        return false;
    }
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        for (ExprPtr a : fc->args) if (depends_on_v(a, v)) return true;
        return false;
    }
    return false;
}

// Cerca DiracDelta(arg) nella expr (top-level o dentro Product).  Se trovato,
// restituisce (delta_factor_index, arg_of_delta).  delta_factor_index = -1
// significa expr === DiracDelta(arg) direttamente.
struct DeltaMatch {
    ExprPtr arg;
    std::vector<ExprPtr> other_factors;
    bool found{false};
};

[[nodiscard]] DeltaMatch find_dirac_in_product(ExprPtr expr) {
    DeltaMatch m;
    if (const auto* fc = expr_cast<FuncCall>(expr);
        fc && fc->func_id == BuiltinOp::DiracDelta && fc->args.size() == 1U) {
        m.found = true;
        m.arg = fc->args[0];
        return m;
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr f : prod->factors) {
            if (const auto* fc = expr_cast<FuncCall>(f);
                fc && fc->func_id == BuiltinOp::DiracDelta
                && fc->args.size() == 1U) {
                if (m.found) {
                    // Più di un delta: scope cap.10+, non gestiamo.
                    m.found = false;
                    return m;
                }
                m.found = true;
                m.arg = fc->args[0];
            } else {
                m.other_factors.push_back(f);
            }
        }
    } else if (const auto* bin = expr_cast<Binary>(expr); bin && bin->op == BinaryOp::Mul) {
        for (ExprPtr f : {bin->left, bin->right}) {
            if (const auto* fc = expr_cast<FuncCall>(f);
                fc && fc->func_id == BuiltinOp::DiracDelta
                && fc->args.size() == 1U) {
                if (m.found) { m.found = false; return m; }
                m.found = true;
                m.arg = fc->args[0];
            } else {
                m.other_factors.push_back(f);
            }
        }
    }
    return m;
}

}  // namespace

Result<ExprPtr> try_integrate_dirac_sifting(
    ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<ExprPtr>(
            "calculus", "try_integrate_dirac_sifting",
            msg,
            cas::error::reason_codes::LAPLACE_UNKNOWN_FORM,
            "Dirac sifting: arg non lineare o Dirac multipli — scope cap.10+",
            "F0.8");
    };

    auto m = find_dirac_in_product(expr);
    if (!m.found) return fail_unimpl("integrand non contiene δ(arg) riconoscibile");

    ExprPtr arg = m.arg;
    if (!depends_on_v(arg, var))
        return fail_unimpl("δ(arg) costante in var — integrale degenere");

    // Verifica arg lineare in var: c = d(arg)/d(var) indip da var.
    auto c_res = diff(arg, var, 1U, ctx);
    if (c_res.is_error()) return fail_unimpl("impossibile derivare arg(t)");
    ExprPtr c = c_res.value();
    if (depends_on_v(c, var))
        return fail_unimpl("δ(arg) con arg non lineare in var");

    // Risolve arg(var = a) = 0  ⇒  a = (arg(0)) / (-c)  (linear arg = c·var + b,
    // dove b = arg − c·var, b indipendente da var).
    ExprPtr cv = arena.make<Binary>(BinaryOp::Mul, c, arena.make<Symbol>(var));
    ExprPtr b = arena.make<Binary>(BinaryOp::Sub, arg, cv);
    auto b_tog = algebra::together(b, ctx);
    if (b_tog.is_error()) return fail_unimpl("impossibile estrarre b");
    auto b_simp = ctx.simplify(b_tog.value());
    if (b_simp.is_error()) return fail_unimpl("impossibile semplificare b");
    if (depends_on_v(b_simp.value(), var))
        return fail_unimpl("arg non riconducibile a c·var + b");

    // a = -b/c.
    ExprPtr neg_b = arena.make<Unary>(UnaryOp::Neg, b_simp.value());
    ExprPtr a_val = arena.make<Binary>(BinaryOp::Div, neg_b, c);
    auto a_simp = ctx.simplify(a_val);
    if (a_simp.is_error()) return fail_unimpl("impossibile semplificare a");
    a_val = a_simp.value();

    // Costruisci f(var) = prodotto degli altri fattori; valuta in var = a.
    ExprPtr f_of_var;
    if (m.other_factors.empty()) {
        f_of_var = arena.make<IntegerLit>(BigInt(1));
    } else if (m.other_factors.size() == 1U) {
        f_of_var = m.other_factors[0];
    } else {
        f_of_var = arena.make<Product>(std::vector<ExprPtr>(m.other_factors));
    }
    auto f_at_a = ctx.substitute(f_of_var, var, a_val);
    if (f_at_a.is_error()) return fail_unimpl("substitute(f, var, a) failed");

    // Divide per |c|: rappresentato come abs(c) tramite FuncCall(Abs, c).
    // Se c è IntegerLit positivo, abs = c.  Caso generale rimane symbolic abs.
    bool c_positive = false;
    if (const auto* il = expr_cast<IntegerLit>(c); il && !il->value.is_negative()
        && !il->value.is_zero()) c_positive = true;
    if (const auto* rl = expr_cast<RationalLit>(c);
        rl && !rl->numerator.is_negative() && !rl->numerator.is_zero()) c_positive = true;

    ExprPtr abs_c = c_positive
        ? c
        : arena.make<FuncCall>(BuiltinOp::Abs, std::vector<ExprPtr>{c});
    ExprPtr result = arena.make<Binary>(BinaryOp::Div, f_at_a.value(), abs_c);
    auto res_simp = ctx.simplify(result);
    return res_simp.is_ok() ? res_simp : ok(result);
}

}  // namespace cas::calculus
