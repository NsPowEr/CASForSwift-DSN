// A26 / HC-A26-PRIMITIVE-PARAMQ-RATIONAL — rational limited integration Q(x).
//
// Vertical split from risch_parametric.cpp (anti-monolith, CLAUDE.md §ANTI-MONOLITH).
// Contains:  solve_param_limited_integration_rational_q
//
// Solve  y' = Σ c_i g_i  for c_i ∈ Q and y ∈ Q(x), given rational g_i (the
// f = 0 case of the parametric Risch DE; Bronstein §7.2/§7.3).  This is exactly
// the forcing the primitive (log) tower descent produces, where the base
// solver solve_risch_de_parametric_q (Q[x] only) bails on a rational g.
//
// Method: integrate each g_i with the engine's complete rational-function
// integrator, split ∫g_i = R_i (rational) + Σ_k a_{ik} A_k where each A_k is a
// log/arctan atom with rational coefficient a_{ik} (the residues).  Σ c_i g_i
// integrates rationally iff Σ_i c_i a_{ik} = 0 for every atom A_k; the solution
// is then y = Σ c_i R_i.  The c-basis is the null space of the atom-coefficient
// matrix.  SOUND BY CONSTRUCTION: every candidate (y, c) is verified by
// back-substitution D(y) ≡ Σ c_i g_i, and unverifiable ones are dropped — so an
// incomplete atom split yields incompleteness (diagnostic), never a wrong
// answer (REGOLA ZERO).  Algebraic (non-rational) residues are not handled and
// produce a clean Unimplemented.

#include "calculus_internal.hpp"
#include "risch_parametric_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

using detail::as_rational;
using detail::null_space_basis;
using detail::row_echelon;

// e is a rational function of var (built from var, constants, +,−,×,÷ and
// integer powers).  A var-dependent FuncCall or non-integer power is not.
[[nodiscard]] bool is_rational_in_var(ExprPtr e, const Symbol& var) {
    if (!depends_on(e, var)) return true;
    if (expr_cast<Symbol>(e)) return true;
    if (const auto* u = expr_cast<Unary>(e); u && u->op == UnaryOp::Neg)
        return is_rational_in_var(u->operand, var);
    if (const auto* b = expr_cast<Binary>(e)) {
        switch (b->op) {
            case BinaryOp::Add: case BinaryOp::Sub:
            case BinaryOp::Mul: case BinaryOp::Div:
                return is_rational_in_var(b->left, var) &&
                       is_rational_in_var(b->right, var);
            case BinaryOp::Pow:
                return is_rational_in_var(b->left, var) &&
                       expr_cast<IntegerLit>(b->right) != nullptr;
            default: return false;
        }
    }
    if (const auto* s = expr_cast<Sum>(e)) {
        for (ExprPtr t : s->terms) if (!is_rational_in_var(t, var)) return false;
        return true;
    }
    if (const auto* p = expr_cast<Product>(e)) {
        for (ExprPtr fct : p->factors) if (!is_rational_in_var(fct, var)) return false;
        return true;
    }
    return false;
}

// Flatten a multiplicative term into factors (Neg → explicit −1 factor).
void flatten_factors(ExprPtr e, std::vector<ExprPtr>& out, AstArena& arena) {
    if (const auto* p = expr_cast<Product>(e)) {
        for (ExprPtr f : p->factors) flatten_factors(f, out, arena);
        return;
    }
    if (const auto* b = expr_cast<Binary>(e); b && b->op == BinaryOp::Mul) {
        flatten_factors(b->left, out, arena);
        flatten_factors(b->right, out, arena);
        return;
    }
    if (const auto* u = expr_cast<Unary>(e); u && u->op == UnaryOp::Neg) {
        out.push_back(arena.make<IntegerLit>(BigInt(-1)));
        flatten_factors(u->operand, out, arena);
        return;
    }
    out.push_back(e);
}

// Flatten an additive expression into terms.
void collect_additive_terms(ExprPtr e, std::vector<ExprPtr>& out, AstArena& arena) {
    if (const auto* s = expr_cast<Sum>(e)) {
        for (ExprPtr t : s->terms) collect_additive_terms(t, out, arena);
        return;
    }
    if (const auto* b = expr_cast<Binary>(e)) {
        if (b->op == BinaryOp::Add) {
            collect_additive_terms(b->left, out, arena);
            collect_additive_terms(b->right, out, arena);
            return;
        }
        if (b->op == BinaryOp::Sub) {
            collect_additive_terms(b->left, out, arena);
            collect_additive_terms(arena.make<Unary>(UnaryOp::Neg, b->right), out, arena);
            return;
        }
    }
    if (const auto* u = expr_cast<Unary>(e); u && u->op == UnaryOp::Neg) {
        std::vector<ExprPtr> inner;
        collect_additive_terms(u->operand, inner, arena);
        for (ExprPtr t : inner) out.push_back(arena.make<Unary>(UnaryOp::Neg, t));
        return;
    }
    out.push_back(e);
}

// A transcendental term: (rational constant) · (single log/arctan atom).
[[nodiscard]] std::optional<std::pair<Rational, ExprPtr>>
extract_coeff_atom(ExprPtr term, const Symbol& var, symbolic::CASContext& ctx) {
    std::vector<ExprPtr> factors;
    flatten_factors(term, factors, ctx.arena());
    Rational coeff(BigInt(1));
    std::vector<ExprPtr> atoms;
    for (ExprPtr f : factors) {
        if (!depends_on(f, var)) {
            auto r = as_rational(f);
            if (!r) { if (auto s = ctx.simplify(f); s.is_ok()) r = as_rational(s.value()); }
            if (!r) return std::nullopt;
            coeff = coeff * *r;
        } else {
            atoms.push_back(f);
        }
    }
    if (atoms.size() != 1U) return std::nullopt;
    if (!expr_cast<FuncCall>(atoms[0])) return std::nullopt;
    return std::make_pair(coeff, atoms[0]);
}

[[nodiscard]] ExprPtr rational_to_expr(const Rational& c, AstArena& arena) {
    return (c.denominator() == BigInt(1))
        ? static_cast<ExprPtr>(arena.make<IntegerLit>(c.numerator()))
        : static_cast<ExprPtr>(arena.make<RationalLit>(c.numerator(), c.denominator()));
}

// Σ_i c_i · terms[i]  (skips zero coefficients).
[[nodiscard]] ExprPtr scaled_sum(
    const std::vector<Rational>& c, const std::vector<ExprPtr>& terms, AstArena& arena) {
    std::vector<ExprPtr> acc;
    for (std::size_t i = 0; i < terms.size(); ++i) {
        if (c[i].numerator().is_zero()) continue;
        acc.push_back(arena.make<Binary>(BinaryOp::Mul, rational_to_expr(c[i], arena), terms[i]));
    }
    if (acc.empty()) return arena.make<IntegerLit>(BigInt(0));
    if (acc.size() == 1U) return acc.front();
    return arena.make<Sum>(std::move(acc));
}

[[nodiscard]] bool simplifies_to_zero(ExprPtr e, symbolic::CASContext& ctx) {
    auto tog = algebra::together(e, ctx);
    ExprPtr x = tog.is_ok() ? tog.value() : e;
    auto s = ctx.simplify(x);
    if (s.is_error()) return false;
    if (const auto* il = expr_cast<IntegerLit>(s.value())) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(s.value())) return rl->numerator.is_zero();
    return false;
}

}  // namespace

Result<std::vector<ParametricRischDeQSolution>>
solve_param_limited_integration_rational_q(
    const std::vector<ExprPtr>& g_vec, const Symbol& var, symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();
    const std::size_t m = g_vec.size();

    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<std::vector<ParametricRischDeQSolution>>(
            "calculus", "solve_param_limited_integration_rational_q", msg,
            cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
            "Rational limited integration: only rational residues handled; general "
            "f != 0 rational ParamRischDE pending (HC-A26-PRIMITIVE-PARAMQ-RATIONAL)",
            "HC-A26-PRIMITIVE-PARAMQ-RATIONAL");
    };
    if (m == 0U) return fail_unimpl("empty forcing vector");

    std::vector<ExprPtr> R(m);
    std::vector<ExprPtr> atom_keys;
    std::vector<std::vector<Rational>> atom_coeff;

    for (std::size_t i = 0; i < m; ++i) {
        auto Gi = integrate(g_vec[i], var, ctx);
        if (Gi.is_error())
            return fail_unimpl("engine could not integrate a rational forcing g_i");
        ExprPtr Gi_e = Gi.value();
        if (auto s = ctx.simplify(Gi_e); s.is_ok()) Gi_e = s.value();

        std::vector<ExprPtr> terms;
        collect_additive_terms(Gi_e, terms, arena);
        std::vector<ExprPtr> rat_terms;
        for (ExprPtr t : terms) {
            if (is_rational_in_var(t, var)) { rat_terms.push_back(t); continue; }
            auto ca = extract_coeff_atom(t, var, ctx);
            if (!ca)
                return fail_unimpl("antiderivative not in rational + Σ c·log/arctan form");
            std::size_t idx = atom_keys.size();
            for (std::size_t k = 0; k < atom_keys.size(); ++k)
                if (structural_equal(atom_keys[k], ca->second)) { idx = k; break; }
            if (idx == atom_keys.size()) {
                atom_keys.push_back(ca->second);
                atom_coeff.emplace_back(m, Rational(BigInt(0)));
            }
            atom_coeff[idx][i] = atom_coeff[idx][i] + ca->first;
        }
        R[i] = rat_terms.empty() ? arena.make<IntegerLit>(BigInt(0))
             : (rat_terms.size() == 1U ? rat_terms.front()
                                       : arena.make<Sum>(std::move(rat_terms)));
    }

    const std::size_t K = atom_keys.size();
    std::vector<std::vector<Rational>> basis;
    if (K == 0U) {
        for (std::size_t i = 0; i < m; ++i) {
            std::vector<Rational> e_i(m, Rational(BigInt(0)));
            e_i[i] = Rational(BigInt(1));
            basis.push_back(std::move(e_i));
        }
    } else {
        std::vector<std::vector<Rational>> M = atom_coeff;
        auto pivots = row_echelon(M, m);
        basis = null_space_basis(M, pivots, m);
    }

    std::vector<ParametricRischDeQSolution> out;
    for (auto& v : basis) {
        ExprPtr y = scaled_sum(v, R, arena);
        if (auto tog = algebra::together(y, ctx); tog.is_ok()) {
            if (auto s = ctx.simplify(tog.value()); s.is_ok()) y = s.value();
        }
        auto dy = diff(y, var, 1U, ctx);
        if (dy.is_error()) continue;
        ExprPtr rhs = scaled_sum(v, g_vec, arena);
        ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, dy.value(), rhs);
        if (!simplifies_to_zero(delta, ctx)) continue;
        out.push_back({y, v});
    }

    return ok(std::move(out));
}

}  // namespace cas::calculus
