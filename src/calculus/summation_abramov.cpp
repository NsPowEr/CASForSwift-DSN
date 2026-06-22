// F5.7 — Abramov rational summation helpers.
//
// Sub-block 0 (Gosper):  definite sum via Newton-Leibniz finite calculus.
// Sub-block 1 (Polygamma): A/(linear(k))^m atoms via ψ^(m-1) antidifference.
// Sub-block 2 (Abramov-Full): partial-fraction decomposition + per-atom routing.
// B6-bis (Quadratic-RootOf): Q-irreducible quadratic atoms (B₁k+B₀)/Q(k) via
//   C·ψ(k−α) + D·ψ(k−β) where α,β = RootOf(Q̃,t,0/1) and C,D are residues.

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include "summation_internal.hpp"
#include "summation_abramov_internal.hpp"
#include "../symbolic/summation_gosper.hpp"
#include <optional>
#include <vector>

namespace cas::calculus {

using namespace abramov_detail;

Result<ExprPtr> try_gosper_definite(
    ExprPtr term, const Symbol& var,
    ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx) {
    auto antidiff = symbolic::gosper_sum(term, var, ctx);
    if (antidiff.is_error()) return fail<ExprPtr>(antidiff.error());
    if (!antidiff.value().has_value()) {
        return fail<ExprPtr>(CASError{
            .kind    = CASErrorKind::Unimplemented,
            .message = "Gosper: term is not Gosper-summable",
        });
    }
    AstArena& arena = ctx.arena();
    ExprPtr S = antidiff.value().value();
    ExprPtr upper_plus_one = arena.make<Binary>(BinaryOp::Add, upper,
        arena.make<IntegerLit>(BigInt(1)));
    auto S_upper = ctx.substitute(S, var, upper_plus_one);
    if (S_upper.is_error()) return fail<ExprPtr>(S_upper.error());
    auto S_lower = ctx.substitute(S, var, lower);
    if (S_lower.is_error()) return fail<ExprPtr>(S_lower.error());
    ExprPtr diff = arena.make<Binary>(BinaryOp::Sub,
        S_upper.value(), S_lower.value());
    return ctx.simplify(diff);
}

Result<ExprPtr> try_polygamma_definite(
    ExprPtr term, const Symbol& var,
    ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx) {
    auto antidiff_opt = try_polygamma_antidiff(term, var, ctx);
    if (!antidiff_opt.has_value()) {
        return fail<ExprPtr>(CASError{
            .kind    = CASErrorKind::Unimplemented,
            .message = "Polygamma summation: term is not of the form "
                       "A/(linear(k))^m",
        });
    }
    AstArena& arena = ctx.arena();
    ExprPtr S = antidiff_opt.value();
    ExprPtr upper_plus_one = arena.make<Binary>(BinaryOp::Add, upper,
        arena.make<IntegerLit>(BigInt(1)));
    auto S_upper = ctx.substitute(S, var, upper_plus_one);
    if (S_upper.is_error()) return fail<ExprPtr>(S_upper.error());
    auto S_lower = ctx.substitute(S, var, lower);
    if (S_lower.is_error()) return fail<ExprPtr>(S_lower.error());
    ExprPtr diff = arena.make<Binary>(BinaryOp::Sub,
        S_upper.value(), S_lower.value());
    return ctx.simplify(diff);
}

Result<ExprPtr> try_abramov_definite(
    ExprPtr term, const Symbol& var,
    ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx) {
    if (!has_rational_dependency(term, var.name)) {
        return fail<ExprPtr>(CASError{
            .kind    = CASErrorKind::Unimplemented,
            .message = "Abramov: term has no rational dependency on the "
                       "summation variable",
        });
    }

    auto together_res = algebra::together(term, ctx);
    if (together_res.is_error()) return fail<ExprPtr>(together_res.error());
    auto parts = algebra::apart_num_den(together_res.value(), ctx);
    if (parts.is_error()) return fail<ExprPtr>(parts.error());

    auto deg_N_res = algebra::polynomial_degree(parts.value().numerator, var, ctx);
    auto deg_D_res = algebra::polynomial_degree(parts.value().denominator, var, ctx);
    if (deg_N_res.is_error() || deg_D_res.is_error()) {
        return fail<ExprPtr>(CASError{
            .kind    = CASErrorKind::Unimplemented,
            .message = "Abramov: rational term must be polynomial-over-"
                       "polynomial in the summation variable",
        });
    }
    const std::size_t deg_N = deg_N_res.value();
    const std::size_t deg_D = deg_D_res.value();
    AstArena& arena = ctx.arena();
    ExprPtr N = parts.value().numerator;
    ExprPtr D = parts.value().denominator;

    ExprPtr polynomial_part  = arena.make<IntegerLit>(BigInt(0));
    ExprPtr proper_remainder = term;

    if (deg_N >= deg_D) {
        auto divmod = algebra::polynomial_divmod(N, D, var, ctx);
        if (divmod.is_error()) return fail<ExprPtr>(divmod.error());
        polynomial_part  = divmod.value().quotient;
        proper_remainder = arena.make<Binary>(BinaryOp::Div,
            divmod.value().remainder, D);
        auto simp_r = ctx.simplify(proper_remainder);
        if (simp_r.is_error()) return fail<ExprPtr>(simp_r.error());
        proper_remainder = simp_r.value();
    }

    auto fractions = algebra::partial_fractions(proper_remainder, var, ctx);
    if (fractions.is_error()) return fail<ExprPtr>(fractions.error());

    ExprPtr upper_plus_one = arena.make<Binary>(BinaryOp::Add, upper,
        arena.make<IntegerLit>(BigInt(1)));

    std::vector<ExprPtr> total_terms;

    bool poly_part_is_zero = false;
    if (const auto* il = expr_cast<IntegerLit>(polynomial_part))
        poly_part_is_zero = il->value.is_zero();
    if (!poly_part_is_zero) {
        auto poly_def = try_gosper_definite(polynomial_part, var, lower, upper, ctx);
        if (poly_def.is_ok()) {
            total_terms.push_back(poly_def.value());
        } else {
            return poly_def;
        }
    }

    for (ExprPtr atom : fractions.value()) {
        auto atom_simp = ctx.simplify(atom);
        if (atom_simp.is_error()) return fail<ExprPtr>(atom_simp.error());

        // Q-linear atoms: A/(c₁k+c₀)^m → polygamma antidifference.
        auto antidiff_opt = try_polygamma_antidiff(atom_simp.value(), var, ctx);

        // Q-irreducible quadratic atoms: (B₁k+B₀)/Q(k)^m → RootOf polygamma.
        if (!antidiff_opt.has_value())
            antidiff_opt = try_quadratic_atom_antidiff(atom_simp.value(), var, ctx);

        if (!antidiff_opt.has_value()) {
            return fail<ExprPtr>(CASError{
                .kind    = CASErrorKind::Unimplemented,
                .message = "Abramov: partial-fraction atom not of supported form. "
                           "Handled: Q-linear A/(ck+d)^m and Q-irreducible "
                           "quadratic (B₁k+B₀)/Q(k)^m. "
                           "Unhandled: Q-irreducible deg≥3.",
            });
        }

        ExprPtr S = antidiff_opt.value();
        auto S_upper = ctx.substitute(S, var, upper_plus_one);
        if (S_upper.is_error()) return fail<ExprPtr>(S_upper.error());
        auto S_lower = ctx.substitute(S, var, lower);
        if (S_lower.is_error()) return fail<ExprPtr>(S_lower.error());
        total_terms.push_back(
            arena.make<Binary>(BinaryOp::Sub,
                S_upper.value(), S_lower.value()));
    }

    if (total_terms.empty())
        return ok(arena.make<IntegerLit>(BigInt(0)));

    ExprPtr combined = (total_terms.size() == 1U)
        ? total_terms.front()
        : arena.make<Sum>(std::move(total_terms));
    return ctx.simplify(combined);
}

} // namespace cas::calculus
