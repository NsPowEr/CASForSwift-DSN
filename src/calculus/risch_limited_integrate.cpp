// risch_limited_integrate.cpp — A38: the Limited Integration Problem over a
// differential tower K = Q(x, t_1, ..., t_n), Bronstein "Symbolic Integration
// I" §7.2.
//
// Problem (7.30):  given f, w_1, ..., w_m ∈ K, decide whether there exist
// constants c_1, ..., c_m ∈ Const(K) = Q and v ∈ K with
//
//      f = D(v) + c_1·w_1 + ... + c_m·w_m .
//
// Bronstein states two equivalent routes (§7.2, opening discussion).  We take
// the second one, which he calls "applicable for arbitrary w_i's": (7.30) IS a
// parametric Risch differential equation for v, so it is solved by the §7.1
// algorithm — already built, tower-recursive and sound (A1 + A26,
// solve_risch_de_parametric_field).  Following (7.35)/(7.36) verbatim, we
// introduce the extra unknown constant c_0 on f and solve the homogeneous
// parametric equation
//
//      D(v) + 0·v = c_0·f + Σ_i c_i·(−w_i)                              (7.36)
//
// then impose the affine constraint c_0 = 1.  The §7.1 solver returns a basis
// {(y_j, c_j)} of the (linear) solution space; the constraint c_0 = 1 is a
// single linear equation on the basis coordinates, so a solution with c_0 = 1
// exists iff some basis vector has c_{j,0} ≠ 0, and dividing that vector by
// c_{j,0} yields one.  If every basis vector has c_{j,0} = 0 the answer is a
// legitimate NEGATIVE ("no such v and c_i exist"), reported as a distinct
// diagnostic rather than a silent failure.
//
// LimitedIntegrateReduce (the a, b, h, N pre-reduction of §7.2) is, in
// Bronstein's own words, "a simplified version of the algorithm of Sect. 7.1
// that takes advantage of [no cancellation]": it sharpens the denominator and
// degree bounds, i.e. it is a performance refinement of the very path taken
// here, not an independent correctness requirement.  It is deliberately not
// implemented yet — the general §7.1 descent already bounds denominators and
// degrees soundly.
//
// SOUND BY CONSTRUCTION (REGOLA ZERO): the returned (v, c_i) is re-verified by
// exact back-substitution in the tower, f − D(v) − Σ c_i·w_i ≡ 0; a candidate
// that does not verify is dropped.  At worst incomplete (diagnostic), never a
// wrong answer.

#include "calculus_internal.hpp"
#include "risch_parametric_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

using detail::rational_to_expr;

// e ≡ 0 in the tower?  A bare simplify() does not always collapse a rational
// expression to the literal 0, so clear denominators first and test the
// numerator (same robustness pattern as solve_param_limited_integration_rational_q).
[[nodiscard]] bool tower_is_zero(ExprPtr e, symbolic::CASContext& ctx) {
    auto tog = algebra::together(e, ctx);
    ExprPtr x = tog.is_ok() ? tog.value() : e;
    auto s = ctx.simplify(x);
    ExprPtr z = s.is_ok() ? s.value() : x;
    if (const auto* il = expr_cast<IntegerLit>(z)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(z)) return rl->numerator.is_zero();
    auto parts = algebra::apart_num_den(z, ctx);
    if (parts.is_error()) return false;
    ExprPtr num = parts.value().numerator;
    auto num_exp = algebra::expand(num, ctx);
    if (num_exp.is_ok()) num = num_exp.value();
    auto num_simp = ctx.simplify(num);
    if (num_simp.is_ok()) num = num_simp.value();
    if (const auto* il = expr_cast<IntegerLit>(num)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(num)) return rl->numerator.is_zero();
    return false;
}

}  // namespace

Result<LimitedIntegrationFieldSolution> limited_integrate_field(
    ExprPtr f,
    const std::vector<ExprPtr>& w_vec,
    const DifferentialField& field,
    symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();
    const std::size_t m = w_vec.size();

    auto diagnostic = [&](const char* msg, const char* code, const char* hint) {
        return make_unimplemented<LimitedIntegrationFieldSolution>(
            "calculus", "limited_integrate_field", msg, code, hint, "A38");
    };

    // (7.36): unknowns (c_0, c_1, ..., c_m) on the forcings (f, −w_1, ..., −w_m).
    std::vector<ExprPtr> g_vec;
    g_vec.reserve(m + 1U);
    g_vec.push_back(f);
    for (ExprPtr w : w_vec) {
        ExprPtr neg = arena.make<Unary>(UnaryOp::Neg, w);
        if (auto s = ctx.simplify(neg); s.is_ok()) neg = s.value();
        g_vec.push_back(neg);
    }

    ExprPtr zero = arena.make<IntegerLit>(BigInt(0));
    auto sols = solve_risch_de_parametric_field(
        zero, g_vec, field.extensions().size(), field, ctx);
    if (sols.is_error())
        return fail<LimitedIntegrationFieldSolution>(sols.error());

    // Impose c_0 = 1 (7.36).  c_0 is a linear coordinate on the solution space,
    // so a solution with c_0 = 1 exists iff some basis vector has c_0 ≠ 0.
    for (const auto& sol : sols.value()) {
        if (sol.c.empty() || sol.c[0].numerator().is_zero()) continue;
        const Rational inv = Rational(BigInt(1)) / sol.c[0];

        ExprPtr v = arena.make<Binary>(BinaryOp::Mul, rational_to_expr(inv, arena), sol.y);
        if (auto tog = algebra::together(v, ctx); tog.is_ok()) v = tog.value();
        if (auto s = ctx.simplify(v); s.is_ok()) v = s.value();

        std::vector<Rational> c;
        c.reserve(m);
        for (std::size_t i = 0; i < m; ++i) c.push_back(sol.c[i + 1U] * inv);

        // Back-substitution check: f − D(v) − Σ c_i·w_i ≡ 0 in the tower.
        auto Dv = field.derive_in_generators(v, ctx);
        if (Dv.is_error()) continue;
        std::vector<ExprPtr> residue_terms;
        residue_terms.push_back(f);
        residue_terms.push_back(arena.make<Unary>(UnaryOp::Neg, Dv.value()));
        for (std::size_t i = 0; i < m; ++i) {
            if (c[i].numerator().is_zero()) continue;
            residue_terms.push_back(arena.make<Unary>(UnaryOp::Neg,
                arena.make<Binary>(BinaryOp::Mul, rational_to_expr(c[i], arena), w_vec[i])));
        }
        ExprPtr residue = residue_terms.size() == 1U
            ? residue_terms.front()
            : static_cast<ExprPtr>(arena.make<Sum>(std::move(residue_terms)));
        if (!tower_is_zero(residue, ctx)) continue;

        return ok(LimitedIntegrationFieldSolution{v, std::move(c)});
    }

    if (sols.value().empty()) {
        return diagnostic(
            "parametric Risch DE (7.36) returned no solution basis",
            cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
            "Limited integration (Bronstein §7.2): the §7.1 tower descent found "
            "no solution of D(v) = Σ c_i·g_i for this forcing");
    }
    return diagnostic(
        "no solution of (7.36) has c_0 = 1",
        cas::error::reason_codes::RISCH_NO_MATCH,
        "Limited integration (Bronstein §7.2): f is provably not of the form "
        "D(v) + Σ c_i·w_i over this tower (every solution of the associated "
        "parametric Risch DE has c_0 = 0), or the verified candidate failed "
        "exact back-substitution");
}

}  // namespace cas::calculus
