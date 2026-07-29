// A1 — Parametric PolyRischDE, non-cancellation case (Bronstein, "Symbolic
// Integration I", §7.1, algorithm ParamPolyRischDENoCancel1, "When deg(b) is
// Large Enough").  This is the df>0 branch of solve_risch_de_parametric_field
// (risch_rde_bronstein.cpp), previously Unimplemented (task A1, ledger
// HC-F8-PENDING-17), now reachable + tested via the A26 tower harness.
//
// Setting.  After the denominator-clearing step in the field solver we must
// solve the polynomial parametric equation over a differential field K[t]:
//
//     D(q) + b·q = Σ_i c_i·g_i        (b = f_new, g_i = g_new_i)         (7.16)
//
// for q ∈ K[t] (deg ≤ N) and constants c_i ∈ Const(K) = Q, where deg_t(b) =
// d_b > 0.  For the log (δ(t)=0) and exp (δ(t)=1) monomials handled here,
// d_b > 0 ⇒ d_b > max(0, δ(t)−1), so the leading terms of D(q) and b·q never
// sum to zero — the *non-cancellation* regime (Lemma 6.5.1(i)).
//
// Algorithm (verbatim §7.1).  Peel the top coefficient top-down: at degree n,
//     s_i = coefficient(g_i, t^{n+d_b}) / lc(b) ∈ K,     h_i += s_i·t^n,
//     g_i ← g_i − D(s_i·t^n) − b·s_i·t^n
// for n = N..0.  The leftover constraint Σ_i c_i·g_i = 0 (g_i ∈ K[t] now) is a
// homogeneous linear system on the constants c_i — Bronstein's ConstantSystem.
// Here K = Q(x): each residual t-coefficient is a rational function of x, so
// clearing denominators and equating x-power coefficients gives Q-linear rows.
// Its null space parametrises q = Σ_i c_i·h_i.  The reduction to a constant
// system is done for any tower K = k(t_1,…,t_j) by constant_system_nullspace
// (Bronstein §7.1 ConstantSystem / Lemma 7.1.2), not just K = Q(x).
//
// SOUNDNESS (REGOLA ZERO).  Every candidate (q, c) is verified by field
// back-substitution D(q) + b·q ≡ Σ_i c_i·g_i in K[t]; unverified candidates are
// dropped.  At worst incomplete, never a wrong answer.

#include "calculus_internal.hpp"
#include "risch_parametric_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace cas::calculus {

namespace {

using detail::as_rational;

[[nodiscard]] ExprPtr int_lit(std::int64_t v, AstArena& arena) {
    return arena.make<IntegerLit>(BigInt(v));
}

// t^n as an expression.
[[nodiscard]] ExprPtr t_power(const Symbol& t, std::size_t n, AstArena& arena) {
    if (n == 0U) return int_lit(1, arena);
    ExprPtr ts = arena.make<Symbol>(t);
    if (n == 1U) return ts;
    return arena.make<Binary>(BinaryOp::Pow, ts,
                              int_lit(static_cast<std::int64_t>(n), arena));
}

// Coefficient of t^k in a PolyExpr (0 if absent).
[[nodiscard]] ExprPtr poly_coeff(const algebra::PolyExpr& p, std::size_t k, AstArena& arena) {
    const auto& cs = p.coefficients();
    if (k < cs.size() && cs[k]) return cs[k];
    return int_lit(0, arena);
}

// Q-coefficient vector of a polynomial-in-var expr, nullopt if not a polynomial
// over Q in var (e.g. a lower-tower generator θ still present).
[[nodiscard]] std::optional<std::vector<Rational>>
qx_coeffs(ExprPtr e, const Symbol& var, symbolic::CASContext& ctx) {
    auto pr = algebra::parse_polynomial(e, var, ctx);
    if (pr.is_error()) return std::nullopt;
    std::vector<Rational> out;
    out.reserve(pr.value().coefficients().size());
    for (ExprPtr c : pr.value().coefficients()) {
        auto r = as_rational(c);
        if (!r) { if (auto s = ctx.simplify(c); s.is_ok()) r = as_rational(s.value()); }
        if (!r) return std::nullopt;
        out.push_back(*r);
    }
    return out;
}

[[nodiscard]] ExprPtr rat_lit(const Rational& c, AstArena& arena) {
    return (c.denominator() == BigInt(1))
        ? static_cast<ExprPtr>(arena.make<IntegerLit>(c.numerator()))
        : static_cast<ExprPtr>(arena.make<RationalLit>(c.numerator(), c.denominator()));
}

}  // namespace

Result<std::vector<ParametricRischDeQSolution>>
solve_param_poly_risch_de_nocancel1(
    ExprPtr f_new, const std::vector<ExprPtr>& g_new_vec, int N,
    const Symbol& t, const DifferentialField& field, symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();
    const std::size_t m = g_new_vec.size();
    const Symbol& base = field.base_var();

    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<std::vector<ParametricRischDeQSolution>>(
            "calculus", "solve_param_poly_risch_de_nocancel1", msg,
            cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
            "Parametric PolyRischDE non-cancellation (Bronstein Symbolic "
            "Integration I §7.1, ParamPolyRischDENoCancel1)",
            "HC-F8-PENDING-17");
    };
    if (m == 0U) return fail_unimpl("empty forcing vector");
    if (N < 0) N = 0;

    // b = f_new as a polynomial in t; d_b = deg_t(b) > 0, b_d = lc(b) ∈ K.
    auto b_poly_res = algebra::parse_polynomial(f_new, t, ctx);
    if (b_poly_res.is_error()) return fail_unimpl("f_new not polynomial in t");
    const algebra::PolyExpr& b_poly = b_poly_res.value();
    const std::size_t d_b = poly_degree(b_poly);
    ExprPtr b_d = leading_coefficient(b_poly);
    if (!b_d || d_b == 0U) return fail_unimpl("deg_t(f_new) must be > 0 for non-cancellation");

    // Peel top coefficients top-down (Bronstein §7.1).  h[i] holds the particular
    // polynomial h_i (coeff vector by t-power); g_expr[i] is the running residual.
    std::vector<ExprPtr> g_expr(g_new_vec);
    std::vector<std::vector<ExprPtr>> h(m);

    for (int n = N; n >= 0; --n) {
        const std::size_t nz = static_cast<std::size_t>(n);
        for (std::size_t i = 0; i < m; ++i) {
            auto gp = algebra::parse_polynomial(g_expr[i], t, ctx);
            if (gp.is_error()) return fail_unimpl("residual g_i not polynomial in t");
            ExprPtr ck = poly_coeff(gp.value(), nz + d_b, arena);
            ExprPtr s_i = arena.make<Binary>(BinaryOp::Div, ck, b_d);
            if (auto s = ctx.simplify(s_i); s.is_ok()) s_i = s.value();

            while (h[i].size() <= nz) h[i].push_back(int_lit(0, arena));
            h[i][nz] = s_i;

            // g_i ← g_i − D(s_i·t^n) − b·s_i·t^n.
            ExprPtr term = arena.make<Binary>(BinaryOp::Mul, s_i, t_power(t, nz, arena));
            auto dterm = field.derive(term, ctx);
            if (dterm.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(dterm.error());
            ExprPtr bterm = arena.make<Binary>(BinaryOp::Mul, f_new, term);
            ExprPtr upd = arena.make<Binary>(BinaryOp::Sub,
                              arena.make<Binary>(BinaryOp::Sub, g_expr[i], dterm.value()), bterm);
            if (auto s = ctx.simplify(upd); s.is_ok()) upd = s.value();
            g_expr[i] = upd;
        }
    }

    // Residual constraint Σ_i c_i·g_i(t) ≡ 0 in K[t].  For each t-power j the
    // K-linear form Σ_i c_i·coeff_j(g_i) must vanish, giving a homogeneous system
    // A·c = 0 with entries in K.  ConstantSystem (Bronstein §7.1, Lemma 7.1.2)
    // reduces it to the system for its constant solutions and returns a Q null-
    // space basis — valid for any tower K = k(t_1,…,t_j), not just K = Q(x).
    std::vector<algebra::PolyExpr> res;
    res.reserve(m);
    std::size_t max_tdeg = 0;
    for (std::size_t i = 0; i < m; ++i) {
        auto rp = algebra::parse_polynomial(g_expr[i], t, ctx);
        if (rp.is_error()) return fail_unimpl("residual not polynomial in t");
        max_tdeg = std::max(max_tdeg, rp.value().coefficients().size());
        res.push_back(rp.value());
    }

    std::vector<std::vector<ExprPtr>> A;
    A.reserve(max_tdeg);
    for (std::size_t j = 0; j < max_tdeg; ++j) {
        std::vector<ExprPtr> row(m);
        for (std::size_t i = 0; i < m; ++i) row[i] = poly_coeff(res[i], j, arena);
        A.push_back(std::move(row));
    }

    auto basis_res = constant_system_nullspace(std::move(A), m, field, ctx);
    if (basis_res.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(basis_res.error());
    const std::vector<std::vector<Rational>>& basis = basis_res.value();

    // Verify a candidate (q, c) by field back-substitution: the residual
    // D(q) + f_new·q − Σ_i c_i·g_new_i must vanish identically in K[t].
    auto verify = [&](ExprPtr q, const std::vector<Rational>& c) -> bool {
        auto dq = field.derive(q, ctx);
        if (dq.is_error()) return false;
        ExprPtr lhs = arena.make<Binary>(BinaryOp::Add, dq.value(),
                          arena.make<Binary>(BinaryOp::Mul, f_new, q));
        ExprPtr delta = lhs;
        for (std::size_t i = 0; i < m; ++i) {
            if (c[i].numerator().is_zero()) continue;
            ExprPtr term = arena.make<Binary>(BinaryOp::Mul, rat_lit(c[i], arena), g_new_vec[i]);
            delta = arena.make<Binary>(BinaryOp::Sub, delta, term);
        }
        auto tog = algebra::together(delta, ctx);
        ExprPtr z = tog.is_ok() ? tog.value() : delta;
        if (auto s = ctx.simplify(z); s.is_ok()) z = s.value();
        if (const auto* il = expr_cast<IntegerLit>(z)) return il->value.is_zero();
        if (const auto* rl = expr_cast<RationalLit>(z)) return rl->numerator.is_zero();
        // Structural fallback: as a polynomial in t, every coefficient must be 0.
        auto zp = algebra::parse_polynomial(z, t, ctx);
        if (zp.is_error()) return false;
        for (ExprPtr coeff : zp.value().coefficients()) {
            if (!coeff) continue;
            auto ct = algebra::together(coeff, ctx);
            ExprPtr cc = ct.is_ok() ? ct.value() : coeff;
            if (auto s = ctx.simplify(cc); s.is_ok()) cc = s.value();
            auto parts = algebra::apart_num_den(cc, ctx);
            if (parts.is_error()) return false;
            auto nc = qx_coeffs(parts.value().numerator, base, ctx);
            if (!nc) return false;
            for (const auto& v : *nc) if (!v.numerator().is_zero()) return false;
        }
        return true;
    };

    std::vector<ParametricRischDeQSolution> out;
    for (auto& c : basis) {
        // q = Σ_i c_i·h_i  (h_i = Σ_k h[i][k]·t^k).
        std::vector<ExprPtr> terms;
        for (std::size_t i = 0; i < m; ++i) {
            if (c[i].numerator().is_zero()) continue;
            ExprPtr ci = rat_lit(c[i], arena);
            for (std::size_t k = 0; k < h[i].size(); ++k) {
                ExprPtr hk = h[i][k];
                if (const auto* il = expr_cast<IntegerLit>(hk); il && il->value.is_zero()) continue;
                terms.push_back(arena.make<Product>(std::vector<ExprPtr>{
                    ci, hk, t_power(t, k, arena)}));
            }
        }
        ExprPtr q = terms.empty() ? int_lit(0, arena)
                  : (terms.size() == 1U ? terms.front()
                                        : static_cast<ExprPtr>(arena.make<Sum>(std::move(terms))));
        if (auto s = ctx.simplify(q); s.is_ok()) q = s.value();

        if (!verify(q, c)) continue;
        out.push_back({q, c});
    }

    return ok(std::move(out));
}

}  // namespace cas::calculus
