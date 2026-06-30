// A26 / HC-A26-PRIMITIVE-PARAMQ-RATIONAL — parametric Risch DE over Q(x).
// Vertical split from risch_parametric.cpp (anti-monolith, CLAUDE.md §ANTI-MONOLITH).
// Contains:
//   solve_param_limited_integration_rational_q   (f = 0 case, §7.2/§7.3)
//   solve_param_risch_de_rational_q              (f ≠ 0 case, §6.1/§7.1)
// Both are SOUND BY CONSTRUCTION (REGOLA ZERO): every candidate is verified
// exactly (f=0: symbolic back-substitution D(y) ≡ Σ c_i g_i; f≠0: Q-rational
// residual of the cleared identity) and unverified candidates are dropped — at
// worst incomplete (diagnostic), never a wrong answer.  Per-function headers
// below carry the algorithm detail.

#include "calculus_internal.hpp"
#include "risch_parametric_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

using detail::as_rational;
using detail::null_space_basis;
using detail::row_echelon;

// Forward decl (defined below; used by the robust zero-check).
[[nodiscard]] std::optional<std::vector<Rational>>
poly_coeffs_q(ExprPtr e, const Symbol& var, symbolic::CASContext& ctx);

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

// Robust rational-function zero check: e ≡ 0 iff, after clearing denominators,
// its numerator is the zero polynomial in var.  (A bare simplify() does not
// always collapse a rational expression to a literal 0, so we test the
// numerator structurally.)
[[nodiscard]] bool simplifies_to_zero(
    ExprPtr e, const Symbol& var, symbolic::CASContext& ctx) {
    auto tog = algebra::together(e, ctx);
    ExprPtr x = tog.is_ok() ? tog.value() : e;
    auto s = ctx.simplify(x);
    ExprPtr z = s.is_ok() ? s.value() : x;
    if (const auto* il = expr_cast<IntegerLit>(z)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(z)) return rl->numerator.is_zero();
    auto parts = algebra::apart_num_den(z, ctx);
    if (parts.is_error()) return false;
    auto nc = poly_coeffs_q(parts.value().numerator, var, ctx);
    if (!nc) return false;
    for (const auto& c : *nc) if (!c.numerator().is_zero()) return false;
    return true;
}

void strip_trailing_q(std::vector<Rational>& v) {
    while (!v.empty() && v.back().numerator().is_zero()) v.pop_back();
}

// Q-coefficient vector of a polynomial-in-var ExprPtr (nullopt if not a
// polynomial over Q in var).
[[nodiscard]] std::optional<std::vector<Rational>>
poly_coeffs_q(ExprPtr e, const Symbol& var, symbolic::CASContext& ctx) {
    auto pr = algebra::parse_polynomial(e, var, ctx);
    if (pr.is_error()) return std::nullopt;
    std::vector<Rational> out;
    out.reserve(pr.value().size());
    for (ExprPtr c : pr.value().coefficients()) {
        auto r = as_rational(c);
        if (!r) { if (auto s = ctx.simplify(c); s.is_ok()) r = as_rational(s.value()); }
        if (!r) return std::nullopt;
        out.push_back(*r);
    }
    return out;
}

// Σ_k P_k · var^k  as ExprPtr.
[[nodiscard]] ExprPtr build_poly_expr_q(
    const std::vector<Rational>& P, const Symbol& var, AstArena& arena) {
    std::vector<ExprPtr> terms;
    ExprPtr vs = arena.make<Symbol>(var);
    for (std::size_t k = 0; k < P.size(); ++k) {
        if (P[k].numerator().is_zero()) continue;
        ExprPtr ce = rational_to_expr(P[k], arena);
        if (k == 0U) { terms.push_back(ce); continue; }
        ExprPtr xk = (k == 1U) ? vs
            : static_cast<ExprPtr>(arena.make<Binary>(BinaryOp::Pow, vs,
                  arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(k)))));
        const bool is_one = P[k].numerator() == BigInt(1) && P[k].denominator() == BigInt(1);
        terms.push_back(is_one ? xk : arena.make<Binary>(BinaryOp::Mul, ce, xk));
    }
    if (terms.empty()) return arena.make<IntegerLit>(BigInt(0));
    if (terms.size() == 1U) return terms.front();
    return arena.make<Sum>(std::move(terms));
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
        if (!simplifies_to_zero(delta, var, ctx)) continue;
        out.push_back({y, v});
    }

    return ok(std::move(out));
}

// --- f ≠ 0 : full parametric Risch DE  y' + f·y = Σ c_i g_i over Q(x). ---
//
// Ansatz y = P/D with D = lcm(den f, den g_i).  Substituting into the
// LCM-cleared form  D·y' + Fn·y = Σ c_i Gn_i  (Fn = fn·D/fd, Gn_i = gn_i·D/gd_i)
// gives the polynomial identity (Bronstein §6.1, parametrized over c_i):
//
//     D·P' + (Fn − D')·P − Σ_i c_i (Gn_i·D)  =  0
//
// a homogeneous Q-linear system in the coefficients of P (deg ≤ M_bound) and
// the scalars c_1..c_m.  Its null space is a basis of the solution space; each
// basis vector (P, c) yields y = P/D, VERIFIED by back-substitution
// D(y)+f·y ≡ Σ c_i g_i (sound by construction — unverified candidates dropped).
//
// LIMITATION (shared verbatim with the non-parametric solve_risch_de_rational_q):
// the ansatz finds only solutions whose denominator divides D.  Negative-integer
// residues at simple poles of f (WeakNormalizer, Bronstein 6.1.1) can give y a
// higher-order pole than den(f) carries; such solutions are silently absent —
// incompleteness, never a wrong answer.  Tracked as an engine-wide gap.
Result<std::vector<ParametricRischDeQSolution>>
solve_param_risch_de_rational_q(
    ExprPtr f_expr, const std::vector<ExprPtr>& g_vec, const Symbol& var,
    symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();
    const std::size_t m = g_vec.size();

    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<std::vector<ParametricRischDeQSolution>>(
            "calculus", "solve_param_risch_de_rational_q", msg,
            cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
            "Parametric Risch DE over Q(x): den(y) must divide lcm(den f, den g_i); "
            "negative-integer-residue poles (WeakNormalizer, Bronstein 6.1.1) unhandled",
            "HC-A26-PRIMITIVE-PARAMQ-RATIONAL");
    };
    if (m == 0U) return fail_unimpl("empty forcing vector");

    // Split f and each g_i into numerator/denominator.
    auto fp = algebra::apart_num_den(f_expr, ctx);
    if (fp.is_error()) return fail_unimpl("cannot split f into num/den");
    ExprPtr fn = fp.value().numerator, fd = fp.value().denominator;
    std::vector<ExprPtr> gn(m), gd(m);
    for (std::size_t i = 0; i < m; ++i) {
        auto gp = algebra::apart_num_den(g_vec[i], ctx);
        if (gp.is_error()) return fail_unimpl("cannot split g_i into num/den");
        gn[i] = gp.value().numerator;
        gd[i] = gp.value().denominator;
    }

    // D = lcm(fd, gd_0, ..., gd_{m-1}).
    ExprPtr D_expr = fd;
    for (std::size_t i = 0; i < m; ++i) {
        ExprPtr prod = arena.make<Binary>(BinaryOp::Mul, D_expr, gd[i]);
        auto g = algebra::polynomial_gcd(D_expr, gd[i], var, ctx);
        if (g.is_ok()) {
            auto q = ctx.simplify(arena.make<Binary>(BinaryOp::Div, prod, g.value()));
            if (q.is_error()) return fail_unimpl("lcm denominator failed");
            D_expr = q.value();
        } else {
            auto s = ctx.simplify(prod);
            if (s.is_error()) return fail_unimpl("lcm denominator failed");
            D_expr = s.value();
        }
    }

    // Fn = fn·(D/fd) ; Gn_i = gn_i·(D/gd_i)  (polynomial since fd, gd_i | D).
    auto mul_div = [&](ExprPtr num, ExprPtr den) -> Result<ExprPtr> {
        auto q = ctx.simplify(arena.make<Binary>(BinaryOp::Div, D_expr, den));
        if (q.is_error()) return q;
        return ctx.simplify(arena.make<Binary>(BinaryOp::Mul, num, q.value()));
    };
    auto Fn_e = mul_div(fn, fd);
    if (Fn_e.is_error()) return fail_unimpl("Fn build failed");
    std::vector<std::vector<Rational>> Gn_c(m);
    for (std::size_t i = 0; i < m; ++i) {
        auto r = mul_div(gn[i], gd[i]);
        if (r.is_error()) return fail_unimpl("Gn_i build failed");
        auto c = poly_coeffs_q(r.value(), var, ctx);
        if (!c) return fail_unimpl("Gn_i not polynomial over Q");
        Gn_c[i] = std::move(*c);
    }

    auto D_c = poly_coeffs_q(D_expr, var, ctx);
    auto Fn_c = poly_coeffs_q(Fn_e.value(), var, ctx);
    if (!D_c || !Fn_c) return fail_unimpl("D or Fn not polynomial over Q");
    std::vector<Rational> D = std::move(*D_c), Fn = std::move(*Fn_c);
    strip_trailing_q(D);
    strip_trailing_q(Fn);
    if (D.empty()) return fail_unimpl("zero denominator");

    // D'.
    std::vector<Rational> Dp;
    if (D.size() >= 2U) {
        Dp.assign(D.size() - 1U, Rational(BigInt(0)));
        for (std::size_t k = 1; k < D.size(); ++k)
            Dp[k - 1U] = Rational(BigInt(static_cast<std::int64_t>(k))) * D[k];
    }
    // H = Fn − D'.
    std::vector<Rational> H(std::max(Fn.size(), Dp.size()), Rational(BigInt(0)));
    for (std::size_t i = 0; i < Fn.size(); ++i) H[i] = H[i] + Fn[i];
    for (std::size_t i = 0; i < Dp.size(); ++i) H[i] = H[i] - Dp[i];
    strip_trailing_q(H);

    // Gd_i = Gn_i · D.
    std::vector<std::vector<Rational>> Gd(m);
    int deg_Gd_max = -1;
    for (std::size_t i = 0; i < m; ++i) {
        strip_trailing_q(Gn_c[i]);
        if (Gn_c[i].empty()) continue;
        std::vector<Rational> prod(Gn_c[i].size() + D.size() - 1U, Rational(BigInt(0)));
        for (std::size_t a = 0; a < Gn_c[i].size(); ++a)
            for (std::size_t b = 0; b < D.size(); ++b)
                prod[a + b] = prod[a + b] + Gn_c[i][a] * D[b];
        strip_trailing_q(prod);
        deg_Gd_max = std::max(deg_Gd_max, static_cast<int>(prod.size()) - 1);
        Gd[i] = std::move(prod);
    }

    const int deg_D = static_cast<int>(D.size()) - 1;
    const int deg_H = H.empty() ? -1 : static_cast<int>(H.size()) - 1;

    // Degree bound for P (dominant-degree analysis of the identity), mirroring
    // solve_risch_de_rational_q.  Over-bound is safe: back-substitution drops
    // any spurious candidate.
    int M_bound = (deg_H >= deg_D - 1) ? (deg_Gd_max - deg_H)
                                       : (deg_Gd_max - deg_D + 1);
    if (deg_H == deg_D - 1 && !D.empty() && !H.empty()) {
        const auto& D_d = D.back();
        const auto& H_h = H.back();
        if (!D_d.numerator().is_zero()) {
            auto m_hom_rat = - (H_h / D_d);
            if (m_hom_rat.denominator() == BigInt(1) && !m_hom_rat.numerator().is_negative()) {
                int m_hom = static_cast<int>(m_hom_rat.numerator().to_u64());
                if (m_hom > M_bound) {
                    M_bound = m_hom;
                }
            }
        }
    }
    if (M_bound < deg_D) M_bound = deg_D;   // cancellation slack
    if (M_bound < 0) M_bound = 0;
    const int cap = static_cast<int>(ctx.max_risch_rational_ansatz_degree());
    if (M_bound > cap)
        return fail_unimpl("degree bound for P exceeds ctx.max_risch_rational_ansatz_degree");

    const std::size_t n_P = static_cast<std::size_t>(M_bound) + 1U;
    const std::size_t n_unk = n_P + m;
    int max_eq = std::max({deg_Gd_max, deg_D + M_bound - 1,
                           deg_H >= 0 ? deg_H + M_bound : -1, deg_D - 1});
    if (max_eq < 0) max_eq = 0;
    const std::size_t n_eq = static_cast<std::size_t>(max_eq) + 1U;

    // Homogeneous system: columns 0..M_bound = P_k, columns n_P..n_P+m-1 = c_i.
    std::vector<std::vector<Rational>> Mtx(n_eq, std::vector<Rational>(n_unk, Rational(BigInt(0))));
    for (std::size_t j = 0; j < n_eq; ++j) {
        // D·P' : Σ_{p+k-1=j, k≥1} k·D_p·P_k.
        for (std::size_t k = 1; k <= static_cast<std::size_t>(M_bound); ++k) {
            std::size_t pn = j + 1U;
            if (pn < k) continue;
            std::size_t p = pn - k;
            if (p >= D.size()) continue;
            Mtx[j][k] = Mtx[j][k] + Rational(BigInt(static_cast<std::int64_t>(k))) * D[p];
        }
        // (Fn − D')·P = H·P : Σ_{p+k=j} H_p·P_k.
        for (std::size_t k = 0; k <= static_cast<std::size_t>(M_bound); ++k) {
            if (k > j) break;
            std::size_t p = j - k;
            if (p >= H.size()) continue;
            Mtx[j][k] = Mtx[j][k] + H[p];
        }
        // −Σ_i c_i (Gn_i·D).
        for (std::size_t i = 0; i < m; ++i)
            if (j < Gd[i].size()) Mtx[j][n_P + i] = Mtx[j][n_P + i] - Gd[i][j];
    }

    auto pivots = row_echelon(Mtx, n_unk);
    auto basis = null_space_basis(Mtx, pivots, n_unk);

    // Exact verification of a candidate (P, c): the residual polynomial
    // D·P' + H·P − Σ_i c_i (Gn_i·D)  must vanish identically.  This is the very
    // identity the null-space enforces, recomputed in pure Q-rational arithmetic
    // — exact and independent of the symbolic simplifier, so it stays sound even
    // for high-degree / many-forcing systems.  Sound by construction (REGOLA
    // ZERO): a candidate failing this is dropped.
    auto residual_is_zero = [&](const std::vector<Rational>& P,
                                const std::vector<Rational>& c) -> bool {
        std::size_t maxd = D.size() + P.size() + 2U;
        for (std::size_t i = 0; i < m; ++i) maxd = std::max(maxd, Gd[i].size());
        for (std::size_t j = 0; j < maxd; ++j) {
            Rational r(BigInt(0));
            for (std::size_t k = 1; k < P.size(); ++k) {           // D·P'
                if (P[k].numerator().is_zero()) continue;
                std::size_t pn = j + 1U;
                if (pn < k) continue;
                std::size_t p = pn - k;
                if (p >= D.size()) continue;
                r = r + Rational(BigInt(static_cast<std::int64_t>(k))) * D[p] * P[k];
            }
            for (std::size_t k = 0; k < P.size(); ++k) {           // H·P
                if (P[k].numerator().is_zero() || k > j) continue;
                std::size_t p = j - k;
                if (p >= H.size()) continue;
                r = r + H[p] * P[k];
            }
            for (std::size_t i = 0; i < m; ++i) {                  // −Σ c_i Gd_i
                if (c[i].numerator().is_zero() || j >= Gd[i].size()) continue;
                r = r - c[i] * Gd[i][j];
            }
            if (!r.numerator().is_zero()) return false;
        }
        return true;
    };

    std::vector<ParametricRischDeQSolution> out;
    for (auto& v : basis) {
        std::vector<Rational> P(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(n_P));
        std::vector<Rational> c(v.begin() + static_cast<std::ptrdiff_t>(n_P), v.end());
        strip_trailing_q(P);
        if (!residual_is_zero(P, c)) continue;

        // Emit y as a cleared-integer-numerator fraction  (L·P)(x) / (L·D(x)),
        // where L is a common denominator of P's coefficients — a canonical
        // integer-coefficient representation of the (exactly verified) solution.
        BigInt L(1);
        for (const auto& pk : P) L = L * pk.denominator();
        std::vector<Rational> P_int(P.size(), Rational(BigInt(0)));
        for (std::size_t k = 0; k < P.size(); ++k) P_int[k] = P[k] * Rational(L);
        ExprPtr num = build_poly_expr_q(P_int, var, arena);
        ExprPtr den = arena.make<Binary>(BinaryOp::Mul,
                          arena.make<IntegerLit>(L), D_expr);
        ExprPtr y = arena.make<Binary>(BinaryOp::Div, num, den);
        out.push_back({y, std::move(c)});
    }
    return ok(std::move(out));
}

}  // namespace cas::calculus
