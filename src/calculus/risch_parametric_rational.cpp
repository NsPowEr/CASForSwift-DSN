// A26 / HC-A26-PRIMITIVE-PARAMQ-RATIONAL — parametric Risch DE over Q(x),
// f ≠ 0 case (full parametric Risch DE via the P/D ansatz, §6.1/§7.1).
// Anti-monolith split history: vertical split from risch_parametric.cpp; the
// f = 0 case (solve_param_limited_integration_rational_q) moved to
// risch_parametric_limited.cpp (A29, zero logic changes).
// SOUND BY CONSTRUCTION (REGOLA ZERO): every candidate is verified exactly via
// the Q-rational residual of the cleared identity and unverified candidates are
// dropped — at worst incomplete (diagnostic), never a wrong answer.

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

using detail::null_space_basis;
using detail::poly_coeffs_q;
using detail::rational_to_expr;
using detail::row_echelon;

void strip_trailing_q(std::vector<Rational>& v) {
    while (!v.empty() && v.back().numerator().is_zero()) v.pop_back();
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

// WeakNormalizer-style denominator inflation (Bronstein 6.1.1).  At a SIMPLE
// pole of f (squarefree factor s | den(f), multiplicity 1) with positive-integer
// residue n, the homogeneous solution y = (s)^−n has a pole of order n, but
// lcm(den f, den g_i) carries s only to order 1.  Residue n at the roots of
// s_n = gcd(fn − n·fd', s) (Rothstein-Trager), so multiply D by s_n^(n−1) to
// bring those solutions into the ansatz.  Order-≥2 poles of f are already
// covered by lcm and skipped.  Sound regardless: back-substitution rejects any
// non-solution, so an imperfect inflation only affects completeness.
[[nodiscard]] ExprPtr inflate_denominator(
    ExprPtr D_expr, ExprPtr fn, ExprPtr fd, const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    auto sqf = algebra::square_free_factorization(fd, var, ctx);
    auto fdp = diff(fd, var, 1U, ctx);
    if (sqf.is_error() || fdp.is_error()) return D_expr;
    const std::size_t cap = ctx.max_risch_rational_ansatz_degree();
    ExprPtr result = D_expr;
    for (const auto& sf : sqf.value().factors) {
        if (sf.multiplicity != 1U || !depends_on(sf.factor, var)) continue;
        for (std::size_t n = 2; n <= cap; ++n) {
            ExprPtr nfdp = arena.make<Binary>(BinaryOp::Mul,
                arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(n))), fdp.value());
            ExprPtr cand = arena.make<Binary>(BinaryOp::Sub, fn, nfdp);
            ExprPtr v;
            if (auto cs = ctx.simplify(cand); cs.is_ok() &&
                    expr_cast<IntegerLit>(cs.value()) &&
                    expr_cast<IntegerLit>(cs.value())->value.is_zero()) {
                v = sf.factor;                       // residue n on all of s
            } else {
                auto g = algebra::polynomial_gcd(cand, sf.factor, var, ctx);
                if (g.is_error() || !depends_on(g.value(), var)) continue;
                v = g.value();
            }
            for (std::size_t e = 0; e + 1U < n; ++e)
                result = arena.make<Binary>(BinaryOp::Mul, result, v);
        }
    }
    auto s = ctx.simplify(result);
    return s.is_ok() ? s.value() : result;
}

}  // namespace

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
// DENOMINATOR BOUND (WeakNormalizer, Bronstein 6.1.1).  The ansatz finds only
// solutions whose denominator divides D.  At a SIMPLE pole of f with positive-
// integer residue n the homogeneous solution y_h has a pole of order n while
// lcm(den f, den g_i) carries it only to order 1; inflate_denominator() lifts D
// by the Rothstein-Trager factor (residue n at the roots of gcd(fn−n·fd', s)) so
// those solutions enter the ansatz (verified below: WeakNormalizer_* tests).
// Residual incompleteness (never a wrong answer — back-substitution drops any
// non-solution): poles of f of MULTIPLICITY ≥ 2, residues that are not integers
// in [2, cap], and algebraic (non-rational) residues.  The non-parametric
// solve_risch_de_rational_q does NOT yet inflate — it carries the order-1 gap.
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
            "Parametric Risch DE over Q(x): WeakNormalizer (Bronstein 6.1.1) "
            "inflates D for simple positive-integer-residue poles; residual cases "
            "are multiplicity-≥2 poles, residue beyond cap, and algebraic residues",
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

    // Inflate D to cover positive-integer-residue simple poles (WeakNormalizer).
    D_expr = inflate_denominator(D_expr, fn, fd, var, ctx);

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
        // Verify up to the exact degree of the residual D·P' + H·P − Σ c_i Gd_i:
        // deg(D·P')=deg_D+deg_P−1 (size D.size()+P.size()−1), deg(H·P)=deg_H+deg_P
        // (size H.size()+P.size()−1), and each Gd_i term.  Derived bound, no padding.
        std::size_t maxd = 0;
        if (!P.empty()) {
            maxd = std::max(maxd, D.size() + P.size() - 1U);
            if (!H.empty()) maxd = std::max(maxd, H.size() + P.size() - 1U);
        }
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
