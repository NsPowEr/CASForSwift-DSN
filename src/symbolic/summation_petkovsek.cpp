// F5.7 — Petkovšek closed form for a higher-order linear recurrence.
//
// Given  Σ_{i=0}^J p[i](n)·S(n+i)=0  and the J initial values S(n0..n0+J−1),
// expresses S(n) = Σ c_i h_i(n) as a ℚ-linear combination of the recurrence's
// hypergeometric solutions h_i (found by hyper_all_ratios, each reconstructed
// as z^{n−n0}·∏Pochhammer terms).  The constants c_i are fitted from the
// initial values and the result is VERIFIED by extending S through the
// recurrence well beyond the fitted points — so a sequence that does not lie in
// the hypergeometric span yields ok(nullopt), never a wrong closed form.

#include "summation_hyper.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <vector>

namespace cas::symbolic {
namespace {

ExprPtr int_lit(AstArena& a, long long v) { return a.make<IntegerLit>(BigInt(v)); }
ExprPtr rat_lit(AstArena& a, const Rational& r) {
    if (r.denominator() == BigInt(1)) return a.make<IntegerLit>(r.numerator());
    return a.make<RationalLit>(r.numerator(), r.denominator());
}

std::optional<Rational> as_rational(ExprPtr e) {
    if (const auto* il = expr_cast<IntegerLit>(e)) return Rational(il->value);
    if (const auto* rl = expr_cast<RationalLit>(e)) return Rational(rl->numerator, rl->denominator);
    if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
        if (const auto* il2 = expr_cast<IntegerLit>(un->operand)) return Rational(-il2->value);
        if (const auto* rl2 = expr_cast<RationalLit>(un->operand)) return Rational(-rl2->numerator, rl2->denominator);
    }
    return std::nullopt;
}

// Exact value of the definite sum  Σ_{k=k_lo}^{n_val} F(n_val, k)  (upper = n).
std::optional<Rational> eval_definite_sum(
    ExprPtr F, const Symbol& n_sym, const Symbol& k_sym,
    long long n_val, long long k_lo, symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();
    auto Fn = ctx.substitute(F, n_sym, int_lit(a, n_val));
    if (Fn.is_error()) return std::nullopt;
    ExprPtr acc = int_lit(a, 0);
    for (long long kv = k_lo; kv <= n_val; ++kv) {
        auto Fk = ctx.substitute(Fn.value(), k_sym, int_lit(a, kv));
        if (Fk.is_error()) return std::nullopt;
        acc = a.make<Binary>(BinaryOp::Add, acc, Fk.value());
    }
    auto s = ctx.simplify(acc);
    if (s.is_error()) return std::nullopt;
    return as_rational(s.value());
}

// Evaluate a univariate expression at n = v, returning an exact rational.
std::optional<Rational> eval_rational(ExprPtr e, const Symbol& n, long long v,
                                      symbolic::CASContext& ctx) {
    auto sub = ctx.substitute(e, n, int_lit(ctx.arena(), v));
    if (sub.is_error()) return std::nullopt;
    auto s = ctx.simplify(sub.value());
    if (s.is_error()) return std::nullopt;
    return as_rational(s.value());
}

// ∏_{m=from}^{to-1} ρ(m) as an exact rational (1 when from==to).
std::optional<Rational> ratio_product(ExprPtr rho, const Symbol& n,
                                      long long from, long long to,
                                      symbolic::CASContext& ctx) {
    Rational acc(1);
    for (long long m = from; m < to; ++m) {
        auto v = eval_rational(rho, n, m, ctx);
        if (!v) return std::nullopt;
        acc = acc * (*v);
    }
    return acc;
}

// Rational roots (with multiplicity) of a monic-or-not polynomial in n, plus its
// rational leading coefficient.  Fails (nullopt) if the polynomial does not split
// into rational linear factors — i.e. it has a degree-≥1 irreducible factor or an
// irrational root.
struct SplitPoly { Rational lead; std::vector<Rational> roots; };
std::optional<SplitPoly> split_rational(ExprPtr poly, const Symbol& n,
                                        symbolic::CASContext& ctx) {
    auto cs = algebra::univariate_coefficients(poly, n, ctx);
    if (cs.is_error()) return std::nullopt;
    std::vector<Rational> coeffs;
    for (ExprPtr c : cs.value()) {
        auto sc = ctx.simplify(c);
        auto r = as_rational(sc.is_ok() ? sc.value() : c);
        if (!r) return std::nullopt;
        coeffs.push_back(*r);
    }
    while (coeffs.size() > 1U && coeffs.back().numerator().is_zero()) coeffs.pop_back();
    if (coeffs.empty()) return std::nullopt;

    SplitPoly sp;
    sp.lead = coeffs.back();
    const std::size_t deg = coeffs.size() - 1U;
    if (deg == 0U) return sp;  // nonzero constant

    auto roots = algebra::solve_polynomial(poly, n, ctx);
    if (roots.is_error()) return std::nullopt;
    for (ExprPtr r : roots.value()) {
        auto sr = ctx.simplify(r);
        auto rr = as_rational(sr.is_ok() ? sr.value() : r);
        if (!rr) return std::nullopt;        // irrational/algebraic root
        sp.roots.push_back(*rr);
    }
    if (sp.roots.size() != deg) return std::nullopt;  // did not fully split over ℚ
    return sp;
}

// Reconstruct h(n) with h(n0)=1 and h(n+1)/h(n)=ρ(n):
//   h(n) = z^{n−n0} · ∏_i Poch(n0−r_i, n−n0) / ∏_j Poch(n0−s_j, n−n0),
// with z = lc(num)/lc(den), r_i roots of num, s_j roots of den.  nullopt if ρ
// does not factor into rational linear factors.
std::optional<ExprPtr> ratio_to_term(ExprPtr rho, const Symbol& n, long long n0,
                                      symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();
    auto tog = algebra::together(rho, ctx);
    ExprPtr combined = tog.is_ok() ? tog.value() : rho;
    auto parts = algebra::apart_num_den(combined, ctx);
    if (parts.is_error()) return std::nullopt;
    auto sn = split_rational(parts.value().numerator, n, ctx);
    auto sd = split_rational(parts.value().denominator, n, ctx);
    if (!sn || !sd) return std::nullopt;

    Rational z = sn->lead / sd->lead;
    ExprPtr count = a.make<Binary>(BinaryOp::Sub, a.make<Symbol>(n), int_lit(a, n0));

    auto poch = [&](const Rational& root) {
        Rational base = Rational(BigInt(n0)) - root;   // n0 − root
        return a.make<FuncCall>(BuiltinOp::Pochhammer,
            std::vector<ExprPtr>{rat_lit(a, base), count});
    };

    ExprPtr term = a.make<Binary>(BinaryOp::Pow, rat_lit(a, z), count);  // z^{n−n0}
    for (const Rational& r : sn->roots)
        term = a.make<Binary>(BinaryOp::Mul, term, poch(r));
    for (const Rational& s : sd->roots)
        term = a.make<Binary>(BinaryOp::Div, term, poch(s));

    auto simp = ctx.simplify(term);
    return simp.is_ok() ? simp.value() : term;
}

// Solve the m×m rational system  M·c = b  by exact Gaussian elimination.
std::optional<std::vector<Rational>> solve_rational_system(
    std::vector<std::vector<Rational>> M, std::vector<Rational> b) {
    const std::size_t m = M.size();
    if (m == 0U) return std::vector<Rational>{};
    for (std::size_t col = 0U; col < m; ++col) {
        std::size_t piv = col;
        while (piv < m && M[piv][col].numerator().is_zero()) ++piv;
        if (piv == m) return std::nullopt;            // singular
        std::swap(M[piv], M[col]); std::swap(b[piv], b[col]);
        Rational inv = Rational(1) / M[col][col];
        for (std::size_t j = 0U; j < m; ++j) M[col][j] = M[col][j] * inv;
        b[col] = b[col] * inv;
        for (std::size_t r = 0U; r < m; ++r) {
            if (r == col || M[r][col].numerator().is_zero()) continue;
            Rational f = M[r][col];
            for (std::size_t j = 0U; j < m; ++j) M[r][j] = M[r][j] - f * M[col][j];
            b[r] = b[r] - f * b[col];
        }
    }
    return b;
}

// Value of a closed form S(n) at n = v, as an exact rational.
std::optional<Rational> eval_closed_form(ExprPtr S, const Symbol& n, long long v,
                                         symbolic::CASContext& ctx) {
    return eval_rational(S, n, v, ctx);
}

}  // namespace

Result<std::optional<ExprPtr>> solve_recurrence_closed_form(
    const std::vector<ExprPtr>& p, const std::vector<ExprPtr>& init,
    const Symbol& n, long long n0, symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();
    using OptE = std::optional<ExprPtr>;

    if (p.size() < 2U) return ok(OptE{std::nullopt});
    const std::size_t J = p.size() - 1U;
    if (init.size() < J) return ok(OptE{std::nullopt});

    // Hypergeometric solutions of the recurrence.
    std::vector<ExprPtr> ratios;
    auto hr = hyper_all_ratios(p, n, ratios, ctx);
    if (hr.is_error()) return fail<OptE>(hr.error());
    if (ratios.empty()) return ok(OptE{std::nullopt});

    // Reconstruct each solution h_i and keep those expressible over ℚ.
    std::vector<ExprPtr> rho_used, h_used;
    for (ExprPtr rho : ratios) {
        auto h = ratio_to_term(rho, n, n0, ctx);
        if (!h) continue;
        rho_used.push_back(rho);
        h_used.push_back(*h);
    }
    const std::size_t m = rho_used.size();
    if (m == 0U) return ok(OptE{std::nullopt});

    // Exact sequence values S(n0+t): the J supplied initial values, extended via
    // the recurrence  S(n0+J+s) = −(Σ_{i<J} p_i(n0+s)·S(n0+s+i)) / p_J(n0+s).
    std::vector<Rational> S;
    for (std::size_t t = 0U; t < J; ++t) {
        auto v = eval_rational(init[t], n, n0, ctx);  // init is a value; n absent ⇒ identity
        if (!v) return ok(OptE{std::nullopt});
        S.push_back(*v);
    }
    const std::size_t need = m + J + 2U;  // fit m, verify the rest
    for (std::size_t s = 0U; S.size() < need; ++s) {
        auto pJ = eval_rational(p[J], n, n0 + static_cast<long long>(s), ctx);
        if (!pJ || pJ->numerator().is_zero()) return ok(OptE{std::nullopt});
        Rational acc(0);
        bool ok_row = true;
        for (std::size_t i = 0U; i < J; ++i) {
            auto pi = eval_rational(p[i], n, n0 + static_cast<long long>(s), ctx);
            if (!pi) { ok_row = false; break; }
            acc = acc + (*pi) * S[s + i];
        }
        if (!ok_row) return ok(OptE{std::nullopt});
        S.push_back((Rational(-1) * acc) / (*pJ));
    }

    // Numeric value of h_i at n0+t via its term ratio (independent of Pochhammer
    // evaluation): h_i(n0+t) = ∏_{x=n0}^{n0+t−1} ρ_i(x).
    auto h_val = [&](std::size_t i, std::size_t t) -> std::optional<Rational> {
        return ratio_product(rho_used[i], n, n0, n0 + static_cast<long long>(t), ctx);
    };

    // Fit constants: M[t][i] = h_i(n0+t), b[t] = S(n0+t), t = 0..m−1.
    std::vector<std::vector<Rational>> Mtx(m, std::vector<Rational>(m));
    std::vector<Rational> b(m);
    for (std::size_t t = 0U; t < m; ++t) {
        for (std::size_t i = 0U; i < m; ++i) {
            auto hv = h_val(i, t);
            if (!hv) return ok(OptE{std::nullopt});
            Mtx[t][i] = *hv;
        }
        b[t] = S[t];
    }
    auto c = solve_rational_system(Mtx, b);
    if (!c) return ok(OptE{std::nullopt});

    // Verify Σ c_i h_i(n0+t) = S(n0+t) for ALL computed t (beyond the fitted m).
    for (std::size_t t = 0U; t < S.size(); ++t) {
        Rational sum(0);
        for (std::size_t i = 0U; i < m; ++i) {
            auto hv = h_val(i, t);
            if (!hv) return ok(OptE{std::nullopt});
            sum = sum + (*c)[i] * (*hv);
        }
        if (sum != S[t]) return ok(OptE{std::nullopt});  // not in hypergeometric span
    }

    // Build the symbolic closed form  Σ c_i h_i(n).
    ExprPtr result = nullptr;
    for (std::size_t i = 0U; i < m; ++i) {
        if ((*c)[i].numerator().is_zero()) continue;
        ExprPtr term = a.make<Binary>(BinaryOp::Mul, rat_lit(a, (*c)[i]), h_used[i]);
        result = result ? a.make<Binary>(BinaryOp::Add, result, term) : term;
    }
    if (!result) result = int_lit(a, 0);
    auto simp = ctx.simplify(result);
    return ok(OptE{simp.is_ok() ? simp.value() : result});
}

Result<std::optional<ExprPtr>> sum_closed_form_from_recurrence(
    const std::vector<ExprPtr>& p, ExprPtr F,
    const Symbol& n, const Symbol& k, ExprPtr lower, symbolic::CASContext& ctx) {
    using OptE = std::optional<ExprPtr>;
    if (p.size() < 2U) return ok(OptE{std::nullopt});
    const std::size_t J = p.size() - 1U;

    const auto* lo_il = expr_cast<IntegerLit>(lower);
    if (!lo_il) return ok(OptE{std::nullopt});   // direct summation needs integer lower
    const long long lo = lo_il->value.to_double();

    // Initial values S(lo+t) = Σ_{k=lo}^{lo+t} F(lo+t, k), t = 0..J−1.
    std::vector<ExprPtr> init;
    for (std::size_t t = 0U; t < J; ++t) {
        auto v = eval_definite_sum(F, n, k, lo + static_cast<long long>(t), lo, ctx);
        if (!v) return ok(OptE{std::nullopt});
        init.push_back(rat_lit(ctx.arena(), *v));
    }

    auto closed = solve_recurrence_closed_form(p, init, n, lo, ctx);
    if (closed.is_error()) return closed;
    if (!closed.value().has_value()) return ok(OptE{std::nullopt});
    ExprPtr S = *closed.value();

    // Cross-verify against directly-computed sums beyond the fitted points —
    // this is what guards against telescoping boundary-term corrections.
    for (long long extra = 0; extra < 3; ++extra) {
        const long long N = lo + static_cast<long long>(J) + extra;
        auto cf = eval_closed_form(S, n, N, ctx);
        auto ds = eval_definite_sum(F, n, k, N, lo, ctx);
        if (!cf || !ds || *cf != *ds) return ok(OptE{std::nullopt});
    }
    return ok(OptE{S});
}

}  // namespace cas::symbolic
