// polynomial_gcd_brown.cpp — F3.1 GCD dispatch and legacy eval/interp-over-Z.
//
// HONEST naming (T3-Opus A2 fix):
//   * gcd_brown()          → public dispatcher: tries REAL modular (gcd_brown_modular
//                            in polynomial_gcd_brown_modular.cpp) first; on failure
//                            falls back to gcd_eval_interp_z below.
//   * gcd_eval_interp_z()  → Lagrange interpolation directly in Z[x_n] over recursive
//                            calls in Z[x_1..x_{n-1}].  Coefficient growth UNBOUNDED;
//                            kept ONLY as a robust fallback for inputs that defeat the
//                            modular path.  This is NOT Brown's modular algorithm.
//   * gcd_zippel_sparse()  → public dispatcher: tries REAL Prony sparse interpolation
//                            (gcd_zippel_prony) first; falls back to gcd_brown_modular
//                            and then to gcd_eval_interp_z.
//   * gcd_ez()             → cofactor-certified wrapper around the dispatch chain.
//
// Ref: Geddes-Czapor-Labahn §7.4–7.5.  BigInt-only; Result<T>; no throw.

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"
#include "cas/numtheory.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace cas::algebra {

namespace {

using Monomial   = std::vector<unsigned int>;
using SparsePoly = std::map<Monomial, BigInt>;

[[nodiscard]] static Result<IntPoly> mv_to_intpoly_b(const MultivariatePolynomial& poly,
                                                       const Symbol& var) {
    if (poly.is_zero()) return ok(IntPoly{});
    std::map<std::size_t, BigInt> cm;
    std::size_t md = 0U;
    for (const auto& term : poly.terms()) {
        std::size_t d = 0U;
        for (const auto& [sym, exp] : term.factors) {
            if (sym.name == var.name) d = exp;
            else return fail<IntPoly>(make_error(CASErrorKind::Unimplemented,
                "gcd_brown: expected univariate polynomial after specialization"));
        }
        cm[d] += term.coefficient; md = std::max(md, d);
    }
    IntPoly r; r.resize(md + 1U, BigInt(0));
    for (const auto& [d, c] : cm) r[d] = c;
    r.normalize([](const BigInt& v){ return v.is_zero(); });
    return ok(std::move(r));
}

[[nodiscard]] static MultivariatePolynomial intpoly_to_mv_b(const IntPoly& poly,
                                                              const Symbol& var) {
    std::vector<MultivariateTerm> terms;
    for (std::size_t d = 0; d < poly.size(); ++d) {
        if (poly[d].is_zero()) continue;
        std::vector<std::pair<Symbol, unsigned int>> f;
        if (d > 0U) f.emplace_back(var, static_cast<unsigned int>(d));
        terms.push_back(MultivariateTerm{ .coefficient = poly[d], .factors = std::move(f) });
    }
    return MultivariatePolynomial(std::move(terms));
}

static std::vector<Symbol> collect_vars_brown(const MultivariatePolynomial& p,
                                               const MultivariatePolynomial& q) {
    auto vp = p.variables();
    auto vq = q.variables();
    std::vector<Symbol> all = vp;
    for (const auto& s : vq)
        if (std::none_of(all.begin(), all.end(),
                         [&](const Symbol& a){ return a.name == s.name; }))
            all.push_back(s);
    std::sort(all.begin(), all.end(),
              [](const Symbol& a, const Symbol& b){ return a.name < b.name; });
    return all;
}

static std::size_t deg_in_var_b(const MultivariatePolynomial& p, const Symbol& v) {
    std::size_t d = 0;
    for (const auto& t : p.terms())
        for (const auto& [sym, exp] : t.factors)
            if (sym.name == v.name) d = std::max(d, static_cast<std::size_t>(exp));
    return d;
}

static BigInt bigint_gcd_b(BigInt a, BigInt b) {
    a = a.abs(); b = b.abs();
    while (!b.is_zero()) { BigInt t = a % b; a = std::move(b); b = std::move(t); }
    return a;
}

static BigInt poly_content_b(const MultivariatePolynomial& p) {
    BigInt c(0);
    for (const auto& t : p.terms()) c = bigint_gcd_b(c, t.coefficient.abs());
    return c.is_zero() ? BigInt(1) : c;
}

static MultivariatePolynomial make_primitive_b(const MultivariatePolynomial& p) {
    if (p.is_zero()) return p;
    BigInt c = poly_content_b(p);
    std::vector<MultivariateTerm> terms;
    terms.reserve(p.terms().size());
    for (const auto& t : p.terms()) {
        BigInt nc = t.coefficient / c;
        if (!nc.is_zero())
            terms.push_back(MultivariateTerm{ .coefficient = nc, .factors = t.factors });
    }
    if (!terms.empty()) {
        bool all_neg = true;
        for (const auto& t : terms) if (!t.coefficient.is_negative()) { all_neg = false; break; }
        if (all_neg) for (auto& t : terms) t.coefficient = -t.coefficient;
    }
    return MultivariatePolynomial(std::move(terms));
}

static SparsePoly to_sparse_b(const MultivariatePolynomial& p,
                               const std::vector<Symbol>& vars) {
    SparsePoly sp;
    for (const auto& term : p.terms()) {
        Monomial mono(vars.size(), 0U);
        for (const auto& [sym, exp] : term.factors)
            for (std::size_t i = 0; i < vars.size(); ++i)
                if (vars[i].name == sym.name) mono[i] = exp;
        sp[mono] += term.coefficient;
        if (sp[mono].is_zero()) sp.erase(mono);
    }
    return sp;
}

static MultivariatePolynomial from_sparse_b(const SparsePoly& sp,
                                             const std::vector<Symbol>& vars) {
    std::vector<MultivariateTerm> terms;
    for (const auto& [mono, coeff] : sp) {
        if (coeff.is_zero()) continue;
        std::vector<std::pair<Symbol, unsigned int>> factors;
        for (std::size_t i = 0; i < vars.size(); ++i)
            if (mono[i] > 0) factors.emplace_back(vars[i], mono[i]);
        terms.push_back(MultivariateTerm{ .coefficient = coeff, .factors = std::move(factors) });
    }
    return MultivariatePolynomial(std::move(terms));
}

[[nodiscard]] static std::optional<MultivariatePolynomial> exact_div_b(
        const MultivariatePolynomial& p,
        const MultivariatePolynomial& d,
        const std::vector<Symbol>& vars) {
    if (d.is_zero()) return std::nullopt;
    if (p.is_zero()) return MultivariatePolynomial{};
    SparsePoly rem = to_sparse_b(p, vars);
    SparsePoly div = to_sparse_b(d, vars);
    if (div.empty()) return std::nullopt;
    SparsePoly quo;
    const std::size_t budget = (rem.size() + 1U) * (div.size() + 1U) + 8U;
    const auto [dlm, dlc] = *std::prev(div.end());
    std::size_t steps = 0;
    while (!rem.empty()) {
        if (++steps > budget) return std::nullopt;
        const auto [rlm, rlc] = *std::prev(rem.end());
        bool mon_ok = true;
        for (std::size_t i = 0; i < vars.size(); ++i) {
            std::size_t ri = (i < rlm.size()) ? rlm[i] : 0U;
            std::size_t di = (i < dlm.size()) ? dlm[i] : 0U;
            if (ri < di) { mon_ok = false; break; }
        }
        if (!mon_ok) return std::nullopt;
        if ((rlc % dlc) != BigInt(0)) return std::nullopt;
        BigInt tc = rlc / dlc;
        Monomial tm(vars.size(), 0U);
        for (std::size_t i = 0; i < vars.size(); ++i) {
            std::size_t ri = (i < rlm.size()) ? rlm[i] : 0U;
            std::size_t di = (i < dlm.size()) ? dlm[i] : 0U;
            tm[i] = static_cast<unsigned int>(ri - di);
        }
        quo[tm] += tc;
        if (quo[tm].is_zero()) quo.erase(tm);
        for (const auto& [dm, dc] : div) {
            Monomial nm(vars.size(), 0U);
            for (std::size_t i = 0; i < vars.size(); ++i)
                nm[i] = tm[i] + ((i < dm.size()) ? dm[i] : 0U);
            rem[nm] -= tc * dc;
            if (rem[nm].is_zero()) rem.erase(nm);
        }
    }
    return from_sparse_b(quo, vars);
}

[[nodiscard]] static bool divides_b(const MultivariatePolynomial& p,
                                     const MultivariatePolynomial& d,
                                     const std::vector<Symbol>& vars) {
    return exact_div_b(p, d, vars).has_value();
}

[[nodiscard]] Result<MultivariatePolynomial> gcd_eval_interp_z_impl(
        const MultivariatePolynomial& P, const MultivariatePolynomial& Q,
        symbolic::CASContext& ctx, std::size_t depth);

[[nodiscard]] Result<MultivariatePolynomial> gcd_eval_interp_z_impl(
        const MultivariatePolynomial& P, const MultivariatePolynomial& Q,
        symbolic::CASContext& ctx, std::size_t depth) {

    if (P.is_zero()) return ok(make_primitive_b(Q));
    if (Q.is_zero()) return ok(make_primitive_b(P));

    const std::vector<Symbol> vars = collect_vars_brown(P, Q);

    if (vars.empty()) {
        BigInt g = bigint_gcd_b(P.terms()[0].coefficient.abs(), Q.terms()[0].coefficient.abs());
        return ok(MultivariatePolynomial({{.coefficient = g, .factors = {}}}));
    }
    if (vars.size() == 1U) {
        auto pi = mv_to_intpoly_b(P, vars[0]);
        if (pi.is_error()) return fail<MultivariatePolynomial>(pi.error());
        auto qi = mv_to_intpoly_b(Q, vars[0]);
        if (qi.is_error()) return fail<MultivariatePolynomial>(qi.error());
        IntPoly g = gcd_integer_poly_primitive(std::move(pi.value()), std::move(qi.value()));
        if (!g.empty() && g.leading_coeff().is_negative())
            for (auto& c : g.coefficients()) c = -c;
        return ok(intpoly_to_mv_b(g, vars[0]));
    }

    if (depth > ctx.max_gcd_recursion_depth())
        return make_unimplemented<MultivariatePolynomial>("algebra", "gcd_eval_interp_z",
            "depth=" + std::to_string(depth), "GCD_EVAL_INTERP_Z_DEPTH_EXCEEDED",
            "Increase ctx.max_gcd_recursion_depth() (F3.1)", "F3.1");

    BigInt cont_g = bigint_gcd_b(poly_content_b(P), poly_content_b(Q));
    MultivariatePolynomial pp = make_primitive_b(P);
    MultivariatePolynomial pq = make_primitive_b(Q);

    const Symbol& eval_var   = vars.back();
    const std::vector<Symbol> sub_vars(vars.begin(), vars.end() - 1U);
    const std::size_t deg_bound =
        std::min(deg_in_var_b(pp, eval_var), deg_in_var_b(pq, eval_var));

    const double delta = ctx.gcd_error_probability();
    const std::size_t log2_inv = (delta > 0.0 && delta < 1.0)
        ? static_cast<std::size_t>(std::ceil(std::log2(1.0 / delta))) : std::size_t{10U};
    const std::size_t extra    = std::max(std::size_t{2U}, log2_inv);
    const std::size_t need_pts = deg_bound + 1U + extra;

    if (need_pts > ctx.max_gcd_total_calls())
        return make_unimplemented<MultivariatePolynomial>("algebra", "gcd_eval_interp_z",
            "need_pts=" + std::to_string(need_pts), "GCD_EVAL_INTERP_Z_POINTS_BUDGET_EXCEEDED",
            "Increase ctx.max_gcd_total_calls() (F3.1)", "F3.1");

    struct EvalSample { BigInt pt; MultivariatePolynomial gcd; };
    std::vector<EvalSample> lucky;
    lucky.reserve(need_pts);
    std::optional<std::vector<std::size_t>> expected_sub_degs;

    const long long max_pt = static_cast<long long>(need_pts * 4LL + 50LL);
    for (long long pt = 1LL; static_cast<std::size_t>(lucky.size()) < need_pts && pt <= max_pt; ++pt) {
        BigInt val(pt);
        ExprPtr val_expr = ctx.arena().make<IntegerLit>(val);
        auto pp_ev = pp.evaluate_at(eval_var, val_expr);
        if (pp_ev.is_error()) continue;
        auto pq_ev = pq.evaluate_at(eval_var, val_expr);
        if (pq_ev.is_error()) continue;
        auto g_ev = gcd_eval_interp_z_impl(pp_ev.value(), pq_ev.value(), ctx, depth + 1U);
        if (g_ev.is_error()) continue;
        MultivariatePolynomial ge = make_primitive_b(g_ev.value());
        std::vector<std::size_t> sub_degs;
        for (const auto& sv : sub_vars) sub_degs.push_back(deg_in_var_b(ge, sv));
        if (!expected_sub_degs.has_value()) {
            expected_sub_degs = sub_degs;
        } else {
            bool lower = false;
            for (std::size_t i = 0; i < sub_degs.size() && i < expected_sub_degs->size(); ++i)
                if (sub_degs[i] < (*expected_sub_degs)[i]) { lower = true; break; }
            if (lower) { expected_sub_degs = sub_degs; lucky.clear(); }
            if (sub_degs != *expected_sub_degs) continue;
        }
        lucky.push_back({ .pt = val, .gcd = std::move(ge) });
    }

    if (lucky.size() < deg_bound + 1U)
        return make_unimplemented<MultivariatePolynomial>("algebra", "gcd_eval_interp_z",
            "lucky=" + std::to_string(lucky.size()) + " < need=" + std::to_string(deg_bound + 1U),
            "GCD_EVAL_INTERP_Z_NOT_ENOUGH_LUCKY_POINTS",
            "All evaluation points unlucky — adversarial input or large degree (F3.1)", "F3.1");

    std::set<Monomial> all_monomials_set;
    for (const auto& s : lucky) {
        SparsePoly gsp = to_sparse_b(s.gcd, sub_vars);
        for (const auto& [mono, _coeff] : gsp) all_monomials_set.insert(mono);
    }

    const std::size_t n_pts = std::min(lucky.size(), deg_bound + 1U);
    std::vector<BigInt> pts;
    pts.reserve(n_pts);
    for (std::size_t i = 0; i < n_pts; ++i) pts.push_back(lucky[i].pt);

    std::vector<BigInt> denoms(n_pts);
    BigInt lcm_den(1);
    for (std::size_t i = 0; i < n_pts; ++i) {
        BigInt d(1);
        for (std::size_t j = 0; j < n_pts; ++j) if (j != i) d *= (pts[i] - pts[j]);
        if (d.is_zero())
            return make_unimplemented<MultivariatePolynomial>("algebra", "gcd_eval_interp_z",
                "n_pts=" + std::to_string(n_pts), "GCD_EVAL_INTERP_Z_DUPLICATE_EVAL_POINTS",
                "Duplicate evaluation points — internal error (F3.1)", "F3.1");
        denoms[i] = d;
        BigInt ad = d.abs();
        BigInt g2 = bigint_gcd_b(lcm_den.abs(), ad);
        lcm_den = lcm_den * (ad / g2);
    }

    std::vector<std::vector<BigInt>> num_polys(n_pts);
    for (std::size_t i = 0; i < n_pts; ++i) {
        num_polys[i] = { BigInt(1) };
        for (std::size_t j = 0; j < n_pts; ++j) {
            if (j == i) continue;
            std::vector<BigInt> np(num_polys[i].size() + 1U, BigInt(0));
            for (std::size_t k = 0; k < num_polys[i].size(); ++k) {
                np[k + 1U] += num_polys[i][k];
                np[k]      -= num_polys[i][k] * pts[j];
            }
            num_polys[i] = std::move(np);
        }
        BigInt scale = lcm_den / denoms[i];
        for (auto& c : num_polys[i]) c *= scale;
    }

    std::vector<MultivariateTerm> candidate_terms;
    for (const auto& mono : all_monomials_set) {
        std::vector<BigInt> coeffs_scaled(n_pts, BigInt(0));
        for (std::size_t i = 0; i < n_pts; ++i) {
            SparsePoly gsp = to_sparse_b(lucky[i].gcd, sub_vars);
            auto it = gsp.find(mono);
            BigInt val_i = (it != gsp.end()) ? it->second : BigInt(0);
            if (val_i.is_zero()) continue;
            for (std::size_t k = 0; k < num_polys[i].size() && k < coeffs_scaled.size(); ++k)
                coeffs_scaled[k] += num_polys[i][k] * val_i;
        }
        for (std::size_t k = 0; k < coeffs_scaled.size(); ++k) {
            if (coeffs_scaled[k].is_zero()) continue;
            if ((coeffs_scaled[k] % lcm_den) != BigInt(0)) continue;
            BigInt int_coeff = coeffs_scaled[k] / lcm_den;
            if (int_coeff.is_zero()) continue;
            std::vector<std::pair<Symbol, unsigned int>> factors;
            for (std::size_t vi = 0; vi < sub_vars.size(); ++vi)
                if (vi < mono.size() && mono[vi] > 0)
                    factors.emplace_back(sub_vars[vi], mono[vi]);
            if (k > 0) factors.emplace_back(eval_var, static_cast<unsigned int>(k));
            candidate_terms.push_back(MultivariateTerm{ .coefficient = int_coeff,
                                                         .factors = std::move(factors) });
        }
    }

    MultivariatePolynomial candidate = make_primitive_b(
        MultivariatePolynomial(std::move(candidate_terms)));

    if (candidate.is_zero())
        return ok(MultivariatePolynomial({{.coefficient = BigInt(1), .factors = {}}}));

    if (!divides_b(pp, candidate, vars) || !divides_b(pq, candidate, vars))
        return make_unimplemented<MultivariatePolynomial>("algebra", "gcd_eval_interp_z",
            "vars=" + std::to_string(vars.size()), "GCD_EVAL_INTERP_Z_CERTIFICATE_FAILED",
            "Interpolated candidate fails divisibility — caller should fall back to GCDHEU (F3.1)",
            "F3.1");

    if (cont_g != BigInt(1) && !cont_g.is_zero()) {
        std::vector<MultivariateTerm> terms;
        for (const auto& t : candidate.terms())
            terms.push_back(MultivariateTerm{
                .coefficient = t.coefficient * cont_g, .factors = t.factors });
        return ok(MultivariatePolynomial(std::move(terms)));
    }
    return ok(candidate);
}

}  // namespace

// Honest legacy entry: Lagrange interpolation over Z (coefficient growth UNBOUNDED).
// NOT Brown's modular GCD.  Kept as a robust fallback path.
Result<MultivariatePolynomial> gcd_eval_interp_z(
        const MultivariatePolynomial& P, const MultivariatePolynomial& Q,
        symbolic::CASContext& ctx) {
    return gcd_eval_interp_z_impl(P, Q, ctx, 0U);
}

// Public dispatcher: tries REAL Brown's modular GCD first (no Z coefficient growth),
// falls back to the legacy eval/interp-Z path on failure.  Both certify divisibility.
Result<MultivariatePolynomial> gcd_brown(
        const MultivariatePolynomial& P, const MultivariatePolynomial& Q,
        symbolic::CASContext& ctx) {
    auto modular = gcd_brown_modular(P, Q, ctx);
    if (modular.is_ok()) return modular;
    return gcd_eval_interp_z_impl(P, Q, ctx, 0U);
}

// Public sparse dispatcher: tries REAL Zippel Prony first; on failure cascades to
// Brown's modular, then to the legacy eval/interp-Z path.  Every transition emits
// a diagnostic via the per-stage Result<T> error path.
Result<MultivariatePolynomial> gcd_zippel_sparse(
        const MultivariatePolynomial& P, const MultivariatePolynomial& Q,
        symbolic::CASContext& ctx) {
    bool run_zippel = false;
    std::vector<Symbol> vars = P.variables();
    for (const auto& s : Q.variables()) {
        if (std::none_of(vars.begin(), vars.end(), [&](const Symbol& v) { return v.name == s.name; })) {
            vars.push_back(s);
        }
    }
    if (vars.size() >= 2U) {
        std::size_t n = vars.size();
        std::size_t d = std::max(P.total_degree(), Q.total_degree());
        auto bin_res = numtheory::binomial(BigInt(n + d), BigInt(d));
        double M = 1.0;
        if (bin_res.is_ok()) {
            M = bin_res.value().to_double();
        }
        double actual_terms = static_cast<double>(std::max<std::size_t>(1U, std::max(P.terms().size(), Q.terms().size())));
        double density_ratio = M / actual_terms;
        if (density_ratio >= ctx.zippel_density_threshold()) {
            run_zippel = true;
        }
    }

    if (run_zippel) {
        auto prony = gcd_zippel_prony(P, Q, ctx);
        if (prony.is_ok()) return prony;
    }
    auto modular = gcd_brown_modular(P, Q, ctx);
    if (modular.is_ok()) return modular;
    return gcd_eval_interp_z_impl(P, Q, ctx, 0U);
}

Result<GcdWithCofactors> gcd_ez(
        const MultivariatePolynomial& P, const MultivariatePolynomial& Q,
        symbolic::CASContext& ctx) {
    if (P.is_zero())
        return ok(GcdWithCofactors{ .gcd = Q, .cofactor_p = MultivariatePolynomial{},
            .cofactor_q = MultivariatePolynomial({{.coefficient = BigInt(1), .factors = {}}}), });
    if (Q.is_zero())
        return ok(GcdWithCofactors{ .gcd = P,
            .cofactor_p = MultivariatePolynomial({{.coefficient = BigInt(1), .factors = {}}}),
            .cofactor_q = MultivariatePolynomial{}, });

    const std::vector<Symbol> vars = collect_vars_brown(P, Q);
    // Prefer REAL Brown's modular GCD; on failure fall back to legacy paths.
    auto g_res = gcd_brown_modular(P, Q, ctx);
    if (g_res.is_error()) g_res = gcd_eval_interp_z_impl(P, Q, ctx, 0U);
    if (g_res.is_error()) g_res = gcd_multivariate_eval_interp(P, Q, ctx);
    if (g_res.is_error()) return fail<GcdWithCofactors>(g_res.error());
    MultivariatePolynomial g = std::move(g_res.value());

    auto cofp = exact_div_b(P, g, vars);
    if (!cofp.has_value())
        return make_unimplemented<GcdWithCofactors>("algebra", "gcd_ez",
            "vars=" + std::to_string(vars.size()), "GCD_EZ_COFACTOR_P_NOT_EXACT",
            "P / gcd(P,Q) is not exact — GCD may be incorrect (F3.1)", "F3.1");
    auto cofq = exact_div_b(Q, g, vars);
    if (!cofq.has_value())
        return make_unimplemented<GcdWithCofactors>("algebra", "gcd_ez",
            "vars=" + std::to_string(vars.size()), "GCD_EZ_COFACTOR_Q_NOT_EXACT",
            "Q / gcd(P,Q) is not exact — GCD may be incorrect (F3.1)", "F3.1");

    auto sp_gp = to_sparse_b(g * cofp.value(), vars);
    auto sp_p  = to_sparse_b(P, vars);
    auto sp_gq = to_sparse_b(g * cofq.value(), vars);
    auto sp_q  = to_sparse_b(Q, vars);

    if (sp_gp != sp_p || sp_gq != sp_q)
        return make_unimplemented<GcdWithCofactors>("algebra", "gcd_ez",
            "cert_p=" + std::to_string(sp_gp == sp_p) + ",cert_q=" + std::to_string(sp_gq == sp_q),
            "GCD_EZ_CERTIFICATE_FAILED",
            "g * cofactor ≠ input — GCD or cofactor computation error (F3.1)", "F3.1");

    return ok(GcdWithCofactors{
        .gcd        = std::move(g),
        .cofactor_p = std::move(cofp.value()),
        .cofactor_q = std::move(cofq.value()),
    });
}

}  // namespace cas::algebra
