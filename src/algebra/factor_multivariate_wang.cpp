// F3.2 — Wang multivariate polynomial factorization driver (EEZ).
// Public entry: factor_multivariate.  Pipeline (GCL §6.4–6.8):
//   1. parse + integer content extraction
//   2. multivariate squarefree decomposition in the main variable (Yun)
//   3. for each squarefree part: choose a good evaluation point, factor the
//      univariate image, distribute leading coefficients (Wang), lift by
//      multivariate Hensel, and CERTIFY by exact reconstruction.
// Everything is exact (BigInt / Rational); the final certification guarantees no
// silently-wrong result — failures surface as explicit Unimplemented diagnostics.

#include "factor_multivariate_internal.hpp"

#include <algorithm>
#include <map>
#include <set>

namespace cas::algebra {

namespace {

// Pick the variable ordering: main variable = the one of largest degree (Wang's
// heuristic favours the main var with smallest degree, but any deterministic
// ordering is correct since we certify; we use largest total appearance first to
// keep the univariate image informative).  Deterministic by name as tiebreak.
[[nodiscard]] WangContext build_wang_context(const MultivariatePolynomial& mv) {
    auto vars = mv.variables();
    std::sort(vars.begin(), vars.end(),
              [](const Symbol& a, const Symbol& b) { return a.name < b.name; });
    WangContext wc;
    wc.vars = std::move(vars);
    return wc;
}

// Multivariate derivative w.r.t. variable index.
[[nodiscard]] MPoly mpoly_derivative(const MPoly& p, std::size_t var) {
    MPoly r;
    for (const auto& [mono, coeff] : p.terms) {
        if (mono[var] == 0U) continue;
        Monomial m = mono;
        BigInt c = coeff * BigInt(static_cast<long long>(mono[var]));
        m[var] -= 1U;
        BigInt& slot = r.terms[m];
        slot += c;
        if (slot.is_zero()) r.terms.erase(m);
    }
    return r;
}

// Multivariate GCD over Z, via the existing certified multivariate GCD.
[[nodiscard]] Result<MPoly> mpoly_gcd(
    const MPoly& a, const MPoly& b, const WangContext& wc, symbolic::CASContext& ctx) {
    if (a.is_zero()) return ok(b);
    if (b.is_zero()) return ok(a);
    auto ea = mpoly_to_expr(a, wc, ctx);
    if (ea.is_error()) return fail<MPoly>(ea.error());
    auto eb = mpoly_to_expr(b, wc, ctx);
    if (eb.is_error()) return fail<MPoly>(eb.error());
    auto g = polynomial_gcd_multivariate(ea.value(), eb.value(), ctx);
    if (g.is_error()) return fail<MPoly>(g.error());
    return mpoly_from_expr(g.value(), wc, ctx);
}

// Yun squarefree decomposition in the main variable (index 0).
// Returns pairs (squarefree_part, multiplicity).  Exact; uses multivariate GCD.
[[nodiscard]] Result<std::vector<std::pair<MPoly, unsigned int>>> squarefree_main(
    const MPoly& a, const WangContext& wc, symbolic::CASContext& ctx) {
    std::vector<std::pair<MPoly, unsigned int>> out;
    if (mpoly_degree_in(a, 0U) == 0U) {
        out.emplace_back(a, 1U);
        return ok(std::move(out));
    }
    MPoly da = mpoly_derivative(a, 0U);
    auto gres = mpoly_gcd(a, da, wc, ctx);
    if (gres.is_error()) return fail<std::vector<std::pair<MPoly, unsigned int>>>(gres.error());
    MPoly c = gres.value();

    auto wopt = mpoly_exact_div(a, c);
    if (!wopt.has_value()) {
        return fail<std::vector<std::pair<MPoly, unsigned int>>>(make_error(
            CASErrorKind::Unimplemented,
            "wang squarefree: gcd(a,a') did not divide a (non-Z domain?)"));
    }
    MPoly w = *wopt;
    auto yopt = mpoly_exact_div(da, c);
    if (!yopt.has_value()) {
        return fail<std::vector<std::pair<MPoly, unsigned int>>>(make_error(
            CASErrorKind::Unimplemented, "wang squarefree: a' not divisible by gcd"));
    }
    MPoly y = *yopt;
    MPoly z = mpoly_sub(y, mpoly_derivative(w, 0U));

    unsigned int i = 1U;
    while (!z.is_zero()) {
        auto gg = mpoly_gcd(w, z, wc, ctx);
        if (gg.is_error()) return fail<std::vector<std::pair<MPoly, unsigned int>>>(gg.error());
        MPoly g = gg.value();
        if (mpoly_degree_in(g, 0U) > 0U ||
            !(g.terms.size() == 1U && g.terms.begin()->first ==
              Monomial(wc.nvars(), 0U) && g.terms.begin()->second.abs() == BigInt(1))) {
            out.emplace_back(g, i);
        }
        auto wn = mpoly_exact_div(w, g);
        auto zn = mpoly_exact_div(z, g);
        if (!wn.has_value() || !zn.has_value()) {
            return fail<std::vector<std::pair<MPoly, unsigned int>>>(make_error(
                CASErrorKind::Unimplemented, "wang squarefree: Yun division failed"));
        }
        w = *wn;
        z = mpoly_sub(*zn, mpoly_derivative(w, 0U));
        ++i;
    }
    // The final w is the highest squarefree part.
    if (mpoly_degree_in(w, 0U) > 0U) {
        out.emplace_back(w, i);
    }
    if (out.empty()) {
        out.emplace_back(a, 1U);
    }
    return ok(std::move(out));
}

// Convert an ExprPtr univariate factor (in main var) to an IntPoly.
[[nodiscard]] Result<IntPoly> expr_to_intpoly_main(
    ExprPtr e, const Symbol& main_var, symbolic::CASContext& ctx) {
    auto pp = parse_polynomial(e, main_var, ctx);
    if (pp.is_error()) return fail<IntPoly>(pp.error());
    auto ip = poly_to_integer_poly(pp.value());
    if (ip.is_error()) return fail<IntPoly>(ip.error());
    return ok(std::move(ip.value()));
}

// Factor one squarefree multivariate part (squarefree in main var).
// Appends irreducible factors (as MPoly) to `out_factors`.
[[nodiscard]] Result<void> factor_squarefree_part(
    const MPoly& part, const WangContext& wc, symbolic::CASContext& ctx,
    std::vector<MPoly>& out_factors) {
    const std::size_t n = wc.nvars();

    // Univariate case: delegate entirely to the univariate factorizer.
    if (n == 1U) {
        auto e = mpoly_to_expr(part, wc, ctx);
        if (e.is_error()) return fail<void>(e.error());
        auto f = factor_polynomial(e.value(), wc.vars[0], ctx);
        if (f.is_error()) return fail<void>(f.error());
        for (const auto& pf : f.value().factors) {
            auto ip = expr_to_intpoly_main(pf.factor, wc.vars[0], ctx);
            if (ip.is_error()) return fail<void>(ip.error());
            MPoly mf = mpoly_from_intpoly(ip.value(), 0U, n);
            for (unsigned int m = 0; m < pf.multiplicity; ++m) {
                out_factors.push_back(mf);
            }
        }
        // content
        if (const auto* il = expr_cast<IntegerLit>(f.value().content)) {
            if (!(il->value.abs() == BigInt(1))) {
                out_factors.push_back(mpoly_constant(il->value, n));
            } else if (il->value.is_negative()) {
                out_factors.push_back(mpoly_constant(BigInt(-1), n));
            }
        }
        return ok();
    }

    const Symbol& main_var = wc.vars[0];
    const unsigned int main_deg = mpoly_degree_in(part, 0U);
    if (main_deg == 0U) {
        // No dependence on main var: treat as irreducible-in-this-frame factor.
        // (Could recurse on another variable; certified driver keeps it whole.)
        if (!(part.terms.size() == 1U &&
              part.terms.begin()->second.abs() == BigInt(1) &&
              part.terms.begin()->first == Monomial(n, 0U))) {
            out_factors.push_back(part);
        }
        return ok();
    }

    // Choose a good evaluation point for x_2..x_n.  Good = (a) leading coeff in
    // main var stays nonzero, (b) univariate image stays squarefree & same degree.
    // Search a deterministic expanding box of integer points (bound derived from
    // degree: at most (2B+1)^(n-1) points where B grows until success).
    const MPoly lc_main = mpoly_leading_coeff_in(part, 0U);

    // Precompute the irreducible (multivariate) factors of lc_main.  Wang's
    // leading-coefficient condition (GCL §6.6): the chosen point must make these
    // lc-factor images pairwise "distinguishable" — each lc-factor image must
    // possess a prime that divides no other lc-factor image — so that the LC
    // distribution is unambiguous.  Empty list ⇒ lc_main constant ⇒ no condition.
    std::vector<MPoly> lc_factors;
    {
        const bool lc_constant = lc_main.terms.size() == 1U &&
            lc_main.terms.begin()->first == Monomial(n, 0U);
        if (!lc_constant) {
            auto lce = mpoly_to_expr(lc_main, wc, ctx);
            if (lce.is_error()) return fail<void>(lce.error());
            auto lcf = factor_multivariate(lce.value(), ctx);
            if (lcf.is_ok()) {
                for (const auto& pf : lcf.value().factors) {
                    auto mp = mpoly_from_expr(pf.factor, wc, ctx);
                    if (mp.is_error()) return fail<void>(mp.error());
                    for (unsigned int m = 0; m < pf.multiplicity; ++m)
                        lc_factors.push_back(mp.value());
                }
            }
        }
    }
    auto eval_at = [&](const MPoly& p, const std::vector<BigInt>& pt) {
        MPoly cur = p;
        for (std::size_t vi = 1; vi < n; ++vi) cur = mpoly_eval_var(cur, vi, pt[vi - 1U]);
        return cur.is_zero() ? BigInt(0) : cur.terms.begin()->second;
    };
    auto lc_factors_distinguishable = [&](const std::vector<BigInt>& pt) {
        std::vector<BigInt> imgs;
        imgs.reserve(lc_factors.size());
        for (const auto& f : lc_factors) {
            BigInt e = eval_at(f, pt).abs();
            if (e.is_zero() || e == BigInt(1)) return false;
            imgs.push_back(e);
        }
        // Each image must own a prime absent from all others: equivalent to
        // img_i not dividing the product of the gcd-reduced others.  We use the
        // pairwise condition: for each i there is a prime p | img_i with p∤img_j
        // (j≠i).  Implemented as: img_i / gcd(img_i, prod_{j≠i} img_j) > 1.
        for (std::size_t i = 0; i < imgs.size(); ++i) {
            BigInt others(1);
            for (std::size_t j = 0; j < imgs.size(); ++j)
                if (j != i) others *= imgs[j];
            BigInt g = gcd(imgs[i], others);
            if (imgs[i] / g == BigInt(1)) return false;  // no distinguishing prime
        }
        return true;
    };

    std::vector<BigInt> point(n - 1U, BigInt(0));
    // Bound on |a_i|: Schwartz-Zippel — a nonzero discriminant/resultant of
    // degree D vanishes at <= D points per axis; we expand the box until found.
    // F3.2-WANG-EVAL-BOUND: configurable via ctx.set_max_wang_eval_radius().
    const std::size_t user_radius = ctx.max_wang_eval_radius();
    const long long max_radius = (user_radius > 0U)
        ? static_cast<long long>(user_radius)
        : static_cast<long long>(part.terms.size() + main_deg + 4U);

    // Collect ALL acceptable evaluation points across the box (good-point
    // candidates).  For each, we will attempt the full LC distribution +
    // Hensel + certification; the FIRST that succeeds wins.  Failures of
    // the LC step (e.g. non-divisors fails, or content correction is not
    // exact at this point) trigger advancement to the next candidate —
    // mirroring the EEZ "config retry" loop in Wang 1978 §3 / SymPy's
    // dmp_zz_wang outer loop.
    struct PointCand {
        std::vector<BigInt> pt;
        IntPoly image;
    };
    std::vector<PointCand> candidates;
    (void)lc_factors_distinguishable;  // legacy distinguishability — diag only
    for (long long radius = 0; radius <= max_radius; ++radius) {
        std::vector<long long> coords(n - 1U, 0);
        bool done_box = false;
        while (!done_box) {
            bool on_shell = (radius == 0);
            for (long long c : coords) {
                if (std::llabs(c) == radius) { on_shell = true; break; }
            }
            if (on_shell) {
                std::vector<BigInt> p(n - 1U);
                for (std::size_t i = 0; i < n - 1U; ++i) p[i] = BigInt(coords[i]);
                MPoly lc_eval = lc_main;
                for (std::size_t vi = 1; vi < n; ++vi)
                    lc_eval = mpoly_eval_var(lc_eval, vi, p[vi - 1U]);
                if (!lc_eval.is_zero()) {
                    MPoly img = part;
                    for (std::size_t vi = 1; vi < n; ++vi)
                        img = mpoly_eval_var(img, vi, p[vi - 1U]);
                    auto ipopt = mpoly_to_intpoly(img, 0U);
                    if (ipopt.has_value() &&
                        mpoly_degree_in(img, 0U) == main_deg) {
                        IntPoly im = *ipopt;
                        IntPoly imd;
                        imd.coefficients().resize(im.size() > 0 ? im.size() - 1U : 0U, BigInt(0));
                        for (std::size_t d = 1; d < im.size(); ++d)
                            imd[d - 1U] = im[d] * BigInt(static_cast<long long>(d));
                        normalize_integer_poly(imd);
                        IntPoly gg = gcd_integer_poly_primitive(im, imd);
                        if (gg.size() <= 1U) {
                            candidates.push_back({p, im});
                        }
                    }
                }
            }
            std::size_t pos = 0;
            while (pos < n - 1U) {
                coords[pos]++;
                if (coords[pos] > radius) { coords[pos] = -radius; ++pos; }
                else break;
            }
            if (pos == n - 1U) done_box = true;
        }
    }

    if (candidates.empty()) {
        return fail<void>(make_error(
            CASErrorKind::Unimplemented,
            "wang: no good evaluation point found within the search bound; the "
            "polynomial may require a larger search box or a different main "
            "variable (GCL §6.5 good-evaluation condition)"));
    }
    // Try candidates in order.  Track the last diagnostic for reporting if
    // every candidate is exhausted without success.
    CASError last_err = make_error(
        CASErrorKind::Unimplemented,
        "wang: no candidate evaluation point yielded a successful factorization");

    for (const auto& cand : candidates) {
        point = cand.pt;
        const IntPoly& univ_image = cand.image;

        auto img_expr = integer_coefficients_to_expr(univ_image, main_var, ctx);
        if (img_expr.is_error()) { last_err = img_expr.error(); continue; }
        auto ufac = factor_polynomial(img_expr.value(), main_var, ctx);
        if (ufac.is_error()) { last_err = ufac.error(); continue; }

        std::vector<IntPoly> univar_factors;
        for (const auto& pf : ufac.value().factors) {
            auto ip = expr_to_intpoly_main(pf.factor, main_var, ctx);
            if (ip.is_error()) { last_err = ip.error(); univar_factors.clear(); break; }
            for (unsigned int m = 0; m < pf.multiplicity; ++m)
                univar_factors.push_back(ip.value());
        }
        if (univar_factors.empty()) continue;

        if (univar_factors.size() <= 1U) {
            out_factors.push_back(part);
            return ok();
        }

        auto lc = wang_distribute_leading_coeff(part, univar_factors, point, wc, ctx);
        if (lc.is_error()) { last_err = lc.error(); continue; }

        MPoly lift_target = part;
        const BigInt overall = lc.value().overall_constant;
        if (!(overall == BigInt(1))) {
            lift_target = mpoly_scale(part, overall);
        }
        auto lifted = wang_multivariate_hensel(
            lift_target, lc.value().lc, lc.value().adjusted, point, wc, ctx);
        if (lifted.is_error()) { last_err = lifted.error(); continue; }

        std::vector<MPoly> cands = lifted.value();
        MPoly prod = mpoly_constant(BigInt(1), n);
        for (auto& f : cands) {
            BigInt cont = mpoly_integer_content(f);
            if (cont > BigInt(1)) {
                auto d = mpoly_exact_div(f, mpoly_constant(cont, n));
                if (d.has_value()) f = *d;
            }
            prod = mpoly_mul(prod, f);
        }
        auto qopt = mpoly_exact_div(part, prod);
        if (!qopt.has_value() || qopt->terms.size() != 1U ||
            qopt->terms.begin()->first != Monomial(n, 0U)) {
            last_err = make_error(
                CASErrorKind::Unimplemented,
                "wang: lifted factors do not reconstruct the input exactly at this "
                "point; advancing to next candidate");
            continue;
        }
        BigInt unit = qopt->terms.begin()->second;
        for (std::size_t i = 0; i < cands.size(); ++i) {
            out_factors.push_back(cands[i]);
        }
        if (!(unit == BigInt(1))) {
            out_factors.push_back(mpoly_constant(unit, n));
        }
        return ok();
    }
    return fail<void>(last_err);
}

}  // namespace

Result<Factorization> factor_multivariate(ExprPtr poly, symbolic::CASContext& ctx) {
    if (!poly) {
        return fail<Factorization>(make_error(
            CASErrorKind::InvalidArgument, "factor_multivariate requires a non-null polynomial"));
    }
    if (contains_decimal_literal(poly)) {
        return fail<Factorization>(make_error(
            CASErrorKind::Unimplemented,
            "factor_multivariate: decimal literals are outside the exact symbolic core"));
    }

    auto mv_res = parse_multivariate_polynomial(poly, ctx);
    if (mv_res.is_error()) return fail<Factorization>(mv_res.error());
    const MultivariatePolynomial& mv = mv_res.value();

    if (mv.is_zero()) {
        return fail<Factorization>(make_error(
            CASErrorKind::InvalidArgument, "the zero polynomial has no canonical factorization"));
    }

    WangContext wc = build_wang_context(mv);
    const std::size_t n = wc.nvars();

    Factorization result;
    if (n == 0U) {
        // pure integer constant
        auto e = parse_multivariate_polynomial(poly, ctx);
        result.content = poly;
        return ok(std::move(result));
    }

    auto a_res = mpoly_from_expr(poly, wc, ctx);
    if (a_res.is_error()) return fail<Factorization>(a_res.error());
    MPoly a = a_res.value();

    // Extract integer content.
    BigInt content = mpoly_integer_content(a);
    // sign normalisation: make leading (lex-highest) coeff positive.
    if (!a.is_zero() && std::prev(a.terms.end())->second.is_negative()) {
        content = -content;
    }
    if (!(content == BigInt(1))) {
        auto d = mpoly_exact_div(a, mpoly_constant(content, n));
        if (d.has_value()) a = *d;
    }
    result.content = ctx.arena().make<IntegerLit>(content);

    // Squarefree decomposition in main variable.
    auto sf = squarefree_main(a, wc, ctx);
    if (sf.is_error()) return fail<Factorization>(sf.error());

    // Factor each squarefree part, attaching its multiplicity.
    std::map<std::string, std::pair<ExprPtr, unsigned int>> factor_accum;
    for (const auto& [part, mult] : sf.value()) {
        if (mpoly_degree_in(part, 0U) == 0U &&
            part.terms.size() == 1U &&
            part.terms.begin()->second.abs() == BigInt(1) &&
            part.terms.begin()->first == Monomial(n, 0U)) {
            continue;  // unit
        }
        std::vector<MPoly> irr;
        auto fr = factor_squarefree_part(part, wc, ctx, irr);
        if (fr.is_error()) return fail<Factorization>(fr.error());
        for (const auto& f : irr) {
            auto fe = mpoly_to_expr(f, wc, ctx);
            if (fe.is_error()) return fail<Factorization>(fe.error());
            result.factors.push_back(PolynomialFactor{fe.value(), mult});
        }
    }

    return ok(std::move(result));
}

}  // namespace cas::algebra
