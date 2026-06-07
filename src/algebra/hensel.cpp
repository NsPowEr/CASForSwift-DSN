#include "cas/algebra/hensel.hpp"
#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "factor_multivariate_internal.hpp"

namespace cas::algebra {

MultivariatePolynomial operator-(const MultivariatePolynomial& a, const MultivariatePolynomial& b) {
    std::vector<MultivariateTerm> terms = a.terms();
    for (const auto& tb : b.terms()) {
        terms.push_back({-tb.coefficient, tb.factors});
    }
    return MultivariatePolynomial(terms);
}

[[nodiscard]] Result<std::pair<MultivariatePolynomial, MultivariatePolynomial>> linear_hensel_step(
    const MultivariatePolynomial& f, 
    const MultivariatePolynomial& g0, 
    const MultivariatePolynomial& h0, 
    const BezoutCoeffs& bezout,
    const Ideal& mod_ideal,
    symbolic::CASContext& ctx
) {
    auto vars = f.variables();
    Symbol main_var("");
    for (const auto& v : vars) {
        bool found = false;
        for (const auto& ev : mod_ideal.point.vars) {
            if (v.name == ev.name) { found = true; break; }
        }
        if (!found) { main_var = v; break; }
    }
    if (main_var.name.empty() && !vars.empty()) {
        main_var = vars[0];
    } else if (main_var.name.empty()) {
        main_var = Symbol("x");
    }

    auto e = operator-(f, g0 * h0);
    
    auto e_expr = multivariate_to_expr(e, ctx);
    if (e_expr.is_error()) return fail<std::pair<MultivariatePolynomial, MultivariatePolynomial>>(e_expr.error());
    
    auto h0_expr = multivariate_to_expr(h0, ctx);
    if (h0_expr.is_error()) return fail<std::pair<MultivariatePolynomial, MultivariatePolynomial>>(h0_expr.error());
    
    auto g0_expr = multivariate_to_expr(g0, ctx);
    if (g0_expr.is_error()) return fail<std::pair<MultivariatePolynomial, MultivariatePolynomial>>(g0_expr.error());
    
    auto s_expr = multivariate_to_expr(bezout.s, ctx);
    if (s_expr.is_error()) return fail<std::pair<MultivariatePolynomial, MultivariatePolynomial>>(s_expr.error());
    
    auto t_expr = multivariate_to_expr(bezout.t, ctx);
    if (t_expr.is_error()) return fail<std::pair<MultivariatePolynomial, MultivariatePolynomial>>(t_expr.error());

    auto e_s = multiply_exprs(e_expr.value(), s_expr.value(), ctx);
    if (e_s.is_error()) return fail<std::pair<MultivariatePolynomial, MultivariatePolynomial>>(e_s.error());
    
    auto e_t = multiply_exprs(e_expr.value(), t_expr.value(), ctx);
    if (e_t.is_error()) return fail<std::pair<MultivariatePolynomial, MultivariatePolynomial>>(e_t.error());

    auto dm_h = polynomial_divmod(e_s.value(), h0_expr.value(), main_var, ctx);
    if (dm_h.is_error()) return fail<std::pair<MultivariatePolynomial, MultivariatePolynomial>>(dm_h.error());
    
    auto dm_g = polynomial_divmod(e_t.value(), g0_expr.value(), main_var, ctx);
    if (dm_g.is_error()) return fail<std::pair<MultivariatePolynomial, MultivariatePolynomial>>(dm_g.error());

    auto ph_poly = parse_multivariate_polynomial(dm_h.value().remainder, ctx);
    if (ph_poly.is_error()) return fail<std::pair<MultivariatePolynomial, MultivariatePolynomial>>(ph_poly.error());
    
    auto pg_poly = parse_multivariate_polynomial(dm_g.value().remainder, ctx);
    if (pg_poly.is_error()) return fail<std::pair<MultivariatePolynomial, MultivariatePolynomial>>(pg_poly.error());

    return ok(std::make_pair(g0 + pg_poly.value(), h0 + ph_poly.value()));
}

[[nodiscard]] Result<std::vector<MultivariatePolynomial>> multivariate_hensel_lift(
    const MultivariatePolynomial& f,
    const std::vector<MultivariatePolynomial>& mod_factors,
    const EvaluationPoint& eval_pt,
    const BigInt& prime,
    symbolic::CASContext& ctx
) {
    (void)prime;

    // 1. Convert everything to MPoly and use wang_multivariate_hensel.
    // We need WangContext.
    WangContext wc;
    auto vars = f.variables();
    // main var is the one not in eval_pt
    Symbol main_var("");
    for (const auto& v : vars) {
        bool found = false;
        for (const auto& ev : eval_pt.vars) {
            if (v.name == ev.name) { found = true; break; }
        }
        if (!found) { main_var = v; break; }
    }
    if (main_var.name.empty() && !vars.empty()) main_var = vars[0];
    
    wc.vars.push_back(main_var);
    for (const auto& ev : eval_pt.vars) {
        wc.vars.push_back(ev);
    }
    
    auto f_expr = multivariate_to_expr(f, ctx);
    if (f_expr.is_error()) return fail<std::vector<MultivariatePolynomial>>(f_expr.error());
    auto f_mpoly = mpoly_from_expr(f_expr.value(), wc, ctx);
    if (f_mpoly.is_error()) return fail<std::vector<MultivariatePolynomial>>(f_mpoly.error());

    // mod_factors -> univariate factors? 
    // They should be converted to IntPoly for univariate_factors
    std::vector<IntPoly> univar_factors;
    for (const auto& mf : mod_factors) {
        auto mf_expr = multivariate_to_expr(mf, ctx);
        if (mf_expr.is_error()) return fail<std::vector<MultivariatePolynomial>>(mf_expr.error());
        auto coeffs = univariate_coefficients(mf_expr.value(), main_var, ctx);
        if (coeffs.is_error()) return fail<std::vector<MultivariatePolynomial>>(coeffs.error());
        
        std::vector<BigInt> int_coeffs;
        for (const auto& c : coeffs.value()) {
            auto ic = expr_to_integer_coefficient(c);
            if (ic.is_error()) return fail<std::vector<MultivariatePolynomial>>(ic.error());
            int_coeffs.push_back(ic.value());
        }
        univar_factors.push_back(IntPoly(int_coeffs));
    }

    // Since we need lifted_lc, we compute it via distribute_lc
    auto ldist = wang_distribute_leading_coeff(f_mpoly.value(), univar_factors, eval_pt.values, wc, ctx);
    if (ldist.is_error()) return fail<std::vector<MultivariatePolynomial>>(ldist.error());

    auto result_mpoly = wang_multivariate_hensel(f_mpoly.value(), ldist.value().lc, ldist.value().adjusted, eval_pt.values, wc, ctx);
    if (result_mpoly.is_error()) return fail<std::vector<MultivariatePolynomial>>(result_mpoly.error());

    std::vector<MultivariatePolynomial> final_res;
    for (const auto& p : result_mpoly.value()) {
        auto e = mpoly_to_expr(p, wc, ctx);
        if (e.is_error()) return fail<std::vector<MultivariatePolynomial>>(e.error());
        auto mp = parse_multivariate_polynomial(e.value(), ctx);
        if (mp.is_error()) return fail<std::vector<MultivariatePolynomial>>(mp.error());
        final_res.push_back(mp.value());
    }

    return ok(std::move(final_res));
}

} // namespace cas::algebra
