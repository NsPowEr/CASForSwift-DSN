// L3-06: factorization of f(x) in Q[x] over a two-level algebraic tower
// Q(alpha_1, alpha_2) via composite Trager shift + iterated absolute
// resultant.
//
// Algorithm:
//
//   1. Build minimal polynomials in two fresh symbols y1, y2:
//        m1(y1)         in Q[y1]
//        m2(y1, y2)     in Q[y1, y2]
//   2. Form the symbolic shift  s = s1*y1 + s2*y2  and the shifted
//      polynomial f_shifted(x, y1, y2) = f(x - s) in Q[x, y1, y2].
//   3. Compute the absolute norm in two steps:
//        N1(x, y1) = Res_{y2} ( m2(y1, y2), f_shifted(x, y1, y2) )
//        N (x)     = Res_{y1} ( m1(y1),     N1(x, y1) )
//      N(x) lies in Q[x].  Denominators are cleared via lcm-scaling
//      before factor_over_integers (no Categoria-4 bail-out on type).
//   4. Iterate (s1, s2) until N(x) is square-free in Q[x] (rigorous
//      gcd(N, N') == 1 check).  Bound: ctx.max_trager_tower_shift_attempts
//      or the discriminant collision default.  Timeout / InternalError
//      errors propagate immediately, never swallowed.
//   5. Factor N(x) = prod N_i(x) over Q.  Each rational factor N_i lifts
//      to a factor f_i(x) of f(x) over Q(alpha_1, alpha_2) via
//        f_i(x) = gcd_{Q(alpha_1, alpha_2)[x]} ( f(x), N_i(x + s1*alpha_1 + s2*alpha_2) )
//      computed with AlgebraicElement<AlgebraicNumber> arithmetic.
//   6. Post-condition: prod_i deg(f_i) == deg(f).  Mismatch -> InternalError.
//      Leading coefficient of f is preserved in out.content so that
//      content * prod(factors) reconstructs the input polynomial exactly.

#include "cas/algebra.hpp"
#include "cas/algebraic_tower_bridge.hpp"
#include "cas/symbolic.hpp"

#include "factorization_tower_internal.hpp"
#include "polynomial_internal.hpp"

#include <cstddef>
#include <string>
#include <utility>

namespace cas {
namespace algebra {

namespace fti = factorization_tower_internal;

Result<Factorization> factor_polynomial_tower(
    ExprPtr poly,
    const Symbol& var,
    const TowerGenerators& gens,
    symbolic::CASContext& ctx) {
    if (!poly) {
        return fail<Factorization>(fti::tower_error(
            CASErrorKind::InvalidArgument,
            "factor_polynomial_tower: null polynomial"));
    }
    if (gens.min_poly_1.size() < 2U || gens.min_poly_2.size() < 2U) {
        return fail<Factorization>(fti::tower_error(
            CASErrorKind::InvalidArgument,
            "factor_polynomial_tower: tower generators must have minimal polynomials of positive degree"));
    }
    if (!fti::tower_min_polys_well_formed(gens)) {
        return fail<Factorization>(fti::tower_error(
            CASErrorKind::InvalidArgument,
            "factor_polynomial_tower: tower minimal polynomials must have non-zero leading coefficient"));
    }

    auto parsed = parse_polynomial(poly, var, ctx);
    if (parsed.is_error()) return fail<Factorization>(parsed.error());
    auto f_rat = poly_to_rational_poly(parsed.value());
    if (f_rat.is_error()) {
        return fail<Factorization>(fti::tower_error(
            CASErrorKind::Unimplemented,
            "factor_polynomial_tower: input polynomial must lie in Q[var]"));
    }
    RatPoly f_poly = std::move(f_rat.value());
    normalize_rational_coefficients(f_poly);
    if (f_poly.is_zero()) {
        return fail<Factorization>(fti::tower_error(
            CASErrorKind::InvalidArgument,
            "factor_polynomial_tower: zero polynomial has no canonical factorization"));
    }
    const std::size_t deg_f = f_poly.degree();
    const std::size_t deg_m1 = gens.min_poly_1.size() - 1U;
    const std::size_t deg_m2 = gens.min_poly_2.size() - 1U;
    if (deg_f == 0U) {
        Factorization trivial;
        trivial.content = fti::rational_lit(ctx.arena(), f_poly.constant_term());
        return ok(std::move(trivial));
    }

    Symbol y1_sym = ctx.make_fresh_symbol("__l306_y1");
    Symbol y2_sym = ctx.make_fresh_symbol("__l306_y2");

    auto m1_expr_res = fti::rational_coeffs_to_expr(gens.min_poly_1, y1_sym, ctx);
    if (m1_expr_res.is_error()) return fail<Factorization>(m1_expr_res.error());
    auto m2_expr_res = fti::lift_outer_min_poly(gens.min_poly_2, y1_sym, y2_sym, ctx);
    if (m2_expr_res.is_error()) return fail<Factorization>(m2_expr_res.error());

    PolyExpr f_poly_expr;
    f_poly_expr.reserve(f_poly.size());
    for (const Rational& r : f_poly.coefficients()) {
        f_poly_expr.push_back(r.numerator().is_zero() ? ExprPtr{} : fti::rational_lit(ctx.arena(), r));
    }
    normalize_poly(f_poly_expr);
    auto f_expr_res = polynomial_to_expr(f_poly_expr, var, ctx);
    if (f_expr_res.is_error()) return fail<Factorization>(f_expr_res.error());

    const std::size_t user_bound = ctx.max_trager_tower_shift_attempts();
    const std::size_t default_bound = fti::compute_default_shift_bound(deg_f, deg_m1, deg_m2);
    const std::size_t max_attempts = (user_bound > 0U) ? user_bound : default_bound;

    std::size_t attempts = 0U;
    for (std::size_t weight = 0U; attempts < max_attempts; ++weight) {
        for (std::size_t s1 = 0U; s1 <= weight && attempts < max_attempts; ++s1) {
            const std::size_t s2 = weight - s1;
            ++attempts;

            ExprPtr shift_minus = fti::build_shift_minus(var, y1_sym, y2_sym, s1, s2, ctx.arena());
            auto shifted_f_raw = ctx.substitute(f_expr_res.value(), var, shift_minus);
            if (shifted_f_raw.is_error()) {
                if (fti::is_fatal_inner_error(shifted_f_raw.error())) return fail<Factorization>(shifted_f_raw.error());
                continue;
            }
            auto shifted_f_exp = expand(shifted_f_raw.value(), ctx);
            if (shifted_f_exp.is_error()) {
                if (fti::is_fatal_inner_error(shifted_f_exp.error())) return fail<Factorization>(shifted_f_exp.error());
                continue;
            }

            auto n1_res = polynomial_resultant(m2_expr_res.value(), shifted_f_exp.value(), y2_sym, ctx);
            if (n1_res.is_error()) {
                if (fti::is_fatal_inner_error(n1_res.error())) return fail<Factorization>(n1_res.error());
                continue;
            }
            auto n1_exp = expand(n1_res.value(), ctx);
            if (n1_exp.is_error()) {
                if (fti::is_fatal_inner_error(n1_exp.error())) return fail<Factorization>(n1_exp.error());
                continue;
            }

            auto norm_res = polynomial_resultant(m1_expr_res.value(), n1_exp.value(), y1_sym, ctx);
            if (norm_res.is_error()) {
                if (fti::is_fatal_inner_error(norm_res.error())) return fail<Factorization>(norm_res.error());
                continue;
            }
            auto norm_exp = expand(norm_res.value(), ctx);
            if (norm_exp.is_error()) {
                if (fti::is_fatal_inner_error(norm_exp.error())) return fail<Factorization>(norm_exp.error());
                continue;
            }

            auto sf = fti::is_square_free_over_q(norm_exp.value(), var, ctx);
            if (sf.is_error()) {
                if (fti::is_fatal_inner_error(sf.error())) return fail<Factorization>(sf.error());
                continue;
            }
            if (!sf.value()) continue;

            auto norm_parsed = parse_polynomial(norm_exp.value(), var, ctx);
            if (norm_parsed.is_error()) {
                if (fti::is_fatal_inner_error(norm_parsed.error())) return fail<Factorization>(norm_parsed.error());
                continue;
            }
            auto norm_rat = poly_to_rational_poly(norm_parsed.value());
            if (norm_rat.is_error()) {
                if (fti::is_fatal_inner_error(norm_rat.error())) return fail<Factorization>(norm_rat.error());
                continue;
            }
            auto integer_norm = fti::clear_denominators_to_integer_poly(norm_rat.value(), var, ctx);
            if (integer_norm.is_error()) return fail<Factorization>(integer_norm.error());

            auto norm_factorization = factor_over_integers(integer_norm.value(), var, ctx);
            if (norm_factorization.is_error()) {
                if (fti::is_fatal_inner_error(norm_factorization.error())) return fail<Factorization>(norm_factorization.error());
                continue;
            }

            auto f_tower_raw = fti::rational_poly_to_tower_coeffs(f_poly, gens);
            auto f_tower = fti::monic_tower(std::move(f_tower_raw));
            if (f_tower.is_error()) return fail<Factorization>(f_tower.error());

            // Preserve leading coefficient of f as factorization content so
            // that content * prod(factors) == f (B-L3-06-CRITICO fix).
            Factorization out;
            const Rational& f_lead = f_poly.leading_coeff();
            out.content = fti::rational_lit(ctx.arena(), f_lead);
            std::size_t cumulative_degree = 0U;

            for (const auto& pf : norm_factorization.value().factors) {
                auto ni_parsed = parse_polynomial(pf.factor, var, ctx);
                if (ni_parsed.is_error()) {
                    if (fti::is_fatal_inner_error(ni_parsed.error())) return fail<Factorization>(ni_parsed.error());
                    continue;
                }
                auto ni_rat = poly_to_rational_poly(ni_parsed.value());
                if (ni_rat.is_error()) return fail<Factorization>(ni_rat.error());

                auto gi_tower_res = fti::shift_rational_factor_in_tower(ni_rat.value(), s1, s2, gens);
                if (gi_tower_res.is_error()) return fail<Factorization>(gi_tower_res.error());
                auto gi_tower = fti::monic_tower(std::move(gi_tower_res.value()));
                if (gi_tower.is_error()) return fail<Factorization>(gi_tower.error());

                auto gcd_res = tower_detail::poly_extended_gcd<AlgebraicTowerTwoLevel>(
                    f_tower.value(), gi_tower.value());
                if (gcd_res.is_error()) return fail<Factorization>(gcd_res.error());
                auto [g_coeffs, s_coeffs, t_coeffs] = std::move(gcd_res.value());
                (void)s_coeffs;
                (void)t_coeffs;
                tower_detail::strip_trailing(g_coeffs);
                if (g_coeffs.size() <= 1U) continue;

                auto monic_g = fti::monic_tower(std::move(g_coeffs));
                if (monic_g.is_error()) return fail<Factorization>(monic_g.error());
                cumulative_degree += monic_g.value().size() - 1U;

                auto factor_expr = fti::tower_coeffs_to_expr(monic_g.value(), var, gens, ctx);
                if (factor_expr.is_error()) return fail<Factorization>(factor_expr.error());
                out.factors.push_back(PolynomialFactor{factor_expr.value(), 1U});
            }

            if (out.factors.empty()) {
                auto f_back = polynomial_to_expr(f_poly_expr, var, ctx);
                if (f_back.is_error()) return fail<Factorization>(f_back.error());
                // Irreducible monic part: divide f by leading coefficient.
                // f / leading = monic representative.
                Rational inv_lead = Rational(BigInt(1)) / f_lead;
                PolyExpr monic_poly;
                monic_poly.reserve(f_poly.size());
                for (const Rational& r : f_poly.coefficients()) {
                    const Rational scaled = r * inv_lead;
                    monic_poly.push_back(scaled.numerator().is_zero()
                        ? ExprPtr{}
                        : fti::rational_lit(ctx.arena(), scaled));
                }
                normalize_poly(monic_poly);
                auto monic_expr = polynomial_to_expr(monic_poly, var, ctx);
                if (monic_expr.is_error()) return fail<Factorization>(monic_expr.error());
                out.factors.push_back(PolynomialFactor{monic_expr.value(), 1U});
                return ok(std::move(out));
            }

            if (cumulative_degree != deg_f) {
                return fail<Factorization>(fti::tower_error(
                    CASErrorKind::InternalError,
                    "factor_polynomial_tower: cumulative lifted degree " +
                    std::to_string(cumulative_degree) +
                    " != deg(f) " + std::to_string(deg_f) +
                    "; Trager invariant violated"));
            }
            return ok(std::move(out));
        }
    }

    return fail<Factorization>(fti::tower_error(
        CASErrorKind::Unimplemented,
        "factor_polynomial_tower: no square-free composite Trager shift found within ctx.max_trager_tower_shift_attempts"));
}

}  // namespace algebra
}  // namespace cas
