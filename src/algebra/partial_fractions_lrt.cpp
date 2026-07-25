#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/ast_debug.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include "partial_fractions_rioboo.hpp"

#include <algorithm>
#include <vector>

namespace cas::algebra {

namespace {

[[nodiscard]] CASError make_error_local(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

// A45: the multiplication helper used to read `if (!a || is_zero(a)) return b;`,
// i.e. it treated a *zero* operand as the multiplicative identity — 0·b returned
// b, which injected a spurious log term into every Rioboo conversion whose R(z)
// has no linear coefficient. A46 moved it to partial_fractions_logtoreal.cpp so
// both Rioboo paths share one definition; `rioboo_mul` is that function.
[[nodiscard]] ExprPtr mul_expr_lrt(AstArena& arena, ExprPtr a, ExprPtr b) {
    return rioboo_mul(arena, a, b);
}

[[nodiscard]] Result<PolyExpr> poly_derivative_lrt(const PolyExpr& poly, symbolic::CASContext& ctx) {
    if (poly.empty()) return ok(PolyExpr{});
    PolyExpr result;
    result.resize(poly.size() - 1, poly_make_integer(ctx.arena(), 0));
    for (std::size_t i = 1; i < poly.size(); ++i) {
        if (!poly[i]) continue;
        ExprPtr i_expr = poly_make_integer(ctx.arena(), static_cast<long long>(i));
        auto term_res = poly_simplify_expr(ctx.arena().make<Binary>(BinaryOp::Mul, i_expr, poly[i]), ctx);
        if (term_res.is_error()) return fail<PolyExpr>(term_res.error());
        result[i - 1] = term_res.value();
    }
    normalize_poly(result);
    return ok(result);
}

// A45: the local subresultant PRS that used to live here (a `g`/`h` reduced-PRS
// recursion) has been removed. It computed the chain correctly but took the last
// degree-0 element as the resultant, which is only valid when deg(R_{k-1}) = 1;
// on a *defective* last step it silently returned a resultant of too low a
// degree (measured: 4 of the 21 cases in scripts/a45_prs_simulation.py, all four
// in the golden corpus). resultant_generic (polynomial_resultant_generic.hpp)
// implements the same Collins/Brown recursion of Bronstein's spec but carries
// the gamma factor across steps, which applies the tau_k correction of Theorem
// 1.5.3 — so integrate_rational_lrt now sources both the resultant and the PRS
// from there, and there is only one such recursion in the engine.

// Mulders' correction, Bronstein Symbolic_Integration_I.md:1829-1830 (the step
// missing from the original publication of Theorem 2.5.1, see footnote at line
// 1836). Requires gcd(lc_x(S_i), Q_i) = 1; when it is not, S_i must be divided
// by gcd(A_j, Q_i)^j for each factor A_j of the squarefree decomposition of
// lc_x(S_i), otherwise S_i(alpha, x) can vanish at a root alpha of Q_i.
[[nodiscard]] Result<PolyExpr> mulders_correction(
    PolyExpr S_i, const RatPoly& Q_i, const Symbol& z_var, symbolic::CASContext& ctx) {
    if (is_zero_poly(S_i) || Q_i.is_zero()) return ok(std::move(S_i));

    auto lc_poly = parse_polynomial(leading_coefficient(S_i), z_var, ctx);
    if (lc_poly.is_error()) return ok(std::move(S_i));
    auto lc_rat = poly_to_rational_poly(lc_poly.value());
    if (lc_rat.is_error()) return ok(std::move(S_i));

    auto [common, s_cofactor, t_cofactor] = extended_gcd_rational_poly(lc_rat.value(), Q_i);
    (void)s_cofactor;
    (void)t_cofactor;
    // Coprime leading coefficient: the spec's requirement already holds and no
    // division is needed. This is the overwhelmingly common case.
    if (common.is_zero() || common.degree() == 0U) return ok(std::move(S_i));

    // Non-coprime: divide S_i by gcd(A_j, Q_i)^j over the squarefree
    // decomposition of lc_x(S_i), coefficient by coefficient in K[z].
    auto lc_expr = polynomial_to_expr(lc_poly.value(), z_var, ctx);
    if (lc_expr.is_error()) return fail<PolyExpr>(lc_expr.error());
    auto lc_sf = square_free_factorization(lc_expr.value(), z_var, ctx);
    if (lc_sf.is_error()) return fail<PolyExpr>(lc_sf.error());

    for (const auto& a_j : lc_sf.value().factors) {
        auto a_poly = parse_polynomial(a_j.factor, z_var, ctx);
        if (a_poly.is_error()) continue;
        auto a_rat = poly_to_rational_poly(a_poly.value());
        if (a_rat.is_error()) continue;
        auto [g_j, sj, tj] = extended_gcd_rational_poly(a_rat.value(), Q_i);
        (void)sj;
        (void)tj;
        if (g_j.is_zero() || g_j.degree() == 0U) continue;

        for (unsigned int power = 0; power < a_j.multiplicity; ++power) {
            PolyExpr divided;
            divided.reserve(S_i.size());
            bool exact = true;
            for (const auto& coeff : S_i.coefficients()) {
                if (!coeff) {
                    divided.push_back(poly_make_integer(ctx.arena(), 0));
                    continue;
                }
                auto c_poly = parse_polynomial(coeff, z_var, ctx);
                if (c_poly.is_error()) { exact = false; break; }
                auto c_rat = poly_to_rational_poly(c_poly.value());
                if (c_rat.is_error()) { exact = false; break; }
                auto [quot, rem] = div_rem_rational_poly(c_rat.value(), g_j);
                if (!rem.is_zero()) { exact = false; break; }
                PolyExpr quot_expr;
                quot_expr.reserve(quot.size());
                for (const auto& qc : quot.coefficients()) quot_expr.push_back(make_rational_expr(ctx.arena(), qc));
                auto as_expr = polynomial_to_expr(quot_expr, z_var, ctx);
                if (as_expr.is_error()) { exact = false; break; }
                divided.push_back(as_expr.value());
            }
            // The spec guarantees exactness; a non-exact division means the
            // input violated a precondition, so stop correcting rather than
            // silently returning a wrong S_i.
            if (!exact) return ok(std::move(S_i));
            normalize_poly(divided);
            S_i = std::move(divided);
        }
    }
    return ok(std::move(S_i));
}

[[nodiscard]] Result<ExprPtr> rioboo_conversion(ExprPtr R_i_z, PolyExpr G_z_x, const Symbol& var, const Symbol& z_var, symbolic::CASContext& ctx) {
    if (!R_i_z) return ok(ExprPtr{});
    auto r_poly_res = parse_polynomial(R_i_z, z_var, ctx);
    if (r_poly_res.is_error()) return fail<ExprPtr>(r_poly_res.error());
    PolyExpr r_poly = r_poly_res.value();
    
    if (poly_degree(r_poly) == 1) {
        auto coeffs = poly_to_rational_poly(r_poly);
        if (coeffs.is_error()) return fail<ExprPtr>(coeffs.error());
        Rational c = -coeffs.value().constant_term() / coeffs.value().leading_coeff();
        
        PolyExpr G_c = G_z_x;
        for (auto& coeff : G_c.coefficients()) {
            if (coeff) {
                auto sub = ctx.substitute(coeff, z_var, make_rational_expr(ctx.arena(), c));
                if (sub.is_error()) return fail<ExprPtr>(sub.error());
                coeff = sub.value();
            } else {
                coeff = poly_make_integer(ctx.arena(), 0);
            }
        }
        auto G_expr = polynomial_to_expr(G_c, var, ctx);
        if (G_expr.is_error()) return fail<ExprPtr>(G_expr.error());
        
        ExprPtr ln_G = ctx.arena().make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{G_expr.value()});
        return ok(mul_expr_lrt(ctx.arena(), make_rational_expr(ctx.arena(), c), ln_G));
    }
    
    if (poly_degree(r_poly) == 2) {
        auto r_coeffs_res = poly_to_rational_poly(r_poly);
        if (r_coeffs_res.is_error()) return fail<ExprPtr>(r_coeffs_res.error());
        auto r_coeffs = r_coeffs_res.value();

        Rational lc = r_coeffs.leading_coeff();
        ExprPtr a_expr = make_rational_expr(ctx.arena(), r_coeffs[1] / lc);
        ExprPtr b_expr = make_rational_expr(ctx.arena(), r_coeffs[0] / lc);
        return rioboo_quadratic_real_form(a_expr, b_expr, G_z_x, var, z_var, ctx);
    }

    // A46: un quartico irriducibile su Q i cui residui vivono in un'estensione
    // quadratica di Q si spezza in due quadratici REALI, e ognuno passa per la
    // stessa forma chiusa di Rioboo del caso razionale (spec §2.8). Fuori da
    // quella classe si cade sulla somma formale `RootSum` qui sotto.
    if (poly_degree(r_poly) == 4) {
        auto quadratics = real_quadratic_factors_of_quartic(r_poly, z_var, ctx);
        if (quadratics.is_ok() && quadratics.value().size() == 2U) {
            std::vector<ExprPtr> terms;
            bool complete = true;
            for (const auto& [a_expr, b_expr] : quadratics.value()) {
                auto part = rioboo_quadratic_real_form(a_expr, b_expr, G_z_x, var, z_var, ctx);
                if (part.is_error() || !part.value()) { complete = false; break; }
                terms.push_back(part.value());
            }
            if (complete && terms.size() == 2U) {
                return ctx.simplify(ctx.arena().make<Sum>(std::move(terms)));
            }
        }
    }


    auto G_expr_res = polynomial_to_expr(G_z_x, var, ctx);
    if (G_expr_res.is_ok()) {
        ExprPtr G_expr = G_expr_res.value();
        ExprPtr ln_G = ctx.arena().make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{G_expr});
        ExprPtr z_sym = ctx.arena().make<Symbol>(z_var.name);
        ExprPtr term = mul_expr_lrt(ctx.arena(), z_sym, ln_G);
        
        return ok(ctx.arena().make<FuncCall>(BuiltinOp::RootSum, std::vector<ExprPtr>{R_i_z, z_sym, term}));
    }

    auto roots_res = solve_polynomial(R_i_z, z_var, ctx);
    if (roots_res.is_error()) return fail<ExprPtr>(roots_res.error());
    std::vector<ExprPtr> terms;
    for (auto root : roots_res.value()) {
        PolyExpr G_root = G_z_x;
        for (auto& coeff : G_root.coefficients()) {
            if (coeff) {
                auto sub = ctx.substitute(coeff, z_var, root);
                if (sub.is_error()) return fail<ExprPtr>(sub.error());
                coeff = sub.value();
            }
        }
        auto G_expr = polynomial_to_expr(G_root, var, ctx);
        if (G_expr.is_error()) return fail<ExprPtr>(G_expr.error());
        ExprPtr ln_G = ctx.arena().make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{G_expr.value()});
        terms.push_back(mul_expr_lrt(ctx.arena(), root, ln_G));
    }
    return ok(ctx.arena().make<Sum>(std::move(terms)));
}

} // namespace

Result<ExprPtr> integrate_rational_lrt(ExprPtr P_expr, ExprPtr Q_expr, const Symbol& var, symbolic::CASContext& ctx) {
    auto P_poly_res = parse_polynomial(P_expr, var, ctx);
    auto Q_poly_res = parse_polynomial(Q_expr, var, ctx);
    if (P_poly_res.is_error() || Q_poly_res.is_error()) {
        return fail<ExprPtr>(make_error_local(CASErrorKind::InvalidArgument, "Failed to parse polynomials"));
    }
    
    PolyExpr P = P_poly_res.value();
    PolyExpr Q = Q_poly_res.value();

    // Short-circuit: 0/Q integrates to 0. Without this guard, the resultant
    // Res(Q, P - z·Q') reduces to ±z·Res(Q, Q')^… ≠ 0, factors give a non-
    // empty root set, and rioboo_conversion emits a spurious log(Q) term.
    // (Hermite reduction frequently leaves remaining_P = 0 after fully
    //  absorbing the integrand into the rational part — e.g. ∫2/(t+1)² dt.)
    if (is_zero_poly(P)) return ok(make_integer(ctx.arena(), 0));

    auto Q_prime_res = poly_derivative_lrt(Q, ctx);
    if (Q_prime_res.is_error()) return fail<ExprPtr>(Q_prime_res.error());
    PolyExpr Q_prime = Q_prime_res.value();
    
    // Divieto hardcode cat. 7: the parameter t of the Rothstein-Trager
    // resultant must not collide with a user symbol, so it is minted fresh
    // rather than being the literal "__z" this code used to hardcode.
    Symbol z_var = ctx.make_fresh_symbol("lrt_t");
    ExprPtr z_sym = ctx.arena().make<Symbol>(z_var.name);
    PolyExpr z_poly_const = poly_make_monomial(z_sym, 0);

    auto zQprime_res = poly_multiply(Q_prime, z_poly_const, ctx);
    if (zQprime_res.is_error()) return fail<ExprPtr>(zQprime_res.error());

    auto target_poly_res = poly_subtract(P, zQprime_res.value(), ctx);
    if (target_poly_res.is_error()) return fail<ExprPtr>(target_poly_res.error());
    PolyExpr target_poly = target_poly_res.value();

    // (R, (R_0, ..., R_k)) <- SubResultant_x(D, A - t·dD/dx)   [spec line 1823]
    //
    // Divieto hardcode cat. 9: the chain honours a configured budget instead of
    // only discovering it has overrun once finished. An explicit hard deadline
    // wins; otherwise the per-operation ctx.timeout() applies. resultant_generic
    // checks this between steps and reports Unimplemented with a diagnostic.
    ResultantDeadline deadline = std::nullopt;
    if (ctx.hard_deadline() != std::chrono::steady_clock::time_point::max()) {
        deadline = ctx.hard_deadline();
    } else if (ctx.timeout().count() > 0) {
        deadline = std::chrono::steady_clock::now() + ctx.timeout();
    }

    std::vector<std::vector<ExprPtr>> chain;
    auto res_res = resultant_generic<ExprPtr>(
        Q.coefficients(), target_poly.coefficients(), &ctx, deadline, &chain);
    if (res_res.is_error()) return fail<ExprPtr>(res_res.error());

    auto R_z_res = poly_simplify_expr(res_res.value(), ctx);
    if (R_z_res.is_error()) return fail<ExprPtr>(R_z_res.error());
    ExprPtr R_z = R_z_res.value();
    if (is_zero_poly(PolyExpr({R_z}))) return ok(make_integer(ctx.arena(), 0));

    // (Q_1, ..., Q_n) <- SquareFree(R)                          [spec line 1824]
    // The multiplicity i is what selects the subresultant below; a full
    // irreducible factorization would throw that information away.
    auto sf_res = square_free_factorization(R_z, z_var, ctx);
    if (sf_res.is_error()) return fail<ExprPtr>(sf_res.error());

    const std::size_t deg_Q = poly_degree(Q);
    std::vector<ExprPtr> integral_terms;
    for (const auto& sf_factor : sf_res.value().factors) {
        // for i such that deg_t(Q_i) > 0                        [spec line 1825]
        auto q_i_res = parse_polynomial(sf_factor.factor, z_var, ctx);
        if (q_i_res.is_error() || poly_degree(q_i_res.value()) == 0U) continue;
        auto q_i_rat = poly_to_rational_poly(q_i_res.value());

        const std::size_t i = static_cast<std::size_t>(sf_factor.multiplicity);
        PolyExpr S_i;
        if (i == deg_Q) {
            S_i = Q;                                          // [spec line 1826]
        } else {
            // S_i <- R_m where deg_x(R_m) = i                 [spec line 1828]
            for (const auto& step : chain) {
                PolyExpr candidate{step};
                normalize_poly(candidate);
                if (!is_zero_poly(candidate) && poly_degree(candidate) == i) {
                    S_i = std::move(candidate);
                    break;
                }
            }
        }
        if (is_zero_poly(S_i)) continue;

        if (q_i_rat.is_ok()) {
            auto corrected = mulders_correction(std::move(S_i), q_i_rat.value(), z_var, ctx);
            if (corrected.is_error()) return fail<ExprPtr>(corrected.error());
            S_i = std::move(corrected.value());
        }

        // Presentation only: the roots of Q_i split over its irreducible
        // factors, so summing per factor is the same formal sum as spec line
        // 1831 — and it keeps the closed forms (log for linear, log+arctan for
        // quadratic) that rioboo_conversion already emits, instead of
        // collapsing everything into a RootSum. Explicitly allowed by the spec
        // discussion at line 2378.
        auto q_i_factored = factor_over_integers(sf_factor.factor, z_var, ctx);
        std::vector<ExprPtr> root_polys;
        if (q_i_factored.is_ok()) {
            for (const auto& f : q_i_factored.value().factors) root_polys.push_back(f.factor);
        }
        if (root_polys.empty()) root_polys.push_back(sf_factor.factor);

        for (const auto& root_poly : root_polys) {
            auto conversion = rioboo_conversion(root_poly, S_i, var, z_var, ctx);
            if (conversion.is_ok() && conversion.value()) {
                integral_terms.push_back(conversion.value());
            }
        }
    }

    if (integral_terms.empty()) return ok(make_integer(ctx.arena(), 0));
    if (integral_terms.size() == 1U) return ok(integral_terms.front());
    
    auto sum_node = ctx.arena().make<Sum>(std::move(integral_terms));
    return ctx.simplify(sum_node);
}

} // namespace cas::algebra
