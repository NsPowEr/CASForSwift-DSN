#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

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

[[nodiscard]] ExprPtr mul_expr_lrt(AstArena& arena, ExprPtr a, ExprPtr b) {
    if (!a || is_zero_poly(PolyExpr({a}))) return b;
    if (!b || is_zero_poly(PolyExpr({b}))) return a;
    return arena.make<Binary>(BinaryOp::Mul, a, b);
}

[[nodiscard]] ExprPtr div_expr_lrt(AstArena& arena, ExprPtr a, ExprPtr b) {
    return arena.make<Binary>(BinaryOp::Div, a, b);
}

[[nodiscard]] ExprPtr pow_expr_lrt(AstArena& arena, ExprPtr base, ExprPtr exp) {
    return arena.make<Binary>(BinaryOp::Pow, base, exp);
}

[[nodiscard]] Result<PolyExpr> poly_derivative_lrt(const PolyExpr& poly, symbolic::CASContext& ctx) {
    if (poly.empty()) return ok(PolyExpr{});
    PolyExpr result;
    result.resize(poly.size() - 1, nullptr);
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

[[nodiscard]] Result<PolyExpr> pseudo_remainder_lrt(PolyExpr A, const PolyExpr& B, symbolic::CASContext& ctx) {
    if (is_zero_poly(B))
        return fail<PolyExpr>(make_error_local(CASErrorKind::InvalidArgument, "Divisor cannot be zero"));

    std::size_t m = poly_degree(A);
    std::size_t n = poly_degree(B);
    if (m < n) return ok(std::move(A));

    ExprPtr b_n = leading_coefficient(B);
    PolyExpr R = A;

    for (std::size_t step = 0; step <= m - n; ++step) {
        if (is_zero_poly(R)) break;
        std::size_t deg_r = poly_degree(R);

        // Scale R by b_n
        for (auto& coeff : R.coefficients()) {
            if (coeff) {
                auto s = poly_simplify_expr(mul_expr_lrt(ctx.arena(), coeff, b_n), ctx);
                if (s.is_error()) return fail<PolyExpr>(s.error());
                coeff = s.value();
            }
        }

        if (deg_r == m - step) {
            ExprPtr lc_r = leading_coefficient(R);
            PolyExpr term = poly_make_monomial(lc_r, deg_r - n);

            auto sub_term_res = poly_multiply(B, term, ctx);
            if (sub_term_res.is_error()) return fail<PolyExpr>(sub_term_res.error());

            auto sub_res = poly_subtract(R, sub_term_res.value(), ctx);
            if (sub_res.is_error()) return fail<PolyExpr>(sub_res.error());
            R = sub_res.value();
        }
    }
    normalize_poly(R);
    return ok(std::move(R));
}

struct SubresultantStep {
    PolyExpr P;
    std::size_t degree;
};

[[nodiscard]] Result<std::vector<SubresultantStep>> subresultant_prs_lrt(PolyExpr P1, PolyExpr P2, symbolic::CASContext& ctx) {
    std::vector<SubresultantStep> steps;
    
    ExprPtr g = poly_make_integer(ctx.arena(), 1);
    ExprPtr h = poly_make_integer(ctx.arena(), 1);
    
    steps.push_back({P1, poly_degree(P1)});
    steps.push_back({P2, poly_degree(P2)});

    while (true) {
        std::size_t d = poly_degree(P1) - poly_degree(P2);
        
        auto prem_res = pseudo_remainder_lrt(P1, P2, ctx);
        if (prem_res.is_error()) return fail<std::vector<SubresultantStep>>(prem_res.error());
        PolyExpr P3 = prem_res.value();
        
        if (is_zero_poly(P3)) break;

        ExprPtr divisor;
        {
            auto h_pow = poly_simplify_expr(pow_expr_lrt(ctx.arena(), h, poly_make_integer(ctx.arena(), d)), ctx);
            if (h_pow.is_error()) return fail<std::vector<SubresultantStep>>(h_pow.error());
            auto div_res = poly_simplify_expr(mul_expr_lrt(ctx.arena(), g, h_pow.value()), ctx);
            if (div_res.is_error()) return fail<std::vector<SubresultantStep>>(div_res.error());
            divisor = div_res.value();
        }
        
        auto p3_div = poly_divide_by_scalar(P3, divisor, ctx);
        if (p3_div.is_error()) return fail<std::vector<SubresultantStep>>(p3_div.error());
        P3 = p3_div.value();
        
        P1 = P2;
        P2 = P3;
        g = leading_coefficient(P1);
        
        if (d == 1) {
            h = g;
        } else {
            auto g_pow = poly_simplify_expr(pow_expr_lrt(ctx.arena(), g, poly_make_integer(ctx.arena(), d)), ctx);
            auto h_pow_inv = poly_simplify_expr(pow_expr_lrt(ctx.arena(), h, poly_make_integer(ctx.arena(), d - 1)), ctx);
            if (g_pow.is_error() || h_pow_inv.is_error()) return fail<std::vector<SubresultantStep>>(g_pow.error());
            auto h_next = poly_simplify_expr(div_expr_lrt(ctx.arena(), g_pow.value(), h_pow_inv.value()), ctx);
            if (h_next.is_error()) return fail<std::vector<SubresultantStep>>(h_next.error());
            h = h_next.value();
        }
        
        steps.push_back({P2, poly_degree(P2)});
    }
    
    return ok(std::move(steps));
}

[[nodiscard]] Result<ExprPtr> rioboo_conversion(ExprPtr R_i_z, PolyExpr G_z_x, const Symbol& var, const Symbol& z_var, symbolic::CASContext& ctx) {
    auto r_poly_res = parse_polynomial(R_i_z, z_var, ctx);
    if (r_poly_res.is_error()) return fail<ExprPtr>(r_poly_res.error());
    PolyExpr r_poly = r_poly_res.value();
    
    if (poly_degree(r_poly) == 1) {
        // z - c = 0 => z = c
        auto coeffs = poly_to_rational_poly(r_poly);
        if (coeffs.is_error()) return fail<ExprPtr>(coeffs.error());
        Rational c = -coeffs.value().constant_term() / coeffs.value().leading_coeff();
        
        PolyExpr G_c = G_z_x;
        for (auto& coeff : G_c.coefficients()) {
            if (coeff) {
                auto sub = ctx.substitute(coeff, z_var, make_rational_expr(ctx.arena(), c));
                if (sub.is_error()) return fail<ExprPtr>(sub.error());
                coeff = sub.value();
            }
        }
        auto G_expr = polynomial_to_expr(G_c, var, ctx);
        if (G_expr.is_error()) return fail<ExprPtr>(G_expr.error());
        
        ExprPtr ln_G = ctx.arena().make<FuncCall>("ln", std::vector<ExprPtr>{G_expr.value()});
        return ok(mul_expr_lrt(ctx.arena(), make_rational_expr(ctx.arena(), c), ln_G));
    }
    
    if (poly_degree(r_poly) == 2) {
        auto r_coeffs_res = poly_to_rational_poly(r_poly);
        if (r_coeffs_res.is_error()) return fail<ExprPtr>(r_coeffs_res.error());
        auto r_coeffs = r_coeffs_res.value();
        
        Rational lc = r_coeffs.leading_coeff();
        Rational a = r_coeffs[1] / lc;
        Rational b = r_coeffs[0] / lc;
        
        PolyExpr Q1, Q0;
        Q1.resize(G_z_x.size(), nullptr);
        Q0.resize(G_z_x.size(), nullptr);
        
        for (std::size_t i = 0; i < G_z_x.size(); ++i) {
            if (!G_z_x[i]) continue;
            auto g_i_z_res = parse_polynomial(G_z_x[i], z_var, ctx);
            if (g_i_z_res.is_error()) continue;
            auto g_i_z = poly_to_rational_poly(g_i_z_res.value());
            if (g_i_z.is_error()) continue;
            
            // rem = g_i(z) mod (z^2 + az + b)
            RatPoly monic_r({b, a, Rational(1)});
            auto [quot, rem] = div_rem_rational_poly(g_i_z.value(), monic_r);
            if (rem.size() > 0) Q0[i] = make_rational_expr(ctx.arena(), rem[0]);
            if (rem.size() > 1) Q1[i] = make_rational_expr(ctx.arena(), rem[1]);
        }
        normalize_poly(Q1);
        normalize_poly(Q0);
        
        auto Q1_res = polynomial_to_expr(Q1, var, ctx);
        ExprPtr Q1_x = Q1_res.is_ok() ? Q1_res.value() : make_integer(ctx.arena(), 0);
        auto Q0_res = polynomial_to_expr(Q0, var, ctx);
        ExprPtr Q0_x = Q0_res.is_ok() ? Q0_res.value() : make_integer(ctx.arena(), 0);
        
        ExprPtr a_expr = make_rational_expr(ctx.arena(), a);
        ExprPtr b_expr = make_rational_expr(ctx.arena(), b);
        
        // D = 4b - a^2
        ExprPtr disc_expr = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Sub,
            mul_expr_lrt(ctx.arena(), make_integer(ctx.arena(), 4), b_expr),
            mul_expr_lrt(ctx.arena(), a_expr, a_expr))).value();
        
        ExprPtr sqrt_disc = ctx.arena().make<FuncCall>("sqrt", std::vector<ExprPtr>{disc_expr});
        
        // Norm = Q0^2 - a Q1 Q0 + b Q1^2
        ExprPtr norm = ctx.simplify(ctx.arena().make<Sum>(std::vector<ExprPtr>{
            mul_expr_lrt(ctx.arena(), Q0_x, Q0_x),
            mul_expr_lrt(ctx.arena(), ctx.arena().make<Unary>(UnaryOp::Neg, a_expr), mul_expr_lrt(ctx.arena(), Q1_x, Q0_x)),
            mul_expr_lrt(ctx.arena(), b_expr, mul_expr_lrt(ctx.arena(), Q1_x, Q1_x))
        })).value();
        
        ExprPtr ln_term = mul_expr_lrt(ctx.arena(), make_rational_expr(ctx.arena(), -a/Rational(2)),
            ctx.arena().make<FuncCall>("ln", std::vector<ExprPtr>{norm}));
            
        ExprPtr atan_arg = div_expr_lrt(ctx.arena(),
            mul_expr_lrt(ctx.arena(), sqrt_disc, Q1_x),
            ctx.arena().make<Binary>(BinaryOp::Sub, 
                mul_expr_lrt(ctx.arena(), make_integer(ctx.arena(), 2), Q0_x),
                mul_expr_lrt(ctx.arena(), a_expr, Q1_x)));
        
        ExprPtr atan_term = mul_expr_lrt(ctx.arena(), ctx.arena().make<Unary>(UnaryOp::Neg, sqrt_disc),
            ctx.arena().make<FuncCall>("arctan", std::vector<ExprPtr>{atan_arg}));
            
        return ok(ctx.arena().make<Sum>(std::vector<ExprPtr>{ln_term, atan_term}));
    }
    
    // Fallback: sum over roots
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
        ExprPtr ln_G = ctx.arena().make<FuncCall>("ln", std::vector<ExprPtr>{G_expr.value()});
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
    
    auto Q_prime_res = poly_derivative_lrt(Q, ctx);
    if (Q_prime_res.is_error()) return fail<ExprPtr>(Q_prime_res.error());
    PolyExpr Q_prime = Q_prime_res.value();
    
    Symbol z_var("__z");
    ExprPtr z_sym = ctx.arena().make<Symbol>(z_var.name);
    PolyExpr z_poly_const = poly_make_monomial(z_sym, 0);
    
    auto zQprime_res = poly_multiply(Q_prime, z_poly_const, ctx);
    if (zQprime_res.is_error()) return fail<ExprPtr>(zQprime_res.error());
    
    auto target_poly_res = poly_subtract(P, zQprime_res.value(), ctx);
    if (target_poly_res.is_error()) return fail<ExprPtr>(target_poly_res.error());
    PolyExpr target_poly = target_poly_res.value();
    
    auto prs_res = subresultant_prs_lrt(Q, target_poly, ctx);
    if (prs_res.is_error()) return fail<ExprPtr>(prs_res.error());
    const auto& prs = prs_res.value();
    
    PolyExpr R_poly;
    for (auto it = prs.rbegin(); it != prs.rend(); ++it) {
        if (it->degree == 0) {
            R_poly = it->P;
            break;
        }
    }
    
    if (is_zero_poly(R_poly)) return ok(make_integer(ctx.arena(), 0));
    
    ExprPtr R_z = R_poly[0];
    auto r_fact_res = factor_over_integers(R_z, z_var, ctx);
    if (r_fact_res.is_error()) return fail<ExprPtr>(r_fact_res.error());
    
    std::vector<ExprPtr> integral_terms;
    for (const auto& fact : r_fact_res.value().factors) {
        std::size_t d = 0;
        PolyExpr Gd;
        
        for (const auto& step : prs) {
            bool divides_all = true;
            for (const auto& coeff : step.P.coefficients()) {
                if (!coeff) continue;
                auto c_z_res = parse_polynomial(coeff, z_var, ctx);
                if (c_z_res.is_error()) { divides_all = false; break; }
                auto f_z_res = parse_polynomial(fact.factor, z_var, ctx);
                if (f_z_res.is_error()) { divides_all = false; break; }
                
                auto div = divide_poly_with_remainder(c_z_res.value(), f_z_res.value(), ctx);
                if (div.is_error() || !is_zero_poly(div.value().remainder)) {
                    divides_all = false;
                    break;
                }
            }
            if (!divides_all) {
                d = step.degree;
                Gd = step.P;
                break;
            }
        }
        
        auto conversion = rioboo_conversion(fact.factor, Gd, var, z_var, ctx);
        if (conversion.is_ok()) {
            integral_terms.push_back(conversion.value());
        }
    }
    
    auto final_res = ctx.simplify(ctx.arena().make<Sum>(std::move(integral_terms)));
    return ok(final_res.is_ok() ? final_res.value() : make_integer(ctx.arena(), 0));
}

} // namespace cas::algebra
