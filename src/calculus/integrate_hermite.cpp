#include "integrate_engine.hpp"
#include "cas/algebra.hpp"
#include "cas/differential_algebra.hpp"
#include "../algebra/algebra_internal.hpp"
#include "../algebra/polynomial_internal.hpp"

namespace cas::calculus {

using namespace cas::algebra;

namespace {

// Helper to convert RatPoly to ExprPtr
[[nodiscard]] Result<ExprPtr> ratpoly_to_expr(const RatPoly& poly, const Symbol& var, symbolic::CASContext& ctx) {
    if (poly.empty()) return ok(make_integer(ctx.arena(), 0));
    std::vector<ExprPtr> terms;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        if (poly[i].numerator().is_zero()) continue;
        ExprPtr coeff = make_rational_expr(ctx.arena(), poly[i]);
        if (i == 0) {
            terms.push_back(coeff);
        } else {
            ExprPtr x_pow = (i == 1) 
                ? static_cast<ExprPtr>(ctx.arena().make<Symbol>(var.name))
                : static_cast<ExprPtr>(ctx.arena().make<Binary>(BinaryOp::Pow, 
                    ctx.arena().make<Symbol>(var.name), 
                    make_integer(ctx.arena(), static_cast<long long>(i))));
            terms.push_back(multiply_exprs(coeff, x_pow, ctx).value());
        }
    }
    if (terms.empty()) return ok(make_integer(ctx.arena(), 0));
    if (terms.size() == 1) return ok(terms[0]);
    return ok(ctx.arena().make<Sum>(std::move(terms)));
}

// Derivative of RatPoly
[[nodiscard]] RatPoly ratpoly_derivative(const RatPoly& p) {
    if (p.size() <= 1) return {};
    RatPoly res;
    res.reserve(p.size() - 1);
    for (std::size_t i = 1; i < p.size(); ++i) {
        res.push_back(p[i] * Rational(BigInt(static_cast<long long>(i))));
    }
    return res;
}

} // namespace

Result<HermiteReduction> hermite_reduction_exact(
    ExprPtr P_expr,
    ExprPtr Q_expr,
    const Symbol& var,
    symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();

    auto P_poly_res = parse_polynomial(P_expr, var, ctx);
    if (P_poly_res.is_error()) return fail<HermiteReduction>(P_poly_res.error());
    auto Q_poly_res = parse_polynomial(Q_expr, var, ctx);
    if (Q_poly_res.is_error()) return fail<HermiteReduction>(Q_poly_res.error());

    auto P_rat_res = poly_to_rational_poly(P_poly_res.value());
    if (P_rat_res.is_error()) return fail<HermiteReduction>(P_rat_res.error());
    auto Q_rat_res = poly_to_rational_poly(Q_poly_res.value());
    if (Q_rat_res.is_error()) return fail<HermiteReduction>(Q_rat_res.error());

    RatPoly P = P_rat_res.value();
    RatPoly Q = Q_rat_res.value();
    normalize_rational_coefficients(P);
    normalize_rational_coefficients(Q);

    if (Q.empty()) return fail<HermiteReduction>(make_error(CASErrorKind::Undefined, "Division by zero in Hermite reduction"));

    // Square-free factorization of Q
    BigInt common_den(1);
    for (const auto& c : Q.coefficients()) {
        common_den = (common_den * c.denominator()) / gcd(common_den, c.denominator());
    }
    IntPoly Q_int;
    for (const auto& c : Q.coefficients()) {
        Q_int.push_back((c * Rational(common_den)).numerator());
    }
    normalize_integer_poly(Q_int);
    BigInt cont = integer_content(Q_int);
    divide_integer_coefficients_by_scalar(Q_int, cont);

    auto sqf_res = square_free_factorize_integer_poly(Q_int, ctx);
    if (sqf_res.is_error()) return fail<HermiteReduction>(sqf_res.error());
    const auto& sqf = sqf_res.value();

    ExprPtr rational_part = make_integer(arena, 0);

    for (const auto& factor_info : sqf) {
        if (factor_info.multiplicity < 2) continue;

        RatPoly V_base;
        for (const auto& c : factor_info.factor.coefficients()) {
            V_base.push_back(Rational(c));
        }
        normalize_rational_coefficients(V_base);

        for (unsigned int k = factor_info.multiplicity; k >= 2; --k) {
            RatPoly Vk = RatPoly({Rational(BigInt(1))});
            for (unsigned int i = 0; i < k; ++i) Vk = mul_rational_poly(Vk, V_base);
            
            auto div_U_res = div_rem_rational_poly(Q, Vk);
            RatPoly U = div_U_res.first;
            
            RatPoly V_prime = ratpoly_derivative(V_base);
            
            RatPoly UV_prime = mul_rational_poly(U, V_prime);
            RatPoly coeff_poly = RatPoly({Rational(BigInt(-static_cast<long long>(k - 1)))});
            RatPoly target = mul_rational_poly(coeff_poly, UV_prime);
            
            auto [g, s, t] = extended_gcd_rational_poly(target, V_base);
            Rational g0 = g.empty() ? Rational(BigInt(0)) : g[0];
            if (g0.numerator().is_zero()) {
                return fail<HermiteReduction>(make_error(CASErrorKind::InternalError, "GCD is zero in Hermite reduction"));
            }
            for (auto& c : s.coefficients()) {
                auto div_res = checked_divide(c, g0);
                if (div_res.is_ok()) c = div_res.value();
            }
            normalize_rational_coefficients(s);

            auto [q_A, rem_A] = div_rem_rational_poly(mul_rational_poly(P, s), V_base);
            RatPoly A = rem_A;

            auto A_expr_res = ratpoly_to_expr(A, var, ctx);
            if (A_expr_res.is_error()) return fail<HermiteReduction>(A_expr_res.error());
            ExprPtr A_expr = A_expr_res.value();

            RatPoly Vk_minus_1 = RatPoly({Rational(BigInt(1))});
            for (unsigned int i = 0; i < k - 1; ++i) Vk_minus_1 = mul_rational_poly(Vk_minus_1, V_base);
            auto Vk_m1_expr_res = ratpoly_to_expr(Vk_minus_1, var, ctx);
            if (Vk_m1_expr_res.is_error()) return fail<HermiteReduction>(Vk_m1_expr_res.error());
            ExprPtr Vk_m1_expr = Vk_m1_expr_res.value();
            
            auto term = divide_exprs(A_expr, Vk_m1_expr, ctx).value();
            rational_part = add_exprs(rational_part, term, ctx).value();

            RatPoly A_prime = ratpoly_derivative(A);
            RatPoly A_prime_U_V = mul_rational_poly(mul_rational_poly(A_prime, U), V_base);
            RatPoly nmo_A_U_Vp = mul_rational_poly(mul_rational_poly(coeff_poly, mul_rational_poly(A, U)), V_prime);
            
            RatPoly P_new_full = sub_rational_poly(sub_rational_poly(P, A_prime_U_V), nmo_A_U_Vp);
            auto [P_new, rem_P] = div_rem_rational_poly(P_new_full, V_base);
            P = P_new;
            
            auto [Q_new, rem_Q] = div_rem_rational_poly(Q, V_base);
            Q = Q_new;
        }
    }

    auto final_P_expr = ratpoly_to_expr(P, var, ctx).value();
    auto final_Q_expr = ratpoly_to_expr(Q, var, ctx).value();

    return ok(HermiteReduction{
        .rational_part = simplify_expr(rational_part, ctx).value(),
        .remaining_P = final_P_expr,
        .remaining_Q = final_Q_expr
    });
}

} // namespace cas::calculus
