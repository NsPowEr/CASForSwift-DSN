#include "cas/algebra.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <vector>

namespace cas::algebra {

[[nodiscard]] ExprPtr make_rational_expr(AstArena& arena, const Rational& value) {
    if (value.denominator() == BigInt(1)) {
        return arena.make<IntegerLit>(value.numerator());
    }
    return arena.make<RationalLit>(value.numerator(), value.denominator());
}

[[nodiscard]] static Result<ExprPtr> build_factor_power_expr(
    ExprPtr factor,
    unsigned int power,
    symbolic::CASContext& ctx) {
    if (power == 1U) {
        return ok(factor);
    }
    return pow_expr(factor, static_cast<std::size_t>(power), ctx);
}

// Build polynomial numerator A_0 + A_1*x + ... + A_{d-1}*x^{d-1} from coefficient vector.
[[nodiscard]] static Result<ExprPtr> build_polynomial_numerator(
    const std::vector<Rational>& coeffs,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    std::vector<ExprPtr> sum_parts;
    for (std::size_t k = 0U; k < coeffs.size(); ++k) {
        if (coeffs[k].numerator().is_zero()) continue;
        ExprPtr term_expr = make_rational_expr(ctx.arena(), coeffs[k]);
        if (k > 0U) {
            ExprPtr x_sym = ctx.arena().make<Symbol>(var.name);
            ExprPtr x_pow = (k == 1U)
                ? x_sym
                : static_cast<ExprPtr>(ctx.arena().make<Binary>(
                    BinaryOp::Pow, x_sym,
                    ctx.arena().make<IntegerLit>(BigInt(static_cast<long long>(k)))));
            auto mult = multiply_exprs(term_expr, x_pow, ctx);
            if (mult.is_error()) return mult;
            term_expr = mult.value();
        }
        sum_parts.push_back(term_expr);
    }
    if (sum_parts.empty()) return ok(make_integer(ctx.arena(), 0));
    if (sum_parts.size() == 1U) return ok(sum_parts[0]);
    return simplify_expr(ctx.arena().make<Sum>(std::move(sum_parts)), ctx);
}

Result<std::vector<ExprPtr>> partial_fractions(
    ExprPtr rational_expr,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    if (!rational_expr) {
        return fail<std::vector<ExprPtr>>(make_error(
            CASErrorKind::InvalidArgument,
            "partial_fractions richiede un'espressione razionale non nulla"));
    }

    auto parts = apart_num_den(rational_expr, ctx);
    if (parts.is_error()) {
        return fail<std::vector<ExprPtr>>(parts.error());
    }

    auto numerator_poly = parse_polynomial(parts.value().numerator, var, ctx);
    if (numerator_poly.is_error()) {
        return fail<std::vector<ExprPtr>>(numerator_poly.error());
    }
    auto denominator_poly = parse_polynomial(parts.value().denominator, var, ctx);
    if (denominator_poly.is_error()) {
        return fail<std::vector<ExprPtr>>(denominator_poly.error());
    }

    if (is_zero_poly(denominator_poly.value())) {
        return fail<std::vector<ExprPtr>>(make_error(
            CASErrorKind::Undefined,
            "partial_fractions ha ricevuto un denominatore nullo"));
    }
    if (is_zero_poly(numerator_poly.value())) {
        return ok(std::vector<ExprPtr>{make_integer(ctx.arena(), 0)});
    }
    if (poly_degree(numerator_poly.value()) >= poly_degree(denominator_poly.value())) {
        return fail<std::vector<ExprPtr>>(make_error(
            CASErrorKind::InvalidArgument,
            "partial_fractions richiede una funzione razionale propria"));
    }

    auto denominator_factors = factor_over_integers(parts.value().denominator, var, ctx);
    if (denominator_factors.is_error()) {
        return fail<std::vector<ExprPtr>>(denominator_factors.error());
    }

    auto numerator_rat_res = poly_to_rational_poly(numerator_poly.value());
    if (numerator_rat_res.is_error()) {
        return fail<std::vector<ExprPtr>>(numerator_rat_res.error());
    }
    RatPoly N = numerator_rat_res.value();
    normalize_rational_coefficients(N);

    auto denominator_rat_res = poly_to_rational_poly(denominator_poly.value());
    if (denominator_rat_res.is_error()) {
        return fail<std::vector<ExprPtr>>(denominator_rat_res.error());
    }
    RatPoly D = denominator_rat_res.value();
    normalize_rational_coefficients(D);

    struct FactorData {
        RatPoly d;
        unsigned int k;
        RatPoly V; // d^k
        ExprPtr d_expr;
    };
    std::vector<FactorData> factors;

    RatPoly P_V;
    P_V.push_back(Rational(BigInt(1)));

    for (const PolynomialFactor& pf : denominator_factors.value().factors) {
        auto factor_poly = parse_polynomial(pf.factor, var, ctx);
        if (factor_poly.is_error()) return fail<std::vector<ExprPtr>>(factor_poly.error());
        auto factor_rat = poly_to_rational_poly(factor_poly.value());
        if (factor_rat.is_error()) return fail<std::vector<ExprPtr>>(factor_rat.error());
        
        RatPoly d = factor_rat.value();
        normalize_rational_coefficients(d);
        if (d.size() <= 1) continue; // Constant factor
        
        RatPoly V;
        V.push_back(Rational(BigInt(1)));
        for (unsigned int p = 0; p < pf.multiplicity; ++p) {
            V = mul_rational_poly(V, d);
        }
        
        P_V = mul_rational_poly(P_V, V);
        factors.push_back({d, pf.multiplicity, V, pf.factor});
    }

    if (factors.empty()) {
        return fail<std::vector<ExprPtr>>(make_error(
            CASErrorKind::Unimplemented,
            "partial_fractions richiede almeno un fattore nel denominatore"));
    }

    Rational lc_D = D.leading_coeff();
    Rational lc_PV = P_V.leading_coeff();
    auto C_res = checked_divide(lc_D, lc_PV);
    if (C_res.is_error()) return fail<std::vector<ExprPtr>>(C_res.error());
    Rational C = C_res.value();

    for (auto& coeff : N.coefficients()) {
        auto div_res = checked_divide(coeff, C);
        if (div_res.is_error()) return fail<std::vector<ExprPtr>>(div_res.error());
        coeff = div_res.value();
    }
    normalize_rational_coefficients(N);

    std::vector<ExprPtr> terms;

    for (std::size_t i = 0; i < factors.size(); ++i) {
        RatPoly U = factors[i].V;
        RatPoly W;
        W.push_back(Rational(BigInt(1)));
        for (std::size_t j = i + 1; j < factors.size(); ++j) {
            W = mul_rational_poly(W, factors[j].V);
        }

        RatPoly A;
        if (W.size() <= 1) {
            A = N;
        } else {
            auto [g, s, t] = extended_gcd_rational_poly(U, W);
            Rational g0 = g.empty() ? Rational(BigInt(0)) : g[0];
            if (g0.numerator().is_zero()) {
                return fail<std::vector<ExprPtr>>(make_error(
                    CASErrorKind::InternalError, "GCD is zero in partial fractions"));
            }
            
            for (auto& c : s.coefficients()) {
                auto res = checked_divide(c, g0);
                if (res.is_ok()) c = res.value();
            }
            for (auto& c : t.coefficients()) {
                auto res = checked_divide(c, g0);
                if (res.is_ok()) c = res.value();
            }
            normalize_rational_coefficients(s);
            normalize_rational_coefficients(t);

            auto [q, rem_A] = div_rem_rational_poly(mul_rational_poly(N, t), U);
            A = rem_A;

            RatPoly AW = mul_rational_poly(A, W);
            RatPoly N_minus_AW = sub_rational_poly(N, AW);
            auto [B, rem_B] = div_rem_rational_poly(N_minus_AW, U);
            N = B;
        }

        RatPoly current_A = A;
        for (unsigned int j = factors[i].k; j >= 1; --j) {
            auto [q_A, rem_A] = div_rem_rational_poly(current_A, factors[i].d);
            current_A = q_A;
            if (!rem_A.empty()) {
                auto num_expr = build_polynomial_numerator(rem_A.coefficients(), var, ctx);
                if (num_expr.is_error()) return fail<std::vector<ExprPtr>>(num_expr.error());
                
                auto denom_expr = build_factor_power_expr(factors[i].d_expr, j, ctx);
                if (denom_expr.is_error()) return fail<std::vector<ExprPtr>>(denom_expr.error());
                
                auto term = divide_exprs(num_expr.value(), denom_expr.value(), ctx);
                if (term.is_error()) return fail<std::vector<ExprPtr>>(term.error());
                terms.push_back(term.value());
            }
        }
    }

    if (terms.empty()) {
        terms.push_back(make_integer(ctx.arena(), 0));
    }
    return ok(std::move(terms));
}

} // namespace cas::algebra
