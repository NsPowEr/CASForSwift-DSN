#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "polynomial_internal.hpp"
#include "algebra_internal.hpp"
#include <algorithm>
#include "polynomial_gcd_multivariate_helpers.hpp"

namespace cas {
namespace algebra {

namespace {

[[nodiscard]] Result<PolyExpr> poly_derivative_expr(const PolyExpr& poly, symbolic::CASContext& ctx) {
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
    return ok(std::move(result));
}

} // namespace

[[nodiscard]] Result<SquareFreeFactorization> square_free_factorization(ExprPtr poly, const Symbol& var, symbolic::CASContext& ctx) {
    SquareFreeFactorization result;
    result.content = poly_make_integer(ctx.arena(), 1);

    auto parse_res = parse_polynomial(poly, var, ctx);
    if (parse_res.is_error()) return fail<SquareFreeFactorization>(parse_res.error());
    PolyExpr A_expr = parse_res.value();

    if (is_zero_poly(A_expr)) {
        result.content = poly_make_integer(ctx.arena(), 0);
        return ok(std::move(result));
    }

    if (poly_degree(A_expr) == 0) {
        result.content = A_expr.coefficients()[0];
        return ok(std::move(result));
    }

    // Try to convert to IntPoly to use efficient Subresultant GCD
    auto A_int_res = poly_to_integer_poly(A_expr);
    if (A_int_res.is_ok()) {
        IntPoly A = std::move(A_int_res.value());
        BigInt cont = integer_content(A);
        A = primitive_integer_poly(std::move(A));
        result.content = ctx.arena().make<IntegerLit>(cont);

        IntPoly A_prime;
        if (A.size() > 1) {
            A_prime.resize(A.size() - 1);
            for (std::size_t i = 1; i < A.size(); ++i) {
                A_prime[i - 1] = A[i] * BigInt(static_cast<long long>(i));
            }
        }
        normalize_integer_poly(A_prime);

        IntPoly G = gcd_integer_poly_dispatch(A, A_prime, ctx).gcd;
        
        // C1 = A / G, D1 = A' / G - C1'
        // In Z[x], we use exact division for A/G and A'/G
        // Since G = gcd(A, A'), they are guaranteed to divide.
        // We use pseudo_remainder or a dedicated exact division.
        auto exact_div = [](const IntPoly& num, const IntPoly& den) -> IntPoly {
            if (den.is_zero() || num.is_zero()) return num;
            if (den.degree() == 0) {
                if (den.constant_term().is_zero()) return num;  // guard against /0
                IntPoly res = num;
                divide_integer_coefficients_by_scalar(res, den.constant_term());
                return res;
            }
            // For square-free, den is a factor of num, so pseudo_remainder should result in 0
            // but we want the quotient.
            // We can use a simple division here because we know it's exact.
            if (num.degree() < den.degree()) return IntPoly{};  // safety: underflow guard
            IntPoly q;
            IntPoly r = num;
            q.resize(num.degree() - den.degree() + 1);
            BigInt lc_den = den.leading_coeff();
            for (int i = static_cast<int>(num.degree() - den.degree()); i >= 0; --i) {
                q[i] = r[i + den.degree()] / lc_den;
                for (std::size_t j = 0; j < den.size(); ++j) {
                    r[i + j] -= q[i] * den[j];
                }
            }
            normalize_integer_poly(q);
            return q;
        };

        IntPoly C = exact_div(A, G);
        IntPoly D = exact_div(A_prime, G);
        
        IntPoly C_prime;
        if (C.size() > 1) {
            C_prime.resize(C.size() - 1);
            for (std::size_t i = 1; i < C.size(); ++i) {
                C_prime[i - 1] = C[i] * BigInt(static_cast<long long>(i));
            }
        }
        normalize_integer_poly(C_prime);
        
        // D = D - C'
        if (D.size() < C_prime.size()) D.resize(C_prime.size(), BigInt(0));
        for (std::size_t i = 0; i < C_prime.size(); ++i) D[i] -= C_prime[i];
        normalize_integer_poly(D);

        unsigned int i = 1;
        while (!C.is_zero() && C.degree() > 0) {
            IntPoly factor_poly = gcd_integer_poly_dispatch(C, D, ctx).gcd;
            if (factor_poly.degree() > 0) {
                auto factor_expr = integer_coefficients_to_expr(factor_poly, var, ctx);
                if (factor_expr.is_ok()) {
                    result.factors.push_back({factor_expr.value(), i});
                }
            }

            C = exact_div(C, factor_poly);
            D = exact_div(D, factor_poly);
            
            IntPoly next_C_prime;
            if (C.size() > 1) {
                next_C_prime.resize(C.size() - 1);
                for (std::size_t i = 1; i < C.size(); ++i) {
                    next_C_prime[i - 1] = C[i] * BigInt(static_cast<long long>(i));
                }
            }
            normalize_integer_poly(next_C_prime);
            
            // D = D - C'
            if (D.size() < next_C_prime.size()) D.resize(next_C_prime.size(), BigInt(0));
            for (std::size_t j = 0; j < next_C_prime.size(); ++j) D[j] -= next_C_prime[j];
            normalize_integer_poly(D);
            i++;
        }
        return ok(std::move(result));
    }

    // Fallback to PolyExpr if not integer coefficients
    auto A_prime_res = poly_derivative_expr(A_expr, ctx);
    if (A_prime_res.is_error()) return fail<SquareFreeFactorization>(A_prime_res.error());
    PolyExpr A_prime = A_prime_res.value();

    auto A_expr_val = polynomial_to_expr(A_expr, var, ctx);
    auto A_prime_expr_val = polynomial_to_expr(A_prime, var, ctx);
    if (A_expr_val.is_error() || A_prime_expr_val.is_error()) return fail<SquareFreeFactorization>(make_error(CASErrorKind::InternalError, "Conversion error"));

    auto G_expr_res = polynomial_gcd(A_expr_val.value(), A_prime_expr_val.value(), var, ctx);
    if (G_expr_res.is_error()) return fail<SquareFreeFactorization>(G_expr_res.error());
    auto G_poly_res = parse_polynomial(G_expr_res.value(), var, ctx);
    if (G_poly_res.is_error()) return fail<SquareFreeFactorization>(G_poly_res.error());
    PolyExpr G = G_poly_res.value();

    auto C_div = divide_poly_with_remainder(A_expr, G, ctx);
    if (C_div.is_error()) return fail<SquareFreeFactorization>(C_div.error());
    PolyExpr C = C_div.value().quotient;

    auto D_div = divide_poly_with_remainder(A_prime, G, ctx);
    if (D_div.is_error()) return fail<SquareFreeFactorization>(D_div.error());
    PolyExpr D = D_div.value().quotient;

    auto C_prime_res = poly_derivative_expr(C, ctx);
    if (C_prime_res.is_error()) return fail<SquareFreeFactorization>(C_prime_res.error());
    auto D_minus_C_prime = poly_subtract(D, C_prime_res.value(), ctx);
    if (D_minus_C_prime.is_error()) return fail<SquareFreeFactorization>(D_minus_C_prime.error());
    D = D_minus_C_prime.value();

    unsigned int i = 1;
    while (poly_degree(C) > 0) {
        auto C_expr = polynomial_to_expr(C, var, ctx).value();
        auto D_expr = polynomial_to_expr(D, var, ctx).value();
        auto factor_expr_res = polynomial_gcd(C_expr, D_expr, var, ctx);
        if (factor_expr_res.is_error()) return fail<SquareFreeFactorization>(factor_expr_res.error());
        ExprPtr factor_expr = factor_expr_res.value();

        auto factor_poly_res = parse_polynomial(factor_expr, var, ctx);
        if (factor_poly_res.is_error()) return fail<SquareFreeFactorization>(factor_poly_res.error());
        PolyExpr factor_poly = factor_poly_res.value();

        if (poly_degree(factor_poly) > 0) {
            result.factors.push_back({factor_expr, i});
        }

        auto C_next_div = divide_poly_with_remainder(C, factor_poly, ctx);
        if (C_next_div.is_error()) return fail<SquareFreeFactorization>(C_next_div.error());
        C = C_next_div.value().quotient;

        auto D_next_div = divide_poly_with_remainder(D, factor_poly, ctx);
        if (D_next_div.is_error()) return fail<SquareFreeFactorization>(D_next_div.error());
        D = D_next_div.value().quotient;

        auto next_C_prime_res = poly_derivative_expr(C, ctx);
        if (next_C_prime_res.is_error()) return fail<SquareFreeFactorization>(next_C_prime_res.error());
        auto next_D_res = poly_subtract(D, next_C_prime_res.value(), ctx);
        if (next_D_res.is_error()) return fail<SquareFreeFactorization>(next_D_res.error());
        D = next_D_res.value();
        i++;
    }
    
    return ok(std::move(result));
}

[[nodiscard]] Result<std::vector<MultivariateSquareFreeFactor>> square_free_factorize_multivariate(
    const MultivariatePolynomial& poly, symbolic::CASContext& ctx) {
    
    std::vector<MultivariateSquareFreeFactor> factors;
    if (poly.is_zero()) {
        return ok(std::move(factors));
    }
    
    // Extract integer content
    BigInt int_cont = poly.integer_content();
    MultivariatePolynomial P = poly;
    if (int_cont > BigInt(1) || int_cont < BigInt(-1)) {
        factors.push_back({MultivariatePolynomial{{MultivariateTerm{int_cont, {}}}}, 1U});
        // Divide P by int_cont
        std::vector<MultivariateTerm> new_terms;
        new_terms.reserve(P.terms().size());
        for (const auto& term : P.terms()) {
            new_terms.push_back({term.coefficient / int_cont, term.factors});
        }
        P = MultivariatePolynomial(std::move(new_terms));
    } else if (int_cont == BigInt(-1)) {
        factors.push_back({MultivariatePolynomial{{MultivariateTerm{BigInt(-1), {}}}}, 1U});
        std::vector<MultivariateTerm> new_terms;
        new_terms.reserve(P.terms().size());
        for (const auto& term : P.terms()) {
            new_terms.push_back({-term.coefficient, term.factors});
        }
        P = MultivariatePolynomial(std::move(new_terms));
    }
    
    auto vars = P.variables();
    if (vars.empty()) {
        return ok(std::move(factors));
    }
    
    const Symbol& main_var = vars[0];
    std::size_t deg = degree_in_var(P, main_var);
    
    // Compute content with respect to main_var
    MultivariatePolynomial C;
    for (std::size_t i = 0; i <= deg; ++i) {
        MultivariatePolynomial coeff_poly = coefficient_poly_in_var(P, main_var, i);
        if (coeff_poly.is_zero()) continue;
        
        if (C.terms().empty()) {
            C = coeff_poly;
        } else {
            auto gcd_res = gcd_brown(C, coeff_poly, ctx);
            if (gcd_res.is_error()) return fail<std::vector<MultivariateSquareFreeFactor>>(gcd_res.error());
            C = gcd_res.value();
            if (is_unit_polynomial(C)) break;
        }
    }
    
    if (!C.terms().empty() && !is_unit_polynomial(C)) {
        auto content_factors_res = square_free_factorize_multivariate(C, ctx);
        if (content_factors_res.is_error()) return fail<std::vector<MultivariateSquareFreeFactor>>(content_factors_res.error());
        for (const auto& f : content_factors_res.value()) {
            factors.push_back(f);
        }
        
        auto P_div = exact_quotient(P, C, vars, ctx);
        if (P_div.is_error()) return fail<std::vector<MultivariateSquareFreeFactor>>(P_div.error());
        P = P_div.value().value_or(MultivariatePolynomial{});
    }
    
    if (is_unit_polynomial(P)) {
        return ok(std::move(factors));
    }
    
    // Yun's algorithm on primitive part P
    MultivariatePolynomial A = P;
    MultivariatePolynomial A_prime = A.derivative(main_var);
    
    if (A_prime.is_zero()) {
        // Since P is primitive, if A_prime is zero, it shouldn't have variables
        // if it did, it would be caught by content
        return fail<std::vector<MultivariateSquareFreeFactor>>(make_error(CASErrorKind::InternalError, "Primitive part derivative is zero"));
    }
    
    auto G_res = gcd_brown(A, A_prime, ctx);
    if (G_res.is_error()) return fail<std::vector<MultivariateSquareFreeFactor>>(G_res.error());
    MultivariatePolynomial G = G_res.value();
    
    auto C_yun_div = exact_quotient(A, G, vars, ctx);
    if (C_yun_div.is_error()) return fail<std::vector<MultivariateSquareFreeFactor>>(C_yun_div.error());
    MultivariatePolynomial C_yun = C_yun_div.value().value_or(MultivariatePolynomial{});
    
    auto D_yun_div = exact_quotient(A_prime, G, vars, ctx);
    if (D_yun_div.is_error()) return fail<std::vector<MultivariateSquareFreeFactor>>(D_yun_div.error());
    MultivariatePolynomial D_yun = D_yun_div.value().value_or(MultivariatePolynomial{});
    
    MultivariatePolynomial C_prime = C_yun.derivative(main_var);
    MultivariatePolynomial D_minus_C_prime = D_yun + (MultivariatePolynomial{{MultivariateTerm{BigInt(-1), {}}}} * C_prime);
    D_yun = D_minus_C_prime;
    
    unsigned int i = 1;
    while (!C_yun.is_zero() && C_yun.total_degree() > 0 && !is_unit_polynomial(C_yun)) {
        auto factor_res = gcd_brown(C_yun, D_yun, ctx);
        if (factor_res.is_error()) return fail<std::vector<MultivariateSquareFreeFactor>>(factor_res.error());
        MultivariatePolynomial factor_poly = factor_res.value();
        
        if (!is_unit_polynomial(factor_poly)) {
            factors.push_back({factor_poly, i});
        }
        
        auto C_next_div = exact_quotient(C_yun, factor_poly, vars, ctx);
        if (C_next_div.is_error()) return fail<std::vector<MultivariateSquareFreeFactor>>(C_next_div.error());
        C_yun = C_next_div.value().value_or(MultivariatePolynomial{});
        
        auto D_next_div = exact_quotient(D_yun, factor_poly, vars, ctx);
        if (D_next_div.is_error()) return fail<std::vector<MultivariateSquareFreeFactor>>(D_next_div.error());
        D_yun = D_next_div.value().value_or(MultivariatePolynomial{});
        
        MultivariatePolynomial next_C_prime = C_yun.derivative(main_var);
        MultivariatePolynomial next_D = D_yun + (MultivariatePolynomial{{MultivariateTerm{BigInt(-1), {}}}} * next_C_prime);
        D_yun = next_D;
        i++;
    }
    
    return ok(std::move(factors));
}

}  // namespace algebra
}  // namespace cas
