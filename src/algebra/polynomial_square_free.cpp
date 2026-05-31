#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "polynomial_internal.hpp"
#include "algebra_internal.hpp"
#include <algorithm>

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

}  // namespace algebra
}  // namespace cas
