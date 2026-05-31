#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/error_helpers.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <algorithm>
#include <map>
#include <vector>

namespace cas::algebra {

Result<std::vector<ExprPtr>> univariate_coefficients(
    ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    auto poly = parse_multivariate_polynomial(expr, ctx);
    if (poly.is_error()) return fail<std::vector<ExprPtr>>(poly.error());
    return poly.value().to_univariate_coefficients(var, ctx);
}

Result<PolynomialBezout> polynomial_bezout(
    ExprPtr a, ExprPtr b, const Symbol& var, symbolic::CASContext& ctx) {
    auto pa = parse_polynomial(a, var, ctx);
    if (pa.is_error()) return fail<PolynomialBezout>(pa.error());
    auto pb = parse_polynomial(b, var, ctx);
    if (pb.is_error()) return fail<PolynomialBezout>(pb.error());
    
    // Univariate polynomial extended GCD in Q[x]
    auto xgcd = poly_extended_gcd(pa.value(), pb.value(), ctx);
    if (xgcd.is_error()) return fail<PolynomialBezout>(xgcd.error());
    
    auto g_expr = polynomial_to_expr(xgcd.value().gcd, var, ctx);
    if (g_expr.is_error()) return fail<PolynomialBezout>(g_expr.error());
    auto s_expr = polynomial_to_expr(xgcd.value().s, var, ctx);
    if (s_expr.is_error()) return fail<PolynomialBezout>(s_expr.error());
    auto t_expr = polynomial_to_expr(xgcd.value().t, var, ctx);
    if (t_expr.is_error()) return fail<PolynomialBezout>(t_expr.error());
    
    return ok(PolynomialBezout{g_expr.value(), s_expr.value(), t_expr.value()});
}

Result<ExprPtr> polynomial_exact_divide(
    ExprPtr a, ExprPtr b, const Symbol& var, symbolic::CASContext& ctx) {
    auto pa = parse_polynomial(a, var, ctx);
    if (pa.is_error()) return fail<ExprPtr>(pa.error());
    auto pb = parse_polynomial(b, var, ctx);
    if (pb.is_error()) return fail<ExprPtr>(pb.error());
    auto dm = divide_poly_with_remainder(pa.value(), pb.value(), ctx);
    if (dm.is_error()) return fail<ExprPtr>(dm.error());
    if (!is_zero_poly(dm.value().remainder)) {
        return fail<ExprPtr>(CASError{CASErrorKind::InvalidArgument,
            "polynomial_exact_divide: division is not exact", std::nullopt});
    }
    return polynomial_to_expr(dm.value().quotient, var, ctx);
}

Result<PolynomialDivMod> polynomial_divmod(
    ExprPtr a, ExprPtr b, const Symbol& var, symbolic::CASContext& ctx) {
    auto pa = parse_polynomial(a, var, ctx);
    if (pa.is_error()) return fail<PolynomialDivMod>(pa.error());
    auto pb = parse_polynomial(b, var, ctx);
    if (pb.is_error()) return fail<PolynomialDivMod>(pb.error());
    auto dm = divide_poly_with_remainder(pa.value(), pb.value(), ctx);
    if (dm.is_error()) return fail<PolynomialDivMod>(dm.error());
    
    auto q_expr = polynomial_to_expr(dm.value().quotient, var, ctx);
    if (q_expr.is_error()) return fail<PolynomialDivMod>(q_expr.error());
    auto r_expr = polynomial_to_expr(dm.value().remainder, var, ctx);
    if (r_expr.is_error()) return fail<PolynomialDivMod>(r_expr.error());
    
    return ok(PolynomialDivMod{q_expr.value(), r_expr.value()});
}

Result<std::size_t> polynomial_degree(
    ExprPtr a, const Symbol& var, symbolic::CASContext& ctx) {
    auto pa = parse_polynomial(a, var, ctx);
    if (pa.is_error()) return fail<std::size_t>(pa.error());
    if (pa.value().is_zero()) return ok(static_cast<std::size_t>(0));
    return ok(pa.value().degree());
}

} // namespace cas::algebra
