#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/symbolic.hpp"

#include "polynomial_internal.hpp"

namespace cas {
namespace algebra {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string msg) {
    return CASError{.kind = kind, .message = std::move(msg), .hint = std::nullopt};
}

}  // namespace

[[nodiscard]] Result<ExprPtr> polynomial_resultant(ExprPtr p, ExprPtr q, const Symbol& var, symbolic::CASContext& ctx) {
    auto res_f = parse_polynomial(p, var, ctx);
    if (res_f.is_error()) return fail<ExprPtr>(res_f.error());
    auto res_g = parse_polynomial(q, var, ctx);
    if (res_g.is_error()) return fail<ExprPtr>(res_g.error());

    auto res = resultant_generic<ExprPtr>(
        res_f.value().coefficients(),
        res_g.value().coefficients(),
        &ctx);
    if (res.is_error()) return fail<ExprPtr>(res.error());

    auto normalized = ctx.simplify(res.value());
    return normalized.is_ok() ? normalized : ok(res.value());
}

[[nodiscard]] Result<ExprPtr> polynomial_discriminant(ExprPtr p, const Symbol& var, symbolic::CASContext& ctx) {
    if (!p) return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Discriminant requires non-null polynomial"));

    auto res_f = parse_polynomial(p, var, ctx);
    if (res_f.is_error()) return fail<ExprPtr>(res_f.error());
    PolyExpr f = std::move(res_f.value());

    if (is_zero_poly(f)) return ok(poly_make_integer(ctx.arena(), 0));
    std::size_t n = poly_degree(f);
    if (n < 2) return ok(poly_make_integer(ctx.arena(), 0));

    ExprPtr an = leading_coefficient(f);

    auto deriv_res = cas::calculus::diff(p, var, 1U, ctx);
    if (deriv_res.is_error()) return fail<ExprPtr>(deriv_res.error());

    auto res_val = polynomial_resultant(p, deriv_res.value(), var, ctx);
    if (res_val.is_error()) return res_val;

    const int sign = ((n % 4U == 2U) || (n % 4U == 3U)) ? -1 : 1;

    auto signed_res = poly_simplify_expr(
        ctx.arena().make<Binary>(BinaryOp::Mul, poly_make_integer(ctx.arena(), sign), res_val.value()),
        ctx);
    if (signed_res.is_error()) return fail<ExprPtr>(signed_res.error());

    auto disc = poly_simplify_expr(
        ctx.arena().make<Binary>(BinaryOp::Div, signed_res.value(), an),
        ctx);
    if (disc.is_error()) return fail<ExprPtr>(disc.error());

    return disc;
}

}  // namespace algebra
}  // namespace cas
