#include "cas/calculus.hpp"
#include "integrate_engine.hpp"

namespace cas::calculus {

Result<ExprPtr> integrate(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    auto primitive = integrate_detail::integrate_indefinite_impl(expr, var, ctx);
    if (primitive.is_error()) {
        return primitive;
    }
    return symbolic::materialize_expr(primitive.value(), ctx.arena());
}

Result<ExprPtr> definite_integral(ExprPtr expr, const Symbol& var, ExprPtr lower, ExprPtr upper, symbolic::CASContext& ctx) {
    auto primitive = integrate(expr, var, ctx);
    if (primitive.is_error()) {
        return primitive;
    }

    auto lower_value = ctx.substitute(primitive.value(), var, lower);
    if (lower_value.is_error()) {
        return lower_value;
    }

    auto upper_value = ctx.substitute(primitive.value(), var, upper);
    if (upper_value.is_error()) {
        return upper_value;
    }

    return ctx.simplify(integrate_detail::make_sum(ctx.arena(), {
        upper_value.value(),
        integrate_detail::make_unary(ctx.arena(), UnaryOp::Neg, lower_value.value()),
    }));
}

}  // namespace cas::calculus
