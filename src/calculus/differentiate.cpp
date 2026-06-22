#include "cas/calculus.hpp"
#include "calculus_internal.hpp"
#include "differentiate_internal.hpp"
#include "cas/error.hpp"

#include <string>
#include <utility>
#include <vector>

namespace cas::calculus {

Result<ExprPtr> Differentiator::differentiate(ExprPtr expr, const Symbol& var, unsigned int order) {
    if (order == 0U) {
        return ok(expr);
    }

    ExprPtr current = expr;
    for (unsigned int iteration = 0; iteration < order; ++iteration) {
        auto step = differentiate_once(current, var);
        if (step.is_error()) {
            return step;
        }

        auto simplified = context_.simplify(step.value());
        if (simplified.is_error()) {
            return simplified;
        }
        current = simplified.value();
    }

    return ok(current);
}

Result<ExprPtr> diff(ExprPtr expr, const Symbol& var, unsigned int order, symbolic::CASContext& ctx) {
    if (ctx.is_caching_enabled()) {
        auto key = symbolic::CASContext::DiffKey{expr, var.name, order};
        if (auto cached = ctx.diff_cache_.get(key)) {
            return ok(*cached);
        }
    }

    auto differentiated = Differentiator(ctx).differentiate(expr, var, order);
    if (differentiated.is_error()) {
        return differentiated;
    }

    auto result = differentiated.value();
    if (ctx.is_caching_enabled()) {
        auto key = symbolic::CASContext::DiffKey{expr, var.name, order};
        ctx.diff_cache_.put(key, result);
    }

    return symbolic::materialize_expr(result, ctx.arena());
}

Result<ExprPtr> partial_diff(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    return diff(expr, var, 1U, ctx);
}

Result<ExprPtr> implicit_diff(
    ExprPtr relation,
    const Symbol& dependent,
    const Symbol& independent,
    symbolic::CASContext& ctx) {
    if (!relation) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::InvalidArgument,
            "Cannot compute implicit differentiation for a null relation"));
    }
    if (dependent.name == independent.name) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::InvalidArgument,
            "Dependent and independent variables must be distinct in implicit differentiation"));
    }
    if (!depends_on(relation, dependent)) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::InvalidArgument,
            "Implicit differentiation requires the dependent variable to appear in the relation"));
    }

    auto fx = partial_diff(relation, independent, ctx);
    if (fx.is_error()) {
        return fx;
    }

    auto fy = partial_diff(relation, dependent, ctx);
    if (fy.is_error()) {
        return fy;
    }

    auto simplified_fy = ctx.simplify(fy.value());
    if (simplified_fy.is_error()) {
        return simplified_fy;
    }
    if (is_exact_zero(simplified_fy.value())) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::Undefined,
            "Implicit differentiation requires a nonzero partial derivative with respect to the dependent variable"));
    }

    auto result = ctx.simplify(make_binary(
        ctx.arena(),
        BinaryOp::Div,
        make_unary(ctx.arena(), UnaryOp::Neg, fx.value()),
        simplified_fy.value()));
    if (result.is_error()) {
        return result;
    }
    return result;
}

Result<std::vector<ExprPtr>> gradient(ExprPtr expr, const std::vector<Symbol>& vars, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<std::vector<ExprPtr>>(make_error(CASErrorKind::InvalidArgument, "Cannot compute the gradient of a null expression"));
    }

    std::vector<ExprPtr> components;
    components.reserve(vars.size());
    for (const Symbol& var : vars) {
        auto derivative = partial_diff(expr, var, ctx);
        if (derivative.is_error()) {
            return fail<std::vector<ExprPtr>>(derivative.error());
        }
        components.push_back(derivative.value());
    }

    return ok(std::move(components));
}

Result<ExprPtr> jacobian(const std::vector<ExprPtr>& funcs, const std::vector<Symbol>& vars, symbolic::CASContext& ctx) {
    std::vector<ExprPtr> elements;
    elements.reserve(funcs.size() * vars.size());

    for (ExprPtr function : funcs) {
        if (!function) {
            return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot compute the Jacobian of a null function"));
        }

        for (const Symbol& var : vars) {
            auto derivative = partial_diff(function, var, ctx);
            if (derivative.is_error()) {
                return derivative;
            }
            elements.push_back(derivative.value());
        }
    }

    return ok(make_matrix(ctx.arena(), funcs.size(), vars.size(), std::move(elements)));
}

Result<ExprPtr> hessian(ExprPtr expr, const std::vector<Symbol>& vars, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot compute the Hessian of a null expression"));
    }

    auto first_order = gradient(expr, vars, ctx);
    if (first_order.is_error()) {
        return fail<ExprPtr>(first_order.error());
    }

    std::vector<ExprPtr> elements;
    elements.reserve(vars.size() * vars.size());
    for (ExprPtr component : first_order.value()) {
        for (const Symbol& var : vars) {
            auto derivative = partial_diff(component, var, ctx);
            if (derivative.is_error()) {
                return derivative;
            }
            elements.push_back(derivative.value());
        }
    }

    return ok(make_matrix(ctx.arena(), vars.size(), vars.size(), std::move(elements)));
}

}  // namespace cas::calculus
