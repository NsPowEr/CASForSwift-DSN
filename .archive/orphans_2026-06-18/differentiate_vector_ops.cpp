#include "cas/calculus.hpp"
#include "calculus_internal.hpp"
#include "cas/error.hpp"

#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

[[nodiscard]] ExprPtr make_integer(AstArena& arena, long long value) {
    return arena.make<IntegerLit>(BigInt(value));
}

[[nodiscard]] ExprPtr make_unary(AstArena& arena, UnaryOp op, ExprPtr operand) {
    return arena.make<Unary>(op, operand);
}

[[nodiscard]] ExprPtr make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs) {
    return arena.make<Binary>(op, lhs, rhs);
}

[[nodiscard]] ExprPtr make_matrix(AstArena& arena, std::size_t rows, std::size_t cols, std::vector<ExprPtr> elements) {
    return arena.make<Matrix>(rows, cols, std::move(elements));
}

[[nodiscard]] bool is_exact_zero(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value.is_zero();
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator.is_zero();
    }
    return false;
}

}  // namespace

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
