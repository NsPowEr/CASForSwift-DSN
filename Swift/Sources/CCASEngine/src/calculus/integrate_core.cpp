#include "integrate_engine.hpp"

#include "cas/algebra.hpp"

#include <utility>
#include <vector>

namespace cas::calculus::integrate_detail {

Integrator::Integrator(symbolic::CASContext& context) noexcept : context_(context), arena_(context.arena()) {}

Result<ExprPtr> Integrator::integrate(ExprPtr expr, const Symbol& var) {
    if (expr_is<FuncCall>(expr)) {
        auto direct_unsimplified = integrate_once(expr, var);
        if (direct_unsimplified.is_ok()) {
            return context_.simplify(direct_unsimplified.value());
        }
    }

    if (algebra::partial_fractions(expr, var, context_).is_ok()) {
        auto direct_result = integrate_once(expr, var);
        if (direct_result.is_ok()) {
            return context_.simplify(direct_result.value());
        }
    }

    auto simplified = context_.simplify(expr);
    if (simplified.is_error()) {
        return simplified;
    }

    auto result = integrate_once(simplified.value(), var);
    if (result.is_error()) {
        return result;
    }

    return context_.simplify(result.value());
}

Result<bool> Integrator::expressions_match_after_simplify(ExprPtr lhs, ExprPtr rhs) {
    auto simplified_lhs = context_.simplify(lhs);
    if (simplified_lhs.is_error()) {
        return fail<bool>(simplified_lhs.error());
    }
    auto simplified_rhs = context_.simplify(rhs);
    if (simplified_rhs.is_error()) {
        return fail<bool>(simplified_rhs.error());
    }
    return ok(structural_equal(simplified_lhs.value(), simplified_rhs.value()));
}

Result<ExprPtr> Integrator::integrate_once(ExprPtr expr, const Symbol& var) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot integrate a null expression"));
    }

    if (expr_is<IntegerLit>(expr) || expr_is<RationalLit>(expr) || expr_is<Constant>(expr)) {
        return ok(make_product(arena_, {expr, arena_.make<Symbol>(var)}));
    }
    if (expr_is<DecimalLit>(expr)) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Decimal literals are not supported in symbolic integration"));
    }
    if (const auto* symbol = expr_cast<Symbol>(expr)) {
        if (symbol->name == var.name) {
            return ok(make_product(arena_, {
                make_rational(arena_, 1, 2),
                make_binary(arena_, BinaryOp::Pow, expr, make_integer(arena_, 2)),
            }));
        }
        return ok(make_product(arena_, {expr, arena_.make<Symbol>(var)}));
    }
    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op == UnaryOp::Neg) {
            if (matches_reciprocal_sqrt_one_minus_square(unary->operand, var)) {
                return ok(make_function(arena_, "arccos", {arena_.make<Symbol>(var)}));
            }
            auto inner = integrate_once(unary->operand, var);
            if (inner.is_error()) {
                return inner;
            }
            return ok(make_unary(arena_, UnaryOp::Neg, inner.value()));
        }
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Factorial integration is not implemented"));
    }
    if (const auto* binary = expr_cast<Binary>(expr)) {
        return integrate_binary(*binary, var);
    }
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        return integrate_function(*call, var);
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> terms;
        terms.reserve(sum->terms.size());
        for (ExprPtr term : sum->terms) {
            auto integrated = integrate_once(term, var);
            if (integrated.is_error()) {
                return integrated;
            }
            terms.push_back(integrated.value());
        }
        return ok(make_sum(arena_, std::move(terms)));
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        return integrate_product(*product, var);
    }

    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Symbolic integration is not implemented for this expression kind"));
}

Result<ExprPtr> Integrator::integrate_binary(const Binary& binary, const Symbol& var) {
    switch (binary.op) {
    case BinaryOp::Add: {
        auto lhs = integrate_once(binary.left, var);
        if (lhs.is_error()) {
            return lhs;
        }
        auto rhs = integrate_once(binary.right, var);
        if (rhs.is_error()) {
            return rhs;
        }
        return ok(make_sum(arena_, {lhs.value(), rhs.value()}));
    }
    case BinaryOp::Sub: {
        auto lhs = integrate_once(binary.left, var);
        if (lhs.is_error()) {
            return lhs;
        }
        auto rhs = integrate_once(binary.right, var);
        if (rhs.is_error()) {
            return rhs;
        }
        return ok(make_sum(arena_, {lhs.value(), make_unary(arena_, UnaryOp::Neg, rhs.value())}));
    }
    case BinaryOp::Mul:
        return integrate_product(Product({binary.left, binary.right}), var);
    case BinaryOp::Div:
        if (is_one(binary.left) && is_same_symbol(binary.right, var)) {
            return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {arena_.make<Symbol>(var)})}));
        }
        if (auto numerator = exact_scalar_from_expr(binary.left);
            numerator.has_value() && numerator->numerator() != BigInt(0)) {
            if (auto affine = extract_affine_argument(binary.right, var);
                affine.has_value() && affine->coefficient.numerator() != BigInt(0)) {
                return ok(make_product(arena_, {
                    make_rational(arena_, (*numerator) / affine->coefficient),
                    make_function(arena_, "ln", {make_function(arena_, "abs", {binary.right})}),
                }));
            }
        }
        if (is_one(binary.left) && matches_one_plus_square(binary.right, var)) {
            return ok(make_function(arena_, "arctan", {arena_.make<Symbol>(var)}));
        }
        {
            ExprPtr arctan_base{};
            if (is_one(binary.left) && matches_square_plus_constant_square(binary.right, var, arctan_base)) {
                ExprPtr x = arena_.make<Symbol>(var);
                return ok(make_binary(arena_, BinaryOp::Div,
                    make_function(arena_, "arctan", {make_binary(arena_, BinaryOp::Div, x, arctan_base)}),
                    arctan_base));
            }
        }
        if (is_one(binary.left)) {
            const auto* sqrt_call = expr_cast<FuncCall>(binary.right);
            if (sqrt_call != nullptr && sqrt_call->name == "sqrt" && sqrt_call->args.size() == 1U &&
                matches_one_minus_square(sqrt_call->args.front(), var)) {
                return ok(make_function(arena_, "arcsin", {arena_.make<Symbol>(var)}));
            }
            ExprPtr arcsin_base{};
            if (sqrt_call != nullptr && sqrt_call->name == "sqrt" && sqrt_call->args.size() == 1U &&
                matches_constant_square_minus_variable_square(sqrt_call->args.front(), var, arcsin_base)) {
                ExprPtr x = arena_.make<Symbol>(var);
                return ok(make_function(arena_, "arcsin", {make_binary(arena_, BinaryOp::Div, x, arcsin_base)}));
            }
            ExprPtr constant_base{};
            if (sqrt_call != nullptr && sqrt_call->name == "sqrt" && sqrt_call->args.size() == 1U &&
                matches_square_plus_constant_square(sqrt_call->args.front(), var, constant_base)) {
                ExprPtr x = arena_.make<Symbol>(var);
                return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {make_sum(arena_, {x, binary.right})})}));
            }
            if (matches_square_minus_constant_square(binary.right, var, constant_base)) {
                ExprPtr x = arena_.make<Symbol>(var);
                return ok(make_product(arena_, {
                    make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), make_product(arena_, {make_integer(arena_, 2), constant_base})),
                    make_function(arena_, "ln", {make_function(arena_, "abs", {make_binary(arena_, BinaryOp::Div, make_sum(arena_, {x, make_unary(arena_, UnaryOp::Neg, constant_base)}), make_sum(arena_, {x, constant_base}))})}),
                }));
            }
        }
        if (is_negative_one(binary.left)) {
            const auto* sqrt_call = expr_cast<FuncCall>(binary.right);
            if (sqrt_call != nullptr && sqrt_call->name == "sqrt" && sqrt_call->args.size() == 1U &&
                matches_one_minus_square(sqrt_call->args.front(), var)) {
                return ok(make_function(arena_, "arccos", {arena_.make<Symbol>(var)}));
            }
        }
        // x/sqrt(Q): closed-form via u = Q, du = Q' dx, result = sqrt(Q) or -sqrt(Q)
        if (is_same_symbol(binary.left, var)) {
            const auto* sqrt_call = expr_cast<FuncCall>(binary.right);
            if (sqrt_call != nullptr && sqrt_call->name == "sqrt" && sqrt_call->args.size() == 1U) {
                ExprPtr radicand = sqrt_call->args.front();
                ExprPtr cbase{};
                // x/sqrt(1-x^2) = -sqrt(1-x^2)
                if (matches_one_minus_square(radicand, var)) {
                    return ok(make_unary(arena_, UnaryOp::Neg, binary.right));
                }
                // x/sqrt(a^2-x^2) = -sqrt(a^2-x^2)
                if (matches_constant_square_minus_variable_square(radicand, var, cbase)) {
                    return ok(make_unary(arena_, UnaryOp::Neg, binary.right));
                }
                // x/sqrt(x^2+1) = sqrt(x^2+1)
                if (matches_one_plus_square(radicand, var)) {
                    return ok(binary.right);
                }
                // x/sqrt(x^2+a^2) = sqrt(x^2+a^2)
                if (matches_square_plus_constant_square(radicand, var, cbase)) {
                    return ok(binary.right);
                }
                ExprPtr cbase2{};
                // x/sqrt(x^2-a^2) = sqrt(x^2-a^2)
                if (matches_square_minus_constant_square(radicand, var, cbase2)) {
                    return ok(binary.right);
                }
            }
        }
        if (!depends_on(binary.right, var)) {
            auto res = integrate_once(binary.left, var);
            if (res.is_ok()) {
                return ok(make_binary(arena_, BinaryOp::Div, res.value(), binary.right));
            }
        }
        if (auto quadratic_integral = integrate_linear_over_quadratic(binary, var);
            quadratic_integral.is_ok()) {
            return quadratic_integral;
        }
        if (auto rational_integral = integrate_via_partial_fractions(make_binary(arena_, BinaryOp::Div, binary.left, binary.right), var);
            rational_integral.is_ok()) {
            return rational_integral;
        }
        return integrate_once(make_product(arena_, {binary.left, make_binary(arena_, BinaryOp::Pow, binary.right, make_integer(arena_, -1))}), var);
    case BinaryOp::Pow:
        return integrate_power(binary, var);
    case BinaryOp::Mod:
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Modulo integration is not implemented"));
    }

    return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unknown binary operator"));
}

Result<ExprPtr> integrate_indefinite_impl(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    return Integrator(ctx).integrate(expr, var);
}

}  // namespace cas::calculus::integrate_detail
