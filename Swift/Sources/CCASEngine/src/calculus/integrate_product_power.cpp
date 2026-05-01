#include "integrate_engine.hpp"

#include <utility>
#include <vector>

namespace cas::calculus::integrate_detail {

Result<ExprPtr> Integrator::integrate_product(const Product& product, const Symbol& var) {
    auto substitution = try_u_substitution_for_product(product, var);
    if (substitution.is_ok()) {
        return substitution;
    }

    std::vector<ExprPtr> constant_factors;
    std::vector<ExprPtr> variable_factors;
    constant_factors.reserve(product.factors.size());
    variable_factors.reserve(product.factors.size());

    for (ExprPtr factor : product.factors) {
        if (depends_on(factor, var)) {
            variable_factors.push_back(factor);
        } else {
            constant_factors.push_back(factor);
        }
    }

    if (variable_factors.empty()) {
        constant_factors.push_back(arena_.make<Symbol>(var));
        return ok(make_product(arena_, std::move(constant_factors)));
    }

    if (constant_factors.size() == 1U && is_negative_one(constant_factors.front()) &&
        variable_factors.size() == 1U && matches_reciprocal_sqrt_one_minus_square(variable_factors.front(), var)) {
        return ok(make_function(arena_, "arccos", {arena_.make<Symbol>(var)}));
    }

    if (variable_factors.size() >= 2U) {
        if (auto rational_integral = integrate_via_partial_fractions(arena_.make<Product>(product.factors), var);
            rational_integral.is_ok()) {
            if (constant_factors.empty()) {
                return rational_integral;
            }
            std::vector<ExprPtr> result_factors = constant_factors;
            result_factors.push_back(rational_integral.value());
            return ok(make_product(arena_, std::move(result_factors)));
        }

        auto ibp_res = integrate_by_parts(arena_.make<Product>(variable_factors), var, context_);
        if (ibp_res.is_ok()) {
            if (constant_factors.empty()) {
                return ibp_res;
            }
            std::vector<ExprPtr> result_factors = constant_factors;
            result_factors.push_back(ibp_res.value());
            return ok(make_product(arena_, std::move(result_factors)));
        }
    }

    if (variable_factors.size() == 1U) {
        auto inner = integrate_once(variable_factors.front(), var);
        if (inner.is_error()) {
            auto ibp_res = integrate_by_parts(make_product(arena_, {make_integer(arena_, 1), variable_factors.front()}), var, context_);
            if (ibp_res.is_ok()) {
                if (constant_factors.empty()) {
                    return ibp_res;
                }
                std::vector<ExprPtr> result_factors = constant_factors;
                result_factors.push_back(ibp_res.value());
                return ok(make_product(arena_, std::move(result_factors)));
            }
            return inner;
        }

        if (constant_factors.empty()) {
            return inner;
        }

        constant_factors.push_back(inner.value());
        return ok(make_product(arena_, std::move(constant_factors)));
    }

    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Integration by parts or substitution is not implemented for this product"));
}

Result<ExprPtr> Integrator::integrate_power(const Binary& power, const Symbol& var) {
    if (is_same_symbol(power.right, var) && !depends_on(power.left, var)) {
        if (expr_is<DecimalLit>(power.left)) {
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Decimal literals are not supported in symbolic integration"));
        }
        if (is_one(power.left)) {
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "The integral of 1^x is undefined in this symbolic form"));
        }
        return ok(make_binary(arena_, BinaryOp::Div, make_binary(arena_, BinaryOp::Pow, power.left, arena_.make<Symbol>(var)), make_function(arena_, "ln", {power.left})));
    }

    if (is_rational_value(power.right, -1, 1)) {
        if (const auto* sqrt_call = expr_cast<FuncCall>(power.left);
            sqrt_call != nullptr && sqrt_call->name == "sqrt" && sqrt_call->args.size() == 1U) {
            if (matches_one_minus_square(sqrt_call->args.front(), var)) {
                return ok(make_function(arena_, "arcsin", {arena_.make<Symbol>(var)}));
            }
            ExprPtr constant_base{};
            if (matches_square_plus_constant_square(sqrt_call->args.front(), var, constant_base)) {
                ExprPtr x = arena_.make<Symbol>(var);
                return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {make_sum(arena_, {x, power.left})})}));
            }
        }

        if (matches_one_plus_square(power.left, var)) {
            return ok(make_function(arena_, "arctan", {arena_.make<Symbol>(var)}));
        }

        if (auto affine = extract_affine_argument(power.left, var);
            affine.has_value() && affine->coefficient.numerator() != BigInt(0)) {
            ExprPtr primitive = make_function(arena_, "ln", {make_function(arena_, "abs", {power.left})});
            if (affine->coefficient == Rational(BigInt(1))) {
                return ok(primitive);
            }
            return ok(make_product(arena_, {make_rational(arena_, Rational(BigInt(1)) / affine->coefficient), primitive}));
        }

        ExprPtr constant_base{};
        if (matches_square_minus_constant_square(power.left, var, constant_base)) {
            ExprPtr x = arena_.make<Symbol>(var);
            return ok(make_product(arena_, {
                make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), make_product(arena_, {make_integer(arena_, 2), constant_base})),
                make_function(arena_, "ln", {make_function(arena_, "abs", {make_binary(arena_, BinaryOp::Div, make_sum(arena_, {x, make_unary(arena_, UnaryOp::Neg, constant_base)}), make_sum(arena_, {x, constant_base}))})}),
            }));
        }

        if (auto rational_integral = integrate_via_partial_fractions(make_binary(arena_, BinaryOp::Pow, power.left, power.right), var);
            rational_integral.is_ok()) {
            return rational_integral;
        }
    }

    if (is_rational_value(power.right, -1, 2)) {
        if (matches_one_minus_square(power.left, var)) {
            return ok(make_function(arena_, "arcsin", {arena_.make<Symbol>(var)}));
        }
        ExprPtr constant_base{};
        if (matches_square_plus_constant_square(power.left, var, constant_base)) {
            ExprPtr x = arena_.make<Symbol>(var);
            return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {make_sum(arena_, {x, make_function(arena_, "sqrt", {power.left})})})}));
        }
    }

    if (const auto* integer = expr_cast<IntegerLit>(power.right)) {
        if (auto affine = extract_affine_argument(power.left, var);
            affine.has_value() && affine->coefficient.numerator() != BigInt(0)) {
            const BigInt next = integer->value + BigInt(1);
            if (!next.is_zero()) {
                return ok(make_product(arena_, {make_rational(arena_, Rational(BigInt(1), next) / affine->coefficient), make_binary(arena_, BinaryOp::Pow, power.left, arena_.make<IntegerLit>(next))}));
            }
        }
    }

    if (const auto* call = expr_cast<FuncCall>(power.left); call != nullptr && expr_is<IntegerLit>(power.right)) {
        const auto& exponent = expr_ref<IntegerLit>(power.right);
        const std::string name = canonical_function_name(call->name);
        if (exponent.value == BigInt(2) && call->args.size() == 1U && is_same_symbol(call->args.front(), var)) {
            if (name == "sec") {
                return ok(make_function(arena_, "tan", {arena_.make<Symbol>(var)}));
            }
            if (name == "csc") {
                return ok(make_unary(arena_, UnaryOp::Neg, make_function(arena_, "cot", {arena_.make<Symbol>(var)})));
            }
        }
    }

    if (!is_same_symbol(power.left, var) || depends_on(power.right, var)) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Only powers of the integration variable with constant exponent are implemented"));
    }

    if (is_negative_one(power.right)) {
        return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {arena_.make<Symbol>(var)})}));
    }

    if (const auto* integer = expr_cast<IntegerLit>(power.right)) {
        const BigInt next = integer->value + BigInt(1);
        return ok(make_product(arena_, {arena_.make<RationalLit>(BigInt(1), next), make_binary(arena_, BinaryOp::Pow, arena_.make<Symbol>(var), arena_.make<IntegerLit>(next))}));
    }

    ExprPtr exponent_plus_one = make_sum(arena_, {power.right, make_integer(arena_, 1)});
    return ok(make_binary(arena_, BinaryOp::Div, make_binary(arena_, BinaryOp::Pow, arena_.make<Symbol>(var), exponent_plus_one), exponent_plus_one));
}

}  // namespace cas::calculus::integrate_detail
