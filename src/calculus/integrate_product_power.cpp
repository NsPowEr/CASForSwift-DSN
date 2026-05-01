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

    for (std::size_t index = 0; index < variable_factors.size(); ++index) {
        const auto* power = expr_cast<Binary>(variable_factors[index]);
        if (power == nullptr || power->op != BinaryOp::Pow || !is_rational_value(power->right, -1, 1)) {
            continue;
        }

        std::vector<ExprPtr> numerator_factors = constant_factors;
        for (std::size_t factor_index = 0; factor_index < variable_factors.size(); ++factor_index) {
            if (factor_index != index) {
                numerator_factors.push_back(variable_factors[factor_index]);
            }
        }
        ExprPtr numerator = numerator_factors.empty()
            ? make_integer(arena_, 1)
            : make_product(arena_, std::move(numerator_factors));
        Binary quotient{BinaryOp::Div, numerator, power->left};
        if (auto quadratic_integral = integrate_linear_over_quadratic(quotient, var);
            quadratic_integral.is_ok()) {
            return quadratic_integral;
        }
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
            sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt && sqrt_call->args.size() == 1U) {
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
        BuiltinOp func_id = call->func_id;
        if (exponent.value == BigInt(2) && call->args.size() == 1U && is_same_symbol(call->args.front(), var)) {
            if (func_id == BuiltinOp::Sec) {
                return ok(make_function(arena_, "tan", {arena_.make<Symbol>(var)}));
            }
            if (func_id == BuiltinOp::Csc) {
                return ok(make_unary(arena_, UnaryOp::Neg, make_function(arena_, "cot", {arena_.make<Symbol>(var)})));
            }
        }
    }

    if (const auto* call = expr_cast<FuncCall>(power.left); call != nullptr && call->func_id == BuiltinOp::Ln && call->args.size() == 1U && is_same_symbol(call->args.front(), var) && expr_is<IntegerLit>(power.right)) {
        const auto& exponent = expr_ref<IntegerLit>(power.right);
        if (exponent.value > BigInt(0)) {
            // int ln(x)^n dx = x*ln(x)^n - n * int ln(x)^(n-1) dx
            ExprPtr x = arena_.make<Symbol>(var);
            ExprPtr term1 = make_product(arena_, {x, arena_.make<Binary>(BinaryOp::Pow, power.left, power.right)});
            
            ExprPtr next_exp = arena_.make<IntegerLit>(exponent.value - BigInt(1));
            ExprPtr next_pow = (exponent.value == BigInt(1)) ? make_integer(arena_, 1) : arena_.make<Binary>(BinaryOp::Pow, power.left, next_exp);
            
            auto next_integral = integrate_once(next_pow, var);
            if (next_integral.is_ok()) {
                ExprPtr term2 = make_product(arena_, {arena_.make<IntegerLit>(exponent.value), next_integral.value()});
                return ok(make_sum(arena_, {term1, make_unary(arena_, UnaryOp::Neg, term2)}));
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
