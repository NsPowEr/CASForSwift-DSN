#include "integrate_engine.hpp"

#include <string>
#include <utility>

namespace cas::calculus::integrate_detail {

Result<ExprPtr> Integrator::integrate_function_direct(const std::string& name, ExprPtr argument) {
    BuiltinOp func_id = get_builtin_op(name);
    if (func_id == BuiltinOp::Sin) {
        return ok(make_product(arena_, {make_integer(arena_, -1), make_function(arena_, "cos", {argument})}));
    }
    if (func_id == BuiltinOp::Cos) {
        return ok(make_function(arena_, "sin", {argument}));
    }
    if (func_id == BuiltinOp::Tan) {
        return ok(make_product(arena_, {make_integer(arena_, -1), make_function(arena_, "ln", {make_function(arena_, "abs", {make_function(arena_, "cos", {argument})})})}));
    }
    if (func_id == BuiltinOp::Cot) {
        return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {make_function(arena_, "sin", {argument})})}));
    }
    if (func_id == BuiltinOp::Sec) {
        return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {make_sum(arena_, {make_function(arena_, "sec", {argument}), make_function(arena_, "tan", {argument})})})}));
    }
    if (func_id == BuiltinOp::Csc) {
        return ok(make_product(arena_, {make_integer(arena_, -1), make_function(arena_, "ln", {make_function(arena_, "abs", {make_sum(arena_, {make_function(arena_, "csc", {argument}), make_function(arena_, "cot", {argument})})})})}));
    }
    if (func_id == BuiltinOp::Exp) {
        return ok(make_function(arena_, "exp", {argument}));
    }
    if (func_id == BuiltinOp::Sinh) {
        return ok(make_function(arena_, "cosh", {argument}));
    }
    if (func_id == BuiltinOp::Cosh) {
        return ok(make_function(arena_, "sinh", {argument}));
    }
    if (func_id == BuiltinOp::Ln) {
        ExprPtr x = argument;
        return ok(make_sum(arena_, {make_product(arena_, {x, make_function(arena_, "ln", {x})}), make_unary(arena_, UnaryOp::Neg, x)}));
    }
    if (func_id == BuiltinOp::Atan) {
        ExprPtr x = argument;
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "arctan", {x})}),
            make_product(arena_, {
                make_rational(arena_, -1, 2),
                make_function(arena_, "ln", {make_sum(arena_, {
                    make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2)),
                    make_integer(arena_, 1),
                })}),
            }),
        }));
    }
    if (func_id == BuiltinOp::Sqrt) {
        return integrate_power_direct(argument, make_rational(arena_, 1, 2), Symbol("_u_"));
    }

    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "This elementary function integral is not implemented"));
}

Result<ExprPtr> Integrator::integrate_power_direct(ExprPtr base, ExprPtr exponent, const Symbol& var) {
    if (depends_on(exponent, var)) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Only powers with constant exponent support direct substitution"));
    }

    if (is_negative_one(exponent)) {
        return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {base})}));
    }

    if (const auto* integer = expr_cast<IntegerLit>(exponent)) {
        const BigInt next = integer->value + BigInt(1);
        return ok(make_product(arena_, {arena_.make<RationalLit>(BigInt(1), next), make_binary(arena_, BinaryOp::Pow, base, arena_.make<IntegerLit>(next))}));
    }

    ExprPtr exponent_plus_one = make_sum(arena_, {exponent, make_integer(arena_, 1)});
    return ok(make_binary(arena_, BinaryOp::Div, make_binary(arena_, BinaryOp::Pow, base, exponent_plus_one), exponent_plus_one));
}

Result<ExprPtr> Integrator::integrate_function(const FuncCall& call, const Symbol& var) {
    if (call.args.size() != 1U) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Only unary function integration is implemented"));
    }

    ExprPtr argument = call.args.front();
    BuiltinOp func_id = call.func_id;
    if (func_id == BuiltinOp::Sqrt) {
        return integrate_sqrt_quadratic(argument, var);
    }

    auto affine = extract_affine_argument(argument, var);
    if (!affine.has_value() || affine->coefficient.numerator() == BigInt(0)) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Function integration currently supports only direct or affine exact arguments"));
    }

    auto primitive = integrate_function_direct(call.name, argument);
    if (primitive.is_error()) {
        return primitive;
    }

    if (affine->coefficient == Rational(BigInt(1))) {
        return primitive;
    }

    return ok(make_product(arena_, {make_rational(arena_, Rational(BigInt(1)) / affine->coefficient), primitive.value()}));
}

}  // namespace cas::calculus::integrate_detail
