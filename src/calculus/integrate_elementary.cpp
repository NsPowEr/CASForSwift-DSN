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
    if (func_id == BuiltinOp::Tanh) {
        // ∫tanh(x) dx = ln(cosh(x))
        return ok(make_function(arena_, "ln",
            {make_function(arena_, "cosh", {argument})}));
    }
    if (func_id == BuiltinOp::Ln || func_id == BuiltinOp::Log) {
        // Both Ln and Log are natural log in this engine (see
        // differentiate.cpp same fix). ∫ln(x) dx = x·ln(x) - x.
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
    if (func_id == BuiltinOp::Asin) {
        // ∫asin(x) dx = x·asin(x) + sqrt(1 - x²)
        ExprPtr x = argument;
        ExprPtr one_minus_x2 = make_sum(arena_, {
            make_integer(arena_, 1),
            make_unary(arena_, UnaryOp::Neg,
                make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2))),
        });
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "arcsin", {x})}),
            make_function(arena_, "sqrt", {one_minus_x2}),
        }));
    }
    if (func_id == BuiltinOp::Acos) {
        // ∫acos(x) dx = x·acos(x) - sqrt(1 - x²)
        ExprPtr x = argument;
        ExprPtr one_minus_x2 = make_sum(arena_, {
            make_integer(arena_, 1),
            make_unary(arena_, UnaryOp::Neg,
                make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2))),
        });
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "arccos", {x})}),
            make_unary(arena_, UnaryOp::Neg,
                make_function(arena_, "sqrt", {one_minus_x2})),
        }));
    }
    // F7.5.B1: inverse hyperbolic standalone integrals.
    // These functions are not in BuiltinOp (would touch 76 switch
    // statements under -Wswitch -Werror, deferred to Fase 8); they are
    // parsed as FuncCall(name, …) with BuiltinOp::Custom. Match by
    // canonical name (and 'arc' aliases Maxima sometimes emits).
    auto matches_name = [&](std::initializer_list<const char*> names) {
        for (auto* n : names) if (name == n) return true;
        return false;
    };
    if (matches_name({"asinh", "arcsinh"})) {
        // ∫asinh(x) dx = x·asinh(x) - sqrt(x² + 1)
        ExprPtr x = argument;
        ExprPtr x2_plus_1 = make_sum(arena_, {
            make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2)),
            make_integer(arena_, 1),
        });
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "asinh", {x})}),
            make_unary(arena_, UnaryOp::Neg, make_function(arena_, "sqrt", {x2_plus_1})),
        }));
    }
    if (matches_name({"acosh", "arccosh"})) {
        // ∫acosh(x) dx = x·acosh(x) - sqrt(x² - 1)
        ExprPtr x = argument;
        ExprPtr x2_minus_1 = make_sum(arena_, {
            make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2)),
            make_integer(arena_, -1),
        });
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "acosh", {x})}),
            make_unary(arena_, UnaryOp::Neg, make_function(arena_, "sqrt", {x2_minus_1})),
        }));
    }
    if (matches_name({"atanh", "arctanh"})) {
        // ∫atanh(x) dx = x·atanh(x) + ½·ln(1 - x²)
        ExprPtr x = argument;
        ExprPtr one_minus_x2 = make_sum(arena_, {
            make_integer(arena_, 1),
            make_unary(arena_, UnaryOp::Neg,
                make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2))),
        });
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "atanh", {x})}),
            make_product(arena_, {
                make_rational(arena_, 1, 2),
                make_function(arena_, "ln", {one_minus_x2}),
            }),
        }));
    }
    if (matches_name({"acoth", "arccoth"})) {
        // ∫acoth(x) dx = x·acoth(x) + ½·ln(x² - 1)
        ExprPtr x = argument;
        ExprPtr x2_minus_1 = make_sum(arena_, {
            make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2)),
            make_integer(arena_, -1),
        });
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "acoth", {x})}),
            make_product(arena_, {
                make_rational(arena_, 1, 2),
                make_function(arena_, "ln", {x2_minus_1}),
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
        // F7.5.B1: try quadratic argument first; fall back to power_direct
        // for sqrt(affine in x) — `sqrt(x)` = x^(1/2), affine path scales.
        auto quad = integrate_sqrt_quadratic(argument, var);
        if (quad.is_ok()) return quad;
        auto affine = extract_affine_argument(argument, var);
        if (affine.has_value() && affine->coefficient.numerator() != BigInt(0)) {
            // ∫sqrt(a·x + b) dx = (1/a) · (2/3) · (a·x + b)^(3/2)
            auto primitive = integrate_power_direct(argument,
                make_rational(arena_, 1, 2), var);
            if (primitive.is_ok()) {
                if (affine->coefficient == Rational(BigInt(1))) {
                    return primitive;
                }
                return ok(make_product(arena_, {
                    make_rational(arena_, Rational(BigInt(1)) / affine->coefficient),
                    primitive.value(),
                }));
            }
        }
        return quad;  // propagate original failure
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
