#include "differentiate_internal.hpp"
#include "calculus_internal.hpp"

#include "cas/meijerg.hpp"

namespace cas::calculus {

Result<ExprPtr> Differentiator::differentiate_unary(const Unary& unary, const Symbol& var) {
    if (unary.op == UnaryOp::Neg) {
        auto operand = differentiate_once(unary.operand, var);
        if (operand.is_error()) {
            return operand;
        }
        return ok(make_unary(arena_, UnaryOp::Neg, operand.value()));
    }

    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Factorial differentiation is not implemented"));
}

Result<ExprPtr> Differentiator::differentiate_binary(const Binary& binary, const Symbol& var) {
    switch (binary.op) {
    case BinaryOp::Add: {
        auto lhs = differentiate_once(binary.left, var);
        if (lhs.is_error()) {
            return lhs;
        }
        auto rhs = differentiate_once(binary.right, var);
        if (rhs.is_error()) {
            return rhs;
        }
        return ok(make_sum(arena_, {lhs.value(), rhs.value()}));
    }
    case BinaryOp::Sub: {
        auto lhs = differentiate_once(binary.left, var);
        if (lhs.is_error()) {
            return lhs;
        }
        auto rhs = differentiate_once(binary.right, var);
        if (rhs.is_error()) {
            return rhs;
        }
        return ok(make_sum(arena_, {lhs.value(), make_unary(arena_, UnaryOp::Neg, rhs.value())}));
    }
    case BinaryOp::Mul: {
        auto lhs = differentiate_once(binary.left, var);
        if (lhs.is_error()) {
            return lhs;
        }
        auto rhs = differentiate_once(binary.right, var);
        if (rhs.is_error()) {
            return rhs;
        }
        return ok(make_sum(arena_, {
            make_product(arena_, {lhs.value(), binary.right}),
            make_product(arena_, {binary.left, rhs.value()}),
        }));
    }
    case BinaryOp::Div: {
        auto lhs = differentiate_once(binary.left, var);
        if (lhs.is_error()) {
            return lhs;
        }
        auto rhs = differentiate_once(binary.right, var);
        if (rhs.is_error()) {
            return rhs;
        }
        return ok(make_binary(
            arena_,
            BinaryOp::Div,
            make_sum(arena_, {
                make_product(arena_, {lhs.value(), binary.right}),
                make_unary(arena_, UnaryOp::Neg, make_product(arena_, {binary.left, rhs.value()})),
            }),
            make_power(arena_, binary.right, make_integer(arena_, 2))));
    }
    case BinaryOp::Pow:
        return differentiate_power(binary, var);
    case BinaryOp::Mod:
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Modulo differentiation is not implemented"));
    case BinaryOp::Equal:
    case BinaryOp::Less:
    case BinaryOp::Greater:
    case BinaryOp::LessEqual:
    case BinaryOp::GreaterEqual:
        return ok(arena_.make<IntegerLit>(BigInt(0)));
    }

    return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unknown binary operator"));
}

Result<ExprPtr> Differentiator::differentiate_sum(const Sum& sum, const Symbol& var) {
    std::vector<ExprPtr> derivatives;
    derivatives.reserve(sum.terms.size());
    for (ExprPtr term : sum.terms) {
        auto derivative = differentiate_once(term, var);
        if (derivative.is_error()) {
            return derivative;
        }
        derivatives.push_back(derivative.value());
    }
    return ok(make_sum(arena_, std::move(derivatives)));
}

Result<ExprPtr> Differentiator::differentiate_product(const Product& product, const Symbol& var) {
    if (product.factors.empty()) {
        return ok(make_integer(arena_, 0));
    }

    std::vector<ExprPtr> sum_terms;
    sum_terms.reserve(product.factors.size());
    for (std::size_t index = 0; index < product.factors.size(); ++index) {
        auto derivative = differentiate_once(product.factors[index], var);
        if (derivative.is_error()) {
            return derivative;
        }

        std::vector<ExprPtr> factors = product.factors;
        factors[index] = derivative.value();
        sum_terms.push_back(make_product(arena_, std::move(factors)));
    }

    return ok(make_sum(arena_, std::move(sum_terms)));
}

Result<ExprPtr> Differentiator::differentiate_power(const Binary& power, const Symbol& var) {
    auto base_derivative = differentiate_once(power.left, var);
    if (base_derivative.is_error()) {
        return base_derivative;
    }
    auto exponent_derivative = differentiate_once(power.right, var);
    if (exponent_derivative.is_error()) {
        return exponent_derivative;
    }

    if (!depends_on(power.right, var)) {
        return ok(make_product(arena_, {
            power.right,
            make_power(arena_, power.left, make_sum(arena_, {power.right, make_integer(arena_, -1)})),
            base_derivative.value(),
        }));
    }

    return ok(make_product(arena_, {
        make_power(arena_, power.left, power.right),
        make_sum(arena_, {
            make_product(arena_, {
                exponent_derivative.value(),
                make_function(arena_, "ln", {power.left}),
            }),
            make_binary(
                arena_,
                BinaryOp::Div,
                make_product(arena_, {power.right, base_derivative.value()}),
                power.left),
        }),
    }));
}

Result<ExprPtr> Differentiator::differentiate_function(const FuncCall& call, const Symbol& var) {
    if (call.args.size() == 1U) {
        BuiltinOp func_id = call.func_id;
        ExprPtr argument = call.args.front();

        // Try transcendental rules first (asin, acos, atan, sinh, cosh, tanh)
        auto transc_res = differentiate_transcendental(func_id, argument, var, context_);
        if (transc_res.is_ok()) {
            return transc_res;
        }
        if (transc_res.error().kind != CASErrorKind::Unimplemented) {
            return transc_res;
        }

        // d/dx[ln(abs(u))] = u'/u — must intercept before argument_derivative fails on abs(u)
        if (func_id == BuiltinOp::Ln || func_id == BuiltinOp::Log) {
            if (const auto* abs_call = expr_cast<FuncCall>(argument);
                abs_call != nullptr && abs_call->func_id == BuiltinOp::Abs && abs_call->args.size() == 1U) {
                ExprPtr inner = abs_call->args[0];
                auto inner_deriv = differentiate_once(inner, var);
                if (inner_deriv.is_ok()) {
                    return ok(make_binary(arena_, BinaryOp::Div, inner_deriv.value(), inner));
                }
            }
        }

        auto argument_derivative = differentiate_once(argument, var);
        if (argument_derivative.is_error()) {
            return argument_derivative;
        }

        ExprPtr outer;
        if (func_id == BuiltinOp::Sin) {
            outer = make_function(arena_, "cos", {argument});
        } else if (func_id == BuiltinOp::Cos) {
            outer = make_unary(arena_, UnaryOp::Neg, make_function(arena_, "sin", {argument}));
        } else if (func_id == BuiltinOp::Tan) {
            outer = make_binary(
                arena_,
                BinaryOp::Div,
                make_integer(arena_, 1),
                make_power(arena_, make_function(arena_, "cos", {argument}), make_integer(arena_, 2)));
        } else if (func_id == BuiltinOp::Cot) {
            outer = make_unary(
                arena_,
                UnaryOp::Neg,
                make_binary(
                    arena_,
                    BinaryOp::Div,
                    make_integer(arena_, 1),
                    make_power(arena_, make_function(arena_, "sin", {argument}), make_integer(arena_, 2))));
        } else if (func_id == BuiltinOp::Sec) {
            outer = make_product(arena_, {
                make_function(arena_, "sec", {argument}),
                make_function(arena_, "tan", {argument}),
            });
        } else if (func_id == BuiltinOp::Csc) {
            outer = make_unary(
                arena_,
                UnaryOp::Neg,
                make_product(arena_, {
                    make_function(arena_, "csc", {argument}),
                    make_function(arena_, "cot", {argument}),
                }));
        } else if (func_id == BuiltinOp::Exp) {
            outer = make_function(arena_, "exp", {argument});
        } else if (func_id == BuiltinOp::Ln || func_id == BuiltinOp::Log) {
            outer = make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), argument);
        } else if (func_id == BuiltinOp::Sqrt) {
            outer = make_binary(
                arena_,
                BinaryOp::Div,
                make_integer(arena_, 1),
                make_product(arena_, {
                    make_integer(arena_, 2),
                    make_function(arena_, "sqrt", {argument}),
                }));
        } else if (func_id == BuiltinOp::Erf) {
            outer = make_product(arena_, {
                make_binary(arena_, BinaryOp::Div, make_integer(arena_, 2), 
                    make_function(arena_, "sqrt", {make_constant(arena_, MathConstant::Pi)})),
                make_function(arena_, "exp", {make_unary(arena_, UnaryOp::Neg, make_power(arena_, argument, make_integer(arena_, 2)))})
            });
        } else if (func_id == BuiltinOp::EllipticK) {
            // dK/dk = E(k)/(k(1-k²)) − K(k)/k
            ExprPtr k_arg = argument;
            ExprPtr E_k = make_function(arena_, "EllipticE", {k_arg});
            ExprPtr K_k = make_function(arena_, "EllipticK", {k_arg});
            ExprPtr k_sq = arena_.make<Binary>(BinaryOp::Pow, k_arg,
                make_integer(arena_, 2));
            ExprPtr one_minus_k_sq = arena_.make<Binary>(BinaryOp::Sub,
                make_integer(arena_, 1), k_sq);
            ExprPtr denom = arena_.make<Binary>(BinaryOp::Mul, k_arg,
                one_minus_k_sq);
            ExprPtr first = arena_.make<Binary>(BinaryOp::Div, E_k, denom);
            ExprPtr second = arena_.make<Binary>(BinaryOp::Div, K_k, k_arg);
            outer = arena_.make<Binary>(BinaryOp::Sub, first, second);
        } else if (func_id == BuiltinOp::EllipticE) {
            // dE/dk = (E(k) − K(k))/k
            ExprPtr k_arg = argument;
            ExprPtr E_k = make_function(arena_, "EllipticE", {k_arg});
            ExprPtr K_k = make_function(arena_, "EllipticK", {k_arg});
            ExprPtr E_minus_K = arena_.make<Binary>(BinaryOp::Sub, E_k, K_k);
            outer = arena_.make<Binary>(BinaryOp::Div, E_minus_K, k_arg);
        } else if (func_id == BuiltinOp::Gamma ) {
            outer = make_product(arena_, {
                make_function(arena_, "Gamma", {argument}),
                make_function(arena_, "polygamma", {make_integer(arena_, 0), argument})
            });
        } else if (func_id == BuiltinOp::Abs) {
            if (const auto* symbol = expr_cast<Symbol>(argument);
                symbol != nullptr && !context_.assumptions().could_be_zero(*symbol)) {
                outer = make_function(arena_, "sign", {argument});
            } else if (!depends_on(argument, var)) {
                outer = make_function(arena_, "sign", {argument});
            } else {
                return fail<ExprPtr>(make_error(
                    CASErrorKind::Unimplemented,
                    "Differentiation of abs requires a nonzero assumption on the direct symbolic argument"));
            }
        } else if (call.name == "asinh" || call.name == "arcsinh") {
            outer = make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1),
                make_function(arena_, "sqrt",
                    {make_sum(arena_, {
                        make_power(arena_, argument, make_integer(arena_, 2)),
                        make_integer(arena_, 1)})}));
        } else if (call.name == "acosh" || call.name == "arccosh") {
            outer = make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1),
                make_function(arena_, "sqrt",
                    {make_sum(arena_, {
                        make_power(arena_, argument, make_integer(arena_, 2)),
                        make_integer(arena_, -1)})}));
        } else if (call.name == "atanh" || call.name == "arctanh"
                   || call.name == "acoth" || call.name == "arccoth") {
            outer = make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1),
                make_sum(arena_, {
                    make_integer(arena_, 1),
                    make_unary(arena_, UnaryOp::Neg,
                        make_power(arena_, argument, make_integer(arena_, 2)))}));
        } else {
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation is not implemented for function '" + call.name + "'"));
        }

        if (!depends_on(argument, var)) {
            return ok(make_integer(arena_, 0));
        }
        if (expr_is<IntegerLit>(argument_derivative.value())) {
            const auto& integer = expr_ref<IntegerLit>(argument_derivative.value());
            if (integer.value == BigInt(1)) {
                return ok(outer);
            }
        }

        return ok(make_product(arena_, {outer, argument_derivative.value()}));
    } else if (call.args.size() == 2U) {
        if (call.func_id == BuiltinOp::BesselZero) {
            return ok(make_integer(arena_, 0));
        }
        if (call.func_id == BuiltinOp::BesselJ || call.func_id == BuiltinOp::BesselY || call.func_id == BuiltinOp::BesselI || call.func_id == BuiltinOp::BesselK) {
            ExprPtr nu = call.args[0];
            ExprPtr x = call.args[1];

            if (depends_on(nu, var)) {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation of Bessel functions w.r.t order is not implemented"));
            }

            auto x_derivative = differentiate_once(x, var);
            if (x_derivative.is_error()) return x_derivative;

            ExprPtr outer;
            if (call.func_id == BuiltinOp::BesselJ) {
                outer = make_product(arena_, {
                    make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), make_integer(arena_, 2)),
                    make_sum(arena_, {
                        make_function(arena_, "BesselJ", {make_sum(arena_, {nu, make_integer(arena_, -1)}), x}),
                        make_unary(arena_, UnaryOp::Neg, make_function(arena_, "BesselJ", {make_sum(arena_, {nu, make_integer(arena_, 1)}), x}))
                    })
                });
            } else if (call.func_id == BuiltinOp::BesselY) {
                outer = make_product(arena_, {
                    make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), make_integer(arena_, 2)),
                    make_sum(arena_, {
                        make_function(arena_, "BesselY", {make_sum(arena_, {nu, make_integer(arena_, -1)}), x}),
                        make_unary(arena_, UnaryOp::Neg, make_function(arena_, "BesselY", {make_sum(arena_, {nu, make_integer(arena_, 1)}), x}))
                    })
                });
            } else if (call.func_id == BuiltinOp::BesselI) {
                outer = make_product(arena_, {
                    make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), make_integer(arena_, 2)),
                    make_sum(arena_, {
                        make_function(arena_, "BesselI", {make_sum(arena_, {nu, make_integer(arena_, -1)}), x}),
                        make_function(arena_, "BesselI", {make_sum(arena_, {nu, make_integer(arena_, 1)}), x})
                    })
                });
            } else if (call.func_id == BuiltinOp::BesselK) {
                outer = make_product(arena_, {
                    make_binary(arena_, BinaryOp::Div, make_unary(arena_, UnaryOp::Neg, make_integer(arena_, 1)), make_integer(arena_, 2)),
                    make_sum(arena_, {
                        make_function(arena_, "BesselK", {make_sum(arena_, {nu, make_integer(arena_, -1)}), x}),
                        make_function(arena_, "BesselK", {make_sum(arena_, {nu, make_integer(arena_, 1)}), x})
                    })
                });
            }

            return ok(make_product(arena_, {outer, x_derivative.value()}));
        }
        if (call.func_id == BuiltinOp::Hypergeometric0F1) {
            ExprPtr b = call.args[0];
            ExprPtr z = call.args[1];
            if (depends_on(b, var)) {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "Differentiation of Hypergeometric0F1 w.r.t parametri non supportata"));
            }
            auto z_d = differentiate_once(z, var);
            if (z_d.is_error()) return z_d;
            ExprPtr b_plus_1 = make_sum(arena_, {b, make_integer(arena_, 1)});
            ExprPtr new_fc = make_function(arena_, "Hypergeometric0F1", {b_plus_1, z});
            ExprPtr div_expr = arena_.make<Binary>(BinaryOp::Div, make_integer(arena_, 1), b);
            ExprPtr result = make_product(arena_, {div_expr, new_fc, z_d.value()});
            return ok(result);
        }
    } else if (call.args.size() == 3U) {
        if (call.func_id == BuiltinOp::Hypergeometric1F1) {
            ExprPtr a = call.args[0];
            ExprPtr b = call.args[1];
            ExprPtr z = call.args[2];
            if (depends_on(a, var) || depends_on(b, var)) {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "Differentiation of Hypergeometric1F1 w.r.t parametri non supportata"));
            }
            auto z_d = differentiate_once(z, var);
            if (z_d.is_error()) return z_d;
            ExprPtr a_plus_1 = make_sum(arena_, {a, make_integer(arena_, 1)});
            ExprPtr b_plus_1 = make_sum(arena_, {b, make_integer(arena_, 1)});
            ExprPtr new_fc = make_function(arena_, "Hypergeometric1F1",
                {a_plus_1, b_plus_1, z});
            ExprPtr a_over_b = arena_.make<Binary>(BinaryOp::Div, a, b);
            return ok(make_product(arena_, {a_over_b, new_fc, z_d.value()}));
        }
    } else if (call.args.size() == 4U) {
        if (call.func_id == BuiltinOp::Hypergeometric2F1) {
            ExprPtr a = call.args[0];
            ExprPtr b = call.args[1];
            ExprPtr c = call.args[2];
            ExprPtr z = call.args[3];
            if (depends_on(a, var) || depends_on(b, var) || depends_on(c, var)) {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "Differentiation of Hypergeometric2F1 w.r.t parametri non supportata"));
            }
            auto z_d = differentiate_once(z, var);
            if (z_d.is_error()) return z_d;
            ExprPtr a1 = make_sum(arena_, {a, make_integer(arena_, 1)});
            ExprPtr b1 = make_sum(arena_, {b, make_integer(arena_, 1)});
            ExprPtr c1 = make_sum(arena_, {c, make_integer(arena_, 1)});
            ExprPtr new_fc = make_function(arena_, "Hypergeometric2F1",
                {a1, b1, c1, z});
            ExprPtr ab = make_product(arena_, {a, b});
            ExprPtr ab_over_c = arena_.make<Binary>(BinaryOp::Div, ab, c);
            return ok(make_product(arena_, {ab_over_c, new_fc, z_d.value()}));
        }
    }

    // Meijer G (variable arity): §6.5 h=+1 theta-shift (DLMF 16.19.5 family,
    // numerically certified — see meijerg_derivative_shift):
    //   d/dx G(u|a;b) = (u'/u) * G^{m,n+1}_{p+1,q+1}(u | 0,a ; b, +1).
    if (call.func_id == BuiltinOp::MeijerG) {
        auto view = symbolic::view_meijerg(call);
        if (view.is_error()) return fail<ExprPtr>(view.error());
        for (ExprPtr p : view.value().a) {
            if (depends_on(p, var)) {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "Differentiation of MeijerG w.r.t. parametri non supportata"));
            }
        }
        for (ExprPtr p : view.value().b) {
            if (depends_on(p, var)) {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "Differentiation of MeijerG w.r.t. parametri non supportata"));
            }
        }
        ExprPtr z = view.value().z;
        auto z_d = differentiate_once(z, var);
        if (z_d.is_error()) return z_d;
        auto shifted = symbolic::meijerg_derivative_shift(context_, call);
        if (shifted.is_error()) return shifted;
        ExprPtr chain = arena_.make<Binary>(BinaryOp::Div, z_d.value(), z);
        return ok(make_product(arena_, {chain, shifted.value()}));
    }

    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation is not implemented for function '" + call.name + "' with " + std::to_string(call.args.size()) + " arguments"));
}

}  // namespace cas::calculus
