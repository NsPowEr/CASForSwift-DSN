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

    if (unary.op == UnaryOp::Factorial) {
        // A44: postfix `u!` and `factorial(u)` denote the same function, so the
        // derivative rule lives in exactly one place (the FuncCall branch in
        // differentiate_function). Rewriting here instead of duplicating the
        // Γ(u+1)·ψ(u+1) formula keeps a single source of truth for the identity.
        return differentiate_once(
            make_function(arena_, std::string(builtin_op_name(BuiltinOp::Factorial)),
                          {unary.operand}),
            var);
    }

    return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unknown unary operator in differentiation"));
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
            // A37: the canonical spelling is `gamma` (builtin_op_name). Emitting
            // "Gamma" produced a FuncCall that builtin_from_name does NOT
            // resolve, so the derivative carried an opaque unknown function and
            // no downstream simplify or equality could recognise it as Γ.
            outer = make_product(arena_, {
                make_function(arena_, std::string(builtin_op_name(BuiltinOp::Gamma)), {argument}),
                make_function(arena_, "polygamma", {make_integer(arena_, 0), argument})
            });
        } else if (func_id == BuiltinOp::ExpIntegralEi) {
            // A43 §3 (verificata mpmath): d/du Ei(u) = e^u/u.
            outer = make_binary(arena_, BinaryOp::Div,
                make_function(arena_, "exp", {argument}), argument);
        } else if (func_id == BuiltinOp::SinIntegral) {
            // d/du Si(u) = sin(u)/u
            outer = make_binary(arena_, BinaryOp::Div,
                make_function(arena_, "sin", {argument}), argument);
        } else if (func_id == BuiltinOp::CosIntegral) {
            // d/du Ci(u) = cos(u)/u
            outer = make_binary(arena_, BinaryOp::Div,
                make_function(arena_, "cos", {argument}), argument);
        } else if (func_id == BuiltinOp::SinhIntegral) {
            // d/du Shi(u) = sinh(u)/u
            outer = make_binary(arena_, BinaryOp::Div,
                make_function(arena_, "sinh", {argument}), argument);
        } else if (func_id == BuiltinOp::CoshIntegral) {
            // d/du Chi(u) = cosh(u)/u
            outer = make_binary(arena_, BinaryOp::Div,
                make_function(arena_, "cosh", {argument}), argument);
        } else if (func_id == BuiltinOp::LogIntegral) {
            // d/du li(u) = 1/ln(u)
            outer = make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1),
                make_function(arena_, "ln", {argument}));
        } else if (func_id == BuiltinOp::Dilog) {
            // d/du Li2(u) = -ln(1-u)/u
            outer = make_unary(arena_, UnaryOp::Neg,
                make_binary(arena_, BinaryOp::Div,
                    make_function(arena_, "ln", {make_sum(arena_, {
                        make_integer(arena_, 1),
                        make_unary(arena_, UnaryOp::Neg, argument)})}),
                    argument));
        } else if (func_id == BuiltinOp::Erfi) {
            // d/du erfi(u) = (2/sqrt(pi))*e^{u^2}
            outer = make_product(arena_, {
                make_binary(arena_, BinaryOp::Div, make_integer(arena_, 2),
                    make_function(arena_, "sqrt", {arena_.make<Constant>(MathConstant::Pi)})),
                make_function(arena_, "exp",
                    {make_power(arena_, argument, make_integer(arena_, 2))})});
        } else if (func_id == BuiltinOp::Factorial) {
            // A44: d/du u! = Γ(u+1)·ψ(u+1). Same identity as Γ shifted by one
            // (u! = Γ(u+1)), so it reuses the canonical `gamma`/`polygamma`
            // spellings for the same reason documented in the Gamma branch:
            // a non-canonical name would not resolve through builtin_from_name
            // and the derivative would carry an opaque unknown function.
            ExprPtr arg_plus_one = make_sum(arena_, {argument, make_integer(arena_, 1)});
            outer = make_product(arena_, {
                make_function(arena_, std::string(builtin_op_name(BuiltinOp::Gamma)), {arg_plus_one}),
                make_function(arena_, "polygamma", {make_integer(arena_, 0), arg_plus_one})
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
    }

    // Arity > 1 (Bessel, incomplete gamma, hypergeometric, Meijer G) lives in
    // differentiate_rules_multiarg.cpp — anti-monolith split, no behaviour change.
    return differentiate_function_multiarg(call, var);
}

}  // namespace cas::calculus
