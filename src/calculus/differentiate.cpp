#include "cas/calculus.hpp"
#include "calculus_internal.hpp"
#include "cas/error.hpp"

#include <string>
#include <utility>
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

[[nodiscard]] ExprPtr make_constant(AstArena& arena, MathConstant value) {
    return arena.make<Constant>(value);
}

[[nodiscard]] ExprPtr make_unary(AstArena& arena, UnaryOp op, ExprPtr operand) {
    return arena.make<Unary>(op, operand);
}

[[nodiscard]] ExprPtr make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs) {
    return arena.make<Binary>(op, lhs, rhs);
}

[[nodiscard]] ExprPtr make_sum(AstArena& arena, std::vector<ExprPtr> terms) {
    if (terms.empty()) {
        return make_integer(arena, 0);
    }
    if (terms.size() == 1U) {
        return terms.front();
    }
    return arena.make<Sum>(std::move(terms));
}

[[nodiscard]] ExprPtr make_product(AstArena& arena, std::vector<ExprPtr> factors) {
    if (factors.empty()) {
        return make_integer(arena, 1);
    }
    if (factors.size() == 1U) {
        return factors.front();
    }
    return arena.make<Product>(std::move(factors));
}

[[nodiscard]] ExprPtr make_power(AstArena& arena, ExprPtr base, ExprPtr exponent) {
    return make_binary(arena, BinaryOp::Pow, base, exponent);
}

[[nodiscard]] ExprPtr make_function(AstArena& arena, std::string name, std::vector<ExprPtr> args) {
    return arena.make<FuncCall>(std::move(name), std::move(args));
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

class Differentiator {
public:
    explicit Differentiator(symbolic::CASContext& context) noexcept : context_(context), arena_(context.arena()) {}

    [[nodiscard]] Result<ExprPtr> differentiate(ExprPtr expr, const Symbol& var, unsigned int order) {
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

private:
    [[nodiscard]] Result<ExprPtr> differentiate_once(ExprPtr expr, const Symbol& var) {
        if (!expr) {
            return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot differentiate a null expression"));
        }

        if (expr_is<IntegerLit>(expr) || expr_is<RationalLit>(expr) || expr_is<Constant>(expr)) {
            return ok(make_integer(arena_, 0));
        }
        if (expr_is<DecimalLit>(expr)) {
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Decimal literals are not supported in symbolic differentiation"));
        }
        if (const auto* symbol = expr_cast<Symbol>(expr)) {
            return ok(make_integer(arena_, symbol->name == var.name ? 1 : 0));
        }
        if (const auto* unary = expr_cast<Unary>(expr)) {
            return differentiate_unary(*unary, var);
        }
        if (const auto* binary = expr_cast<Binary>(expr)) {
            return differentiate_binary(*binary, var);
        }
        if (const auto* call = expr_cast<FuncCall>(expr)) {
            return differentiate_function(*call, var);
        }
        if (const auto* sum = expr_cast<Sum>(expr)) {
            return differentiate_sum(*sum, var);
        }
        if (const auto* product = expr_cast<Product>(expr)) {
            return differentiate_product(*product, var);
        }
        if (const auto* derivative = expr_cast<Derivative>(expr)) {
            return differentiate(derivative->expression, derivative->variable, derivative->order + 1U);
        }

        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation is not implemented for this expression kind"));
    }

    [[nodiscard]] Result<ExprPtr> differentiate_unary(const Unary& unary, const Symbol& var) {
        if (unary.op == UnaryOp::Neg) {
            auto operand = differentiate_once(unary.operand, var);
            if (operand.is_error()) {
                return operand;
            }
            return ok(make_unary(arena_, UnaryOp::Neg, operand.value()));
        }

        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Factorial differentiation is not implemented"));
    }

    [[nodiscard]] Result<ExprPtr> differentiate_binary(const Binary& binary, const Symbol& var) {
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
        }

        return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unknown binary operator"));
    }

    [[nodiscard]] Result<ExprPtr> differentiate_sum(const Sum& sum, const Symbol& var) {
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

    [[nodiscard]] Result<ExprPtr> differentiate_product(const Product& product, const Symbol& var) {
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

    [[nodiscard]] Result<ExprPtr> differentiate_power(const Binary& power, const Symbol& var) {
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

    [[nodiscard]] Result<ExprPtr> differentiate_function(const FuncCall& call, const Symbol& var) {
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
            if (func_id == BuiltinOp::Ln) {
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
            } else if (func_id == BuiltinOp::Ln) {
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
            // Binary functions
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
        }

        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation is not implemented for function '" + call.name + "' with " + std::to_string(call.args.size()) + " arguments"));
    }

    symbolic::CASContext& context_;
    AstArena& arena_;
};

}  // namespace

Result<ExprPtr> diff(ExprPtr expr, const Symbol& var, unsigned int order, symbolic::CASContext& ctx) {
    auto differentiated = Differentiator(ctx).differentiate(expr, var, order);
    if (differentiated.is_error()) {
        return differentiated;
    }
    return symbolic::materialize_expr(differentiated.value(), ctx.arena());
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
