#include "calculus_internal.hpp"
#include "cas/error.hpp"

#include <string>
#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] ExprPtr make_integer(AstArena& arena, long long value) {
    return arena.make<IntegerLit>(BigInt(value));
}

[[nodiscard]] ExprPtr make_unary(AstArena& arena, UnaryOp op, ExprPtr operand) {
    return arena.make<Unary>(op, operand);
}

[[nodiscard]] ExprPtr make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs) {
    return arena.make<Binary>(op, lhs, rhs);
}

[[nodiscard]] ExprPtr make_function(AstArena& arena, BuiltinOp op, std::vector<ExprPtr> args) {
    return arena.make<FuncCall>(op, std::move(args));
}

[[nodiscard]] ExprPtr make_product(AstArena& arena, std::vector<ExprPtr> factors) {
    return arena.make<Product>(std::move(factors));
}

[[nodiscard]] ExprPtr make_power(AstArena& arena, ExprPtr base, ExprPtr exponent) {
    return make_binary(arena, BinaryOp::Pow, base, exponent);
}

[[nodiscard]] ExprPtr make_sum(AstArena& arena, std::vector<ExprPtr> terms) {
    return arena.make<Sum>(std::move(terms));
}

// Removed canonical_function_name

}  // namespace

Result<ExprPtr> differentiate_transcendental(
    BuiltinOp func_id,
    ExprPtr argument,
    const Symbol& var,
    symbolic::CASContext& context) {
    
    AstArena& arena = context.arena();

    ExprPtr outer;
    if (func_id == BuiltinOp::Asin) {
        outer = make_binary(arena, BinaryOp::Div, make_integer(arena, 1),
            make_function(arena, BuiltinOp::Sqrt, {
                make_sum(arena, {
                    make_integer(arena, 1),
                    make_unary(arena, UnaryOp::Neg, make_power(arena, argument, make_integer(arena, 2)))
                })
            }));
    } else if (func_id == BuiltinOp::Acos) {
        outer = make_unary(arena, UnaryOp::Neg,
            make_binary(arena, BinaryOp::Div, make_integer(arena, 1),
                make_function(arena, BuiltinOp::Sqrt, {
                    make_sum(arena, {
                        make_integer(arena, 1),
                        make_unary(arena, UnaryOp::Neg, make_power(arena, argument, make_integer(arena, 2)))
                    })
                })));
    } else if (func_id == BuiltinOp::Atan) {
        outer = make_binary(arena, BinaryOp::Div, make_integer(arena, 1),
            make_sum(arena, {
                make_integer(arena, 1),
                make_power(arena, argument, make_integer(arena, 2))
            }));
    } else if (func_id == BuiltinOp::Sinh) {
        outer = make_function(arena, BuiltinOp::Cosh, {argument});
    } else if (func_id == BuiltinOp::Cosh) {
        outer = make_function(arena, BuiltinOp::Sinh, {argument});
    } else if (func_id == BuiltinOp::Tanh) {
        outer = make_binary(arena, BinaryOp::Div, make_integer(arena, 1),
            make_power(arena, make_function(arena, BuiltinOp::Cosh, {argument}), make_integer(arena, 2)));
    } else {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Function is not a transcendental function handled by diff_rules",
        });
    }

    // Chain rule: d/dx f(g(x)) = f'(g(x)) * g'(x)
    auto inner_derivative_res = diff(argument, var, 1U, context);
    if (inner_derivative_res.is_error()) {
        return inner_derivative_res;
    }
    ExprPtr inner_derivative = inner_derivative_res.value();

    // Optimization / Structural Sharing check
    if (expr_is<IntegerLit>(inner_derivative)) {
        const auto& val = expr_ref<IntegerLit>(inner_derivative).value;
        if (val == BigInt(1)) {
            return ok(outer);
        }
        if (val.is_zero()) {
            return ok(make_integer(arena, 0));
        }
    }

    return ok(make_product(arena, {outer, inner_derivative}));
}

}  // namespace cas::calculus
