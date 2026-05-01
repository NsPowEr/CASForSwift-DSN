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

[[nodiscard]] ExprPtr make_function(AstArena& arena, std::string name, std::vector<ExprPtr> args) {
    return arena.make<FuncCall>(std::move(name), std::move(args));
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

[[nodiscard]] std::string canonical_function_name(const std::string& name) {
    if (name == "asin") {
        return "arcsin";
    }
    if (name == "acos") {
        return "arccos";
    }
    if (name == "atan") {
        return "arctan";
    }
    return name;
}

}  // namespace

Result<ExprPtr> differentiate_transcendental(
    const std::string& name,
    ExprPtr argument,
    const Symbol& var,
    symbolic::CASContext& context) {
    
    AstArena& arena = context.arena();
    const std::string cname = canonical_function_name(name);

    ExprPtr outer;
    if (cname == "arcsin") {
        outer = make_binary(arena, BinaryOp::Div, make_integer(arena, 1),
            make_function(arena, "sqrt", {
                make_sum(arena, {
                    make_integer(arena, 1),
                    make_unary(arena, UnaryOp::Neg, make_power(arena, argument, make_integer(arena, 2)))
                })
            }));
    } else if (cname == "arccos") {
        outer = make_unary(arena, UnaryOp::Neg,
            make_binary(arena, BinaryOp::Div, make_integer(arena, 1),
                make_function(arena, "sqrt", {
                    make_sum(arena, {
                        make_integer(arena, 1),
                        make_unary(arena, UnaryOp::Neg, make_power(arena, argument, make_integer(arena, 2)))
                    })
                })));
    } else if (cname == "arctan") {
        outer = make_binary(arena, BinaryOp::Div, make_integer(arena, 1),
            make_sum(arena, {
                make_integer(arena, 1),
                make_power(arena, argument, make_integer(arena, 2))
            }));
    } else if (cname == "sinh") {
        outer = make_function(arena, "cosh", {argument});
    } else if (cname == "cosh") {
        outer = make_function(arena, "sinh", {argument});
    } else if (cname == "tanh") {
        outer = make_binary(arena, BinaryOp::Div, make_integer(arena, 1),
            make_power(arena, make_function(arena, "cosh", {argument}), make_integer(arena, 2)));
    } else {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Function '" + name + "' is not a transcendental function handled by diff_rules",
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
