// builtin_rewrite_trig.cpp — Trigonometric parity and identity rewrite rules.
#include "builtin_rewrite_internal.hpp"

#include <vector>
#include <string>

namespace cas::symbolic {

bool is_odd_parity_function(BuiltinOp func_id) {
    return func_id == BuiltinOp::Sin ||
        func_id == BuiltinOp::Tan ||
        func_id == BuiltinOp::Cot ||
        func_id == BuiltinOp::Csc ||
        func_id == BuiltinOp::Sinh ||
        func_id == BuiltinOp::Tanh ||
        func_id == BuiltinOp::Coth;
}

bool is_even_parity_function(BuiltinOp func_id) {
    return func_id == BuiltinOp::Cos ||
        func_id == BuiltinOp::Sec ||
        func_id == BuiltinOp::Cosh;
}

bool is_parity_rewrite_function(BuiltinOp func_id) {
    return is_odd_parity_function(func_id) || is_even_parity_function(func_id);
}

SquareFunctionKind square_function_kind(ExprPtr expr) {
    const auto* power = expr_cast<Binary>(expr);
    if (power == nullptr || power->op != BinaryOp::Pow) {
        return SquareFunctionKind::None;
    }

    const auto* exponent = expr_cast<IntegerLit>(power->right);
    if (exponent == nullptr || exponent->value != BigInt(2)) {
        return SquareFunctionKind::None;
    }

    const auto* call = expr_cast<FuncCall>(power->left);
    if (call == nullptr || call->args.size() != 1U) {
        return SquareFunctionKind::None;
    }

    if (call->func_id == BuiltinOp::Sin) {
        return SquareFunctionKind::Sin;
    }
    if (call->func_id == BuiltinOp::Cos) {
        return SquareFunctionKind::Cos;
    }
    return SquareFunctionKind::None;
}

void add_trig_rules(std::vector<RewriteRule>& rules, AstArena& arena) {
    const ExprPtr wildcard = arena.make<Symbol>("x_");
    const ExprPtr zero = arena.make<IntegerLit>(BigInt(0));
    const ExprPtr one = arena.make<IntegerLit>(BigInt(1));
    const ExprPtr exponent = arena.make<IntegerLit>(BigInt(2));

    const ExprPtr sin_square = arena.make<Binary>(
        BinaryOp::Pow,
        arena.make<FuncCall>("sin", std::vector<ExprPtr>{wildcard}),
        exponent);
    const ExprPtr cos_square = arena.make<Binary>(
        BinaryOp::Pow,
        arena.make<FuncCall>("cos", std::vector<ExprPtr>{wildcard}),
        exponent);

    rules.push_back(RewriteRule{
        .pattern = arena.make<Sum>(std::vector<ExprPtr>{sin_square, cos_square}),
        .replacement = one,
        .condition = {},
    });
    rules.push_back(RewriteRule{
        .pattern = arena.make<FuncCall>("sin", std::vector<ExprPtr>{zero}),
        .replacement = zero,
        .condition = {},
    });
    rules.push_back(RewriteRule{
        .pattern = arena.make<FuncCall>("cos", std::vector<ExprPtr>{zero}),
        .replacement = one,
        .condition = {},
    });
    rules.push_back(RewriteRule{
        .pattern = arena.make<FuncCall>("tan", std::vector<ExprPtr>{zero}),
        .replacement = zero,
        .condition = {},
    });

    auto add_parity_rule = [&](std::string_view name, bool is_odd) {
        const ExprPtr x = arena.make<Symbol>("x_");
        const ExprPtr neg_x = arena.make<Unary>(UnaryOp::Neg, x);
        if (is_odd) {
            rules.push_back(RewriteRule{
                .pattern = arena.make<FuncCall>(std::string(name), std::vector<ExprPtr>{neg_x}),
                .replacement = arena.make<Unary>(UnaryOp::Neg, 
                    arena.make<FuncCall>(std::string(name), std::vector<ExprPtr>{x})),
                .condition = {},
            });
        } else {
            rules.push_back(RewriteRule{
                .pattern = arena.make<FuncCall>(std::string(name), std::vector<ExprPtr>{neg_x}),
                .replacement = arena.make<FuncCall>(std::string(name), std::vector<ExprPtr>{x}),
                .condition = {},
            });
        }
    };

    add_parity_rule("sin", true);
    add_parity_rule("cos", false);
    add_parity_rule("tan", true);
    add_parity_rule("cot", true);
    add_parity_rule("sec", false);
    add_parity_rule("csc", true);
    add_parity_rule("sinh", true);
    add_parity_rule("cosh", false);
    add_parity_rule("tanh", true);
    add_parity_rule("coth", true);
}

Result<ExprPtr> try_rewrite_trig(ExprPtr expr, AstArena& arena) {
    if (const auto* call = expr_cast<FuncCall>(expr); call != nullptr &&
        call->args.size() == 1U) {
        if (call->func_id == BuiltinOp::Tan) {
            return ok(arena.make<Binary>(
                BinaryOp::Div,
                arena.make<FuncCall>("sin", std::vector<ExprPtr>{call->args.front()}),
                arena.make<FuncCall>("cos", std::vector<ExprPtr>{call->args.front()})));
        }
    }
    return ok(expr);
}

}  // namespace cas::symbolic
