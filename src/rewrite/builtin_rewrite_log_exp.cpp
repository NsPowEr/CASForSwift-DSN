// builtin_rewrite_log_exp.cpp — Logarithmic and exponential rewrite rules.
#include "builtin_rewrite_internal.hpp"

#include <vector>
#include <string>

namespace cas::symbolic {

void add_log_exp_rules(std::vector<RewriteRule>& rules, AstArena& arena) {
    const ExprPtr wildcard = arena.make<Symbol>("x_");
    const ExprPtr zero = arena.make<IntegerLit>(BigInt(0));
    const ExprPtr one = arena.make<IntegerLit>(BigInt(1));
    const ExprPtr e = arena.make<Constant>(MathConstant::E);

    rules.push_back(RewriteRule{
        .pattern = arena.make<FuncCall>("exp", std::vector<ExprPtr>{zero}),
        .replacement = one,
        .condition = {},
    });
    rules.push_back(RewriteRule{
        .pattern = arena.make<FuncCall>("exp", std::vector<ExprPtr>{one}),
        .replacement = e,
        .condition = {},
    });
    rules.push_back(RewriteRule{
        .pattern = arena.make<FuncCall>("ln", std::vector<ExprPtr>{one}),
        .replacement = zero,
        .condition = {},
    });
    rules.push_back(RewriteRule{
        .pattern = arena.make<FuncCall>("ln", std::vector<ExprPtr>{e}),
        .replacement = one,
        .condition = {},
    });
    rules.push_back(RewriteRule{
        .pattern = arena.make<FuncCall>(
            "ln",
            std::vector<ExprPtr>{
                arena.make<Binary>(
                    BinaryOp::Pow,
                    e,
                    wildcard),
            }),
        .replacement = wildcard,
        .condition = {},
    });
    rules.push_back(RewriteRule{
        .pattern = arena.make<FuncCall>(
            "ln",
            std::vector<ExprPtr>{
                arena.make<FuncCall>("exp", std::vector<ExprPtr>{wildcard}),
            }),
        .replacement = wildcard,
        .condition = {},
    });
    rules.push_back(RewriteRule{
        .pattern = arena.make<Binary>(
            BinaryOp::Pow,
            e,
            arena.make<FuncCall>("ln", std::vector<ExprPtr>{wildcard})),
        .replacement = wildcard,
        .condition = [](const MatchMap& m) {
            return exact_expr_is_positive(m.at("x_"));
        },
    });
    rules.push_back(RewriteRule{
        .pattern = arena.make<FuncCall>(
            "exp",
            std::vector<ExprPtr>{
                arena.make<FuncCall>("ln", std::vector<ExprPtr>{wildcard}),
            }),
        .replacement = wildcard,
        .condition = [](const MatchMap& m) {
            return exact_expr_is_positive(m.at("x_"));
        },
    });
}

Result<ExprPtr> try_rewrite_log_exp(ExprPtr expr, AstArena& arena, const Assumptions* assumptions) {
    if (const auto* power = expr_cast<Binary>(expr);
        power != nullptr && power->op == BinaryOp::Pow &&
        expr_cast<Constant>(power->left) != nullptr &&
        expr_ref<Constant>(power->left).value == MathConstant::E) {
        const auto* ln_call = expr_cast<FuncCall>(power->right);
        if (ln_call != nullptr &&
            ln_call->func_id == BuiltinOp::Ln &&
            ln_call->args.size() == 1U &&
            expr_is_positive_under_assumptions(ln_call->args.front(), assumptions)) {
            return ok(ln_call->args.front());
        }
    }

    if (const auto* call = expr_cast<FuncCall>(expr); call != nullptr &&
        call->func_id == BuiltinOp::Ln && call->args.size() == 1U) {
        if (const auto* quotient = expr_cast<Binary>(call->args.front());
            quotient != nullptr && quotient->op == BinaryOp::Div &&
            expr_is_positive_under_assumptions(quotient->left, assumptions) &&
            expr_is_positive_under_assumptions(quotient->right, assumptions)) {
            return ok(arena.make<Sum>(std::vector<ExprPtr>{
                arena.make<FuncCall>("ln", std::vector<ExprPtr>{quotient->left}),
                arena.make<Unary>(
                    UnaryOp::Neg,
                    arena.make<FuncCall>("ln", std::vector<ExprPtr>{quotient->right})),
            }));
        }

        if (const auto* power = expr_cast<Binary>(call->args.front());
            power != nullptr && power->op == BinaryOp::Pow &&
            !(expr_cast<Constant>(power->left) != nullptr &&
              expr_ref<Constant>(power->left).value == MathConstant::E) &&
            expr_is_positive_under_assumptions(power->left, assumptions)) {
            return ok(arena.make<Product>(std::vector<ExprPtr>{
                power->right,
                arena.make<FuncCall>("ln", std::vector<ExprPtr>{power->left}),
            }));
        }

        if (const auto* product = expr_cast<Product>(call->args.front())) {
            std::vector<ExprPtr> terms;
            terms.reserve(product->factors.size());
            for (ExprPtr factor : product->factors) {
                if (!expr_is_positive_under_assumptions(factor, assumptions)) {
                    terms.clear();
                    break;
                }
                terms.push_back(arena.make<FuncCall>("ln", std::vector<ExprPtr>{factor}));
            }
            if (!terms.empty()) {
                return ok(arena.make<Sum>(std::move(terms)));
            }
        }

        if (const auto* sqrt_call = expr_cast<FuncCall>(call->args.front());
            sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt && sqrt_call->args.size() == 1U &&
            expr_is_positive_under_assumptions(sqrt_call->args.front(), assumptions)) {
            return ok(arena.make<Product>(std::vector<ExprPtr>{
                arena.make<RationalLit>(BigInt(1), BigInt(2)),
                arena.make<FuncCall>("ln", std::vector<ExprPtr>{sqrt_call->args.front()}),
            }));
        }
    }

    if (const auto* call = expr_cast<FuncCall>(expr); call != nullptr &&
        call->func_id == BuiltinOp::Exp && call->args.size() == 1U) {
        if (const auto* sum = expr_cast<Sum>(call->args.front()); sum != nullptr) {
            std::vector<ExprPtr> factors;
            factors.reserve(sum->terms.size());
            for (ExprPtr term : sum->terms) {
                factors.push_back(arena.make<FuncCall>("exp", std::vector<ExprPtr>{term}));
            }
            return ok(arena.make<Product>(std::move(factors)));
        }
        if (const auto* binary = expr_cast<Binary>(call->args.front());
            binary != nullptr && binary->op == BinaryOp::Add) {
            return ok(arena.make<Product>(std::vector<ExprPtr>{
                arena.make<FuncCall>("exp", std::vector<ExprPtr>{binary->left}),
                arena.make<FuncCall>("exp", std::vector<ExprPtr>{binary->right}),
            }));
        }
    }

    return ok(expr);
}

}  // namespace cas::symbolic
