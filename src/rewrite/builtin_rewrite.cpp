#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "cas/algebra.hpp"
#include "../symbolic/symbolic_internal.hpp"
#include "builtin_rewrite_internal.hpp"

#include <string>
#include <vector>

namespace cas::symbolic {


bool exact_expr_is_positive(ExprPtr expr) {
    const auto exact = exact_scalar_from_expr(expr);
    return exact.has_value() &&
        !exact->numerator().is_zero() &&
        !exact->numerator().is_negative();
}

bool exact_expr_is_nonnegative(ExprPtr expr) {
    const auto exact = exact_scalar_from_expr(expr);
    return exact.has_value() && !exact->numerator().is_negative();
}

bool exact_expr_is_negative(ExprPtr expr) {
    const auto exact = exact_scalar_from_expr(expr);
    return exact.has_value() && exact->numerator().is_negative();
}

bool expr_is_positive_under_assumptions(ExprPtr expr, const Assumptions* assumptions) {
    if (!expr) {
        return false;
    }

    if (exact_expr_is_positive(expr)) {
        return true;
    }

    if (const auto* symbol = expr_cast<Symbol>(expr)) {
        return assumptions != nullptr && assumptions->is_positive(*symbol);
    }

    if (const auto* constant = expr_cast<Constant>(expr)) {
        switch (constant->value) {
            case MathConstant::Pi:
            case MathConstant::E:
            case MathConstant::EulerGamma:
            case MathConstant::Infinity:
                return true;
            case MathConstant::NegInfinity:
            case MathConstant::ComplexInfinity:
            case MathConstant::Indeterminate:
            case MathConstant::I:
            case MathConstant::NaN:
                return false;
        }
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Pow) {
            return expr_is_positive_under_assumptions(binary->left, assumptions);
        }
        if (binary->op == BinaryOp::Div) {
            return expr_is_positive_under_assumptions(binary->left, assumptions) &&
                   expr_is_positive_under_assumptions(binary->right, assumptions);
        }
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr factor : product->factors) {
            if (!expr_is_positive_under_assumptions(factor, assumptions)) {
                return false;
            }
        }
        return !product->factors.empty();
    }

    return false;
}

bool expr_is_nonnegative_under_assumptions(ExprPtr expr, const Assumptions* assumptions) {
    if (!expr) {
        return false;
    }

    if (exact_expr_is_nonnegative(expr)) {
        return true;
    }

    if (const auto* symbol = expr_cast<Symbol>(expr)) {
        return assumptions != nullptr && assumptions->is_positive(*symbol);
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Pow) {
            return expr_is_positive_under_assumptions(binary->left, assumptions);
        }
        if (binary->op == BinaryOp::Div) {
            return expr_is_nonnegative_under_assumptions(binary->left, assumptions) &&
                   expr_is_positive_under_assumptions(binary->right, assumptions);
        }
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr factor : product->factors) {
            if (!expr_is_nonnegative_under_assumptions(factor, assumptions)) {
                return false;
            }
        }
        return !product->factors.empty();
    }

    return false;
}

bool may_match_builtin_rewrite(ExprPtr expr) {
    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Div) {
            return true;
        }
        return binary->op == BinaryOp::Pow &&
            expr_cast<Constant>(binary->left) != nullptr &&
            expr_ref<Constant>(binary->left).value == MathConstant::E &&
            expr_cast<FuncCall>(binary->right) != nullptr &&
            expr_ref<FuncCall>(binary->right).func_id == BuiltinOp::Ln &&
            expr_ref<FuncCall>(binary->right).args.size() == 1U;
    }

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        return call->args.size() == 1U &&
            (is_parity_rewrite_function(call->func_id) ||
             call->func_id == BuiltinOp::Exp ||
             call->func_id == BuiltinOp::Ln ||
             call->func_id == BuiltinOp::Sqrt);
    }

    const auto* sum = expr_cast<Sum>(expr);
    if (sum == nullptr || sum->terms.size() < 2U) {
        return false;
    }

    bool saw_sin_square = false;
    bool saw_cos_square = false;
    for (ExprPtr term : sum->terms) {
        switch (square_function_kind(term)) {
            case SquareFunctionKind::Sin:
                saw_sin_square = true;
                break;
            case SquareFunctionKind::Cos:
                saw_cos_square = true;
                break;
            case SquareFunctionKind::None:
                break;
        }
        if (saw_sin_square && saw_cos_square) {
            return true;
        }
    }

    return false;
}

class BuiltinRewriteProvider final : public RewriteProvider {
public:
    BuiltinRewriteProvider() {
        add_trig_rules(rules_, rules_arena_);
        add_log_exp_rules(rules_, rules_arena_);
    }

    [[nodiscard]] Result<ExprPtr> try_rewrite(
        ExprPtr expr,
        AstArena& arena,
        const Assumptions* assumptions,
        CASContext* context = nullptr) const override {
        if (!may_match_builtin_rewrite(expr)) {
            return ok(expr);
        }

        // 1. Algebraic & Rational
        auto alg_res = try_rewrite_algebraic(expr, arena, assumptions, context);
        if (alg_res.is_error() || alg_res.value() != expr) {
            return alg_res;
        }

        // 2. Log & Exp
        auto log_res = try_rewrite_log_exp(expr, arena, assumptions);
        if (log_res.is_error() || log_res.value() != expr) {
            return log_res;
        }

        // 3. Trig & Parity
        auto trig_res = try_rewrite_trig(expr, arena);
        if (trig_res.is_error() || trig_res.value() != expr) {
            return trig_res;
        }

        // 4. Fallback: general rule set mapping (pattern rules)
        return apply_rule_set(expr, rules_, arena);
    }

private:
    AstArena rules_arena_;
    std::vector<RewriteRule> rules_;
};

const RewriteProvider& default_rewrite_provider() {
    static const BuiltinRewriteProvider provider;
    return provider;
}

}  // namespace cas::symbolic
