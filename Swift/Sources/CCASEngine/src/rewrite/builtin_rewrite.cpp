#include "cas/symbolic.hpp"
#include "cas/rational.hpp"

#include <string>
#include <vector>

namespace cas::symbolic {
namespace {

[[nodiscard]] std::optional<Rational> exact_scalar_from_expr(ExprPtr expr) {
    if (!expr) {
        return std::nullopt;
    }

    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return Rational(integer->value);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return Rational(rational->numerator, rational->denominator);
    }
    return std::nullopt;
}

[[nodiscard]] bool exact_expr_is_positive(ExprPtr expr) {
    const auto exact = exact_scalar_from_expr(expr);
    return exact.has_value() &&
        !exact->numerator().is_zero() &&
        !exact->numerator().is_negative();
}

[[nodiscard]] bool exact_expr_is_nonnegative(ExprPtr expr) {
    const auto exact = exact_scalar_from_expr(expr);
    return exact.has_value() && !exact->numerator().is_negative();
}

[[nodiscard]] bool exact_expr_is_negative(ExprPtr expr) {
    const auto exact = exact_scalar_from_expr(expr);
    return exact.has_value() && exact->numerator().is_negative();
}

[[nodiscard]] bool is_odd_parity_function(const std::string& name) {
    return name == "sin" ||
        name == "tan" ||
        name == "cot" ||
        name == "csc" ||
        name == "sinh" ||
        name == "tanh" ||
        name == "coth";
}

[[nodiscard]] bool is_even_parity_function(const std::string& name) {
    return name == "cos" ||
        name == "sec" ||
        name == "cosh";
}

[[nodiscard]] bool is_parity_rewrite_function(const std::string& name) {
    return is_odd_parity_function(name) || is_even_parity_function(name);
}

[[nodiscard]] bool expr_is_positive_under_assumptions(ExprPtr expr, const Assumptions* assumptions) {
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
            case MathConstant::Infinity:
                return true;
            case MathConstant::I:
            case MathConstant::NaN:
                return false;
        }
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Pow) {
            // base^exponent is positive if base is positive
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

[[nodiscard]] bool expr_is_nonnegative_under_assumptions(ExprPtr expr, const Assumptions* assumptions) {
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

enum class SquareFunctionKind : std::uint8_t {
    None,
    Sin,
    Cos,
};

[[nodiscard]] SquareFunctionKind square_function_kind(ExprPtr expr) {
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

    if (call->name == "sin") {
        return SquareFunctionKind::Sin;
    }
    if (call->name == "cos") {
        return SquareFunctionKind::Cos;
    }
    return SquareFunctionKind::None;
}

[[nodiscard]] ExprPtr square_function_argument(ExprPtr expr, SquareFunctionKind kind) {
    if (square_function_kind(expr) != kind) {
        return ExprPtr{};
    }

    const auto* power = expr_cast<Binary>(expr);
    const auto* call = power != nullptr ? expr_cast<FuncCall>(power->left) : nullptr;
    return call != nullptr && call->args.size() == 1U ? call->args.front() : ExprPtr{};
}

[[nodiscard]] bool may_match_builtin_rewrite(ExprPtr expr) {
    if (const auto* binary = expr_cast<Binary>(expr)) {
        return binary->op == BinaryOp::Pow &&
            expr_cast<Constant>(binary->left) != nullptr &&
            expr_ref<Constant>(binary->left).value == MathConstant::E &&
            expr_cast<FuncCall>(binary->right) != nullptr &&
            expr_ref<FuncCall>(binary->right).name == "ln" &&
            expr_ref<FuncCall>(binary->right).args.size() == 1U;
    }

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        return call->args.size() == 1U &&
            (is_parity_rewrite_function(call->name) ||
             call->name == "exp" ||
             call->name == "ln" ||
             call->name == "sqrt");
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
        const ExprPtr wildcard = rules_arena_.make<Symbol>("x_");
        const ExprPtr zero = rules_arena_.make<IntegerLit>(BigInt(0));
        const ExprPtr one = rules_arena_.make<IntegerLit>(BigInt(1));
        const ExprPtr e = rules_arena_.make<Constant>(MathConstant::E);
        const ExprPtr exponent = rules_arena_.make<IntegerLit>(BigInt(2));

        const ExprPtr sin_square = rules_arena_.make<Binary>(
            BinaryOp::Pow,
            rules_arena_.make<FuncCall>("sin", std::vector<ExprPtr>{wildcard}),
            exponent);
        const ExprPtr cos_square = rules_arena_.make<Binary>(
            BinaryOp::Pow,
            rules_arena_.make<FuncCall>("cos", std::vector<ExprPtr>{wildcard}),
            exponent);

        rules_.push_back(RewriteRule{
            .pattern = rules_arena_.make<Sum>(std::vector<ExprPtr>{sin_square, cos_square}),
            .replacement = rules_arena_.make<IntegerLit>(BigInt(1)),
            .condition = {},
        });
        rules_.push_back(RewriteRule{
            .pattern = rules_arena_.make<FuncCall>("sin", std::vector<ExprPtr>{zero}),
            .replacement = zero,
            .condition = {},
        });
        rules_.push_back(RewriteRule{
            .pattern = rules_arena_.make<FuncCall>("cos", std::vector<ExprPtr>{zero}),
            .replacement = one,
            .condition = {},
        });
        rules_.push_back(RewriteRule{
            .pattern = rules_arena_.make<FuncCall>("exp", std::vector<ExprPtr>{zero}),
            .replacement = one,
            .condition = {},
        });
        rules_.push_back(RewriteRule{
            .pattern = rules_arena_.make<FuncCall>("exp", std::vector<ExprPtr>{one}),
            .replacement = e,
            .condition = {},
        });
        rules_.push_back(RewriteRule{
            .pattern = rules_arena_.make<FuncCall>("ln", std::vector<ExprPtr>{one}),
            .replacement = zero,
            .condition = {},
        });
        rules_.push_back(RewriteRule{
            .pattern = rules_arena_.make<FuncCall>("ln", std::vector<ExprPtr>{e}),
            .replacement = one,
            .condition = {},
        });
        rules_.push_back(RewriteRule{
            .pattern = rules_arena_.make<FuncCall>(
                "ln",
                std::vector<ExprPtr>{
                    rules_arena_.make<Binary>(
                        BinaryOp::Pow,
                        e,
                        wildcard),
                }),
            .replacement = wildcard,
            .condition = {},
        });
    }

    [[nodiscard]] Result<ExprPtr> try_rewrite(
        ExprPtr expr,
        AstArena& arena,
        const Assumptions* assumptions) const override {
        if (!may_match_builtin_rewrite(expr)) {
            return ok(expr);
        }

        if (const auto* power = expr_cast<Binary>(expr);
            power != nullptr && power->op == BinaryOp::Pow &&
            expr_cast<Constant>(power->left) != nullptr &&
            expr_ref<Constant>(power->left).value == MathConstant::E) {
            const auto* ln_call = expr_cast<FuncCall>(power->right);
            if (ln_call != nullptr &&
                ln_call->name == "ln" &&
                ln_call->args.size() == 1U &&
                expr_is_positive_under_assumptions(ln_call->args.front(), assumptions)) {
                return ok(ln_call->args.front());
            }
        }

        if (const auto* sum = expr_cast<Sum>(expr); sum != nullptr) {
            for (std::size_t sin_index = 0; sin_index < sum->terms.size(); ++sin_index) {
                ExprPtr sin_argument = square_function_argument(sum->terms[sin_index], SquareFunctionKind::Sin);
                if (!sin_argument) {
                    continue;
                }

                for (std::size_t cos_index = 0; cos_index < sum->terms.size(); ++cos_index) {
                    if (sin_index == cos_index) {
                        continue;
                    }

                    ExprPtr cos_argument = square_function_argument(sum->terms[cos_index], SquareFunctionKind::Cos);
                    if (!cos_argument || !structural_equal(sin_argument, cos_argument)) {
                        continue;
                    }

                    std::vector<ExprPtr> reduced_terms;
                    reduced_terms.reserve(sum->terms.size() - 1U);
                    for (std::size_t index = 0; index < sum->terms.size(); ++index) {
                        if (index != sin_index && index != cos_index) {
                            reduced_terms.push_back(sum->terms[index]);
                        }
                    }
                    reduced_terms.push_back(arena.make<IntegerLit>(BigInt(1)));
                    return ok(
                        reduced_terms.size() == 1U ? reduced_terms.front() : arena.make<Sum>(std::move(reduced_terms)));
                }
            }
        }

        if (const auto* call = expr_cast<FuncCall>(expr); call != nullptr &&
            call->name == "ln" && call->args.size() == 1U) {
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
                sqrt_call != nullptr && sqrt_call->name == "sqrt" && sqrt_call->args.size() == 1U &&
                expr_is_positive_under_assumptions(sqrt_call->args.front(), assumptions)) {
                return ok(arena.make<Product>(std::vector<ExprPtr>{
                    arena.make<RationalLit>(BigInt(1), BigInt(2)),
                    arena.make<FuncCall>("ln", std::vector<ExprPtr>{sqrt_call->args.front()}),
                }));
            }
        }

        if (const auto* call = expr_cast<FuncCall>(expr); call != nullptr &&
            call->args.size() == 1U) {
            const auto* unary = expr_cast<Unary>(call->args.front());
            if (unary != nullptr && unary->op == UnaryOp::Neg) {
                if (is_odd_parity_function(call->name)) {
                    return ok(arena.make<Unary>(
                        UnaryOp::Neg,
                        arena.make<FuncCall>(call->name, std::vector<ExprPtr>{unary->operand})));
                }
                if (is_even_parity_function(call->name)) {
                    return ok(arena.make<FuncCall>(call->name, std::vector<ExprPtr>{unary->operand}));
                }
            }

            if (call->name == "tan") {
                return ok(arena.make<Binary>(
                    BinaryOp::Div,
                    arena.make<FuncCall>("sin", std::vector<ExprPtr>{call->args.front()}),
                    arena.make<FuncCall>("cos", std::vector<ExprPtr>{call->args.front()})));
            }
        }

        if (const auto* call = expr_cast<FuncCall>(expr); call != nullptr &&
            call->name == "sqrt" && call->args.size() == 1U) {
            if (const auto* quotient = expr_cast<Binary>(call->args.front());
                quotient != nullptr && quotient->op == BinaryOp::Div &&
                expr_is_nonnegative_under_assumptions(quotient->left, assumptions) &&
                expr_is_positive_under_assumptions(quotient->right, assumptions)) {
                return ok(arena.make<Binary>(
                    BinaryOp::Div,
                    arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{quotient->left}),
                    arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{quotient->right})));
            }

            if (const auto* product = expr_cast<Product>(call->args.front())) {
                std::vector<ExprPtr> numerators;
                std::vector<ExprPtr> denominators;
                for (ExprPtr factor : product->factors) {
                    if (const auto* bin = expr_cast<Binary>(factor); 
                        bin && bin->op == BinaryOp::Pow && 
                        exact_expr_is_negative(bin->right)) {
                        denominators.push_back(bin->left);
                    } else {
                        numerators.push_back(factor);
                    }
                }

                if (!denominators.empty()) {
                    ExprPtr num = numerators.empty() ? arena.make<IntegerLit>(BigInt(1)) : 
                                 (numerators.size() == 1 ? numerators[0] : arena.make<Product>(std::move(numerators)));
                    ExprPtr den = denominators.size() == 1 ? denominators[0] : arena.make<Product>(std::move(denominators));
                    return ok(arena.make<Binary>(
                        BinaryOp::Div,
                        arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{num}),
                        arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{den})));
                }

                std::vector<ExprPtr> terms;
                terms.reserve(product->factors.size());
                for (ExprPtr factor : product->factors) {
                    if (!expr_is_nonnegative_under_assumptions(factor, assumptions)) {
                        terms.clear();
                        break;
                    }
                    terms.push_back(arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{factor}));
                }
                if (!terms.empty()) {
                    return ok(arena.make<Product>(std::move(terms)));
                }
            }
        }

        if (const auto* call = expr_cast<FuncCall>(expr); call != nullptr &&
            call->name == "exp" && call->args.size() == 1U) {
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

        return apply_rule_set(expr, rules_, arena);
    }

private:
    AstArena rules_arena_;
    std::vector<RewriteRule> rules_;
};

}  // namespace

const RewriteProvider& default_rewrite_provider() {
    static const BuiltinRewriteProvider provider;
    return provider;
}

}  // namespace cas::symbolic
