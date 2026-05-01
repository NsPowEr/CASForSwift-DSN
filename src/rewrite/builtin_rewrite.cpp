#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "cas/algebra.hpp"
#include "../symbolic/symbolic_internal.hpp"
#include "../algebra/polynomial_internal.hpp"

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

[[nodiscard]] bool is_odd_parity_function(BuiltinOp func_id) {
    return func_id == BuiltinOp::Sin ||
        func_id == BuiltinOp::Tan ||
        func_id == BuiltinOp::Cot ||
        func_id == BuiltinOp::Csc ||
        func_id == BuiltinOp::Sinh ||
        func_id == BuiltinOp::Tanh ||
        func_id == BuiltinOp::Coth;
}

[[nodiscard]] bool is_even_parity_function(BuiltinOp func_id) {
    return func_id == BuiltinOp::Cos ||
        func_id == BuiltinOp::Sec ||
        func_id == BuiltinOp::Cosh;
}

[[nodiscard]] bool is_parity_rewrite_function(BuiltinOp func_id) {
    return is_odd_parity_function(func_id) || is_even_parity_function(func_id);
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

    if (call->func_id == BuiltinOp::Sin) {
        return SquareFunctionKind::Sin;
    }
    if (call->func_id == BuiltinOp::Cos) {
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
            .replacement = one,
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
            .pattern = rules_arena_.make<FuncCall>("tan", std::vector<ExprPtr>{zero}),
            .replacement = zero,
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
        rules_.push_back(RewriteRule{
            .pattern = rules_arena_.make<FuncCall>(
                "ln",
                std::vector<ExprPtr>{
                    rules_arena_.make<FuncCall>("exp", std::vector<ExprPtr>{wildcard}),
                }),
            .replacement = wildcard,
            .condition = {},
        });
        rules_.push_back(RewriteRule{
            .pattern = rules_arena_.make<Binary>(
                BinaryOp::Pow,
                e,
                rules_arena_.make<FuncCall>("ln", std::vector<ExprPtr>{wildcard})),
            .replacement = wildcard,
            .condition = [](const MatchMap& m) {
                return exact_expr_is_positive(m.at("x_"));
            },
        });
        rules_.push_back(RewriteRule{
            .pattern = rules_arena_.make<FuncCall>(
                "exp",
                std::vector<ExprPtr>{
                    rules_arena_.make<FuncCall>("ln", std::vector<ExprPtr>{wildcard}),
                }),
            .replacement = wildcard,
            .condition = [](const MatchMap& m) {
                return exact_expr_is_positive(m.at("x_"));
            },
        });

        auto add_parity_rule = [&](std::string_view name, bool is_odd) {
            const ExprPtr x = rules_arena_.make<Symbol>("x_");
            const ExprPtr neg_x = rules_arena_.make<Unary>(UnaryOp::Neg, x);
            if (is_odd) {
                rules_.push_back(RewriteRule{
                    .pattern = rules_arena_.make<FuncCall>(std::string(name), std::vector<ExprPtr>{neg_x}),
                    .replacement = rules_arena_.make<Unary>(UnaryOp::Neg, 
                        rules_arena_.make<FuncCall>(std::string(name), std::vector<ExprPtr>{x})),
                    .condition = {},
                });
            } else {
                rules_.push_back(RewriteRule{
                    .pattern = rules_arena_.make<FuncCall>(std::string(name), std::vector<ExprPtr>{neg_x}),
                    .replacement = rules_arena_.make<FuncCall>(std::string(name), std::vector<ExprPtr>{x}),
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

    [[nodiscard]] Result<ExprPtr> try_rewrite(
        ExprPtr expr,
        AstArena& arena,
        const Assumptions* assumptions,
        CASContext* context = nullptr) const override {
        if (!may_match_builtin_rewrite(expr)) {
            return ok(expr);
        }

        if (const auto* binary = expr_cast<Binary>(expr); binary != nullptr && binary->op == BinaryOp::Div) {
            if (context != nullptr) {
                // GCD Reduction for Rational Functions (Pilastro 5)
                auto find_var = [](ExprPtr e) -> std::optional<Symbol> {
                    std::optional<Symbol> found;
                    std::function<void(ExprPtr)> collect = [&](ExprPtr node) {
                        if (found.has_value() || !node) return;
                        if (const auto* s = expr_cast<Symbol>(node)) {
                            found = *s;
                            return;
                        }
                        visit_expr(node, [&](const auto& n) {
                            using T = std::decay_t<decltype(n)>;
                            if constexpr (std::is_same_v<T, Unary>) collect(n.operand);
                            else if constexpr (std::is_same_v<T, Binary>) { collect(n.left); collect(n.right); }
                            else if constexpr (std::is_same_v<T, FuncCall>) { for (auto a : n.args) collect(a); }
                            else if constexpr (std::is_same_v<T, Sum>) { for (auto t : n.terms) collect(t); }
                            else if constexpr (std::is_same_v<T, Product>) { for (auto f : n.factors) collect(f); }
                        });
                    };
                    collect(e);
                    return found;
                };

                auto var = find_var(expr);
                if (var.has_value()) {
                    // Tenta di convertire in polinomi interi per divisione sicura senza simplify
                    auto p1_res = algebra::parse_polynomial(binary->left, *var, *context);
                    auto p2_res = algebra::parse_polynomial(binary->right, *var, *context);
                    
                    if (p1_res.is_ok() && p2_res.is_ok()) {
                        auto i1_res = algebra::poly_to_integer_poly(p1_res.value());
                        auto i2_res = algebra::poly_to_integer_poly(p2_res.value());
                        
                        if (i1_res.is_ok() && i2_res.is_ok()) {
                            auto gcd_poly = algebra::gcd_integer_poly_with_subresultant(i1_res.value(), i2_res.value()).gcd;
                            
                            if (gcd_poly.degree() > 0 || (gcd_poly.size() == 1 && gcd_poly[0] != BigInt(1))) {
                                // Divisione intera esatta
                                auto q1 = algebra::pseudo_remainder_integer_poly(i1_res.value(), gcd_poly);
                                // Nota: pseudo_remainder restituisce il resto. Qui vogliamo il quoziente.
                                // Uso la logica di divisione intera esatta.
                                
                                auto exact_div = [](const algebra::IntPoly& num, const algebra::IntPoly& den) -> algebra::IntPoly {
                                    if (den.is_zero()) return num;
                                    algebra::IntPoly q;
                                    algebra::IntPoly r = num;
                                    if (num.degree() < den.degree()) return q;
                                    q.resize(num.degree() - den.degree() + 1);
                                    BigInt lc_den = den.leading_coeff();
                                    for (int i = static_cast<int>(num.degree() - den.degree()); i >= 0; --i) {
                                        q[i] = r[i + den.degree()] / lc_den;
                                        for (std::size_t j = 0; j < den.size(); ++j) {
                                            r[i + j] -= q[i] * den[j];
                                        }
                                    }
                                    algebra::normalize_integer_poly(q);
                                    return q;
                                };

                                auto n_poly = exact_div(i1_res.value(), gcd_poly);
                                auto d_poly = exact_div(i2_res.value(), gcd_poly);
                                
                                auto n_expr = algebra::integer_coefficients_to_expr(n_poly, *var, *context);
                                auto d_expr = algebra::integer_coefficients_to_expr(d_poly, *var, *context);
                                
                                if (n_expr.is_ok() && d_expr.is_ok()) {
                                    if (!d_poly.is_zero() && d_poly.degree() == 0 && d_poly[0] == BigInt(1)) return n_expr;
                                    return ok(arena.make<Binary>(BinaryOp::Div, n_expr.value(), d_expr.value()));
                                }
                            }
                        }
                    }
                }
            }
        }

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

        if (const auto* sum = expr_cast<Sum>(expr); sum != nullptr) {
            // Expansion of sum inside log is not canonical.
            // But log(a*b) -> log(a) + log(b) is expansion and IS canonical.
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
            call->args.size() == 1U) {
            if (call->func_id == BuiltinOp::Tan) {
                return ok(arena.make<Binary>(
                    BinaryOp::Div,
                    arena.make<FuncCall>("sin", std::vector<ExprPtr>{call->args.front()}),
                    arena.make<FuncCall>("cos", std::vector<ExprPtr>{call->args.front()})));
            }
        }

        if (const auto* call = expr_cast<FuncCall>(expr); call != nullptr &&
            call->func_id == BuiltinOp::Sqrt && call->args.size() == 1U) {
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
