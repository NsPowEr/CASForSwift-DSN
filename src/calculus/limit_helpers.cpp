#include "cas/calculus.hpp"
#include "calculus_internal.hpp"
#include "cas/rational.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"

#include <algorithm>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus {

[[nodiscard]] ExprPtr limit_make_integer(AstArena& arena, long long value) {
    return arena.make<IntegerLit>(BigInt(value));
}

[[nodiscard]] ExprPtr limit_make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs) {
    return arena.make<Binary>(op, lhs, rhs);
}

[[nodiscard]] bool limit_is_zero(ExprPtr expr) {
    if (!expr) return false;
    if (const auto* i = expr_cast<IntegerLit>(expr)) return i->value.is_zero();
    if (const auto* r = expr_cast<RationalLit>(expr)) return r->numerator.is_zero();
    return false;
}

[[nodiscard]] bool limit_is_one(ExprPtr expr) {
    if (!expr) return false;
    if (const auto* i = expr_cast<IntegerLit>(expr)) return i->value == BigInt(1);
    if (const auto* r = expr_cast<RationalLit>(expr)) return r->numerator == r->denominator;
    return false;
}

[[nodiscard]] bool limit_is_infinity(ExprPtr expr) {
    if (!expr) return false;
    if (const auto* c = expr_cast<Constant>(expr)) return c->value == MathConstant::Infinity;
    if (const auto* u = expr_cast<Unary>(expr)) return u->op == UnaryOp::Neg && limit_is_infinity(u->operand);
    return false;
}

[[nodiscard]] bool depends_on(ExprPtr expr, const Symbol& var) {
    if (!expr) return false;
    if (const auto* sym = expr_cast<Symbol>(expr)) return sym->name == var.name;
    
    bool dep = false;
    visit_expr(expr, [&](const auto& node) {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, Unary>) { if (depends_on(node.operand, var)) dep = true; }
        else if constexpr (std::is_same_v<T, Binary>) { if (depends_on(node.left, var) || depends_on(node.right, var)) dep = true; }
        else if constexpr (std::is_same_v<T, FuncCall>) { for (auto a : node.args) if (depends_on(a, var)) dep = true; }
        else if constexpr (std::is_same_v<T, Sum>) { for (auto t : node.terms) if (depends_on(t, var)) dep = true; }
        else if constexpr (std::is_same_v<T, Product>) { for (auto f : node.factors) if (depends_on(f, var)) dep = true; }
    });
    return dep;
}

[[nodiscard]] bool is_bounded(ExprPtr expr, const Symbol& var) {
    if (!depends_on(expr, var)) return true;
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if (call->func_id == BuiltinOp::Sin || call->func_id == BuiltinOp::Cos) return true;
        if (call->func_id == BuiltinOp::Atan) return true;
    }
    return false;
}

[[nodiscard]] unsigned int transcendental_tower_depth(ExprPtr expr, const Symbol& var) {
    if (!expr) return 0U;
    if (!depends_on(expr, var)) return 0U;
    if (expr_is<Symbol>(expr)) return 1U;
    if (const auto* unary = expr_cast<Unary>(expr)) {
        return transcendental_tower_depth(unary->operand, var);
    }
    if (const auto* bin = expr_cast<Binary>(expr)) {
        const unsigned int lh = transcendental_tower_depth(bin->left, var);
        const unsigned int rh = transcendental_tower_depth(bin->right, var);
        return std::max(lh, rh);
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        unsigned int h = 0U;
        for (auto t : sum->terms) h = std::max(h, transcendental_tower_depth(t, var));
        return h;
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        unsigned int h = 0U;
        for (auto f : prod->factors) h = std::max(h, transcendental_tower_depth(f, var));
        return h;
    }
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        unsigned int inner = 0U;
        for (auto a : call->args) inner = std::max(inner, transcendental_tower_depth(a, var));
        // Gruntz tower escalation: each exp(.) wrapping a var-dependent
        // sub-expression adds one comparability level.
        switch (call->func_id) {
        case BuiltinOp::Exp:
            return inner + 1U;
        case BuiltinOp::Ln:
        case BuiltinOp::Log:
        case BuiltinOp::Log10:
            // log peels one level off (asymptotically slower than its argument)
            // but the structural depth of *evaluating* this call still counts.
            return std::max(inner, 1U);
        default:
            return std::max(inner, 1U);
        }
    }
    return 1U;
}

[[nodiscard]] static ExprPtr make_power(AstArena& arena, ExprPtr base, ExprPtr exponent) {
    return arena.make<Binary>(BinaryOp::Pow, base, exponent);
}

[[nodiscard]] static ExprPtr make_sum(AstArena& arena, std::vector<ExprPtr> terms) {
    if (terms.empty()) return limit_make_integer(arena, 0);
    if (terms.size() == 1U) return terms.front();
    return arena.make<Sum>(std::move(terms));
}

[[nodiscard]] static ExprPtr make_product(AstArena& arena, std::vector<ExprPtr> factors) {
    std::vector<ExprPtr> filtered;
    filtered.reserve(factors.size());
    for (auto f : factors) {
        if (const auto* il = expr_cast<IntegerLit>(f)) {
            if (il->value == BigInt(1)) continue;
        }
        filtered.push_back(f);
    }
    if (filtered.empty()) return limit_make_integer(arena, 1);
    if (filtered.size() == 1U) return filtered.front();
    return arena.make<Product>(std::move(filtered));
}

[[nodiscard]] bool is_same_symbol(ExprPtr expr, const Symbol& var);
[[nodiscard]] std::optional<Rational> rational_from_expr(ExprPtr expr);

[[nodiscard]] bool is_same_symbol(ExprPtr expr, const Symbol& var) {
    const auto* symbol = expr_cast<Symbol>(expr);
    return symbol != nullptr && symbol->name == var.name;
}

[[nodiscard]] std::optional<Rational> rational_from_expr(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return Rational(integer->value);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return Rational(rational->numerator, rational->denominator);
    }
    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op != UnaryOp::Neg) return std::nullopt;
        auto operand = rational_from_expr(unary->operand);
        if (!operand.has_value()) return std::nullopt;
        return -operand.value();
    }
    return std::nullopt;
}

std::optional<QuotientView> extract_quotient_view(ExprPtr expr, AstArena& arena) {
    if (!expr) return std::nullopt;

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if (call->func_id == BuiltinOp::Exp) {
            auto inner = call->args[0];
            if (const auto* u = expr_cast<Unary>(inner); u && u->op == UnaryOp::Neg) {
                return QuotientView{
                    .numerator = limit_make_integer(arena, 1),
                    .denominator = arena.make<FuncCall>("exp", std::vector<ExprPtr>{u->operand})
                };
            }
        }
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Div) {
            return QuotientView{.numerator = binary->left, .denominator = binary->right};
        }
        if (binary->op == BinaryOp::Mul) {
            auto lhs_q = extract_quotient_view(binary->left, arena);
            auto rhs_q = extract_quotient_view(binary->right, arena);
            if (lhs_q.has_value() || rhs_q.has_value()) {
                return QuotientView{
                    .numerator = make_product(arena, {lhs_q.has_value() ? lhs_q->numerator : binary->left, 
                                                      rhs_q.has_value() ? rhs_q->numerator : binary->right}),
                    .denominator = make_product(arena, {lhs_q.has_value() ? lhs_q->denominator : limit_make_integer(arena, 1),
                                                        rhs_q.has_value() ? rhs_q->denominator : limit_make_integer(arena, 1)})
                };
            }
        }
        if (binary->op == BinaryOp::Pow) {
            auto exponent = rational_from_expr(binary->right);
            if (exponent.has_value() && exponent->is_integer() && exponent->numerator().is_negative()) {
                const BigInt positive_power = -exponent->numerator();
                ExprPtr denominator = positive_power == BigInt(1)
                    ? binary->left
                    : make_power(arena, binary->left, arena.make<IntegerLit>(positive_power));
                return QuotientView{.numerator = limit_make_integer(arena, 1), .denominator = denominator};
            }
        }
    }

    const auto* product = expr_cast<Product>(expr);
    if (product == nullptr) {
        if (const auto* unary = expr_cast<Unary>(expr)) {
            if (unary->op == UnaryOp::Neg) {
                auto inner = extract_quotient_view(unary->operand, arena);
                if (inner.has_value()) {
                    return QuotientView{
                        .numerator = arena.make<Unary>(UnaryOp::Neg, inner->numerator),
                        .denominator = inner->denominator,
                    };
                }
            }
        }

        if (const auto* sum = expr_cast<Sum>(expr)) {
            std::vector<ExprPtr> numerator_terms;
            ExprPtr common_denominator{};
            for (ExprPtr term : sum->terms) {
                auto term_quotient = extract_quotient_view(term, arena);
                if (!term_quotient.has_value()) return std::nullopt;
                if (!common_denominator) common_denominator = term_quotient->denominator;
                else if (symbolic::canonical_compare(common_denominator, term_quotient->denominator) != 0) return std::nullopt;
                numerator_terms.push_back(term_quotient->numerator);
            }
            if (common_denominator) {
                return QuotientView{
                    .numerator = make_sum(arena, std::move(numerator_terms)),
                    .denominator = common_denominator,
                };
            }
        }
        
        return std::nullopt;
    }

    std::vector<ExprPtr> numerator_factors;
    std::vector<ExprPtr> denominator_factors;
    bool has_den = false;
    for (ExprPtr factor : product->factors) {
        auto factor_quotient = extract_quotient_view(factor, arena);
        if (factor_quotient.has_value()) {
            numerator_factors.push_back(factor_quotient->numerator);
            denominator_factors.push_back(factor_quotient->denominator);
            has_den = true;
        } else {
            numerator_factors.push_back(factor);
        }
    }

    if (!has_den) return std::nullopt;

    return QuotientView{
        .numerator = make_product(arena, std::move(numerator_factors)),
        .denominator = make_product(arena, std::move(denominator_factors)),
    };
}

} // namespace cas::calculus
