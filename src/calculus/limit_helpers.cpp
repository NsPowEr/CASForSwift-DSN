#include "cas/calculus.hpp"
#include "calculus_internal.hpp"
#include "cas/rational.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus {

[[nodiscard]] static CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

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

[[nodiscard]] std::optional<int> sign_of_variable_power_at_infinity(ExprPtr expr, const Symbol& var, int point_sign) {
    if (point_sign == 0) return std::nullopt;
    if (is_same_symbol(expr, var)) return point_sign;
    const auto* power = expr_cast<Binary>(expr);
    if (power == nullptr || power->op != BinaryOp::Pow || !is_same_symbol(power->left, var)) return std::nullopt;
    auto exponent = rational_from_expr(power->right);
    if (!exponent.has_value() || !exponent->is_integer() || exponent->numerator().is_negative()) return std::nullopt;
    if (point_sign > 0) return 1;
    return (exponent->numerator() % BigInt(2)).is_zero() ? 1 : -1;
}

std::optional<QuotientView> extract_quotient_view(ExprPtr expr, AstArena& arena) {
    if (!expr) return std::nullopt;

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Div) {
            return QuotientView{.numerator = binary->left, .denominator = binary->right};
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
        
        return QuotientView{.numerator = expr, .denominator = limit_make_integer(arena, 1)};
    }

    std::vector<ExprPtr> numerator_factors;
    std::vector<ExprPtr> denominator_factors;
    for (ExprPtr factor : product->factors) {
        auto factor_quotient = extract_quotient_view(factor, arena);
        if (factor_quotient.has_value()) {
            numerator_factors.push_back(factor_quotient->numerator);
            denominator_factors.push_back(factor_quotient->denominator);
        } else {
            numerator_factors.push_back(factor);
        }
    }

    return QuotientView{
        .numerator = make_product(arena, std::move(numerator_factors)),
        .denominator = make_product(arena, std::move(denominator_factors)),
    };
}

Result<ExprPtr> try_infinite_limit(ExprPtr expr, const Symbol& var, ExprPtr point, AstArena& arena) {
    if (!depends_on(expr, var)) return ok(expr);
    
    const bool is_pos_inf = !expr_is<Unary>(point);

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if (call->func_id == BuiltinOp::Exp) {
            if (is_pos_inf) return ok(arena.make<Constant>(MathConstant::Infinity));
            return ok(limit_make_integer(arena, 0));
        }
        if (call->func_id == BuiltinOp::Ln) {
            auto inner = try_infinite_limit(call->args[0], var, point, arena);
            if (inner.is_ok() && limit_is_infinity(inner.value())) {
                if (expr_is<Unary>(inner.value())) {
                    return fail<ExprPtr>(make_error(CASErrorKind::Undefined, "ln(-inf) is undefined"));
                }
                return ok(arena.make<Constant>(MathConstant::Infinity));
            }
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (auto t : sum->terms) {
            auto lim_t = try_infinite_limit(t, var, point, arena);
            if (lim_t.is_ok() && limit_is_infinity(lim_t.value())) {
                return lim_t;
            }
        }
    }

    if (const auto* sym = expr_cast<Symbol>(expr)) {
        if (sym->name == var.name) return ok(point);
    }

    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Pow) {
            auto base_lim = try_infinite_limit(bin->left, var, point, arena);
            auto exp_lim = try_infinite_limit(bin->right, var, point, arena);
            if (base_lim.is_ok() && limit_is_infinity(base_lim.value()) && !expr_is<Unary>(base_lim.value())) {
                if (exp_lim.is_ok()) {
                    if (limit_is_infinity(exp_lim.value())) {
                        if (!expr_is<Unary>(exp_lim.value())) return ok(base_lim.value()); // inf^inf = inf
                        return ok(limit_make_integer(arena, 0)); // inf^-inf = 0
                    }
                    auto rat = rational_from_expr(exp_lim.value());
                    if (rat.has_value()) {
                        if (rat->numerator() > BigInt(0)) return ok(base_lim.value());
                        if (rat->numerator() < BigInt(0)) return ok(limit_make_integer(arena, 0));
                    }
                }
            }
        }
    }
    
    return fail<ExprPtr>(CASError{CASErrorKind::Unimplemented, "Infinite limit unimplemented for this form", std::nullopt});
}

} // namespace cas::calculus
