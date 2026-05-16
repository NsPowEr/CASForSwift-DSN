#include "calculus_internal.hpp"
#include "cas/rational.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus {

[[nodiscard]] bool is_same_symbol(ExprPtr expr, const Symbol& var);
[[nodiscard]] std::optional<Rational> rational_from_expr(ExprPtr expr);

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] bool is_logarithmic_in_var(ExprPtr expr, const Symbol& var) {
    const auto* call = expr_cast<FuncCall>(expr);
    return call != nullptr &&
        call->func_id == BuiltinOp::Ln &&
        call->args.size() == 1U &&
        depends_on(call->args.front(), var);
}

[[nodiscard]] bool is_positive_power_growth(ExprPtr expr, const Symbol& var) {
    if (is_same_symbol(expr, var)) return true;
    if (const auto* power = expr_cast<Binary>(expr); power != nullptr && power->op == BinaryOp::Pow && is_same_symbol(power->left, var)) {
        auto exponent = rational_from_expr(power->right);
        return exponent.has_value() && exponent->is_integer() && !exponent->numerator().is_negative() && !exponent->numerator().is_zero();
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr factor : product->factors) {
            if (is_positive_power_growth(factor, var)) return true;
        }
    }
    return false;
}

[[nodiscard]] bool is_reciprocal_positive_power_growth(ExprPtr expr, const Symbol& var) {
    const auto* power = expr_cast<Binary>(expr);
    if (power == nullptr || power->op != BinaryOp::Pow) {
        return false;
    }
    auto exponent = rational_from_expr(power->right);
    return exponent.has_value() &&
        exponent->is_integer() &&
        exponent->numerator().is_negative() &&
        is_positive_power_growth(power->left, var);
}

[[nodiscard]] bool is_reciprocal_logarithmic_growth(ExprPtr expr, const Symbol& var) {
    const auto* power = expr_cast<Binary>(expr);
    if (power == nullptr || power->op != BinaryOp::Pow) {
        return false;
    }
    auto exponent = rational_from_expr(power->right);
    return exponent.has_value() &&
        exponent->is_integer() &&
        exponent->numerator().is_negative() &&
        is_logarithmic_in_var(power->left, var);
}

[[nodiscard]] bool is_exponential_limit(ExprPtr expr, const Symbol& var, ExprPtr point, AstArena& arena, bool to_infinity) {
    const auto* call = expr_cast<FuncCall>(expr);
    if (call == nullptr || call->func_id != BuiltinOp::Exp || call->args.size() != 1U) {
        return false;
    }
    auto inner = try_infinite_limit(call->args.front(), var, point, arena);
    if (!inner.is_ok() || !limit_is_infinity(inner.value())) {
        return false;
    }
    const bool inner_negative_infinity = expr_is<Unary>(inner.value());
    return to_infinity ? !inner_negative_infinity : inner_negative_infinity;
}

void collect_multiplicative_factors(ExprPtr expr, std::vector<ExprPtr>& factors) {
    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr factor : product->factors) {
            collect_multiplicative_factors(factor, factors);
        }
        return;
    }
    if (const auto* binary = expr_cast<Binary>(expr); binary != nullptr && binary->op == BinaryOp::Mul) {
        collect_multiplicative_factors(binary->left, factors);
        collect_multiplicative_factors(binary->right, factors);
        return;
    }
    factors.push_back(expr);
}

[[nodiscard]] std::optional<Result<ExprPtr>> try_multiplicative_growth_limit(
    const std::vector<ExprPtr>& factors,
    const Symbol& var,
    ExprPtr point,
    AstArena& arena) {
    const bool at_positive_infinity = !expr_is<Unary>(point);
    bool has_log = false;
    bool has_exp_to_infinity = false;
    bool has_exp_to_zero = false;
    bool has_positive_power = false;
    bool has_reciprocal_power = false;
    bool has_reciprocal_log = false;
    bool all_other_bounded = true;

    for (ExprPtr factor : factors) {
        if (is_logarithmic_in_var(factor, var)) {
            has_log = true;
            continue;
        }
        if (is_exponential_limit(factor, var, point, arena, true)) {
            has_exp_to_infinity = true;
            continue;
        }
        if (is_exponential_limit(factor, var, point, arena, false)) {
            has_exp_to_zero = true;
            continue;
        }
        if (is_positive_power_growth(factor, var)) {
            has_positive_power = true;
            continue;
        }
        if (is_reciprocal_positive_power_growth(factor, var)) {
            has_reciprocal_power = true;
            continue;
        }
        if (is_reciprocal_logarithmic_growth(factor, var)) {
            has_reciprocal_log = true;
            continue;
        }
        if (!is_bounded(factor, var)) {
            all_other_bounded = false;
            break;
        }
    }

    if (at_positive_infinity && has_log && has_reciprocal_power && all_other_bounded) {
        return ok(limit_make_integer(arena, 0));
    }
    if (at_positive_infinity && has_exp_to_infinity && has_reciprocal_power && all_other_bounded) {
        return ok(arena.make<Constant>(MathConstant::Infinity));
    }
    if (at_positive_infinity && has_exp_to_zero && has_positive_power && all_other_bounded) {
        return ok(limit_make_integer(arena, 0));
    }
    if (at_positive_infinity && has_positive_power && has_reciprocal_log && all_other_bounded) {
        return ok(arena.make<Constant>(MathConstant::Infinity));
    }
    return std::nullopt;
}

struct GrowthRank {
    int level; // 0: bounded, 1: log of poly, 2: poly, 3+: nested exponential tower
};

// Dynamic asymptotic-growth rank for x -> +infinity.
//   bounded            ->  level 0
//   log(arg) [arg→∞]   ->  max(rank(arg).level - 1, 0)
//   polynomial in var  ->  level 2 (matches the prior static convention)
//   exp(arg) [arg→∞]   ->  rank(arg).level + 1   <-- key fix for Cat 10
// This is the original Gruntz idea: exp(exp(x)) lives strictly above
// exp(x^N) because the rank of x^N is the same as rank of x (poly = 2),
// while rank of exp(x) is 3, so exp(exp(x)) has level 4 > 3.
GrowthRank get_growth_rank(ExprPtr expr, const Symbol& var, AstArena& /*arena*/) {
    auto compute_rank = [&](ExprPtr e, auto& self) -> GrowthRank {
        if (!depends_on(e, var)) return {0};

        if (const auto* unary = expr_cast<Unary>(e)) {
            if (unary->op == UnaryOp::Neg) return self(unary->operand, self);
        }

        if (const auto* call = expr_cast<FuncCall>(e); call != nullptr && call->args.size() == 1U) {
            if (call->func_id == BuiltinOp::Exp) {
                GrowthRank inner = self(call->args.front(), self);
                return {inner.level + 1};
            }
            if (call->func_id == BuiltinOp::Ln || call->func_id == BuiltinOp::Log) {
                GrowthRank inner = self(call->args.front(), self);
                return {inner.level > 1 ? inner.level - 1 : 0};
            }
        }

        if (is_logarithmic_in_var(e, var)) return {1};
        if (is_positive_power_growth(e, var)) return {2};

        if (const auto* prod = expr_cast<Product>(e)) {
             int max_rank = 0;
             for (auto f : prod->factors) {
                 max_rank = std::max(max_rank, self(f, self).level);
             }
             return {max_rank};
        }

        if (const auto* sum = expr_cast<Sum>(e)) {
             int max_rank = 0;
             for (auto t : sum->terms) {
                 max_rank = std::max(max_rank, self(t, self).level);
             }
             return {max_rank};
        }

        // Pow: rank carried by the base when the exponent is constant in var;
        // otherwise treat as exp(exponent * ln(base)).
        if (const auto* bin = expr_cast<Binary>(e); bin != nullptr && bin->op == BinaryOp::Pow) {
            const bool exp_dep = depends_on(bin->right, var);
            const bool base_dep = depends_on(bin->left, var);
            if (!exp_dep) {
                return self(bin->left, self);
            }
            if (!base_dep) {
                // c^f(x): treat like exp(f(x) * ln(c)).
                GrowthRank e_rank = self(bin->right, self);
                return {e_rank.level + 1};
            }
            // Both depend on var: f(x)^g(x) = exp(g(x) * ln(f(x))).
            GrowthRank g_rank = self(bin->right, self);
            GrowthRank ln_f_rank = self(bin->left, self);
            int eff = std::max(g_rank.level + (ln_f_rank.level > 1 ? ln_f_rank.level - 1 : 0), 0) + 1;
            return {eff};
        }

        return {2};
    };

    return compute_rank(expr, compute_rank);
}

} // namespace

Result<ExprPtr> try_infinite_limit(ExprPtr expr, const Symbol& var, ExprPtr point, AstArena& arena) {
    if (!depends_on(expr, var)) return ok(expr);
    
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if (call->func_id == BuiltinOp::Exp) {
            auto inner = try_infinite_limit(call->args[0], var, point, arena);
            if (inner.is_ok() && limit_is_infinity(inner.value())) {
                if (expr_is<Unary>(inner.value())) {
                    return ok(limit_make_integer(arena, 0));
                }
                return ok(arena.make<Constant>(MathConstant::Infinity));
            }
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

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op == UnaryOp::Neg) {
            auto inner = try_infinite_limit(unary->operand, var, point, arena);
            if (inner.is_ok() && limit_is_infinity(inner.value())) {
                if (expr_is<Unary>(inner.value())) return ok(arena.make<Constant>(MathConstant::Infinity));
                return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
            }
        }
    }

    std::vector<ExprPtr> multiplicative_factors;
    collect_multiplicative_factors(expr, multiplicative_factors);
    if (multiplicative_factors.size() > 1U) {
        auto growth_limit = try_multiplicative_growth_limit(multiplicative_factors, var, point, arena);
        if (growth_limit.has_value()) {
            return growth_limit.value();
        }
    }

    if (const auto* prod = expr_cast<Product>(expr)) {
        bool has_pos_inf = false;
        bool has_neg_inf = false;
        bool has_zero = false;
        int sign = 1;
        for (auto f : prod->factors) {
            auto lim_f = try_infinite_limit(f, var, point, arena);
            if (lim_f.is_error()) return lim_f;
            
            if (limit_is_infinity(lim_f.value())) {
                if (expr_is<Unary>(lim_f.value())) { has_neg_inf = true; sign = -sign; }
                else { has_pos_inf = true; }
            } else if (auto rat = rational_from_expr(lim_f.value())) {
                if (rat->numerator().is_negative()) sign = -sign;
                else if (rat->numerator().is_zero()) has_zero = true;
            } else {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Finite non-rational limit in Product"));
            }
        }
        if (has_zero && (has_pos_inf || has_neg_inf)) {
            return fail<ExprPtr>(make_error(
                CASErrorKind::Unimplemented,
                "0*inf indeterminate form in try_infinite_limit"));
        }
        if (has_zero) return ok(limit_make_integer(arena, 0));
        if (has_pos_inf || has_neg_inf) {
            if (sign < 0) return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
            return ok(arena.make<Constant>(MathConstant::Infinity));
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        int max_level = -1;
        ExprPtr max_lim = nullptr;
        bool indeterminate = false;

        for (auto t : sum->terms) {
            auto lim_t = try_infinite_limit(t, var, point, arena);
            if (lim_t.is_error()) return lim_t;
            
            if (limit_is_infinity(lim_t.value())) {
                int level = get_growth_rank(t, var, arena).level;
                if (level > max_level) {
                    max_level = level;
                    max_lim = lim_t.value();
                    indeterminate = false;
                } else if (level == max_level) {
                    if (expr_is<Unary>(lim_t.value()) != expr_is<Unary>(max_lim)) {
                        indeterminate = true;
                    }
                }
            }
        }
        if (indeterminate) {
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "inf - inf indeterminate form in try_infinite_limit (Sum)"));
        }
        if (max_lim) return ok(max_lim);
    }

    if (const auto* sym = expr_cast<Symbol>(expr)) {
        if (sym->name == var.name) return ok(point);
    }

    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Add || bin->op == BinaryOp::Sub) {
            auto left_limit = try_infinite_limit(bin->left, var, point, arena);
            auto right_limit = try_infinite_limit(bin->right, var, point, arena);
            
            if (left_limit.is_ok() && limit_is_infinity(left_limit.value())) {
                if (right_limit.is_ok() && limit_is_infinity(right_limit.value())) {
                    const bool left_neg = expr_is<Unary>(left_limit.value());
                    const bool right_neg = expr_is<Unary>(right_limit.value());
                    const bool is_sub = (bin->op == BinaryOp::Sub);
                    const bool right_effectively_neg = (right_neg ^ is_sub);

                    int left_level = get_growth_rank(bin->left, var, arena).level;
                    int right_level = get_growth_rank(bin->right, var, arena).level;

                    if (left_level > right_level) return left_limit;
                    if (right_level > left_level) {
                        if (right_effectively_neg) return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
                        return ok(arena.make<Constant>(MathConstant::Infinity));
                    }

                    if (left_neg == right_effectively_neg) {
                        return left_limit;
                    }
                    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "inf - inf indeterminate form in try_infinite_limit (Binary)"));
                } else {
                    return left_limit;
                }
            } else if (right_limit.is_ok() && limit_is_infinity(right_limit.value())) {
                if (bin->op == BinaryOp::Sub) {
                    if (expr_is<Unary>(right_limit.value())) return ok(arena.make<Constant>(MathConstant::Infinity));
                    return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
                }
                return right_limit;
            }
        }

        if (bin->op == BinaryOp::Mul) {
            auto left_limit = try_infinite_limit(bin->left, var, point, arena);
            auto right_limit = try_infinite_limit(bin->right, var, point, arena);
            if (left_limit.is_ok() && right_limit.is_ok()) {
                const bool left_zero = limit_is_zero(left_limit.value());
                const bool right_zero = limit_is_zero(right_limit.value());
                const bool left_inf = limit_is_infinity(left_limit.value());
                const bool right_inf = limit_is_infinity(right_limit.value());
                if (left_zero && right_zero) return ok(limit_make_integer(arena, 0));
                if ((left_inf || right_inf) && !left_zero && !right_zero) {
                    int sign = 1;
                    if ((left_inf && expr_is<Unary>(left_limit.value())) ||
                        (right_inf && expr_is<Unary>(right_limit.value()))) {
                        sign = -sign;
                    }
                    if (sign < 0) {
                        return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
                    }
                    return ok(arena.make<Constant>(MathConstant::Infinity));
                }
            }
        }
        if (bin->op == BinaryOp::Div) {
            // Handle exp(A) / exp(B) -> exp(A - B)
            const auto* num_call = expr_cast<FuncCall>(bin->left);
            const auto* den_call = expr_cast<FuncCall>(bin->right);
            if (num_call && den_call && num_call->func_id == BuiltinOp::Exp && den_call->func_id == BuiltinOp::Exp) {
                ExprPtr diff = arena.make<Binary>(BinaryOp::Sub, num_call->args[0], den_call->args[0]);
                auto diff_limit = try_infinite_limit(diff, var, point, arena);
                if (diff_limit.is_ok() && limit_is_infinity(diff_limit.value())) {
                    if (expr_is<Unary>(diff_limit.value())) return ok(limit_make_integer(arena, 0));
                    return ok(arena.make<Constant>(MathConstant::Infinity));
                }
            }

            auto numerator_limit = try_infinite_limit(bin->left, var, point, arena);
            auto denominator_limit = try_infinite_limit(bin->right, var, point, arena);
            if (denominator_limit.is_ok() && limit_is_infinity(denominator_limit.value())) {
                if (numerator_limit.is_ok() &&
                    limit_is_infinity(numerator_limit.value()) &&
                    is_exponential_limit(bin->left, var, point, arena, true) &&
                    is_positive_power_growth(bin->right, var)) {
                    return ok(arena.make<Constant>(MathConstant::Infinity));
                }
                if (numerator_limit.is_ok() &&
                    limit_is_infinity(numerator_limit.value()) &&
                    is_positive_power_growth(bin->left, var) &&
                    is_exponential_limit(bin->right, var, point, arena, true)) {
                    return ok(limit_make_integer(arena, 0));
                }
                if (numerator_limit.is_ok() && !limit_is_infinity(numerator_limit.value())) {
                    return ok(limit_make_integer(arena, 0));
                }
                if (numerator_limit.is_ok() &&
                    limit_is_infinity(numerator_limit.value()) &&
                    is_logarithmic_in_var(bin->left, var) &&
                    is_positive_power_growth(bin->right, var)) {
                    return ok(limit_make_integer(arena, 0));
                }
            }
            if (denominator_limit.is_ok() && limit_is_zero(denominator_limit.value()) &&
                is_exponential_limit(bin->right, var, point, arena, false) &&
                numerator_limit.is_ok() && limit_is_infinity(numerator_limit.value())) {
                if (expr_is<Unary>(numerator_limit.value())) {
                    return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
                }
                return ok(arena.make<Constant>(MathConstant::Infinity));
            }
        }
        if (bin->op == BinaryOp::Pow) {
            auto base_lim = try_infinite_limit(bin->left, var, point, arena);
            auto exp_lim = try_infinite_limit(bin->right, var, point, arena);
            if (base_lim.is_ok() && limit_is_zero(base_lim.value())) {
                auto rat = exp_lim.is_ok() ? rational_from_expr(exp_lim.value()) : std::optional<Rational>{};
                if (rat.has_value() && rat->is_integer()) {
                    if (rat->numerator() > BigInt(0)) return ok(limit_make_integer(arena, 0));
                    if (rat->numerator() < BigInt(0)) return ok(arena.make<Constant>(MathConstant::Infinity));
                    return ok(limit_make_integer(arena, 1));
                }
            }
            if (base_lim.is_ok() && limit_is_infinity(base_lim.value())) {
                if (exp_lim.is_ok()) {
                    if (limit_is_infinity(exp_lim.value())) {
                        if (!expr_is<Unary>(exp_lim.value())) return ok(arena.make<Constant>(MathConstant::Infinity));
                        return ok(limit_make_integer(arena, 0));
                    }
                    auto rat = rational_from_expr(exp_lim.value());
                    if (rat.has_value() && rat->is_integer()) {
                        if (rat->numerator() > BigInt(0)) {
                            const bool negative_base = expr_is<Unary>(base_lim.value());
                            const bool odd_integer = !(rat->numerator() % BigInt(2)).is_zero();
                            if (negative_base && odd_integer) {
                                return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
                            }
                            return ok(arena.make<Constant>(MathConstant::Infinity));
                        }
                        if (rat->numerator() < BigInt(0)) return ok(limit_make_integer(arena, 0));
                        return ok(limit_make_integer(arena, 1));
                    }
                }
            }
        }
    }
    
    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Infinite limit unimplemented for this form"));
}

} // namespace cas::calculus
