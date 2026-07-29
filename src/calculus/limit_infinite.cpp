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
    // Ln and Log share one growth class: log_b(x) = ln(x)/ln(b) differs only by a
    // constant factor, so both are logarithmic for asymptotic comparison.
    return call != nullptr &&
        (call->func_id == BuiltinOp::Ln || call->func_id == BuiltinOp::Log) &&
        call->args.size() == 1U &&
        depends_on(call->args.front(), var);
}

// A factor that grows like x^p with p > 0 (rational) as x → +∞. The asymptotic
// growth class depends only on the *sign* of the exponent, not its integrality:
// x^(1/2) dominates every logarithm exactly as x^1 does. `sqrt(u)` is the
// canonical form of u^(1/2), so it is a positive power of u^(1/2)'s base.
[[nodiscard]] bool is_positive_power_growth(ExprPtr expr, const Symbol& var) {
    if (is_same_symbol(expr, var)) return true;
    if (const auto* power = expr_cast<Binary>(expr); power != nullptr && power->op == BinaryOp::Pow) {
        auto exponent = rational_from_expr(power->right);
        return exponent.has_value()
            && !exponent->numerator().is_negative() && !exponent->numerator().is_zero()
            && is_positive_power_growth(power->left, var);
    }
    if (const auto* call = expr_cast<FuncCall>(expr);
        call != nullptr && call->func_id == BuiltinOp::Sqrt && call->args.size() == 1U) {
        return is_positive_power_growth(call->args.front(), var);
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr factor : product->factors) {
            if (is_positive_power_growth(factor, var)) return true;
        }
    }
    return false;
}

// A factor that decays like x^(-p) with p > 0 (rational) as x → +∞, i.e. the
// reciprocal of a positive power. Covers 1/sqrt(x) = Pow(sqrt(x), -1) and
// x^(-1/2) = Pow(x, -1/2) alike — only the negative exponent sign matters.
[[nodiscard]] bool is_reciprocal_positive_power_growth(ExprPtr expr, const Symbol& var) {
    const auto* power = expr_cast<Binary>(expr);
    if (power == nullptr || power->op != BinaryOp::Pow) {
        return false;
    }
    auto exponent = rational_from_expr(power->right);
    return exponent.has_value() &&
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
        exponent->numerator().is_negative() &&
        is_logarithmic_in_var(power->left, var);
}

[[nodiscard]] bool is_exponential_limit(ExprPtr expr, const Symbol& var, ExprPtr point, symbolic::CASContext& ctx, bool to_infinity) {
    const auto* call = expr_cast<FuncCall>(expr);
    if (call == nullptr || call->func_id != BuiltinOp::Exp || call->args.size() != 1U) {
        return false;
    }
    auto inner = try_infinite_limit(call->args.front(), var, point, ctx);
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

// The additive term of `e` with the fastest growth as x → +∞ (Gruntz §3.5,
// via the shared dynamic `compare_growth`). For a non-additive expression it is
// `e` itself. Used to take the leading-term ratio of a quotient of sums.
[[nodiscard]] ExprPtr dominant_growth_term(
    ExprPtr e, const Symbol& var, symbolic::CASContext& ctx) {
    std::vector<ExprPtr> terms;
    if (const auto* sum = expr_cast<Sum>(e)) {
        terms = sum->terms;
    } else if (const auto* bin = expr_cast<Binary>(e);
               bin != nullptr && (bin->op == BinaryOp::Add || bin->op == BinaryOp::Sub)) {
        terms.push_back(bin->left);
        terms.push_back(bin->op == BinaryOp::Sub
            ? static_cast<ExprPtr>(ctx.arena().make<Unary>(UnaryOp::Neg, bin->right))
            : bin->right);
    } else {
        return e;
    }
    if (terms.empty()) return e;
    ExprPtr best = terms.front();
    for (std::size_t i = 1; i < terms.size(); ++i) {
        if (compare_growth(terms[i], best, var, ctx) > 0) best = terms[i];
    }
    return best;
}

[[nodiscard]] std::optional<Result<ExprPtr>> try_multiplicative_growth_limit(
    const std::vector<ExprPtr>& factors,
    const Symbol& var,
    ExprPtr point,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
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
        if (is_exponential_limit(factor, var, point, ctx, true)) {
            has_exp_to_infinity = true;
            continue;
        }
        if (is_exponential_limit(factor, var, point, ctx, false)) {
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

// The legacy static `GrowthRank`/`get_growth_rank` table has been removed.
// The Sum / Binary dominance decisions now use `compare_growth` from
// `limit_mrv.cpp`, which is the single dynamic Gruntz §3.5 comparison shared
// with the MRV solver (CLAUDE.md Cat 10 closure for F5.2 / B3).

} // namespace

Result<ExprPtr> try_infinite_limit(ExprPtr expr, const Symbol& var, ExprPtr point, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    if (!depends_on(expr, var)) return ok(expr);
    
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if (call->func_id == BuiltinOp::Exp) {
            auto inner = try_infinite_limit(call->args[0], var, point, ctx);
            if (inner.is_ok() && limit_is_infinity(inner.value())) {
                if (expr_is<Unary>(inner.value())) {
                    return ok(limit_make_integer(arena, 0));
                }
                return ok(arena.make<Constant>(MathConstant::Infinity));
            }
        }
        if (call->func_id == BuiltinOp::Ln || call->func_id == BuiltinOp::Log) {
            auto inner = try_infinite_limit(call->args[0], var, point, ctx);
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
            auto inner = try_infinite_limit(unary->operand, var, point, ctx);
            if (inner.is_ok() && limit_is_infinity(inner.value())) {
                if (expr_is<Unary>(inner.value())) return ok(arena.make<Constant>(MathConstant::Infinity));
                return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
            }
        }
    }

    // Leading-term reduction for ∞/∞ quotients — handles both Binary Div and the
    // canonical Product·Pow(-1) form via extract_quotient_view. Replace a sum
    // numerator/denominator by its fastest-growing additive term (Gruntz growth
    // comparison) and recurse on the ratio: (x + log x)/(x - log x) → x/x → 1.
    // This runs before the multiplicative/Product handling, which would otherwise
    // mis-report the 0·∞ form as indeterminate. Guarded so at least one side
    // actually reduces ⇒ the ratio of two single terms is structurally simpler,
    // so the recursion terminates.
    if (auto qv = extract_quotient_view(expr, arena); qv.has_value()) {
        ExprPtr dom_num = dominant_growth_term(qv->numerator, var, ctx);
        ExprPtr dom_den = dominant_growth_term(qv->denominator, var, ctx);
        if (dom_num != qv->numerator || dom_den != qv->denominator) {
            auto num_lim = try_infinite_limit(qv->numerator, var, point, ctx);
            auto den_lim = try_infinite_limit(qv->denominator, var, point, ctx);
            if (num_lim.is_ok() && den_lim.is_ok() &&
                limit_is_infinity(num_lim.value()) && limit_is_infinity(den_lim.value())) {
                ExprPtr ratio = arena.make<Binary>(BinaryOp::Div, dom_num, dom_den);
                auto rs = ctx.simplify(ratio);
                ExprPtr rexpr = rs.is_ok() ? rs.value() : ratio;
                auto rl = try_infinite_limit(rexpr, var, point, ctx);
                if (rl.is_ok()) return rl;
            }
        }
    }

    std::vector<ExprPtr> multiplicative_factors;
    collect_multiplicative_factors(expr, multiplicative_factors);
    if (multiplicative_factors.size() > 1U) {
        auto growth_limit = try_multiplicative_growth_limit(multiplicative_factors, var, point, ctx);
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
            auto lim_f = try_infinite_limit(f, var, point, ctx);
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
        // Gruntz §3.5: pick the term whose growth rank dominates the others.
        // `compare_growth` is the dynamic recursive comparison shared with the
        // MRV solver, so nested exp/log towers are handled at arbitrary depth
        // without a static enum table (CLAUDE.md Cat 10).
        ExprPtr max_term = nullptr;
        ExprPtr max_lim = nullptr;
        bool indeterminate = false;

        for (auto t : sum->terms) {
            auto lim_t = try_infinite_limit(t, var, point, ctx);
            if (lim_t.is_error()) return lim_t;

            if (limit_is_infinity(lim_t.value())) {
                if (max_term == nullptr) {
                    max_term = t;
                    max_lim = lim_t.value();
                    indeterminate = false;
                    continue;
                }
                int cmp = compare_growth(t, max_term, var, ctx);
                if (cmp > 0) {
                    max_term = t;
                    max_lim = lim_t.value();
                    indeterminate = false;
                } else if (cmp == 0) {
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
            auto left_limit = try_infinite_limit(bin->left, var, point, ctx);
            auto right_limit = try_infinite_limit(bin->right, var, point, ctx);
            
            if (left_limit.is_ok() && limit_is_infinity(left_limit.value())) {
                if (right_limit.is_ok() && limit_is_infinity(right_limit.value())) {
                    const bool left_neg = expr_is<Unary>(left_limit.value());
                    const bool right_neg = expr_is<Unary>(right_limit.value());
                    const bool is_sub = (bin->op == BinaryOp::Sub);
                    const bool right_effectively_neg = (right_neg ^ is_sub);

                    // Use the dynamic Gruntz comparison (shared with the MRV
                    // engine) instead of the legacy static rank levels.
                    int cmp = compare_growth(bin->left, bin->right, var, ctx);
                    if (cmp > 0) return left_limit;
                    if (cmp < 0) {
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
            auto left_limit = try_infinite_limit(bin->left, var, point, ctx);
            auto right_limit = try_infinite_limit(bin->right, var, point, ctx);
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
                auto diff_limit = try_infinite_limit(diff, var, point, ctx);
                if (diff_limit.is_ok() && limit_is_infinity(diff_limit.value())) {
                    if (expr_is<Unary>(diff_limit.value())) return ok(limit_make_integer(arena, 0));
                    return ok(arena.make<Constant>(MathConstant::Infinity));
                }
            }

            auto numerator_limit = try_infinite_limit(bin->left, var, point, ctx);
            auto denominator_limit = try_infinite_limit(bin->right, var, point, ctx);
            if (denominator_limit.is_ok() && limit_is_infinity(denominator_limit.value())) {
                if (numerator_limit.is_ok() &&
                    limit_is_infinity(numerator_limit.value()) &&
                    is_exponential_limit(bin->left, var, point, ctx, true) &&
                    is_positive_power_growth(bin->right, var)) {
                    return ok(arena.make<Constant>(MathConstant::Infinity));
                }
                if (numerator_limit.is_ok() &&
                    limit_is_infinity(numerator_limit.value()) &&
                    is_positive_power_growth(bin->left, var) &&
                    is_exponential_limit(bin->right, var, point, ctx, true)) {
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
                is_exponential_limit(bin->right, var, point, ctx, false) &&
                numerator_limit.is_ok() && limit_is_infinity(numerator_limit.value())) {
                if (expr_is<Unary>(numerator_limit.value())) {
                    return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
                }
                return ok(arena.make<Constant>(MathConstant::Infinity));
            }
        }
        if (bin->op == BinaryOp::Pow) {
            auto base_lim = try_infinite_limit(bin->left, var, point, ctx);
            auto exp_lim = try_infinite_limit(bin->right, var, point, ctx);
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
