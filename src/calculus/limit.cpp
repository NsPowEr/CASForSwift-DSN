#include "cas/calculus.hpp"
#include "cas/algebra.hpp"

#include "limit_internal.hpp"
#include "calculus_internal.hpp"
#include "cas/error.hpp"

#include <algorithm>
#include <utility>

namespace cas::calculus {
namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

}  // namespace

Result<ExprPtr> LimitEngine::compute(ExprPtr expr, const Symbol& var, ExprPtr point, LimitDirection dir) {
    const unsigned int tower_h = transcendental_tower_depth(expr, var);
    max_depth_budget_ = std::max<unsigned int>(8U, 2U * tower_h + 4U);
    auto together_res = algebra::together(expr, context_);
    auto simplified_expr = context_.simplify(together_res.is_ok() ? together_res.value() : expr);
    if (simplified_expr.is_error()) {
        return simplified_expr;
    }

    auto simplified_point = context_.simplify(point);
    if (simplified_point.is_error()) {
        return simplified_point;
    }

    const bool point_is_pos_inf = limit_is_infinity(simplified_point.value()) &&
                                  !expr_is<Unary>(simplified_point.value());
    const bool point_is_neg_inf = limit_is_infinity(simplified_point.value()) &&
                                  expr_is<Unary>(simplified_point.value());

    if (point_is_pos_inf || point_is_neg_inf) {
        auto infinite_limit = try_infinite_limit(simplified_expr.value(), var, simplified_point.value(), context_);
        if (infinite_limit.is_ok() || infinite_limit.error().kind != CASErrorKind::Unimplemented) {
            return infinite_limit;
        }
    }

    ExprPtr pre_mrv_expr = simplified_expr.value();
    if (auto cancelled = try_cancel_product_pow_inverse(pre_mrv_expr, context_);
        cancelled.has_value()) {
        auto re_simp = context_.simplify(cancelled.value());
        if (re_simp.is_ok()) pre_mrv_expr = re_simp.value();
    }

    if (point_is_pos_inf || point_is_neg_inf) {
        auto q_inf = extract_quotient_view(pre_mrv_expr, arena_);
        if (q_inf.has_value()) {
            LimitDirection inf_dir = point_is_pos_inf
                ? LimitDirection::Right : LimitDirection::Left;
            auto ll = try_log_log_limit(q_inf.value(), var,
                                             simplified_point.value(), inf_dir, 0U);
            if (ll.has_value()) {
                return ll.value();
            }
        }
    }

    if ((point_is_pos_inf || point_is_neg_inf)) {
        if (auto sum_res = try_limit_sum_termwise(
                pre_mrv_expr, var, simplified_point.value(), dir, context_);
            sum_res.has_value()) {
            return *sum_res;
        }
    }

    if (point_is_pos_inf || point_is_neg_inf) {
        auto mrv_res = compute_limit_mrv(pre_mrv_expr, var, simplified_point.value(), context_);
        if (mrv_res.is_ok()) {
            return mrv_res;
        }
        if (mrv_res.error().kind != CASErrorKind::Unimplemented) {
            return fail<ExprPtr>(mrv_res.error());
        }
    }

    if (point_is_pos_inf || point_is_neg_inf) {
        const std::string t_name = var.name + "_inv";
        const Symbol t_var(t_name);
        ExprPtr t_sym = arena_.make<Symbol>(t_name);

        ExprPtr one_over_t = arena_.make<Binary>(BinaryOp::Div, limit_make_integer(arena_, 1), t_sym);
        ExprPtr x_sub = point_is_pos_inf
            ? one_over_t
            : arena_.make<Unary>(UnaryOp::Neg, one_over_t);
        auto expr_sub = context_.substitute(pre_mrv_expr, var, x_sub);
        if (expr_sub.is_error()) return expr_sub;
        
        auto together_t = algebra::together(expr_sub.value(), context_);
        auto expr_t = context_.simplify(together_t.is_ok() ? together_t.value() : expr_sub.value());
        if (expr_t.is_error()) return expr_t;

        ExprPtr zero_pt = limit_make_integer(arena_, 0);
        auto res = compute_recursive(expr_t.value(), t_var, zero_pt, LimitDirection::Right, 0U);
        if (res.is_ok()) return res;

        const auto kind = res.error().kind;
        if (kind != CASErrorKind::Unimplemented &&
            kind != CASErrorKind::Undefined &&
            kind != CASErrorKind::DivisionByZero) {
            return res;
        }
    }

    return compute_recursive(pre_mrv_expr, var, simplified_point.value(), dir, 0U);
}

Result<ExprPtr> LimitEngine::substitute_and_simplify(ExprPtr expr, const Symbol& var, ExprPtr point) {
    auto substituted = context_.substitute(expr, var, point);
    if (substituted.is_error()) {
        return substituted;
    }
    auto res = context_.simplify(substituted.value());
    return res;
}

Result<ExprPtr> LimitEngine::compute_recursive(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    LimitDirection dir,
    unsigned int depth) {
    if (depth >= max_depth_budget_) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::Unimplemented,
            "Il limite richiede piu' iterazioni di quelle supportate"));
    }

    auto direct = substitute_and_simplify(expr, var, point);

    auto quotient = extract_quotient_view(expr, arena_);
    const bool has_quotient = quotient.has_value();

    if (!limit_is_infinity(point)) {
        auto logarithmic_limit = try_logarithmic_root_limit(expr, var, point, dir, arena_);
        if (logarithmic_limit.has_value()) {
            return logarithmic_limit.value();
        }
    }

    if (limit_is_infinity(point)) {
        auto infinite_limit = try_infinite_limit(expr, var, point, context_);
        if (infinite_limit.is_ok() || infinite_limit.error().kind != CASErrorKind::Unimplemented) {
            return infinite_limit;
        }
    }

    if (direct.is_ok() && (!expr_is<Binary>(direct.value()) || direct.value() == expr)) {
        if (!expr_is<Binary>(expr) || (!limit_is_zero(direct.value()) && !limit_is_infinity(direct.value()))) {
            return direct;
        }
    }

    if (const auto* prod = expr_cast<Product>(expr)) {
        bool has_zero = false;
        bool has_infinity = false;
        bool all_bounded = true;
        std::vector<ExprPtr> zero_factors;
        std::vector<ExprPtr> inf_factors;
        std::vector<ExprPtr> non_zero_factors;

        for (auto factor : prod->factors) {
            auto lim_f = compute_recursive(factor, var, point, dir, depth + 1U);
            if (lim_f.is_ok() && limit_is_zero(lim_f.value())) {
                has_zero = true;
                zero_factors.push_back(factor);
            } else if (lim_f.is_ok() && limit_is_infinity(lim_f.value())) {
                has_infinity = true;
                all_bounded = false;
                inf_factors.push_back(factor);
            } else {
                if (!is_bounded(factor, var)) all_bounded = false;
                non_zero_factors.push_back(factor);
            }
        }
        if (has_zero && all_bounded) {
            return ok(limit_make_integer(arena_, 0));
        }
        if (has_zero && has_infinity) {
            auto get_priority = [](ExprPtr f) -> int {
                if (const auto* call = expr_cast<FuncCall>(f)) {
                    if (call->func_id == BuiltinOp::Ln) return 3;
                    if (call->func_id == BuiltinOp::Exp) return 1;
                }
                if (expr_cast<Symbol>(f)) return 2;
                if (const auto* bin = expr_cast<Binary>(f)) {
                    if (bin->op == BinaryOp::Pow) return 2;
                }
                return 0;
            };

            int max_zero_pri = -1;
            for (auto f : zero_factors) max_zero_pri = std::max(max_zero_pri, get_priority(f));
            int max_inf_pri = -1;
            for (auto f : inf_factors) max_inf_pri = std::max(max_inf_pri, get_priority(f));

            ExprPtr num, den_base;
            if (max_zero_pri >= max_inf_pri) {
                num = zero_factors.size() == 1 ? zero_factors[0] : arena_.make<Product>(std::move(zero_factors));
                for (auto factor : non_zero_factors) inf_factors.push_back(factor);
                den_base = inf_factors.size() == 1 ? inf_factors[0] : arena_.make<Product>(std::move(inf_factors));
            } else {
                for (auto factor : non_zero_factors) inf_factors.push_back(factor);
                num = inf_factors.size() == 1 ? inf_factors[0] : arena_.make<Product>(std::move(inf_factors));
                den_base = zero_factors.size() == 1 ? zero_factors[0] : arena_.make<Product>(std::move(zero_factors));
            }

            ExprPtr den = arena_.make<Binary>(BinaryOp::Pow, den_base, limit_make_integer(arena_, -1));
            auto simplified_rewritten = context_.simplify(arena_.make<Binary>(BinaryOp::Div, num, den));
            if (simplified_rewritten.is_error()) return simplified_rewritten;
            
            auto qv = extract_quotient_view(simplified_rewritten.value(), arena_);
            if (qv.has_value()) {
                return compute_quotient_limit(simplified_rewritten.value(), var, point, dir, depth + 1U, qv.value(), direct);
            }
            return compute_recursive(simplified_rewritten.value(), var, point, dir, depth + 1U);
        }
    }

    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Mul) {
            auto lim_l = compute_recursive(bin->left, var, point, dir, depth + 1U);
            auto lim_r = compute_recursive(bin->right, var, point, dir, depth + 1U);

            if (lim_l.is_ok() && lim_r.is_ok()) {
                bool l_zero = limit_is_zero(lim_l.value());
                bool r_zero = limit_is_zero(lim_r.value());
                bool l_inf = limit_is_infinity(lim_l.value());
                bool r_inf = limit_is_infinity(lim_r.value());

                if ((l_zero && is_bounded(bin->right, var)) || (r_zero && is_bounded(bin->left, var))) {
                    return ok(limit_make_integer(arena_, 0));
                }

                if ((l_zero && r_inf) || (r_zero && l_inf)) {
                    auto get_priority = [](ExprPtr f) -> int {
                        if (const auto* call = expr_cast<FuncCall>(f)) {
                            if (call->func_id == BuiltinOp::Ln) return 3;
                            if (call->func_id == BuiltinOp::Exp) return 1;
                        }
                        if (expr_cast<Symbol>(f)) return 2;
                        if (const auto* b = expr_cast<Binary>(f)) {
                            if (b->op == BinaryOp::Pow) return 2;
                        }
                        return 0;
                    };

                    ExprPtr num, den_base;
                    if (get_priority(bin->left) >= get_priority(bin->right)) {
                        num = bin->left;
                        den_base = bin->right;
                    } else {
                        num = bin->right;
                        den_base = bin->left;
                    }

                    ExprPtr den = arena_.make<Binary>(BinaryOp::Pow, den_base, limit_make_integer(arena_, -1));
                    auto simplified_rewritten = context_.simplify(arena_.make<Binary>(BinaryOp::Div, num, den));
                    if (simplified_rewritten.is_error()) return simplified_rewritten;
                    
                    auto qv = extract_quotient_view(simplified_rewritten.value(), arena_);
                    if (qv.has_value()) {
                        return compute_quotient_limit(simplified_rewritten.value(), var, point, dir, depth + 1U, qv.value(), direct);
                    }
                    return compute_recursive(simplified_rewritten.value(), var, point, dir, depth + 1U);
                }
            }
        }
    }

    if (const auto* power = expr_cast<Binary>(expr); power && power->op == BinaryOp::Pow) {
        auto base_lim = compute_recursive(power->left, var, point, dir, depth + 1U);
        auto exp_lim = compute_recursive(power->right, var, point, dir, depth + 1U);

        if (base_lim.is_ok() && exp_lim.is_ok()) {
            const bool base_is_one = limit_is_one(base_lim.value());
            const bool base_is_zero = limit_is_zero(base_lim.value());
            const bool base_is_inf = limit_is_infinity(base_lim.value());
            const bool exp_is_zero = limit_is_zero(exp_lim.value());
            const bool exp_is_inf = limit_is_infinity(exp_lim.value());

            if ((base_is_one && exp_is_inf) || (base_is_zero && exp_is_zero) || (base_is_inf && exp_is_zero)) {
                auto ln_f = arena_.make<FuncCall>("ln", std::vector<ExprPtr>{power->left});
                auto inner_expr = arena_.make<Binary>(BinaryOp::Mul, power->right, ln_f);
                auto inner_lim = compute_recursive(inner_expr, var, point, dir, depth + 1U);
                if (inner_lim.is_ok()) {
                    return ok(arena_.make<FuncCall>("exp", std::vector<ExprPtr>{inner_lim.value()}));
                }
            }
        }
    }

    if (direct.is_error() && !has_quotient) {
        return direct;
    }

    if (limit_is_zero(point) && dir == LimitDirection::Right) {
        if (const auto* power = expr_cast<Binary>(expr);
            power != nullptr &&
            power->op == BinaryOp::Pow &&
            expr_cast<Symbol>(power->left) != nullptr &&
            expr_cast<Symbol>(power->right) != nullptr &&
            expr_ref<Symbol>(power->left).name == var.name &&
            expr_ref<Symbol>(power->right).name == var.name) {
            return ok(limit_make_integer(arena_, 1));
        }
    }

    if (quotient.has_value()) {
        if (auto log_log = try_log_log_limit(quotient.value(), var, point, dir, depth)) {
            return log_log.value();
        }
    }

    if (!quotient.has_value()) {
        return direct;
    }

    return compute_quotient_limit(expr, var, point, dir, depth, quotient.value(), direct);
}

Result<ExprPtr> limit(ExprPtr expr, const Symbol& var, ExprPtr point, LimitDirection dir, symbolic::CASContext& ctx) {
    auto value = LimitEngine(ctx).compute(expr, var, point, dir);
    if (value.is_error()) {
        return value;
    }
    return symbolic::materialize_expr(value.value(), ctx.arena());
}

}  // namespace cas::calculus
