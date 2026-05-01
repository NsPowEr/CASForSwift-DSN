#include "cas/calculus.hpp"
#include "cas/algebra.hpp"

#include "calculus_internal.hpp"
#include "cas/ast_debug.hpp"
#include "cas/error.hpp"
#include "cas/formatter.hpp"

#include <iostream>
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

class LimitEngine {
public:
    explicit LimitEngine(symbolic::CASContext& context) noexcept : context_(context), arena_(context.arena()) {}

    [[nodiscard]] Result<ExprPtr> compute(ExprPtr expr, const Symbol& var, ExprPtr point, LimitDirection dir) {
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
            auto infinite_limit = try_infinite_limit(simplified_expr.value(), var, simplified_point.value(), arena_);
            if (infinite_limit.is_ok() || infinite_limit.error().kind != CASErrorKind::Unimplemented) {
                return infinite_limit;
            }
        }

        if (point_is_pos_inf || point_is_neg_inf) {
            auto q_inf = extract_quotient_view(simplified_expr.value(), arena_);
            if (q_inf.has_value()) {
                LimitDirection inf_dir = point_is_pos_inf
                    ? LimitDirection::Right : LimitDirection::Left;
                if (auto ll = try_log_log_limit(q_inf.value(), var,
                                                 simplified_point.value(), inf_dir, 0U)) {
                    return ll.value();
                }
            }
        }

        if (point_is_pos_inf) {
            auto mrv_res = compute_limit_mrv(simplified_expr.value(), var, simplified_point.value(), context_);
            if (mrv_res.is_ok()) {
                return mrv_res;
            }
        }

        // Sviluppo Asintotico all'infinito via trasformazione x -> ±1/t
        if (point_is_pos_inf || point_is_neg_inf) {
            const std::string t_name = var.name + "_inv";
            const Symbol t_var(t_name);
            ExprPtr t_sym = arena_.make<Symbol>(t_name);

            // lim_{x -> +∞} f(x) = lim_{t -> 0+} f(1/t)
            // lim_{x -> -∞} f(x) = lim_{t -> 0+} f(-1/t)  preserves domain: -1/t < 0 for t > 0
            ExprPtr one_over_t = arena_.make<Binary>(BinaryOp::Div, limit_make_integer(arena_, 1), t_sym);
            ExprPtr x_sub = point_is_pos_inf
                ? one_over_t
                : arena_.make<Unary>(UnaryOp::Neg, one_over_t);
            auto expr_t = context_.substitute(simplified_expr.value(), var, x_sub);
            if (expr_t.is_error()) return expr_t;

            ExprPtr zero_pt = limit_make_integer(arena_, 0);
            return compute_recursive(expr_t.value(), t_var, zero_pt, LimitDirection::Right, 0U);
        }

        return compute_recursive(simplified_expr.value(), var, simplified_point.value(), dir, 0U);
    }

private:
    [[nodiscard]] Result<ExprPtr> substitute_and_simplify(ExprPtr expr, const Symbol& var, ExprPtr point) {
        auto substituted = context_.substitute(expr, var, point);
        if (substituted.is_error()) {
            return substituted;
        }
        return context_.simplify(substituted.value());
    }

    [[nodiscard]] Result<ExprPtr> compute_recursive(
        ExprPtr expr,
        const Symbol& var,
        ExprPtr point,
        LimitDirection dir,
        unsigned int depth) {
        if (depth > 6U) {
            return fail<ExprPtr>(make_error(
                CASErrorKind::Unimplemented,
                "Il limite richiede piu' iterazioni di quelle supportate"));
        }

        auto direct = substitute_and_simplify(expr, var, point);

        // Squeeze Theorem: lim(f*g) = 0 if lim(f)=0 and g is bounded
        if (const auto* prod = expr_cast<Product>(expr)) {
            bool has_zero = false;
            bool all_bounded = true;
            for (auto factor : prod->factors) {
                auto lim_f = compute_recursive(factor, var, point, dir, depth + 1U);
                if (lim_f.is_ok() && limit_is_zero(lim_f.value())) {
                    has_zero = true;
                } else if (!is_bounded(factor, var)) {
                    all_bounded = false;
                    break;
                }
            }
            if (has_zero && all_bounded) {
                return ok(limit_make_integer(arena_, 0));
            }
        }

        if (const auto* bin = expr_cast<Binary>(expr)) {
            if (bin->op == BinaryOp::Mul) {
                auto lim_l = compute_recursive(bin->left, var, point, dir, depth + 1U);
                if (lim_l.is_ok() && limit_is_zero(lim_l.value()) && is_bounded(bin->right, var)) {
                    return ok(limit_make_integer(arena_, 0));
                }
                auto lim_r = compute_recursive(bin->right, var, point, dir, depth + 1U);
                if (lim_r.is_ok() && limit_is_zero(lim_r.value()) && is_bounded(bin->left, var)) {
                    return ok(limit_make_integer(arena_, 0));
                }
            }
        }

        auto quotient = extract_quotient_view(expr, arena_);
        const bool has_quotient = quotient.has_value();

        if (!limit_is_infinity(point)) {
            auto logarithmic_limit = try_logarithmic_root_limit(expr, var, point, dir, arena_);
            if (logarithmic_limit.has_value()) {
                return logarithmic_limit.value();
            }
        }

        if (limit_is_infinity(point)) {
            auto infinite_limit = try_infinite_limit(expr, var, point, arena_);
            if (infinite_limit.is_ok() || infinite_limit.error().kind != CASErrorKind::Unimplemented) {
                return infinite_limit;
            }
        }

        if (direct.is_error() && !has_quotient) {
            return direct;
        }
        if (direct.is_ok() && (!expr_is<Binary>(direct.value()) || direct.value() == expr)) {
            if (!expr_is<Binary>(expr) || (!limit_is_zero(direct.value()) && !limit_is_infinity(direct.value()))) {
                return direct;
            }
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

        // ln(a)/ln(b) con lim(a/b) finito positivo e a,b→∞ → limite = 1
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

    // Regola: ln(a)/ln(b) → 1 quando lim(a/b) finito positivo e a,b→∞
    [[nodiscard]] std::optional<Result<ExprPtr>> try_log_log_limit(
        const QuotientView& quotient, const Symbol& var, ExprPtr point,
        LimitDirection dir, unsigned int depth) {
        if (depth > 3U) return std::nullopt; // evita ricorsione eccessiva

        const auto* num_call = expr_cast<FuncCall>(quotient.numerator);
        if (!num_call || num_call->func_id != BuiltinOp::Ln || num_call->args.size() != 1U)
            return std::nullopt;

        const auto* den_call = expr_cast<FuncCall>(quotient.denominator);
        if (!den_call || den_call->func_id != BuiltinOp::Ln || den_call->args.size() != 1U)
            return std::nullopt;

        ExprPtr a = num_call->args[0];
        ExprPtr b = den_call->args[0];

        // Verifica a,b→∞ al punto (sostituzione produce errore Undefined = diverge)
        auto a_at_pt = substitute_and_simplify(a, var, point);
        if (a_at_pt.is_ok() && !limit_is_infinity(a_at_pt.value())) return std::nullopt;
        auto b_at_pt = substitute_and_simplify(b, var, point);
        if (b_at_pt.is_ok() && !limit_is_infinity(b_at_pt.value())) return std::nullopt;

        // Calcola lim(a/b) ricorsivamente
        ExprPtr a_over_b = arena_.make<Binary>(BinaryOp::Div, a, b);
        auto L = compute_recursive(a_over_b, var, point, dir, depth + 1U);

        if (!L.is_ok()) return std::nullopt;
        if (limit_is_zero(L.value()) || limit_is_infinity(L.value())) return std::nullopt;

        // L finito positivo → lim ln(a)/ln(b) = 1
        return ok(limit_make_integer(arena_, 1));
    }

    [[nodiscard]] Result<ExprPtr> compute_quotient_limit(
        ExprPtr,
        const Symbol& var,
        ExprPtr point,
        LimitDirection dir,
        unsigned int depth,
        const QuotientView& quotient,
        const Result<ExprPtr>& direct) {
        // Divisione per zero durante sostituzione → espressione diverge al punto → trattata come ∞.
        // Eccezione: ln(neg_arg) — se l'argomento ha segno negativo, è un errore di dominio.
        auto as_infinity_if_undefined = [&](ExprPtr original_expr, Result<ExprPtr> r) -> Result<ExprPtr> {
            if (!r.is_error() || r.error().kind != CASErrorKind::Undefined) return r;
            // ln(Unary(Neg, ...)) → argomento negativo → dominio non valido, NON convertire a ∞
            if (const auto* fc = expr_cast<FuncCall>(original_expr)) {
                if (fc->func_id == BuiltinOp::Ln && fc->args.size() == 1U) {
                    if (expr_is<Unary>(fc->args[0]) &&
                        expr_ref<Unary>(fc->args[0]).op == UnaryOp::Neg) {
                        return r;
                    }
                }
            }
            return ok(arena_.make<Constant>(MathConstant::Infinity));
        };

        auto denominator_at_point = as_infinity_if_undefined(quotient.denominator,
            substitute_and_simplify(quotient.denominator, var, point));

        // Bounded / Infinity -> 0 (Squeeze Theorem variant for quotients)
        if (denominator_at_point.is_ok() && limit_is_infinity(denominator_at_point.value())) {
            // Check numerator. If it's also infinity, it's inf/inf.
            auto numerator_check = as_infinity_if_undefined(quotient.numerator,
                substitute_and_simplify(quotient.numerator, var, point));
            if (!numerator_check.is_ok() || !limit_is_infinity(numerator_check.value())) {
                if (is_bounded(quotient.numerator, var)) {
                    return ok(limit_make_integer(arena_, 0));
                }
            }
        }

        auto numerator_at_point = as_infinity_if_undefined(quotient.numerator,
            substitute_and_simplify(quotient.numerator, var, point));
        if (numerator_at_point.is_error()) {
            return numerator_at_point;
        }

        if (denominator_at_point.is_error()) {
            return denominator_at_point;
        }

        if (!limit_is_infinity(point) && limit_is_zero(denominator_at_point.value())) {
            auto pole_limit = try_polynomial_pole_limit(quotient.numerator, quotient.denominator, var, point, dir, arena_);
            if (pole_limit.has_value()) {
                return pole_limit.value();
            }
        }

        const bool zero_over_zero = limit_is_zero(numerator_at_point.value()) && limit_is_zero(denominator_at_point.value());
        const bool infinity_over_infinity =
            limit_is_infinity(numerator_at_point.value()) && limit_is_infinity(denominator_at_point.value());
        if (!zero_over_zero && !infinity_over_infinity) {
            return direct.is_ok()
                ? direct
                : fail<ExprPtr>(make_error(
                    CASErrorKind::Unimplemented,
                    "Il limite non e' una forma indeterminata supportata"));
        }

        if (zero_over_zero) {
            auto canceled = cancel_common_linear_factor(quotient.numerator, quotient.denominator, var, point, arena_);
            if (canceled.has_value()) {
                auto simplified_canceled = context_.simplify(canceled.value());
                if (simplified_canceled.is_ok()) {
                    auto canceled_limit = compute_recursive(simplified_canceled.value(), var, point, dir, depth + 1U);
                    if (canceled_limit.is_ok()) {
                        return canceled_limit;
                    }
                }
            }
        }

        return compute_lhopital_or_taylor(var, point, dir, depth, quotient);
    }

    [[nodiscard]] Result<ExprPtr> compute_lhopital_or_taylor(
        const Symbol& var,
        ExprPtr point,
        LimitDirection dir,
        unsigned int depth,
        const QuotientView& quotient) {
        // 1. Taylor Expansion around finite point (0/0 or inf/inf transformed)
        if (!limit_is_infinity(point)) {
            for (unsigned int expand_order = 4U; expand_order <= 20U; expand_order += 4U) {
                auto taylor_num = taylor_series(quotient.numerator, var, point, expand_order, context_);
                if (taylor_num.is_error()) break;
                auto taylor_den = taylor_series(quotient.denominator, var, point, expand_order, context_);
                if (taylor_den.is_error()) break;

                // Find valuation k: degree of first non-zero term
                auto find_valuation = [&](ExprPtr poly) -> std::optional<unsigned int> {
                    for (unsigned int k = 0U; k <= expand_order; ++k) {
                        Result<ExprPtr> deriv = (k == 0U) ? ok(poly) : diff(poly, var, k, context_);
                        if (deriv.is_error()) return std::nullopt;
                        auto at_pt = substitute_and_simplify(deriv.value(), var, point);
                        if (at_pt.is_ok() && !limit_is_zero(at_pt.value())) return k;
                    }
                    return std::nullopt;
                };

                auto k_num = find_valuation(taylor_num.value().polynomial);
                auto k_den = find_valuation(taylor_den.value().polynomial);

                if (!k_num || !k_den || *k_num == expand_order || *k_den == expand_order) continue;

                // k_num > k_den: numerator vanishes faster than denominator -> 0
                if (*k_num > *k_den) return ok(limit_make_integer(arena_, 0));
                
                // k_num == k_den: limit is ratio of leading coefficients
                if (*k_num == *k_den) {
                    unsigned int k = *k_num;
                    auto d_num = (k == 0U) ? ok(taylor_num.value().polynomial) : diff(taylor_num.value().polynomial, var, k, context_);
                    auto d_den = (k == 0U) ? ok(taylor_den.value().polynomial) : diff(taylor_den.value().polynomial, var, k, context_);
                    if (d_num.is_error() || d_den.is_error()) break;

                    auto c_num = substitute_and_simplify(d_num.value(), var, point);
                    auto c_den = substitute_and_simplify(d_den.value(), var, point);
                    if (c_num.is_error() || c_den.is_error()) break;

                    auto ratio = context_.simplify(arena_.make<Binary>(BinaryOp::Div, c_num.value(), c_den.value()));
                    if (ratio.is_ok()) return ratio;
                    break;
                }

                // k_num < k_den: tends to infinity. Let L'Hopital or recursion handle signs.
                break;
            }
        }

        // 2. L'Hopital rule as fallback
        auto numerator_derivative = diff(quotient.numerator, var, 1U, context_);
        if (numerator_derivative.is_error()) return numerator_derivative;
        auto denominator_derivative = diff(quotient.denominator, var, 1U, context_);
        if (denominator_derivative.is_error()) return denominator_derivative;

        if (structural_equal(numerator_derivative.value(), denominator_derivative.value())) {
            return ok(limit_make_integer(arena_, 1));
        }

        ExprPtr transformed = limit_make_binary(arena_, BinaryOp::Div, numerator_derivative.value(), denominator_derivative.value());
        auto simplified_transformed = context_.simplify(transformed);
        if (simplified_transformed.is_error()) return simplified_transformed;

        return compute_recursive(simplified_transformed.value(), var, point, dir, depth + 1U);
    }

    symbolic::CASContext& context_;
    AstArena& arena_;
};

}  // namespace

Result<ExprPtr> limit(ExprPtr expr, const Symbol& var, ExprPtr point, LimitDirection dir, symbolic::CASContext& ctx) {
    auto value = LimitEngine(ctx).compute(expr, var, point, dir);
    if (value.is_error()) {
        return value;
    }
    return symbolic::materialize_expr(value.value(), ctx.arena());
}

}  // namespace cas::calculus
