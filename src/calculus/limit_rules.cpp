#include "limit_internal.hpp"
#include "cas/calculus.hpp"
#include "cas/algebra.hpp"
#include "cas/error.hpp"

#include <algorithm>
#include <utility>

namespace cas::calculus {

std::optional<Result<ExprPtr>> LimitEngine::try_log_log_limit(
    const QuotientView& quotient, const Symbol& var, ExprPtr point,
    LimitDirection dir, unsigned int depth) {
    // HPP-022 CLOSED: depth bound configurabile via ctx.max_log_log_limit_depth()
    // (default 3, see CASContextParams). Exceeding returns nullopt → caller
    // tries other strategies (never silent wrong answer).
    if (depth > context_.max_log_log_limit_depth()) return std::nullopt;

    const auto* num_call = expr_cast<FuncCall>(quotient.numerator);
    if (!num_call || num_call->func_id != BuiltinOp::Ln || num_call->args.size() != 1U)
        return std::nullopt;

    const auto* den_call = expr_cast<FuncCall>(quotient.denominator);
    if (!den_call || den_call->func_id != BuiltinOp::Ln || den_call->args.size() != 1U)
        return std::nullopt;

    ExprPtr a = num_call->args[0];
    ExprPtr b = den_call->args[0];

    // Verifica a,b→∞ al punto
    auto a_lim = compute_recursive(a, var, point, dir, depth + 1U);
    if (!a_lim.is_ok() || !limit_is_infinity(a_lim.value())) return std::nullopt;
    auto b_lim = compute_recursive(b, var, point, dir, depth + 1U);
    if (!b_lim.is_ok() || !limit_is_infinity(b_lim.value())) return std::nullopt;

    // Calcola lim(a/b) ricorsivamente
    ExprPtr a_over_b = arena_.make<Binary>(BinaryOp::Div, a, b);
    auto L = compute_recursive(a_over_b, var, point, dir, depth + 1U);

    if (!L.is_ok()) return std::nullopt;
    if (limit_is_zero(L.value()) || limit_is_infinity(L.value())) return std::nullopt;

    // L finito positivo → lim ln(a)/ln(b) = 1
    return ok(limit_make_integer(arena_, 1));
}

Result<ExprPtr> LimitEngine::compute_quotient_limit(
    ExprPtr /*expr*/,
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
        compute_recursive(quotient.denominator, var, point, dir, depth + 1U));

    // Bounded / Infinity -> 0 (Squeeze Theorem variant for quotients)
    if (denominator_at_point.is_ok() && limit_is_infinity(denominator_at_point.value())) {
        // Check numerator. If it's also infinity, it's inf/inf.
        auto numerator_check = as_infinity_if_undefined(quotient.numerator,
            compute_recursive(quotient.numerator, var, point, dir, depth + 1U));
        if (!numerator_check.is_ok() || !limit_is_infinity(numerator_check.value())) {
            if (is_bounded(quotient.numerator, var)) {
                return ok(limit_make_integer(arena_, 0));
            }
        }
    }

    auto numerator_at_point = as_infinity_if_undefined(quotient.numerator,
        compute_recursive(quotient.numerator, var, point, dir, depth + 1U));
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
        if (limit_is_zero(denominator_at_point.value()) && !limit_is_zero(numerator_at_point.value())) {
            // Pole: x/0 -> Inf
            return ok(arena_.make<Constant>(MathConstant::Infinity));
        }
        if (direct.is_ok()) return direct;
        return fail<ExprPtr>(direct.error());
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

Result<ExprPtr> LimitEngine::compute_lhopital_or_taylor(
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

}  // namespace cas::calculus
