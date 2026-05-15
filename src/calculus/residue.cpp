#include "cas/calculus.hpp"
#include "cas/algebra.hpp"
#include "calculus_internal.hpp"

#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] ExprPtr integer(AstArena& arena, long long value) {
    return arena.make<IntegerLit>(BigInt(value));
}

[[nodiscard]] bool is_zero_literal(ExprPtr expr) {
    if (const auto* integer_lit = expr_cast<IntegerLit>(expr)) {
        return integer_lit->value.is_zero();
    }
    if (const auto* rational_lit = expr_cast<RationalLit>(expr)) {
        return rational_lit->numerator.is_zero();
    }
    return false;
}

[[nodiscard]] BigInt factorial(unsigned int n) {
    BigInt result(1);
    for (unsigned int k = 2U; k <= n; ++k) {
        result *= BigInt(static_cast<long long>(k));
    }
    return result;
}

[[nodiscard]] Result<ExprPtr> simplify_or_fail(ExprPtr expr, symbolic::CASContext& ctx) {
    auto simplified = ctx.simplify(expr);
    if (simplified.is_error()) {
        return fail<ExprPtr>(simplified.error());
    }
    return simplified;
}

[[nodiscard]] Result<ExprPtr> add(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    return simplify_or_fail(ctx.arena().make<Binary>(BinaryOp::Add, lhs, rhs), ctx);
}

[[nodiscard]] Result<ExprPtr> sub(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    return simplify_or_fail(ctx.arena().make<Binary>(BinaryOp::Sub, lhs, rhs), ctx);
}

[[nodiscard]] Result<ExprPtr> mul(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    return simplify_or_fail(ctx.arena().make<Binary>(BinaryOp::Mul, lhs, rhs), ctx);
}

[[nodiscard]] Result<ExprPtr> div(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    return simplify_or_fail(ctx.arena().make<Binary>(BinaryOp::Div, lhs, rhs), ctx);
}

[[nodiscard]] Result<ExprPtr> local_taylor_coefficient(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr pole,
    unsigned int order,
    symbolic::CASContext& ctx) {
    ExprPtr derivative = expr;
    if (order > 0U) {
        auto diffed = diff(expr, var, order, ctx);
        if (diffed.is_error()) {
            return fail<ExprPtr>(diffed.error());
        }
        derivative = diffed.value();
    }

    auto substituted = symbolic::substitute(derivative, var, pole, ctx);
    if (substituted.is_error()) {
        return fail<ExprPtr>(substituted.error());
    }

    auto coefficient = simplify_or_fail(substituted.value(), ctx);
    if (coefficient.is_error()) {
        return coefficient;
    }

    if (order <= 1U) {
        return coefficient;
    }
    return div(coefficient.value(), ctx.arena().make<IntegerLit>(factorial(order)), ctx);
}

[[nodiscard]] Result<unsigned int> denominator_zero_order(
    ExprPtr denominator,
    const Symbol& var,
    ExprPtr pole,
    symbolic::CASContext& ctx) {
    const unsigned int max_order = static_cast<unsigned int>(ctx.max_integration_depth());
    for (unsigned int order = 0U; order <= max_order; ++order) {
        auto coefficient = local_taylor_coefficient(denominator, var, pole, order, ctx);
        if (coefficient.is_error()) {
            return fail<unsigned int>(coefficient.error());
        }
        if (!is_zero_literal(coefficient.value())) {
            return ok(order);
        }
    }

    return fail<unsigned int>(CASError{
        .kind = CASErrorKind::Unimplemented,
        .message = "Residue denominator zero order exceeds configured integration depth"});
}

[[nodiscard]] Result<ExprPtr> rational_residue_from_series(
    ExprPtr numerator,
    ExprPtr denominator,
    const Symbol& var,
    ExprPtr pole,
    unsigned int denominator_order,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> quotient_coeffs;
    quotient_coeffs.reserve(denominator_order);

    std::vector<ExprPtr> shifted_denominator_coeffs;
    shifted_denominator_coeffs.reserve(denominator_order);
    for (unsigned int i = 0U; i < denominator_order; ++i) {
        auto coeff = local_taylor_coefficient(denominator, var, pole, denominator_order + i, ctx);
        if (coeff.is_error()) {
            return fail<ExprPtr>(coeff.error());
        }
        shifted_denominator_coeffs.push_back(coeff.value());
    }

    for (unsigned int n = 0U; n < denominator_order; ++n) {
        auto numerator_coeff = local_taylor_coefficient(numerator, var, pole, n, ctx);
        if (numerator_coeff.is_error()) {
            return fail<ExprPtr>(numerator_coeff.error());
        }

        ExprPtr correction = integer(arena, 0);
        for (unsigned int i = 1U; i <= n; ++i) {
            auto term = mul(shifted_denominator_coeffs[i], quotient_coeffs[n - i], ctx);
            if (term.is_error()) {
                return fail<ExprPtr>(term.error());
            }
            auto next_correction = add(correction, term.value(), ctx);
            if (next_correction.is_error()) {
                return fail<ExprPtr>(next_correction.error());
            }
            correction = next_correction.value();
        }

        auto adjusted = sub(numerator_coeff.value(), correction, ctx);
        if (adjusted.is_error()) {
            return fail<ExprPtr>(adjusted.error());
        }

        auto quotient_coeff = div(adjusted.value(), shifted_denominator_coeffs[0], ctx);
        if (quotient_coeff.is_error()) {
            return fail<ExprPtr>(quotient_coeff.error());
        }
        quotient_coeffs.push_back(quotient_coeff.value());
    }

    return ok(quotient_coeffs.back());
}

} // namespace

Result<ExprPtr> residue(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr pole,
    symbolic::CASContext& ctx) {
    auto parts = algebra::apart_num_den(expr, ctx);
    if (parts.is_error()) {
        return fail<ExprPtr>(parts.error());
    }

    auto order = denominator_zero_order(parts.value().denominator, var, pole, ctx);
    if (order.is_error()) {
        return fail<ExprPtr>(order.error());
    }
    if (order.value() == 0U) {
        return ok(integer(ctx.arena(), 0));
    }

    return rational_residue_from_series(
        parts.value().numerator,
        parts.value().denominator,
        var,
        pole,
        order.value(),
        ctx);
}

// Build the full Laurent expansion of N/D around `pole`.
// Lemma: write D(x) = (x - x0)^m * D'(x), with D'(x0) != 0.  Then
// for n = k + m,
//   c_k = ( N_n  -  Σ_{i=1..n} D_{m+i} * c_{k - i} )  /  D_m
// where N_n, D_n are Taylor coefficients of N, D at x0.
//
// We extend the recurrence beyond k = -1 (which residue() already supports)
// up to k = positive_order.
//
// `m` (= denominator_order) is detected by denominator_zero_order.
namespace {

[[nodiscard]] Result<LaurentExpansion> rational_laurent_from_series(
    ExprPtr numerator,
    ExprPtr denominator,
    const Symbol& var,
    ExprPtr pole,
    unsigned int denominator_order,        // m: pole multiplicity
    unsigned int positive_order,           // highest positive power requested
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    const int leading_order = -static_cast<int>(denominator_order);
    const unsigned int n_terms = denominator_order + positive_order + 1U;

    // Pre-compute D_{m + i} for i = 0 .. n_terms (need up to index n_terms - 1 inclusive).
    std::vector<ExprPtr> shifted_D_coeffs;
    shifted_D_coeffs.reserve(n_terms);
    for (unsigned int i = 0U; i < n_terms; ++i) {
        auto coeff = local_taylor_coefficient(denominator, var, pole, denominator_order + i, ctx);
        if (coeff.is_error()) return fail<LaurentExpansion>(coeff.error());
        shifted_D_coeffs.push_back(coeff.value());
    }

    // Recurrence: c_k for k = -m .. positive_order.
    std::vector<ExprPtr> coefficients;
    coefficients.reserve(n_terms);
    for (unsigned int idx = 0U; idx < n_terms; ++idx) {
        // n = idx (because k = leading_order + idx and n = k + m = idx).
        const unsigned int n = idx;
        auto N_n = local_taylor_coefficient(numerator, var, pole, n, ctx);
        if (N_n.is_error()) return fail<LaurentExpansion>(N_n.error());

        ExprPtr correction = integer(arena, 0);
        for (unsigned int i = 1U; i <= n; ++i) {
            // term = D_{m+i} * c_{k - i}  =  D_{m+i} * coefficients[idx - i]
            auto term = mul(shifted_D_coeffs[i], coefficients[idx - i], ctx);
            if (term.is_error()) return fail<LaurentExpansion>(term.error());
            auto next_correction = add(correction, term.value(), ctx);
            if (next_correction.is_error()) return fail<LaurentExpansion>(next_correction.error());
            correction = next_correction.value();
        }

        auto adjusted = sub(N_n.value(), correction, ctx);
        if (adjusted.is_error()) return fail<LaurentExpansion>(adjusted.error());
        auto c_k = div(adjusted.value(), shifted_D_coeffs[0], ctx);
        if (c_k.is_error()) return fail<LaurentExpansion>(c_k.error());
        coefficients.push_back(c_k.value());
    }

    ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, arena.make<Symbol>(var), pole);
    ExprPtr remainder_power = arena.make<Binary>(
        BinaryOp::Pow,
        delta,
        arena.make<IntegerLit>(BigInt(static_cast<long long>(positive_order + 1U))));
    ExprPtr remainder = arena.make<FuncCall>("O", std::vector<ExprPtr>{remainder_power});

    return ok(LaurentExpansion{
        .center = pole,
        .leading_order = leading_order,
        .coefficients = std::move(coefficients),
        .positive_order = positive_order,
        .remainder = remainder,
    });
}

}  // namespace

Result<LaurentExpansion> laurent_series(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr center,
    unsigned int positive_order,
    symbolic::CASContext& ctx) {
    auto parts = algebra::apart_num_den(expr, ctx);
    if (parts.is_error()) return fail<LaurentExpansion>(parts.error());

    auto order = denominator_zero_order(parts.value().denominator, var, center, ctx);
    if (order.is_error()) return fail<LaurentExpansion>(order.error());

    return rational_laurent_from_series(
        parts.value().numerator,
        parts.value().denominator,
        var,
        center,
        order.value(),
        positive_order,
        ctx);
}

} // namespace cas::calculus
