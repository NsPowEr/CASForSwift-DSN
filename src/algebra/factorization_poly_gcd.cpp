// factorization_poly_gcd.cpp — polynomial_gcd dispatcher (univariate + multivariate).
// Extracted from factorization_polynomials.cpp (A1 anti-monolith split, F2 Block A).

#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <vector>

namespace cas::algebra {

Result<ExprPtr> polynomial_gcd(ExprPtr p, ExprPtr q, const Symbol& var, symbolic::CASContext& ctx) {
    const bool owns_operation = !ctx.operation_active_;
    if (owns_operation) {
        ctx.operation_active_ = true;
        ctx.trace_capture_active_ = ctx.trace_enabled_;
        ctx.trace_.clear();
        ctx.ops_count_ = 0;
        ctx.operation_started_at_ = std::chrono::steady_clock::now();
    }

    auto record_trace = [&](symbolic::RuleId rule_id, ExprPtr result) {
        if (!ctx.trace_enabled_) {
            return;
        }
        ctx.trace_.push_back(symbolic::TraceStep{
            .rule_id = rule_id,
            .depth = 0U,
            .target_before = p,
            .target_after = q,
            .root_after = result,
        });
    };

    auto finalize = [&]() {
        if (owns_operation) {
            ctx.operation_active_ = false;
            ctx.trace_capture_active_ = false;
            ctx.ops_count_ = 0;
        }
    };

    auto result = [&]() -> Result<ExprPtr> {
    if (!p || !q) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::InvalidArgument,
            "polynomial_gcd richiede due espressioni polinomiali non nulle"));
    }

    // Try multivariate GCD first if there are multiple variables
    auto p_multi = parse_multivariate_polynomial(p, ctx);
    auto q_multi = parse_multivariate_polynomial(q, ctx);

    if (p_multi.is_ok() && q_multi.is_ok()) {
        auto vars = p_multi.value().variables();
        auto vars_q = q_multi.value().variables();
        vars.insert(vars.end(), vars_q.begin(), vars_q.end());
        std::sort(vars.begin(), vars.end(), [](const auto& a, const auto& b){ return a.name < b.name; });
        vars.erase(std::unique(vars.begin(), vars.end(), [](const auto& a, const auto& b){ return a.name == b.name; }), vars.end());

        auto contains_non_main_var = [&](const MultivariatePolynomial& poly) {
            for (const Symbol& candidate : poly.variables()) {
                if (candidate.name != var.name) {
                    return true;
                }
            }
            return false;
        };

        if (contains_non_main_var(p_multi.value()) && contains_non_main_var(q_multi.value())) {
            auto gcd_res = gcd_multivariate_eval_interp(p_multi.value(), q_multi.value(), ctx);
            if (gcd_res.is_ok()) {
                auto gcd_expr = multivariate_to_expr(gcd_res.value(), ctx);
                if (gcd_expr.is_error()) {
                    return fail<ExprPtr>(gcd_expr.error());
                }
                return ctx.simplify(gcd_expr.value());
            }
        }
    }

    auto left = parse_polynomial(p, var, ctx);
    if (left.is_error()) {
        return fail<ExprPtr>(left.error());
    }
    auto right = parse_polynomial(q, var, ctx);
    if (right.is_error()) {
        return fail<ExprPtr>(right.error());
    }

    auto left_integer = poly_to_integer_poly(left.value());
    auto right_integer = poly_to_integer_poly(right.value());
    if (left_integer.is_ok() && right_integer.is_ok()) {
        const BigInt left_content = integer_content(left_integer.value());
        const BigInt right_content = integer_content(right_integer.value());
        IntegerGcdResult gcd_result =
            gcd_integer_poly_dispatch(left_integer.value(), right_integer.value(), ctx);
        IntPoly gcd_poly = std::move(gcd_result.gcd);
        symbolic::RuleId path_rule = symbolic::RuleId::PolynomialGcdPrimitiveFallback;
        switch (gcd_result.path) {
        case IntegerGcdPath::Subresultant:
            path_rule = symbolic::RuleId::PolynomialGcdSubresultant;
            break;
        case IntegerGcdPath::PrimitiveFallbackPsi:
            path_rule = symbolic::RuleId::PolynomialGcdPrimitiveFallbackPsi;
            break;
        case IntegerGcdPath::PrimitiveFallbackBeta:
            path_rule = symbolic::RuleId::PolynomialGcdPrimitiveFallbackBeta;
            break;
        case IntegerGcdPath::PrimitiveFallback:
            path_rule = symbolic::RuleId::PolynomialGcdPrimitiveFallback;
            break;
        case IntegerGcdPath::HalfGcd:
            path_rule = symbolic::RuleId::PolynomialGcdSubresultant;  // closest analogue
            break;
        case IntegerGcdPath::ModularCrt:
            path_rule = symbolic::RuleId::PolynomialGcdSubresultant;  // B2.1 modular CRT path
            break;
        }

        if (is_zero_integer_poly(gcd_poly)) {
            ExprPtr zero = make_integer(ctx.arena(), 0);
            record_trace(path_rule, zero);
            return ok(zero);
        }

        auto gcd_expr = integer_coefficients_to_expr(gcd_poly, var, ctx);
        if (gcd_expr.is_error()) {
            return fail<ExprPtr>(gcd_expr.error());
        }

        Result<ExprPtr> traced_result = ctx.simplify(gcd_expr.value());
        if (left_content == right_content && traced_result.is_ok()) {
            auto parsed_gcd = parse_polynomial(traced_result.value(), var, ctx);
            if (parsed_gcd.is_ok()) {
                auto monic = normalize_poly_monic(parsed_gcd.value(), ctx);
                if (monic.is_ok()) {
                    traced_result = polynomial_to_expr(monic.value(), var, ctx);
                }
            }
        }
        if (traced_result.is_ok()) {
            record_trace(path_rule, traced_result.value());
        }
        return traced_result;
    }

    if (is_zero_poly(left.value()) && is_zero_poly(right.value())) {
        ExprPtr zero = make_integer(ctx.arena(), 0);
        record_trace(symbolic::RuleId::PolynomialGcdSymbolicEuclidean, zero);
        return ok(zero);
    }
    if (is_zero_poly(left.value())) {
        auto normalized = normalize_poly_monic(right.value(), ctx);
        if (normalized.is_error()) {
            return fail<ExprPtr>(normalized.error());
        }
        auto traced_result = polynomial_to_expr(normalized.value(), var, ctx);
        if (traced_result.is_ok()) {
            record_trace(symbolic::RuleId::PolynomialGcdSymbolicEuclidean, traced_result.value());
        }
        return traced_result;
    }
    if (is_zero_poly(right.value())) {
        auto normalized = normalize_poly_monic(left.value(), ctx);
        if (normalized.is_error()) {
            return fail<ExprPtr>(normalized.error());
        }
        auto traced_result = polynomial_to_expr(normalized.value(), var, ctx);
        if (traced_result.is_ok()) {
            record_trace(symbolic::RuleId::PolynomialGcdSymbolicEuclidean, traced_result.value());
        }
        return traced_result;
    }

    PolyExpr a = std::move(left.value());
    PolyExpr b = std::move(right.value());

    while (!is_zero_poly(b)) {
        auto division = divide_poly_with_remainder(a, b, ctx);
        if (division.is_error()) {
            return fail<ExprPtr>(division.error());
        }
        a = std::move(b);
        b = std::move(division.value().remainder);
    }

    auto normalized = normalize_poly_monic(a, ctx);
    if (normalized.is_error()) {
        return fail<ExprPtr>(normalized.error());
    }
    auto traced_result = polynomial_to_expr(normalized.value(), var, ctx);
    if (traced_result.is_ok()) {
        record_trace(symbolic::RuleId::PolynomialGcdSymbolicEuclidean, traced_result.value());
    }
    return traced_result;
    }();

    finalize();
    return result;
}

} // namespace cas::algebra
