#include "integrate_engine.hpp"

#include "cas/algebra.hpp"

#include <string>
#include <utility>
#include <vector>

namespace cas::calculus::integrate_detail {

Result<ExprPtr> Integrator::try_u_substitution_for_product(const Product& product, const Symbol& var) {
    if (product.factors.size() < 2U) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "u-substitution needs at least two factors"));
    }

    for (std::size_t index = 0; index < product.factors.size(); ++index) {
        const auto* call = expr_cast<FuncCall>(product.factors[index]);
        if (call == nullptr || call->args.size() != 1U) {
            continue;
        }

        std::vector<ExprPtr> derivative_factors;
        derivative_factors.reserve(product.factors.size() - 1U);
        for (std::size_t factor_index = 0; factor_index < product.factors.size(); ++factor_index) {
            if (factor_index != index) {
                derivative_factors.push_back(product.factors[factor_index]);
            }
        }

        const ExprPtr derivative_candidate = make_product(arena_, std::move(derivative_factors));
        auto inner_derivative = diff(call->args.front(), var, 1U, context_);
        if (inner_derivative.is_error()) {
            continue;
        }

        auto matches = expressions_match_after_simplify(derivative_candidate, inner_derivative.value());
        if (matches.is_error()) {
            return fail<ExprPtr>(matches.error());
        }
        if (matches.value()) {
            return integrate_function_direct(canonical_function_name(call->name), call->args.front());
        }

        // Scaled u-sub: derivative_candidate = k * inner_derivative where k is constant
        auto ratio_expr = context_.simplify(arena_.make<Binary>(
            BinaryOp::Div, derivative_candidate, inner_derivative.value()));
        if (ratio_expr.is_ok() && !depends_on(ratio_expr.value(), var)) {
            auto primitive = integrate_function_direct(canonical_function_name(call->name), call->args.front());
            if (primitive.is_ok()) {
                auto scaled = context_.simplify(arena_.make<Product>(
                    std::vector<ExprPtr>{ratio_expr.value(), primitive.value()}));
                if (scaled.is_ok()) {
                    return scaled;
                }
            }
        }
    }

    for (std::size_t index = 0; index < product.factors.size(); ++index) {
        const auto* power = expr_cast<Binary>(product.factors[index]);
        if (power == nullptr || power->op != BinaryOp::Pow) {
            continue;
        }

        std::vector<ExprPtr> derivative_factors;
        derivative_factors.reserve(product.factors.size() - 1U);
        for (std::size_t factor_index = 0; factor_index < product.factors.size(); ++factor_index) {
            if (factor_index != index) {
                derivative_factors.push_back(product.factors[factor_index]);
            }
        }

        const ExprPtr derivative_candidate = make_product(arena_, std::move(derivative_factors));
        auto inner_derivative = diff(power->left, var, 1U, context_);
        if (inner_derivative.is_error()) {
            continue;
        }

        auto matches = expressions_match_after_simplify(derivative_candidate, inner_derivative.value());
        if (matches.is_error()) {
            return fail<ExprPtr>(matches.error());
        }
        if (matches.value()) {
            return integrate_power_direct(power->left, power->right, var);
        }

        // Scaled u-sub: derivative_candidate = k * inner_derivative where k is constant
        auto ratio_expr = context_.simplify(arena_.make<Binary>(
            BinaryOp::Div, derivative_candidate, inner_derivative.value()));
        if (ratio_expr.is_ok() && !depends_on(ratio_expr.value(), var)) {
            auto primitive = integrate_power_direct(power->left, power->right, var);
            if (primitive.is_ok()) {
                auto scaled = context_.simplify(arena_.make<Product>(
                    std::vector<ExprPtr>{ratio_expr.value(), primitive.value()}));
                if (scaled.is_ok()) {
                    return scaled;
                }
            }
        }
    }

    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "No supported u-substitution pattern found"));
}

Result<ExprPtr> Integrator::integrate_via_partial_fractions(ExprPtr expr, const Symbol& var) {
    auto terms = algebra::partial_fractions(expr, var, context_);
    if (terms.is_error()) {
        return fail<ExprPtr>(terms.error());
    }

    if (terms.value().empty()) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Partial fraction decomposition produced no terms"));
    }

    // Safety: If only one term is returned, it means no decomposition happened.
    // Proceeding would lead to infinite recursion if integrate_once calls this again.
    if (terms.value().size() <= 1U) {
         return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Partial fractions did not decompose the expression"));
    }

    std::vector<ExprPtr> primitives;
    primitives.reserve(terms.value().size());
    for (ExprPtr term : terms.value()) {
        auto simplified_term = context_.simplify(term);
        if (simplified_term.is_error()) {
            return fail<ExprPtr>(simplified_term.error());
        }
        if (structural_equal(simplified_term.value(), expr)) {
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Partial fractions simplified back to the original integrand"));
        }
        auto primitive = integrate_once(simplified_term.value(), var);
        if (primitive.is_error()) {
            return primitive;
        }
        primitives.push_back(primitive.value());
    }
    return ok(make_sum(arena_, std::move(primitives)));
}

}  // namespace cas::calculus::integrate_detail
