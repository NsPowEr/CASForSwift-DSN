#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/normal_form.hpp"
#include "algebra_internal.hpp"

namespace cas::symbolic {

Result<bool> mathematically_equal(ExprPtr lhs, ExprPtr rhs, CASContext& context) {
    const bool owns_operation = !context.operation_active_;
    if (owns_operation) {
        context.operation_active_ = true;
        context.trace_capture_active_ = false;
        context.trace_.clear();
    }

    auto finalize = [&]() {
        if (owns_operation) context.operation_active_ = false;
    };

    auto lhs_s = context.simplify(lhs);
    if (lhs_s.is_error()) { finalize(); return fail<bool>(lhs_s.error()); }

    auto rhs_s = context.simplify(rhs);
    if (rhs_s.is_error()) { finalize(); return fail<bool>(rhs_s.error()); }

    // Fast path: structural equality after simplify
    if (structural_equal(lhs_s.value(), rhs_s.value())) {
        finalize();
        return ok(true);
    }

    // Mathematical equality using polynomial normal form: normal_form(expand(lhs - rhs)) == 0
    auto diff_expr = context.arena().make<Binary>(BinaryOp::Sub, lhs_s.value(), rhs_s.value());
    auto normal_diff = polynomial_normal_form(diff_expr, context);
    if (normal_diff.is_ok()) {
        if (expr_is<IntegerLit>(normal_diff.value()) && expr_cast<IntegerLit>(normal_diff.value())->value.is_zero()) {
            finalize();
            return ok(true);
        }
    }

    // Rational equality fallback: check num_L * den_R - num_R * den_L == 0 using normal form
    auto lhs_parts = algebra::split_num_den(lhs_s.value(), context);
    auto rhs_parts = algebra::split_num_den(rhs_s.value(), context);
    if (lhs_parts.is_ok() && rhs_parts.is_ok()) {
        auto cross_l = algebra::multiply_exprs(lhs_parts.value().numerator, rhs_parts.value().denominator, context);
        auto cross_r = algebra::multiply_exprs(rhs_parts.value().numerator, lhs_parts.value().denominator, context);
        if (cross_l.is_ok() && cross_r.is_ok()) {
            auto cross_diff_expr = context.arena().make<Binary>(BinaryOp::Sub, cross_l.value(), cross_r.value());
            auto cross_diff = polynomial_normal_form(cross_diff_expr, context);
            if (cross_diff.is_ok() && expr_is<IntegerLit>(cross_diff.value()) && expr_cast<IntegerLit>(cross_diff.value())->value.is_zero()) {
                finalize();
                return ok(true);
            }
        }
    }

    finalize();
    return ok(false);
}

} // namespace cas::symbolic
