// Integrator::integrate_rational — rational-function integration pipeline:
// polynomial part by explicit division, Hermite reduction on the proper
// part, log part via Lazard-Rioboo-Trager (with the Apostol/Bronstein
// inverse-quadratic-power route past the LRT for ∫c/(x²±a²)^n).
// Split from integrate_core.cpp (anti-monolith, A7 step 5 completion).

#include "integrate_engine.hpp"

#include "cas/algebra.hpp"
#include "cas/differential_algebra.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <utility>
#include <vector>

namespace cas::calculus::integrate_detail {

Result<ExprPtr> Integrator::integrate_rational(ExprPtr expr, const Symbol& var) {
    auto parts = algebra::apart_num_den(expr, context_);
    if (parts.is_error()) return fail<ExprPtr>(parts.error());

    ExprPtr P = parts.value().numerator;
    ExprPtr Q = parts.value().denominator;

    // Handle polynomial part if deg(P) >= deg(Q)
    auto P_poly = algebra::parse_polynomial(P, var, context_);
    auto Q_poly = algebra::parse_polynomial(Q, var, context_);

    ExprPtr poly_integral = make_integer(arena_, 0);

    if (P_poly.is_ok() && Q_poly.is_ok() && !algebra::is_zero_poly(Q_poly.value())) {
        if (algebra::poly_degree(P_poly.value()) >= algebra::poly_degree(Q_poly.value())) {
            auto div_res = algebra::divide_poly_with_remainder(P_poly.value(), Q_poly.value(), context_);
            if (div_res.is_ok()) {
                std::vector<ExprPtr> terms;
                for (std::size_t k = 0; k < div_res.value().quotient.size(); ++k) {
                    ExprPtr coeff = div_res.value().quotient[k];
                    if (!coeff) continue;
                    if (const auto* il = expr_cast<IntegerLit>(coeff); il && il->value.is_zero()) continue;

                    ExprPtr kp1 = make_integer(arena_, static_cast<long long>(k + 1));
                    ExprPtr term;
                    if (k == 0) {
                        term = make_product(arena_, {coeff, arena_.make<Symbol>(var)});
                    } else {
                        term = make_product(arena_, {
                            make_binary(arena_, BinaryOp::Div, coeff, kp1),
                            make_binary(arena_, BinaryOp::Pow, arena_.make<Symbol>(var), kp1)
                        });
                    }
                    terms.push_back(term);
                }
                if (!terms.empty()) {
                    auto poly_sum = make_sum(arena_, std::move(terms));
                    auto simplified = context_.simplify(poly_sum);
                    if (simplified.is_ok()) poly_integral = simplified.value();
                }

                auto rem_expr = algebra::polynomial_to_expr(div_res.value().remainder, var, context_);
                if (rem_expr.is_ok()) P = rem_expr.value();
            }
        }
    }

    DifferentialField field(var);
    auto hermite = hermite_reduce(P, Q, var, field, context_);
    if (hermite.is_error()) return fail<ExprPtr>(hermite.error());
    // F7.5: route ∫c/(x²±a²)^n past the (buggy) LRT via Apostol/Bronstein recursion.
    ExprPtr lrt_value;
    if (ExprPtr rP = hermite.value().remaining_P, rQ = hermite.value().remaining_Q; !depends_on(rP, var)) {
        if (auto qd = algebra::polynomial_degree(rQ, var, context_); qd.is_ok() && qd.value() == 2U)
            if (auto r = integrate_inverse_quadratic_power(rQ, BigInt(-1), var); r.is_ok())
                lrt_value = make_product(arena_, {rP, r.value()});
        if (!lrt_value) if (const auto* dp = expr_cast<Binary>(rQ); dp && dp->op == BinaryOp::Pow)
            if (const auto* ei = expr_cast<IntegerLit>(dp->right); ei && !ei->value.is_negative() && ei->value >= BigInt(2))
                if (auto r = integrate_inverse_quadratic_power(dp->left, -ei->value, var); r.is_ok())
                    lrt_value = make_product(arena_, {rP, r.value()});
    }
    if (!lrt_value) {
        auto lrt_res = algebra::integrate_rational_lrt(hermite.value().remaining_P, hermite.value().remaining_Q, var, context_);
        if (lrt_res.is_error()) return fail<ExprPtr>(lrt_res.error());
        lrt_value = lrt_res.value();
    }
    auto lrt_res = ok(lrt_value);

    std::vector<ExprPtr> total;
    auto is_expr_zero = [](ExprPtr e) {
        if (!e) return true;
        if (const auto* il = expr_cast<IntegerLit>(e)) return il->value.is_zero();
        return false;
    };

    if (!is_expr_zero(poly_integral)) total.push_back(poly_integral);
    if (!is_expr_zero(hermite.value().rational_part)) total.push_back(hermite.value().rational_part);
    if (!is_expr_zero(lrt_res.value())) total.push_back(lrt_res.value());

    if (total.empty()) return ok(make_integer(arena_, 0));
    if (total.size() == 1) return ok(total[0]);
    return ok(make_sum(arena_, std::move(total)));
}

}  // namespace cas::calculus::integrate_detail
