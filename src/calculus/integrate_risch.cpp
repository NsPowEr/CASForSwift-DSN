#include "calculus_internal.hpp"
#include "integrate_engine.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/algebra.hpp"
#include "cas/error.hpp"
#include "../algebra/polynomial_internal.hpp"
#include <vector>

namespace cas::calculus {

Result<ExprPtr> integrate_risch(ExprPtr expr, const Symbol& var, symbolic::CASContext& context) {
    AstArena& arena = context.arena();

    // 0. Simple pattern matching for fundamental transcendental functions
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if (call->args.size() == 1 && integrate_detail::is_same_symbol(call->args[0], var)) {
            if (call->func_id == BuiltinOp::Exp) {
                return ok(expr);
            }
            if (call->func_id == BuiltinOp::Ln) {
                ExprPtr x = arena.make<Symbol>(var);
                return ok(arena.make<Sum>(std::vector<ExprPtr>{
                    arena.make<Binary>(BinaryOp::Mul, x, expr),
                    arena.make<Unary>(UnaryOp::Neg, x)
                }));
            }
        }
    }

    // 1. Build Differential Extension Tower
    auto field_res = DifferentialField::build(expr, var, context);
    if (field_res.is_error()) return fail<ExprPtr>(field_res.error());
    const auto& field = field_res.value();

    // 2. Map expression to the differential field (generators t_1, ..., t_n)
    auto gen_expr_res = field.to_field_generators(expr, context);
    if (gen_expr_res.is_error()) return fail<ExprPtr>(gen_expr_res.error());
    ExprPtr gen_expr = gen_expr_res.value();

    // 3. Decompose into P/Q with respect to the topmost generator t_n
    // If no extensions, t = var
    Symbol t_top = field.extensions().empty() ? var : field.extensions().back().t_var;

    auto rational_res = algebra::apart_num_den(gen_expr, context);
    if (rational_res.is_error()) return fail<ExprPtr>(rational_res.error());
    ExprPtr P = rational_res.value().numerator;
    ExprPtr Q = rational_res.value().denominator;

    // 3a. Safeguard: if P or Q still contains the top extension variable as a symbol,
    //     the Rothstein-Trager algorithm would operate in the wrong ring (base variable
    //     instead of the extension variable). Bail out so the elementary integrator handles it.
    if (!field.extensions().empty()) {
        if (integrate_detail::depends_on(P, t_top) || integrate_detail::depends_on(Q, t_top)) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Risch: expression polynomial in transcendental extension — deferred to elementary integrator",
            });
        }
    }

    // 3b. Polynomial part (only valid for the base case: no transcendental extensions).
    //     For extension towers, the generator t_n is transcendental so ∫ t_n^k dt ≠ t_n^(k+1)/(k+1).
    //     ∫ P/Q = ∫ S + ∫ R/Q (S integrated by power rule in the integration variable x).
    ExprPtr poly_integral_part = arena.make<IntegerLit>(BigInt(0));
    if (field.extensions().empty()) {
        auto P_poly = algebra::parse_polynomial(P, t_top, context);
        auto Q_poly = algebra::parse_polynomial(Q, t_top, context);
        if (P_poly.is_ok() && Q_poly.is_ok() && !algebra::is_zero_poly(Q_poly.value())) {
            if (algebra::poly_degree(P_poly.value()) >= algebra::poly_degree(Q_poly.value())) {
                auto div_res = algebra::divide_poly_with_remainder(P_poly.value(), Q_poly.value(), context);
                if (div_res.is_ok()) {
                    const auto& quot = div_res.value().quotient;
                    const auto& rem  = div_res.value().remainder;
                    // Integrate quotient symbolically via the elementary integrator.
                    // Integrate polynomial quotient term-by-term using power rule:
                    // ∫ a_k * t^k dt = a_k / (k+1) * t^(k+1)
                    std::vector<ExprPtr> int_terms;
                    for (std::size_t k = 0; k < quot.size(); ++k) {
                        ExprPtr coeff = quot[k];
                        if (!coeff) continue;
                        bool coeff_zero = false;
                        if (const auto* il = expr_cast<IntegerLit>(coeff)) coeff_zero = il->value.is_zero();
                        if (coeff_zero) continue;
                        // k==0: ∫ c dt = c*t
                        if (k == 0) {
                            int_terms.push_back(arena.make<Binary>(BinaryOp::Mul, coeff, arena.make<Symbol>(t_top)));
                        } else {
                            // ∫ c*t^k dt = c/(k+1) * t^(k+1)
                            ExprPtr kp1 = arena.make<IntegerLit>(BigInt(static_cast<long long>(k + 1)));
                            ExprPtr t_pow = arena.make<Binary>(BinaryOp::Pow, arena.make<Symbol>(t_top), kp1);
                            ExprPtr term = arena.make<Binary>(BinaryOp::Div,
                                arena.make<Binary>(BinaryOp::Mul, coeff, t_pow),
                                kp1);
                            int_terms.push_back(term);
                        }
                    }
                    if (!int_terms.empty()) {
                        ExprPtr raw = int_terms.size() == 1
                            ? int_terms[0]
                            : arena.make<Sum>(std::move(int_terms));
                        auto simp = context.simplify(raw);
                        if (simp.is_ok()) poly_integral_part = simp.value();
                    }
                    // Replace P with the remainder for the rational part.
                    auto rem_expr = algebra::polynomial_to_expr(rem, t_top, context);
                    if (rem_expr.is_ok()) P = rem_expr.value();
                }
            }
        }
    } // end if (field.extensions().empty())

    // 4. Hermite Reduction: ∫ P/Q dt = A/B + ∫ C/D dt where D is square-free
    auto hermite_res = hermite_reduce(P, Q, t_top, field, context);
    if (hermite_res.is_error()) return fail<ExprPtr>(hermite_res.error());
    
    ExprPtr rational_part = hermite_res.value().rational_part;
    ExprPtr rem_P = hermite_res.value().remaining_P;
    ExprPtr rem_Q = hermite_res.value().remaining_Q;

    // 5. Rothstein-Trager for the remaining square-free part ∫ C/D dt
    auto log_part_res = integrate_rothstein_trager(rem_P, rem_Q, t_top, field, context);
    
    if (log_part_res.is_ok()) {
        // combine rational + log parts and map back from field generators
        ExprPtr total_gen;
        bool rat_zero = expr_is<IntegerLit>(rational_part) && expr_ref<IntegerLit>(rational_part).value.is_zero();
        bool log_zero = expr_is<IntegerLit>(log_part_res.value()) && expr_ref<IntegerLit>(log_part_res.value()).value.is_zero();
        if (rat_zero && log_zero) {
            total_gen = arena.make<IntegerLit>(BigInt(0));
        } else if (rat_zero) {
            total_gen = log_part_res.value();
        } else if (log_zero) {
            total_gen = rational_part;
        } else {
            total_gen = arena.make<Sum>(std::vector<ExprPtr>{rational_part, log_part_res.value()});
        }

        auto back_res = field.from_field_generators(total_gen, context);
        if (back_res.is_error()) return back_res;

        // Add polynomial integral part
        bool poly_zero = expr_is<IntegerLit>(poly_integral_part) && expr_ref<IntegerLit>(poly_integral_part).value.is_zero();
        if (poly_zero) return back_res;
        bool back_zero = expr_is<IntegerLit>(back_res.value()) && expr_ref<IntegerLit>(back_res.value()).value.is_zero();
        if (back_zero) return ok(poly_integral_part);
        return context.simplify(arena.make<Sum>(std::vector<ExprPtr>{poly_integral_part, back_res.value()}));
    }

    // Fallback to simpler cases if full Risch fails
    // (Old pattern matching logic could be reintegrated here if needed)
    return fail<ExprPtr>(CASError{
        .kind = CASErrorKind::Unimplemented,
        .message = "Risch algorithm: integrability could not be decided for this transcendental extension",
    });
}

} // namespace cas::calculus
