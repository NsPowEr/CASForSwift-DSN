#include "calculus_internal.hpp"
#include "integrate_engine.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/algebra.hpp"
#include "cas/error.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace cas::calculus {

namespace {

// Integrate the polynomial-in-t part of a single logarithmic extension
// tower:  ∫ Σ_{k=0..n} a_k(x) * t^k dx,   where t = ln(u(x)),  Dt = u'/u.
//
// Standard ansatz:  the antiderivative is again a polynomial in t,
//   B(t) = Σ_{k=0..n} b_k(x) * t^k,   with
//
//   d/dx B(t)
//     = Σ b_k'(x) * t^k + Σ k * b_k(x) * (u'/u) * t^{k-1}
//     = b_n' * t^n + Σ_{k=0..n-1} [ b_k' + (k+1) * b_{k+1} * (u'/u) ] * t^k.
//
// Matching coefficients with a_k * t^k gives the descending recursion
//   b_n = ∫ a_n dx
//   b_k = ∫ [ a_k - (k+1) * b_{k+1} * (u'/u) ] dx,   for k = n-1, n-2, ..., 0.
//
// Each integration is carried out in the lower field (here Q(x) for a single
// log extension) via the existing integrate() routine.  Failures propagate
// as Unimplemented so that the caller can fall back to other strategies.
[[nodiscard]] Result<ExprPtr> integrate_log_polynomial_part(
    const algebra::PolyExpr& quot,
    ExprPtr u_arg,
    const Symbol& t_top,
    const Symbol& var,
    symbolic::CASContext& context) {
    AstArena& arena = context.arena();

    if (quot.empty()) {
        return ok(arena.make<IntegerLit>(BigInt(0)));
    }
    const std::size_t deg = quot.size() - 1U;

    // u'/u (simplified once up front)
    auto du_res = diff(u_arg, var, 1U, context);
    if (du_res.is_error()) return fail<ExprPtr>(du_res.error());
    ExprPtr du_over_u = arena.make<Binary>(BinaryOp::Div, du_res.value(), u_arg);
    if (auto s = context.simplify(du_over_u); s.is_ok()) du_over_u = s.value();

    std::vector<ExprPtr> b(deg + 1U, ExprPtr{});

    for (std::ptrdiff_t k = static_cast<std::ptrdiff_t>(deg); k >= 0; --k) {
        const std::size_t kz = static_cast<std::size_t>(k);
        ExprPtr a_k = (kz < quot.size()) ? quot[kz] : ExprPtr{};
        if (!a_k) a_k = arena.make<IntegerLit>(BigInt(0));

        // rhs = a_k  -  (k+1) * b_{k+1} * (u'/u)
        ExprPtr rhs = a_k;
        if (kz + 1U <= deg && b[kz + 1U]) {
            ExprPtr kp1 = arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(kz + 1U)));
            ExprPtr correction = arena.make<Product>(std::vector<ExprPtr>{kp1, b[kz + 1U], du_over_u});
            rhs = arena.make<Binary>(BinaryOp::Sub, rhs, correction);
        }
        if (auto s = context.simplify(rhs); s.is_ok()) rhs = s.value();

        auto b_k_res = integrate(rhs, var, context);
        if (b_k_res.is_error()) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Risch log: lower-field integration failed at degree k=" + std::to_string(kz),
                .hint = std::nullopt});
        }
        b[kz] = b_k_res.value();
    }

    // Build B(t) = Σ b_k * t^k.
    std::vector<ExprPtr> terms;
    terms.reserve(b.size());
    ExprPtr t_sym = arena.make<Symbol>(t_top.name);
    for (std::size_t k = 0; k < b.size(); ++k) {
        if (!b[k]) continue;
        if (const auto* il = expr_cast<IntegerLit>(b[k]); il && il->value.is_zero()) continue;
        ExprPtr term;
        if (k == 0U) {
            term = b[k];
        } else if (k == 1U) {
            term = arena.make<Binary>(BinaryOp::Mul, b[k], t_sym);
        } else {
            ExprPtr t_pow = arena.make<Binary>(
                BinaryOp::Pow,
                t_sym,
                arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(k))));
            term = arena.make<Binary>(BinaryOp::Mul, b[k], t_pow);
        }
        terms.push_back(term);
    }
    if (terms.empty()) return ok(arena.make<IntegerLit>(BigInt(0)));
    if (terms.size() == 1U) return ok(terms.front());
    ExprPtr raw = arena.make<Sum>(std::move(terms));
    if (auto s = context.simplify(raw); s.is_ok()) return ok(s.value());
    return ok(raw);
}

}  // namespace

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
    Symbol t_top = field.extensions().empty() ? var : field.extensions().back().t_var;

    auto rational_res = algebra::apart_num_den(gen_expr, context);
    if (rational_res.is_error()) return fail<ExprPtr>(rational_res.error());
    ExprPtr P = rational_res.value().numerator;
    ExprPtr Q = rational_res.value().denominator;

    // 3b. Polynomial part
    ExprPtr poly_integral_part = arena.make<IntegerLit>(BigInt(0));
    auto P_poly = algebra::parse_polynomial(P, t_top, context);
    auto Q_poly = algebra::parse_polynomial(Q, t_top, context);

    if (P_poly.is_ok() && Q_poly.is_ok() && !algebra::is_zero_poly(Q_poly.value())) {
        if (algebra::poly_degree(P_poly.value()) >= algebra::poly_degree(Q_poly.value())) {
            auto div_res = algebra::divide_poly_with_remainder(P_poly.value(), Q_poly.value(), context);
            if (div_res.is_ok()) {
                const auto& quot = div_res.value().quotient;
                const auto& rem  = div_res.value().remainder;
                
                std::vector<ExprPtr> int_terms;
                if (field.extensions().empty()) {
                    // Base case: t_top = x, D(t_top) = 1
                    for (std::size_t k = 0; k < quot.size(); ++k) {
                        ExprPtr coeff = quot[k];
                        if (!coeff || algebra::poly_is_zero_expr(coeff)) continue;
                        if (k == 0) {
                            int_terms.push_back(arena.make<Binary>(BinaryOp::Mul, coeff, arena.make<Symbol>(t_top)));
                        } else {
                            ExprPtr kp1 = arena.make<IntegerLit>(BigInt(static_cast<long long>(k + 1)));
                            ExprPtr t_pow = arena.make<Binary>(BinaryOp::Pow, arena.make<Symbol>(t_top), kp1);
                            int_terms.push_back(arena.make<Binary>(BinaryOp::Div, arena.make<Binary>(BinaryOp::Mul, coeff, t_pow), kp1));
                        }
                    }
                } else {
                    // Transcendental case
                    const auto& ext = field.extensions().back();
                    const bool handle_log_polynomial = (ext.type == ExtensionType::Logarithmic);

                    // Logarithmic extension: descending recursion on B(t),
                    // see integrate_log_polynomial_part.  Performs the full
                    // polynomial-in-t integration in one structured pass.
                    if (handle_log_polynomial) {
                        auto B_res = integrate_log_polynomial_part(quot, ext.argument, t_top, var, context);
                        if (B_res.is_error()) return fail<ExprPtr>(B_res.error());
                        int_terms.push_back(B_res.value());
                    }

                    for (std::size_t k = 0; !handle_log_polynomial && k < quot.size(); ++k) {
                        ExprPtr coeff = quot[k];
                        if (!coeff || algebra::poly_is_zero_expr(coeff)) continue;

                        if (ext.type == ExtensionType::Exponential) {
                            // ∫ a_k * t^k where t = exp(u), Dt = u't
                            if (k == 0) {
                                auto base_int = integrate(coeff, var, context);
                                if (base_int.is_ok()) int_terms.push_back(base_int.value());
                                else return base_int;
                            } else {
                                // solve Dy + k*u'*y = a_k
                                auto du_res = diff(ext.argument, var, 1U, context);
                                if (du_res.is_error()) return fail<ExprPtr>(du_res.error());
                                ExprPtr du = du_res.value();
                                
                                auto solve_risch_de = [&](ExprPtr f, ExprPtr g) -> Result<ExprPtr> {
                                    auto g_poly = algebra::parse_polynomial(g, var, context);
                                    if (g_poly.is_ok()) {
                                        std::vector<ExprPtr> y_coeffs(g_poly.value().size());
                                        for (int i = static_cast<int>(g_poly.value().size()) - 1; i >= 0; --i) {
                                            ExprPtr rhs = g_poly.value()[i];
                                            if (i < static_cast<int>(g_poly.value().size()) - 1) {
                                                ExprPtr next_c = arena.make<Binary>(BinaryOp::Mul,
                                                    arena.make<IntegerLit>(BigInt(i + 1)), y_coeffs[i+1]);
                                                auto rhs_sub = context.simplify(arena.make<Binary>(BinaryOp::Sub, rhs, next_c));
                                                if (rhs_sub.is_error()) return fail<ExprPtr>(rhs_sub.error());
                                                rhs = rhs_sub.value();
                                            }
                                            auto div_res_coeff = context.simplify(arena.make<Binary>(BinaryOp::Div, rhs, f));
                                            if (div_res_coeff.is_error()) return fail<ExprPtr>(div_res_coeff.error());
                                            // Verify result is polynomial: Risch DE has polynomial solution only if each coeff is polynomial
                                            auto poly_check = algebra::parse_polynomial(div_res_coeff.value(), var, context);
                                            if (poly_check.is_error()) {
                                                return fail<ExprPtr>(CASError{CASErrorKind::Unimplemented,
                                                    "Risch DE: no polynomial solution (coefficient is not polynomial)", std::nullopt});
                                            }
                                            y_coeffs[i] = div_res_coeff.value();
                                        }
                                        return algebra::polynomial_to_expr(algebra::PolyExpr(y_coeffs), var, context);
                                    }
                                    return fail<ExprPtr>(CASError{CASErrorKind::Unimplemented, "Risch DE: non-poly g", std::nullopt});
                                };

                                auto y_res = solve_risch_de(arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(k)), du), coeff);
                                if (y_res.is_ok()) {
                                    ExprPtr t_pow = (k == 1) ? arena.make<Symbol>(t_top) : arena.make<Binary>(BinaryOp::Pow, arena.make<Symbol>(t_top), arena.make<IntegerLit>(BigInt(k)));
                                    int_terms.push_back(arena.make<Binary>(BinaryOp::Mul, y_res.value(), t_pow));
                                } else {
                                    return fail<ExprPtr>(CASError{CASErrorKind::Unimplemented, "Risch: could not solve DE for exponential term", std::nullopt});
                                }
                            }
                        } else if (ext.type == ExtensionType::Logarithmic) {
                            if (k == 1 && integrate_detail::is_same_symbol(ext.argument, var) && integrate_detail::is_one(coeff)) {
                                int_terms.push_back(arena.make<Sum>(std::vector<ExprPtr>{
                                    arena.make<Binary>(BinaryOp::Mul, arena.make<Symbol>(var), arena.make<Symbol>(t_top)),
                                    arena.make<Unary>(UnaryOp::Neg, arena.make<Symbol>(var))
                                }));
                            } else {
                                return fail<ExprPtr>(CASError{CASErrorKind::Unimplemented, "Risch: log extension integration not fully implemented", std::nullopt});
                            }
                        }
                    }
                }
                
                if (!int_terms.empty()) {
                    ExprPtr raw = int_terms.size() == 1 ? int_terms[0] : arena.make<Sum>(std::move(int_terms));
                    auto simp = context.simplify(raw);
                    if (simp.is_ok()) poly_integral_part = simp.value();
                }
                auto rem_expr = algebra::polynomial_to_expr(rem, t_top, context);
                if (rem_expr.is_ok()) P = rem_expr.value();
            }
        }
    }

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
    return fail<ExprPtr>(CASError{
        .kind = CASErrorKind::Unimplemented,
        .message = "Risch algorithm: integrability could not be decided for this transcendental extension",
        .hint = std::nullopt
    });
}

} // namespace cas::calculus
