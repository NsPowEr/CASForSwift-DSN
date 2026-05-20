// CAS-L3-10 — ODE solver via Laplace transform.
//
// Solves linear ODE with constant coefficients and initial conditions:
//   a_n y^(n) + a_{n-1} y^(n-1) + ... + a_0 y = f(t),   y(0)=y0, y'(0)=y1, ...
//
// Algorithm:
//   1. Apply Laplace L{} to both sides.
//      L{y^(k)} = s^k · Y(s) - s^(k-1)·y(0) - s^(k-2)·y'(0) - ... - y^(k-1)(0)
//   2. Solve algebraically for Y(s).
//   3. Apply inverse Laplace to recover y(t).
//
// MVP: order 1 and 2 with zero initial conditions (most common textbook case).
// Higher orders or non-zero ICs: pass ICs vector explicitly.

#include "calculus_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"

#include <vector>

namespace cas::calculus {

// Solve ODE: coeffs[i] · y^(i) = f(t) per i=0..n.
// initial_conditions[k] = y^(k)(0), k=0..n-1 (zero if not provided).
// Returns y(t) in closed form, or Unimplemented if Laplace pipeline fails.
[[nodiscard]] Result<ExprPtr> solve_ode_laplace(
    const std::vector<ExprPtr>& coeffs,
    ExprPtr forcing,
    const std::vector<ExprPtr>& initial_conditions,
    const Symbol& t,
    symbolic::CASContext& ctx) {
    if (coeffs.size() < 2U) {
        return fail<ExprPtr>(CASError{
            CASErrorKind::InvalidArgument,
            "solve_ode_laplace: requires order ≥ 1 (coeffs.size() ≥ 2)",
            std::nullopt});
    }
    AstArena& arena = ctx.arena();
    const std::size_t n = coeffs.size() - 1;  // ODE order
    Symbol s = ctx.make_fresh_symbol("s");
    ExprPtr s_expr = arena.make<Symbol>(s);

    // Apply Laplace to f(t): F(s) = L{f(t)}.
    auto F_res = laplace_transform(forcing, t, s, ctx);
    if (F_res.is_error()) return fail<ExprPtr>(F_res.error());
    ExprPtr F = F_res.value();

    // Build L{LHS}:
    //   Σ_{k=0..n} coeffs[k] · (s^k · Y - Σ_{j=0..k-1} s^(k-1-j) · y^(j)(0))
    //   = Y · Σ coeffs[k]·s^k  -  Σ_{k≥1} coeffs[k] · Σ_{j<k} s^(k-1-j) · y^(j)(0)
    //
    // Characteristic polynomial: P(s) = Σ coeffs[k]·s^k
    // Initial-condition polynomial: Q(s) = Σ_{k=1..n} coeffs[k] ·
    //                                       Σ_{j=0..k-1} s^(k-1-j) · IC[j]

    ExprPtr P = arena.make<IntegerLit>(BigInt(0));
    for (std::size_t k = 0; k <= n; ++k) {
        ExprPtr term;
        if (k == 0U) {
            term = coeffs[k];
        } else if (k == 1U) {
            term = arena.make<Product>(std::vector<ExprPtr>{coeffs[k], s_expr});
        } else {
            ExprPtr s_pow_k = arena.make<Binary>(BinaryOp::Pow, s_expr,
                arena.make<IntegerLit>(BigInt(static_cast<long long>(k))));
            term = arena.make<Product>(std::vector<ExprPtr>{coeffs[k], s_pow_k});
        }
        P = arena.make<Binary>(BinaryOp::Add, P, term);
    }

    ExprPtr Q_ic = arena.make<IntegerLit>(BigInt(0));
    for (std::size_t k = 1; k <= n; ++k) {
        for (std::size_t j = 0; j < k; ++j) {
            ExprPtr ic = (j < initial_conditions.size())
                ? initial_conditions[j]
                : static_cast<ExprPtr>(arena.make<IntegerLit>(BigInt(0)));
            // s^(k-1-j)
            const long long exp_val = static_cast<long long>(k - 1 - j);
            ExprPtr s_pow;
            if (exp_val == 0) s_pow = arena.make<IntegerLit>(BigInt(1));
            else if (exp_val == 1) s_pow = s_expr;
            else s_pow = arena.make<Binary>(BinaryOp::Pow, s_expr,
                arena.make<IntegerLit>(BigInt(exp_val)));
            ExprPtr term = arena.make<Product>(std::vector<ExprPtr>{
                coeffs[k], s_pow, ic});
            Q_ic = arena.make<Binary>(BinaryOp::Add, Q_ic, term);
        }
    }

    // Simplify P and Q_ic to clean form before assembly.
    auto P_simp = ctx.simplify(P);
    if (P_simp.is_ok()) P = P_simp.value();
    auto Q_simp = ctx.simplify(Q_ic);
    if (Q_simp.is_ok()) Q_ic = Q_simp.value();

    // P(s) · Y(s) - Q_ic(s) = F(s)  →  Y(s) = (F(s) + Q_ic(s)) / P(s)
    ExprPtr Y_num = arena.make<Binary>(BinaryOp::Add, F, Q_ic);
    ExprPtr Y_s = arena.make<Binary>(BinaryOp::Div, Y_num, P);
    auto Y_tog = algebra::together(Y_s, ctx);
    if (Y_tog.is_ok()) Y_s = Y_tog.value();
    auto Y_simp = ctx.simplify(Y_s);
    if (Y_simp.is_ok()) Y_s = Y_simp.value();

    // Inverse Laplace.
    auto y_res = inverse_laplace_transform(Y_s, s, t, ctx);
    if (y_res.is_error()) {
        return fail<ExprPtr>(CASError{
            CASErrorKind::Unimplemented,
            "solve_ode_laplace: inverse Laplace failed on Y(s) — "
            + y_res.error().message,
            std::nullopt});
    }
    return ctx.simplify(y_res.value());
}

}  // namespace cas::calculus
