// L2-11: Improper integral convergence classification and Cauchy principal
// value for simple-pole rational integrands.
//
// See include/cas/improper_integral.hpp for the public contract.

#include "cas/improper_integral.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/symbolic.hpp"

#include <string>
#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] ExprPtr int_lit(AstArena& arena, long long value) {
    return arena.make<IntegerLit>(BigInt(value));
}

[[nodiscard]] bool is_pos_infinity(ExprPtr expr) {
    const auto* c = expr_cast<Constant>(expr);
    return c != nullptr && c->value == MathConstant::Infinity;
}

[[nodiscard]] bool is_neg_infinity(ExprPtr expr) {
    const auto* u = expr_cast<Unary>(expr);
    return u != nullptr && u->op == UnaryOp::Neg && is_pos_infinity(u->operand);
}

[[nodiscard]] bool is_infinity(ExprPtr expr) {
    return is_pos_infinity(expr) || is_neg_infinity(expr);
}

[[nodiscard]] bool is_zero_literal(ExprPtr expr) {
    if (const auto* lit = expr_cast<IntegerLit>(expr)) {
        return lit->value.is_zero();
    }
    if (const auto* lit = expr_cast<RationalLit>(expr)) {
        return lit->numerator.is_zero();
    }
    return false;
}

// Scan up to `scan_window` positive-order coefficients of a Laurent
// expansion to find the first non-zero one and return the effective leading
// order.  Needed because `laurent_series` reports leading_order = -m where m
// is the pole order of the denominator: when the numerator also vanishes at
// the center, the true effective order is larger.  HC-006 closed: window is
// configurable via ctx.improper_leading_order_scan() (default 8).
[[nodiscard]] Result<int> effective_leading_order(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr center,
    symbolic::CASContext& ctx,
    unsigned int scan_window = 0U) {
    if (scan_window == 0U) {
        scan_window = static_cast<unsigned int>(ctx.improper_leading_order_scan());
    }
    auto laurent = laurent_series(expr, var, center, scan_window, ctx);
    if (laurent.is_error()) {
        return fail<int>(laurent.error());
    }
    const auto& exp = laurent.value();
    for (std::size_t idx = 0U; idx < exp.coefficients.size(); ++idx) {
        auto simplified = ctx.simplify(exp.coefficients[idx]);
        if (simplified.is_error()) {
            return fail<int>(simplified.error());
        }
        if (!is_zero_literal(simplified.value())) {
            return ok(exp.leading_order + static_cast<int>(idx));
        }
    }
    // All scanned coefficients vanished: treat as "very analytic" (>= +infty).
    // Return scan_window + 1 to signal "convergent" for any test.
    return ok(static_cast<int>(scan_window) + 1);
}

// Compute the leading order of `expr` (as a rational function of `var`)
// at the finite point `point`.
[[nodiscard]] Result<int> leading_order_at_finite(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    symbolic::CASContext& ctx) {
    return effective_leading_order(expr, var, point, ctx);
}

// Compute the asymptotic leading order of `expr` as |x| -> infinity.
// Substitute x = 1/u; the order of `expr(1/u)` at u=0 (call it L) yields
// asymptotic order of expr at infinity = -L.
//
// Example: 1/(1+x^2)  with x=1/u  =>  u^2/(u^2+1).  L = +2  =>  order at oo = -2.
// Example: 1/x        with x=1/u  =>  u.            L = +1  =>  order at oo = -1.
[[nodiscard]] Result<int> leading_order_at_infinity(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    // Build 1/u using a fresh symbol distinct from var.
    std::string fresh_name = std::string{"__improper_u_"} + var.name;
    ExprPtr u = arena.make<Symbol>(fresh_name);
    ExprPtr one_over_u = arena.make<Binary>(BinaryOp::Div, int_lit(arena, 1), u);

    auto substituted = symbolic::substitute(expr, var, one_over_u, ctx);
    if (substituted.is_error()) {
        return fail<int>(substituted.error());
    }
    auto simplified = ctx.simplify(substituted.value());
    if (simplified.is_error()) {
        return fail<int>(simplified.error());
    }
    auto together = algebra::together(simplified.value(), ctx);
    ExprPtr ready = together.is_ok() ? together.value() : simplified.value();

    auto leading = effective_leading_order(ready, Symbol(fresh_name), int_lit(arena, 0), ctx);
    if (leading.is_error()) {
        return fail<int>(leading.error());
    }
    return ok(-leading.value());
}

[[nodiscard]] std::string format_endpoint(ExprPtr expr) {
    if (is_pos_infinity(expr)) return "+oo";
    if (is_neg_infinity(expr)) return "-oo";
    return "finite";
}

[[nodiscard]] bool endpoint_convergent(int order, bool at_infinity) {
    // At a finite endpoint: convergent iff leading_order >= 0
    //   (analytic or zero at endpoint - integrable).
    //   Actually for finite endpoint we need integrability of (x-a)^k near 0:
    //     k >= 0 -> integrable;
    //     k = -1 -> log divergence;
    //     k <= -2 -> algebraic divergence.
    // At +/- infinity: convergent iff order <= -2
    //   (∫_1^oo x^k dx converges iff k < -1, i.e. k <= -2 in integer order).
    if (at_infinity) {
        return order <= -2;
    }
    return order >= 0;
}

}  // namespace

Result<ConvergenceReport> classify_improper_convergence(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx) {
    // Normalize integrand to a single rational form.
    auto together = algebra::together(expr, ctx);
    ExprPtr integrand = together.is_ok() ? together.value() : expr;
    auto simplified = ctx.simplify(integrand);
    if (simplified.is_ok()) integrand = simplified.value();

    int order_lower = 0;
    int order_upper = 0;

    auto compute_lower = is_infinity(lower)
        ? leading_order_at_infinity(integrand, var, ctx)
        : leading_order_at_finite(integrand, var, lower, ctx);
    if (compute_lower.is_error()) {
        return fail<ConvergenceReport>(compute_lower.error());
    }
    order_lower = compute_lower.value();

    auto compute_upper = is_infinity(upper)
        ? leading_order_at_infinity(integrand, var, ctx)
        : leading_order_at_finite(integrand, var, upper, ctx);
    if (compute_upper.is_error()) {
        return fail<ConvergenceReport>(compute_upper.error());
    }
    order_upper = compute_upper.value();

    const bool lower_conv = endpoint_convergent(order_lower, is_infinity(lower));
    const bool upper_conv = endpoint_convergent(order_upper, is_infinity(upper));

    ConvergenceReport report;
    report.leading_order_at_lower = order_lower;
    report.leading_order_at_upper = order_upper;

    if (lower_conv && upper_conv) {
        report.status = ConvergenceStatus::Convergent;
        report.diagnostic = "Convergent at both endpoints (lower order="
            + std::to_string(order_lower) + ", upper order="
            + std::to_string(order_upper) + ").";
    } else if (!lower_conv && !upper_conv) {
        report.status = ConvergenceStatus::DivergentAtLowerEnd;
        report.diagnostic = "Divergent at both endpoints (lower at "
            + format_endpoint(lower) + " order=" + std::to_string(order_lower)
            + "; upper at " + format_endpoint(upper) + " order="
            + std::to_string(order_upper) + ").";
    } else if (!lower_conv) {
        report.status = ConvergenceStatus::DivergentAtLowerEnd;
        report.diagnostic = "Divergent at lower endpoint (" + format_endpoint(lower)
            + ", leading order " + std::to_string(order_lower) + ").";
    } else {
        report.status = ConvergenceStatus::DivergentAtUpperEnd;
        report.diagnostic = "Divergent at upper endpoint (" + format_endpoint(upper)
            + ", leading order " + std::to_string(order_upper) + ").";
    }
    return ok(report);
}

Result<ExprPtr> cauchy_principal_value(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    ExprPtr pole,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    // Step 1: Laurent expansion at the pole with positive_order = 0.
    auto laurent = laurent_series(expr, var, pole, 0U, ctx);
    if (laurent.is_error()) {
        return fail<ExprPtr>(laurent.error());
    }
    const auto& exp = laurent.value();

    if (exp.leading_order > 0) {
        // Integrand is analytic and vanishes at `pole`: nothing special.
        return definite_integral(expr, var, lower, upper, ctx);
    }
    if (exp.leading_order >= 0) {
        // Integrand analytic at pole — no special handling needed.
        return definite_integral(expr, var, lower, upper, ctx);
    }

    // General Hadamard finite‑part regularisation for a pole of order m ≥ 1.
    // Let f(x) = Σ_{k=-m}^∞ c_k · (x − p)^k  near p.  Write
    //   f(x) = singular(x) + regular(x)
    //   singular(x) = Σ_{k=-m}^{-1} c_k · (x − p)^k
    //   regular(x)  = f(x) − singular(x)   (analytic at p)
    //
    //   ∫_a^b regular(x) dx = F(b) − F(a)            (standard antiderivative)
    //
    //   For the singular part, evaluated as the limit
    //     lim_{ε→0+} [ ∫_a^{p−ε} singular + ∫_{p+ε}^b singular − S(ε) ]
    //   where S(ε) is the diverging tail in ε:
    //
    //     k = −1:    c_{-1} · (ln|b − p| − ln|a − p|)             (PV log)
    //     k ≤ −2:    c_k/(k+1) · ((b − p)^{k+1} − (a − p)^{k+1})
    //                 (Hadamard finite part — divergent symmetric tails
    //                  −2 c_k ε^{k+1}/(k+1) cancel when k is odd, and are
    //                  subtracted by definition when k is even).
    //
    // Both formulas merge cleanly because for k ≤ −2 the bracket
    // (b−p)^{k+1} − (a−p)^{k+1} is exactly the finite part, regardless of
    // parity.  The k = −1 case is special only because (k+1) = 0 forbids
    // the same antiderivative form and the log appears instead.
    const int m = -exp.leading_order;  // m ≥ 1 by the branch above

    // Step 2: build the singular part Σ_{k=-m}^{-1} c_k · (x − p)^k.
    ExprPtr x_minus_p = arena.make<Binary>(BinaryOp::Sub, arena.make<Symbol>(var), pole);
    std::vector<ExprPtr> singular_terms;
    singular_terms.reserve(static_cast<std::size_t>(m));
    for (int k = -m; k <= -1; ++k) {
        const std::size_t idx = static_cast<std::size_t>(k - exp.leading_order);
        if (idx >= exp.coefficients.size()) continue;
        ExprPtr c_k = exp.coefficients[idx];
        // (x − p)^k for k < 0 ⇒ 1 / (x − p)^{-k}
        ExprPtr power_arg = arena.make<IntegerLit>(BigInt(static_cast<long long>(-k)));
        ExprPtr denom_pow = arena.make<Binary>(BinaryOp::Pow, x_minus_p, power_arg);
        ExprPtr term = arena.make<Binary>(BinaryOp::Div, c_k, denom_pow);
        singular_terms.push_back(term);
    }
    ExprPtr singular_expr = arena.make<IntegerLit>(BigInt(0));
    if (singular_terms.size() == 1U) {
        singular_expr = singular_terms.front();
    } else if (singular_terms.size() > 1U) {
        singular_expr = arena.make<Sum>(std::move(singular_terms));
    }
    auto singular_simp = ctx.simplify(singular_expr);
    if (singular_simp.is_error()) return singular_simp;

    // Step 3: regular(x) = f(x) − singular(x).
    ExprPtr regular_expr = arena.make<Binary>(BinaryOp::Sub, expr, singular_simp.value());
    auto regular_simplified = ctx.simplify(regular_expr);
    if (regular_simplified.is_error()) return regular_simplified;
    regular_expr = regular_simplified.value();

    // Step 4a: integrate the regular part and evaluate at the endpoints.
    auto antideriv = integrate(regular_expr, var, ctx);
    if (antideriv.is_error()) return antideriv;
    auto F_upper = ctx.substitute(antideriv.value(), var, upper);
    if (F_upper.is_error()) return F_upper;
    auto F_lower = ctx.substitute(antideriv.value(), var, lower);
    if (F_lower.is_error()) return F_lower;
    ExprPtr regular_part = arena.make<Binary>(BinaryOp::Sub, F_upper.value(), F_lower.value());

    // Step 4b: singular‑part finite‑part contributions, k = −1..−m.
    ExprPtr upper_minus_p = arena.make<Binary>(BinaryOp::Sub, upper, pole);
    ExprPtr lower_minus_p = arena.make<Binary>(BinaryOp::Sub, lower, pole);

    std::vector<ExprPtr> singular_contributions;
    for (int k = -1; k >= -m; --k) {
        const std::size_t idx = static_cast<std::size_t>(k - exp.leading_order);
        if (idx >= exp.coefficients.size()) continue;
        ExprPtr c_k = exp.coefficients[idx];
        if (k == -1) {
            ExprPtr abs_upper = arena.make<FuncCall>(BuiltinOp::Abs, std::vector<ExprPtr>{upper_minus_p});
            ExprPtr abs_lower = arena.make<FuncCall>(BuiltinOp::Abs, std::vector<ExprPtr>{lower_minus_p});
            ExprPtr ln_upper = arena.make<FuncCall>(BuiltinOp::Log, std::vector<ExprPtr>{abs_upper});
            ExprPtr ln_lower = arena.make<FuncCall>(BuiltinOp::Log, std::vector<ExprPtr>{abs_lower});
            ExprPtr log_diff = arena.make<Binary>(BinaryOp::Sub, ln_upper, ln_lower);
            singular_contributions.push_back(arena.make<Binary>(BinaryOp::Mul, c_k, log_diff));
        } else {
            // c_k/(k+1) · ((b − p)^{k+1} − (a − p)^{k+1}),  k+1 < 0.
            const long long expnt = static_cast<long long>(k + 1);
            ExprPtr expnt_lit = arena.make<IntegerLit>(BigInt(expnt));
            ExprPtr ub_pow = arena.make<Binary>(BinaryOp::Pow, upper_minus_p, expnt_lit);
            ExprPtr lb_pow = arena.make<Binary>(BinaryOp::Pow, lower_minus_p, expnt_lit);
            ExprPtr diff_pow = arena.make<Binary>(BinaryOp::Sub, ub_pow, lb_pow);
            ExprPtr k_plus_1_lit = arena.make<IntegerLit>(BigInt(expnt));
            ExprPtr scaled = arena.make<Binary>(BinaryOp::Div, c_k, k_plus_1_lit);
            singular_contributions.push_back(arena.make<Binary>(BinaryOp::Mul, scaled, diff_pow));
        }
    }
    ExprPtr singular_part = arena.make<IntegerLit>(BigInt(0));
    if (singular_contributions.size() == 1U) {
        singular_part = singular_contributions.front();
    } else if (singular_contributions.size() > 1U) {
        singular_part = arena.make<Sum>(std::move(singular_contributions));
    }

    ExprPtr total = arena.make<Binary>(BinaryOp::Add, regular_part, singular_part);
    return ctx.simplify(total);
}

}  // namespace cas::calculus
