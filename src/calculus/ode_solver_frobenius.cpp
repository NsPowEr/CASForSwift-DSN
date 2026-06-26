// L2-01: Frobenius series solution for a regular singular linear 2nd-order ODE.
#include "ode_solver_frobenius_internal.hpp"

#include "cas/ode.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/bigint.hpp"
#include "cas/calculus.hpp"
#include "cas/error.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus {

CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}


BigInt factorial_big(unsigned int n) {
    BigInt result(1);
    for (unsigned int i = 2; i <= n; ++i) {
        result = result * BigInt(static_cast<long long>(i));
    }
    return result;
}

// Test whether a simplified expression is the literal integer zero.
bool is_literal_zero(ExprPtr e) {
    if (const auto* il = expr_cast<IntegerLit>(e)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(e)) return rl->numerator.is_zero();
    return false;
}

// Detect a NaN / Infinity leaking out of simplification (indicates we tried to
// substitute into something singular at x = 0).
bool contains_undefined_constant(ExprPtr e) {
    if (!e) return false;
    if (const auto* c = expr_cast<Constant>(e)) {
        return c->value == MathConstant::NaN || c->value == MathConstant::Infinity;
    }
    switch (e->kind) {
    case ExprKind::Unary:
        return contains_undefined_constant(expr_ref<Unary>(e).operand);
    case ExprKind::Binary: {
        const auto& b = expr_ref<Binary>(e);
        return contains_undefined_constant(b.left) || contains_undefined_constant(b.right);
    }
    case ExprKind::Sum: {
        for (auto t : expr_ref<Sum>(e).terms)
            if (contains_undefined_constant(t)) return true;
        return false;
    }
    case ExprKind::Product: {
        for (auto f : expr_ref<Product>(e).factors)
            if (contains_undefined_constant(f)) return true;
        return false;
    }
    case ExprKind::FuncCall: {
        for (auto a : expr_ref<FuncCall>(e).args)
            if (contains_undefined_constant(a)) return true;
        return false;
    }
    default:
        return false;
    }
}

// Extract the k-th Taylor coefficient of `expr` around x = 0:
//   c_k = (1/k!) * d^k(expr)/dx^k  evaluated at x = 0
Result<ExprPtr> taylor_coefficient(
    ExprPtr expr,
    const Symbol& x,
    unsigned int k,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    Result<ExprPtr> deriv = (k == 0U) ? ok(expr) : diff(expr, x, k, ctx);
    if (deriv.is_error()) return fail<ExprPtr>(deriv.error());

    auto at_zero = ctx.substitute(deriv.value(), x, make_int(arena, 0));
    if (at_zero.is_error()) return fail<ExprPtr>(at_zero.error());

    auto simplified = ctx.simplify(at_zero.value());
    if (simplified.is_error()) return fail<ExprPtr>(simplified.error());

    if (k == 0U) return ok(simplified.value());

    ExprPtr denom = arena.make<IntegerLit>(factorial_big(k));
    ExprPtr scaled = arena.make<Binary>(BinaryOp::Div, simplified.value(), denom);
    return ctx.simplify(scaled);
}

// Build x^r (handling integer / rational / symbolic r without unnecessary noise).
ExprPtr make_x_to_r(ExprPtr r, const Symbol& x, AstArena& arena) {
    ExprPtr x_sym = arena.make<Symbol>(x.name);
    if (const auto* il = expr_cast<IntegerLit>(r); il != nullptr && il->value.is_zero()) {
        return arena.make<IntegerLit>(BigInt(1));
    }
    if (const auto* il = expr_cast<IntegerLit>(r); il != nullptr && il->value == BigInt(1)) {
        return x_sym;
    }
    return arena.make<Binary>(BinaryOp::Pow, x_sym, r);
}

// I(r) = r(r-1) + p0*r + q0 = r^2 + (p0 - 1)*r + q0
Result<ExprPtr> indicial_value(
    ExprPtr p0,
    ExprPtr q0,
    ExprPtr r_arg,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr r_minus_1 = arena.make<Binary>(BinaryOp::Sub, r_arg, make_int(arena, 1));
    ExprPtr r_times_rm1 = arena.make<Binary>(BinaryOp::Mul, r_arg, r_minus_1);
    ExprPtr p0_r = arena.make<Binary>(BinaryOp::Mul, p0, r_arg);
    ExprPtr sum = arena.make<Sum>(std::vector<ExprPtr>{r_times_rm1, p0_r, q0});
    return ctx.simplify(sum);
}

// Detect r1 - r2 = N with N a strictly positive integer.  Returns N on
// success; nullopt otherwise (non-integer gap, zero gap, negative gap, or
// an unsimplifiable expression).
std::optional<unsigned int> integer_gap(
    ExprPtr r1, ExprPtr r2, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    auto diff_res = ctx.simplify(arena.make<Binary>(BinaryOp::Sub, r1, r2));
    if (diff_res.is_error()) return std::nullopt;
    const auto* il = expr_cast<IntegerLit>(diff_res.value());
    if (il == nullptr) return std::nullopt;
    if (il->value.is_negative() || il->value.is_zero()) return std::nullopt;
    // Bound to a sane range (boundary cast, not symbolic arithmetic).
    std::uint64_t v = il->value.to_u64();
    if (v == 0U || v > 1000U) return std::nullopt;
    return static_cast<unsigned int>(v);
}

Result<ExprPtr> solve_ode_frobenius_at_zero(
    ExprPtr a_2,
    ExprPtr a_1,
    ExprPtr a_0,
    const Symbol& x,
    unsigned int num_terms,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    if (!a_2 || !a_1 || !a_0) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument,
                                        "Frobenius: null coefficient passed."));
    }

    // p(x) = a_1 / a_2,  q(x) = a_0 / a_2
    auto p_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Div, a_1, a_2));
    if (p_raw.is_error()) return fail<ExprPtr>(p_raw.error());
    auto q_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Div, a_0, a_2));
    if (q_raw.is_error()) return fail<ExprPtr>(q_raw.error());

    ExprPtr x_sym = arena.make<Symbol>(x.name);

    // p~(x) = x * p(x),  q~(x) = x^2 * q(x).  These MUST be analytic at 0.
    auto canonicalize = [&](ExprPtr e) -> Result<ExprPtr> {
        auto s1 = ctx.simplify(e);
        if (s1.is_error()) return s1;
        auto t = algebra::together(s1.value(), ctx);
        if (t.is_error()) return t;
        auto ex = algebra::expand(t.value(), ctx);
        if (ex.is_error()) return ex;
        return ctx.simplify(ex.value());
    };

    auto p_tilde = canonicalize(arena.make<Binary>(BinaryOp::Mul, x_sym, p_raw.value()));
    if (p_tilde.is_error()) return fail<ExprPtr>(p_tilde.error());
    ExprPtr x_sq = arena.make<Binary>(BinaryOp::Pow, x_sym, make_int(arena, 2));
    auto q_tilde = canonicalize(arena.make<Binary>(BinaryOp::Mul, x_sq, q_raw.value()));
    if (q_tilde.is_error()) return fail<ExprPtr>(q_tilde.error());

    // Sanity-check regular singular: substitute x=0 — must yield a finite value.
    {
        auto p0_check = ctx.substitute(p_tilde.value(), x, make_int(arena, 0));
        auto q0_check = ctx.substitute(q_tilde.value(), x, make_int(arena, 0));
        if (p0_check.is_error() || q0_check.is_error()) {
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                "Frobenius: could not evaluate x*p or x^2*q at x=0; "
                "regular singular check failed."));
        }
        auto p0_s = ctx.simplify(p0_check.value());
        auto q0_s = ctx.simplify(q0_check.value());
        if (p0_s.is_error() || q0_s.is_error()) {
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                "Frobenius: simplification of x*p|_0 or x^2*q|_0 failed."));
        }
        if (contains_undefined_constant(p0_s.value()) ||
            contains_undefined_constant(q0_s.value())) {
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                "Frobenius: x=0 does not appear to be a regular singular point "
                "(x*p or x^2*q is not analytic at 0)."));
        }
    }

    // Taylor coefficients of p_tilde, q_tilde up to order `num_terms`.
    const unsigned int max_k = num_terms;
    std::vector<ExprPtr> p_coeffs(max_k + 1U), q_coeffs(max_k + 1U);
    for (unsigned int k = 0; k <= max_k; ++k) {
        auto pk = taylor_coefficient(p_tilde.value(), x, k, ctx);
        if (pk.is_error()) return fail<ExprPtr>(pk.error());
        auto qk = taylor_coefficient(q_tilde.value(), x, k, ctx);
        if (qk.is_error()) return fail<ExprPtr>(qk.error());
        if (contains_undefined_constant(pk.value()) || contains_undefined_constant(qk.value())) {
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                "Frobenius: Taylor coefficient of x*p or x^2*q at k=" + std::to_string(k) +
                " is undefined; not analytic at 0."));
        }
        p_coeffs[k] = pk.value();
        q_coeffs[k] = qk.value();
    }

    ExprPtr p0 = p_coeffs[0];
    ExprPtr q0 = q_coeffs[0];

    // Indicial polynomial: r^2 + (p0 - 1) r + q0
    Symbol r_sym("r");
    ExprPtr r_ptr = arena.make<Symbol>("r");
    ExprPtr indicial = arena.make<Sum>(std::vector<ExprPtr>{
        arena.make<Binary>(BinaryOp::Pow, r_ptr, make_int(arena, 2)),
        arena.make<Binary>(BinaryOp::Mul,
            arena.make<Binary>(BinaryOp::Sub, p0, make_int(arena, 1)),
            r_ptr),
        q0,
    });
    auto indicial_simp = ctx.simplify(indicial);
    if (indicial_simp.is_error()) return fail<ExprPtr>(indicial_simp.error());

    auto roots_res = algebra::solve_polynomial(indicial_simp.value(), r_sym, ctx);
    if (roots_res.is_error()) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "Frobenius: could not solve indicial polynomial: " + roots_res.error().message));
    }
    auto roots = roots_res.value();
    if (roots.empty()) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "Frobenius: indicial polynomial returned no roots."));
    }

    // De-duplicate roots structurally
    std::vector<ExprPtr> unique_roots;
    for (auto r : roots) {
        bool dup = false;
        for (auto u : unique_roots) {
            auto eq = symbolic::mathematically_equal(r, u, ctx);
            if (eq.is_ok() && eq.value()) { dup = true; break; }
        }
        if (!dup) unique_roots.push_back(r);
    }

    // Resonance detection
    std::optional<unsigned int> N_gap;
    ExprPtr r_large, r_small;
    if (unique_roots.size() == 2U) {
        auto try_gap = integer_gap(unique_roots[0], unique_roots[1], ctx);
        if (try_gap.has_value()) {
            N_gap = try_gap;
            r_large = unique_roots[0];
            r_small = unique_roots[1];
        } else {
            try_gap = integer_gap(unique_roots[1], unique_roots[0], ctx);
            if (try_gap.has_value()) {
                N_gap = try_gap;
                r_large = unique_roots[1];
                r_small = unique_roots[0];
            }
        }
    }

    std::vector<ExprPtr> series_solutions;
    bool use_log_branch = N_gap.has_value() && *N_gap <= num_terms;
    if (use_log_branch) {
        // y_1 at the larger root via standard recurrence
        auto a_coeffs_res = compute_recurrence(r_large, p_coeffs, q_coeffs,
                                               p0, q0, num_terms, ctx);
        if (a_coeffs_res.is_error()) return fail<ExprPtr>(a_coeffs_res.error());
        auto y1_res = build_series(r_large, a_coeffs_res.value(), x, ctx);
        if (y1_res.is_error()) return fail<ExprPtr>(y1_res.error());
        series_solutions.push_back(y1_res.value());
        // y_2 via log-branch construction at the smaller root.
        auto y2_res = build_log_branch(r_large, r_small, *N_gap,
                                       a_coeffs_res.value(),
                                       p_coeffs, q_coeffs, p0, q0,
                                       num_terms, y1_res.value(), x, ctx);
        if (y2_res.is_error()) return fail<ExprPtr>(y2_res.error());
        series_solutions.push_back(y2_res.value());
    } else if (unique_roots.size() == 1U) {
        // Double indicial root (gap N = 0): the quadratic indicial polynomial has
        // a repeated root, so the two independent solutions are y_1 and
        // y_2 = ln(x)·y_1 + x^{r1}·Σ a_n'(r1) x^n (always logarithmic).
        ExprPtr r1 = unique_roots[0];
        auto a_coeffs_res = compute_recurrence(r1, p_coeffs, q_coeffs, p0, q0, num_terms, ctx);
        if (a_coeffs_res.is_error()) return fail<ExprPtr>(a_coeffs_res.error());
        auto y1_res = build_series(r1, a_coeffs_res.value(), x, ctx);
        if (y1_res.is_error()) return fail<ExprPtr>(y1_res.error());
        series_solutions.push_back(y1_res.value());
        auto y2_res = build_double_root_log_branch(
            r1, p_coeffs, q_coeffs, p0, q0, num_terms, y1_res.value(), x, ctx);
        if (y2_res.is_error()) return fail<ExprPtr>(y2_res.error());
        series_solutions.push_back(y2_res.value());
    } else {
        for (auto r : unique_roots) {
            auto coeffs_res = compute_recurrence(r, p_coeffs, q_coeffs, p0, q0, num_terms, ctx);
            if (coeffs_res.is_error()) return fail<ExprPtr>(coeffs_res.error());
            auto y_res = build_series(r, coeffs_res.value(), x, ctx);
            if (y_res.is_error()) return fail<ExprPtr>(y_res.error());
            series_solutions.push_back(y_res.value());
        }
    }

    Symbol c1_sym = ctx.make_fresh_symbol("C");
    Symbol c2_sym = ctx.make_fresh_symbol("C");
    ExprPtr C1 = arena.make<Symbol>(c1_sym);
    ExprPtr C2 = arena.make<Symbol>(c2_sym);

    if (series_solutions.size() == 1U) {
        ExprPtr total = arena.make<Binary>(BinaryOp::Mul, C1, series_solutions[0]);
        return ctx.simplify(total);
    }

    ExprPtr total = arena.make<Sum>(std::vector<ExprPtr>{
        arena.make<Binary>(BinaryOp::Mul, C1, series_solutions[0]),
        arena.make<Binary>(BinaryOp::Mul, C2, series_solutions[1]),
    });
    return ctx.simplify(total);
}

}  // namespace cas::calculus
