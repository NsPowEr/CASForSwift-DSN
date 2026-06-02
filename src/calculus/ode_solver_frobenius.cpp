// L2-01: Frobenius series solution for a regular singular linear 2nd-order ODE.
//
// Given  a_2(x) y'' + a_1(x) y' + a_0(x) y = 0  with a regular singular point
// at x = 0, this constructs (up to `num_terms` terms) the Frobenius series
// solutions y_r(x) = x^r * (1 + c_1 x + c_2 x^2 + ... + c_N x^N) for each root
// of the indicial polynomial I(r) = r(r-1) + p_0 r + q_0.
//
// Notes / current scope:
//   * We treat the "regular singular" check as: x*p and x^2*q must be analytic
//     at x = 0 — verified by substituting x = 0 and ensuring no NaN / Infinity
//     constant leaks out of the simplifier.
//   * Resonant case (roots differing by a non-negative integer with non-zero
//     forcing) returns Unimplemented with a diagnostic message — logarithmic
//     branch construction is left to a follow-up task.
//   * Constants C_1, C_2 use the literal names "_C1_", "_C2_". A true fresh
//     symbol generator on CASContext is a known global gap (CLAUDE.md cat. 7);
//     these underscored names are unlikely to collide with user variables but
//     are NOT proven unique. To be revisited when ctx.make_fresh_symbol exists.

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

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] ExprPtr make_int(AstArena& arena, long long v) {
    return arena.make<IntegerLit>(BigInt(v));
}

[[nodiscard]] BigInt factorial_big(unsigned int n) {
    BigInt result(1);
    for (unsigned int i = 2; i <= n; ++i) result = result * BigInt(static_cast<long long>(i));
    return result;
}

// Test whether a simplified expression is the literal integer zero.
[[nodiscard]] bool is_literal_zero(ExprPtr e) {
    if (const auto* il = expr_cast<IntegerLit>(e)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(e)) return rl->numerator.is_zero();
    return false;
}

// Detect a NaN / Infinity leaking out of simplification (indicates we tried to
// substitute into something singular at x = 0).
[[nodiscard]] bool contains_undefined_constant(ExprPtr e) {
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
[[nodiscard]] Result<ExprPtr> taylor_coefficient(
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
[[nodiscard]] ExprPtr make_x_to_r(ExprPtr r, const Symbol& x, AstArena& arena) {
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
[[nodiscard]] Result<ExprPtr> indicial_value(
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

// Recurrence:
//   c_0 = 1
//   c_n = -[ sum_{k=1..n} ( (n-k+r) * p_k + q_k ) * c_{n-k} ] / I(n+r)
[[nodiscard]] Result<std::vector<ExprPtr>> compute_recurrence(
    ExprPtr root_r,
    const std::vector<ExprPtr>& p_coeffs,
    const std::vector<ExprPtr>& q_coeffs,
    ExprPtr p0,
    ExprPtr q0,
    unsigned int num_terms,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> c;
    c.reserve(num_terms + 1U);
    c.push_back(arena.make<IntegerLit>(BigInt(1)));

    for (unsigned int n = 1U; n <= num_terms; ++n) {
        std::vector<ExprPtr> rhs_terms;
        rhs_terms.reserve(n);
        for (unsigned int k = 1U; k <= n; ++k) {
            if (k >= p_coeffs.size() || k >= q_coeffs.size()) break;
            ExprPtr p_k = p_coeffs[k];
            ExprPtr q_k = q_coeffs[k];
            ExprPtr n_minus_k_plus_r = arena.make<Sum>(std::vector<ExprPtr>{
                make_int(arena, static_cast<long long>(n - k)),
                root_r,
            });
            ExprPtr bracket = arena.make<Sum>(std::vector<ExprPtr>{
                arena.make<Binary>(BinaryOp::Mul, n_minus_k_plus_r, p_k),
                q_k,
            });
            ExprPtr term = arena.make<Binary>(BinaryOp::Mul, bracket, c[n - k]);
            rhs_terms.push_back(term);
        }
        ExprPtr rhs_sum = rhs_terms.empty()
            ? make_int(arena, 0)
            : arena.make<Sum>(std::move(rhs_terms));
        auto rhs_simp = ctx.simplify(rhs_sum);
        if (rhs_simp.is_error()) return fail<std::vector<ExprPtr>>(rhs_simp.error());

        // Denominator: I(n + r)
        ExprPtr n_plus_r = arena.make<Sum>(std::vector<ExprPtr>{
            make_int(arena, static_cast<long long>(n)),
            root_r,
        });
        auto denom_res = indicial_value(p0, q0, n_plus_r, ctx);
        if (denom_res.is_error()) return fail<std::vector<ExprPtr>>(denom_res.error());
        ExprPtr denom = denom_res.value();

        if (is_literal_zero(denom)) {
            if (is_literal_zero(rhs_simp.value())) {
                // Free parameter — pick 0 by convention.
                c.push_back(make_int(arena, 0));
                continue;
            }
            return fail<std::vector<ExprPtr>>(make_error(
                CASErrorKind::Unimplemented,
                "Frobenius resonance at n=" + std::to_string(n) +
                    ": indicial polynomial vanishes with non-zero RHS — "
                    "logarithmic branch required (not yet implemented)."));
        }

        ExprPtr numerator = arena.make<Unary>(UnaryOp::Neg, rhs_simp.value());
        ExprPtr c_n_raw = arena.make<Binary>(BinaryOp::Div, numerator, denom);
        auto c_n_simp = ctx.simplify(c_n_raw);
        if (c_n_simp.is_error()) return fail<std::vector<ExprPtr>>(c_n_simp.error());
        c.push_back(c_n_simp.value());
    }
    return ok(c);
}

// Detect r1 - r2 = N with N a strictly positive integer.  Returns N on
// success; nullopt otherwise (non-integer gap, zero gap, negative gap, or
// an unsimplifiable expression).
[[nodiscard]] std::optional<unsigned int> integer_gap(
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

// Logarithmic Frobenius branch — Coddington-Levinson §4.8.
//
// Setup.  Indicial roots r_1 > r_2 with N = r_1 - r_2 ∈ Z_{>0}.  The first
// solution y_1 = x^{r_1}·Σ a_n x^n is already computed (`a_coeffs`).
// The second independent solution has the form
//   y_2(x) = c·ln(x)·y_1(x) + x^{r_2}·Σ b_n x^n
// with c possibly zero.  Substituting into  L[y] = x²y'' + xP(x)y' + Q(x)y
// (P = p_tilde, Q = q_tilde) and using L[y_1] = 0 yields
//   L[v] = -c·(2x·y_1' + (P-1)·y_1),   v = x^{r_2}·Σ b_n x^n.
// Equating x^{r_2+n} coefficients gives the recurrence
//   I(r_2 + n)·b_n + S_n  =  - c · h_{n-N}    (with h_m = 0 for m < 0)
// where
//   S_n   = Σ_{k=1}^n  ((r_2 + n - k)·p_k + q_k) · b_{n-k}
//   h_m   = (2(r_1 + m) + p_0 - 1)·a_m + Σ_{k=1}^m p_k · a_{m-k}.
// Note h_0 = I'(r_1) = r_1 - r_2 = N ≠ 0.
//
// Algorithm:
//   * n = 0           : b_0 = 1.
//   * 1 ≤ n < N       : b_n = -S_n / I(r_2 + n).  (denominator non-zero.)
//   * n = N           : I(r_2 + N) = I(r_1) = 0.  Solve c·h_0 = -S_N for c,
//                       i.e. c = -S_N / N.  Then b_N is a free parameter
//                       which we fix to 0 by convention.
//   * n > N           : b_n = -(S_n + c·h_{n-N}) / I(r_2 + n).
//
// If c happens to evaluate to zero the second solution is logarithm-free;
// this is the historical "free parameter" case that the standard recurrence
// already accepted with a zero RHS.
[[nodiscard]] Result<ExprPtr> build_log_branch(
    ExprPtr r1,
    ExprPtr r2,
    unsigned int N,
    const std::vector<ExprPtr>& a_coeffs,  // y_1 series coefficients
    const std::vector<ExprPtr>& p_coeffs,
    const std::vector<ExprPtr>& q_coeffs,
    ExprPtr p0,
    ExprPtr q0,
    unsigned int num_terms,
    ExprPtr y_1_series,
    const Symbol& x,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    if (num_terms < N) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "Frobenius log branch: num_terms (" + std::to_string(num_terms) +
            ") smaller than the integer resonance gap N=" + std::to_string(N) +
            "; cannot resolve the c·ln(x) coupling.  Increase the series order."));
    }

    auto coef_term = [&](ExprPtr root_r, unsigned int n, unsigned int k,
                          ExprPtr b_n_minus_k) -> ExprPtr {
        ExprPtr p_k = p_coeffs[k];
        ExprPtr q_k = q_coeffs[k];
        ExprPtr r_plus_nk = arena.make<Sum>(std::vector<ExprPtr>{
            make_int(arena, static_cast<long long>(n - k)),
            root_r});
        ExprPtr bracket = arena.make<Sum>(std::vector<ExprPtr>{
            arena.make<Binary>(BinaryOp::Mul, r_plus_nk, p_k),
            q_k});
        return arena.make<Binary>(BinaryOp::Mul, bracket, b_n_minus_k);
    };

    // Pre-compute h_m for m = 0 .. num_terms - N.
    std::vector<ExprPtr> h(num_terms - N + 1U);
    for (unsigned int m = 0U; m + N <= num_terms && m < a_coeffs.size(); ++m) {
        // (2(r_1 + m) + p_0 - 1) · a_m
        ExprPtr two_r1_m = arena.make<Binary>(BinaryOp::Mul, make_int(arena, 2),
            arena.make<Sum>(std::vector<ExprPtr>{r1, make_int(arena, static_cast<long long>(m))}));
        ExprPtr lead_factor = arena.make<Sum>(std::vector<ExprPtr>{
            two_r1_m, p0, make_int(arena, -1)});
        ExprPtr lead = arena.make<Binary>(BinaryOp::Mul, lead_factor, a_coeffs[m]);
        // Σ_{k=1}^m p_k · a_{m-k}
        std::vector<ExprPtr> cross_terms;
        cross_terms.push_back(lead);
        for (unsigned int k = 1U; k <= m && k < p_coeffs.size(); ++k) {
            cross_terms.push_back(arena.make<Binary>(BinaryOp::Mul,
                p_coeffs[k], a_coeffs[m - k]));
        }
        ExprPtr h_m = (cross_terms.size() == 1U) ? cross_terms[0]
            : arena.make<Sum>(std::move(cross_terms));
        auto hs = ctx.simplify(h_m);
        if (hs.is_error()) return hs;
        h[m] = hs.value();
    }

    // Recurrence for b_n with c determined at n = N.
    std::vector<ExprPtr> b(num_terms + 1U);
    b[0] = make_int(arena, 1);
    ExprPtr c_log = make_int(arena, 0);  // updated at n = N
    bool c_log_resolved = false;

    auto build_S_n = [&](unsigned int n) -> Result<ExprPtr> {
        std::vector<ExprPtr> terms;
        for (unsigned int k = 1U; k <= n; ++k) {
            if (k >= p_coeffs.size() || k >= q_coeffs.size()) break;
            terms.push_back(coef_term(r2, n, k, b[n - k]));
        }
        ExprPtr S = terms.empty() ? make_int(arena, 0)
            : (terms.size() == 1U ? terms[0]
               : arena.make<Sum>(std::move(terms)));
        return ctx.simplify(S);
    };

    for (unsigned int n = 1U; n <= num_terms; ++n) {
        auto S_res = build_S_n(n);
        if (S_res.is_error()) return S_res;
        ExprPtr S_n = S_res.value();

        if (n == N) {
            // c · h_0 = -S_N  →  c = -S_N / h_0,  h_0 = N (≠ 0).
            ExprPtr c_raw = arena.make<Binary>(BinaryOp::Div,
                arena.make<Unary>(UnaryOp::Neg, S_n), h[0]);
            auto cs = ctx.simplify(c_raw);
            if (cs.is_error()) return cs;
            c_log = cs.value();
            c_log_resolved = true;
            b[n] = make_int(arena, 0);  // free parameter
            continue;
        }

        // RHS contribution from c·h_{n-N} when n > N.
        ExprPtr rhs_correction = make_int(arena, 0);
        if (c_log_resolved && n > N) {
            unsigned int m = n - N;
            if (m < h.size()) {
                rhs_correction = arena.make<Binary>(BinaryOp::Mul, c_log, h[m]);
            }
        }

        ExprPtr S_plus_corr = arena.make<Binary>(BinaryOp::Add, S_n, rhs_correction);
        auto S_total = ctx.simplify(S_plus_corr);
        if (S_total.is_error()) return S_total;

        ExprPtr n_plus_r2 = arena.make<Sum>(std::vector<ExprPtr>{
            make_int(arena, static_cast<long long>(n)), r2});
        auto denom_res = indicial_value(p0, q0, n_plus_r2, ctx);
        if (denom_res.is_error()) return denom_res;
        ExprPtr denom = denom_res.value();
        if (is_literal_zero(denom)) {
            // Secondary resonance at a different gap — beyond this branch.
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                "Frobenius log branch: secondary resonance at n=" +
                std::to_string(n) +
                " encountered while building the b_n series.  Multiple "
                "resonance levels require an extended log-power construction "
                "(not yet implemented)."));
        }
        ExprPtr b_n_raw = arena.make<Binary>(BinaryOp::Div,
            arena.make<Unary>(UnaryOp::Neg, S_total.value()), denom);
        auto b_n_simp = ctx.simplify(b_n_raw);
        if (b_n_simp.is_error()) return b_n_simp;
        b[n] = b_n_simp.value();
    }

    // Assemble y_2 = c · ln(x) · y_1  +  x^{r_2} · Σ b_n x^n.
    ExprPtr x_sym = arena.make<Symbol>(x.name);
    std::vector<ExprPtr> series_terms;
    series_terms.push_back(b[0]);
    for (unsigned int n = 1U; n <= num_terms; ++n) {
        if (is_literal_zero(b[n])) continue;
        ExprPtr xn = (n == 1U) ? x_sym
            : arena.make<Binary>(BinaryOp::Pow, x_sym,
                make_int(arena, static_cast<long long>(n)));
        series_terms.push_back(arena.make<Binary>(BinaryOp::Mul, b[n], xn));
    }
    ExprPtr inner = (series_terms.size() == 1U) ? series_terms[0]
        : arena.make<Sum>(std::move(series_terms));
    auto inner_simp = ctx.simplify(inner);
    if (inner_simp.is_error()) return inner_simp;
    ExprPtr power_part = arena.make<Binary>(BinaryOp::Mul,
        make_x_to_r(r2, x, arena), inner_simp.value());

    ExprPtr ln_x = arena.make<FuncCall>(BuiltinOp::Ln,
        std::vector<ExprPtr>{x_sym});
    ExprPtr log_part = arena.make<Binary>(BinaryOp::Mul, c_log,
        arena.make<Binary>(BinaryOp::Mul, ln_x, y_1_series));
    ExprPtr y_2 = arena.make<Binary>(BinaryOp::Add, log_part, power_part);
    return ctx.simplify(y_2);
}

// Build the Frobenius series y_r(x) = x^r * (1 + c_1 x + c_2 x^2 + ... )
[[nodiscard]] Result<ExprPtr> build_series(
    ExprPtr root_r,
    const std::vector<ExprPtr>& c,
    const Symbol& x,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr x_sym = arena.make<Symbol>(x.name);

    std::vector<ExprPtr> inner_terms;
    inner_terms.reserve(c.size());
    inner_terms.push_back(c[0]);  // by construction c_0 = 1
    for (std::size_t n = 1; n < c.size(); ++n) {
        if (is_literal_zero(c[n])) continue;
        ExprPtr x_power = (n == 1)
            ? x_sym
            : arena.make<Binary>(BinaryOp::Pow, x_sym, make_int(arena, static_cast<long long>(n)));
        inner_terms.push_back(arena.make<Binary>(BinaryOp::Mul, c[n], x_power));
    }
    ExprPtr inner = (inner_terms.size() == 1U)
        ? inner_terms[0]
        : arena.make<Sum>(std::move(inner_terms));
    auto inner_simp = ctx.simplify(inner);
    if (inner_simp.is_error()) return fail<ExprPtr>(inner_simp.error());

    ExprPtr xr = make_x_to_r(root_r, x, arena);
    ExprPtr series = arena.make<Binary>(BinaryOp::Mul, xr, inner_simp.value());
    return ctx.simplify(series);
}

}  // namespace

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
    // For inputs where p(x) = a_1/a_2 still carries an x in the denominator
    // (e.g. p = -4/(3x) coming from a_1 = -4x, a_2 = 3x^2), naive simplify
    // does not always cancel the outer multiplication.  Force canonicalization
    // via algebra::together + algebra::expand so the resulting expression is
    // a true polynomial-or-rational in x with the denominator zero analyzable
    // by substitution.
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

    // De-duplicate roots structurally (we don't try to merge by math-equality
    // here; identical roots fold into one and the second independent solution
    // would normally require a log term — flagged when we hit it).
    std::vector<ExprPtr> unique_roots;
    for (auto r : roots) {
        bool dup = false;
        for (auto u : unique_roots) {
            auto eq = symbolic::mathematically_equal(r, u, ctx);
            if (eq.is_ok() && eq.value()) { dup = true; break; }
        }
        if (!dup) unique_roots.push_back(r);
    }

    // Resonance detection: if exactly two roots and their difference is a
    // strictly positive integer, route the smaller-root branch through the
    // logarithmic Frobenius construction (Coddington-Levinson §4.8).  This
    // closes B2c — previously the resonant case returned Unimplemented at
    // the indicial zero.
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
    // Engage log branch only when the requested series order reaches the
    // resonance step; otherwise the standard recurrence at the smaller root
    // terminates before n = N and yields a valid truncated series.
    bool use_log_branch = N_gap.has_value() && *N_gap <= num_terms;
    if (use_log_branch) {
        // y_1 at the larger root via standard recurrence (always works).
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
    } else {
        for (auto r : unique_roots) {
            auto coeffs_res = compute_recurrence(r, p_coeffs, q_coeffs, p0, q0, num_terms, ctx);
            if (coeffs_res.is_error()) return fail<ExprPtr>(coeffs_res.error());
            auto y_res = build_series(r, coeffs_res.value(), x, ctx);
            if (y_res.is_error()) return fail<ExprPtr>(y_res.error());
            series_solutions.push_back(y_res.value());
        }
    }

    // General solution: C_1 * y_1 + C_2 * y_2  (or single y_1 with C_1 if
    // only one root).  Constants generated through ctx.make_fresh_symbol
    // (HC-004 closed): names are guaranteed not to collide with user
    // symbols, no fixed "_C1_"/"_C2_" literal anywhere.
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
