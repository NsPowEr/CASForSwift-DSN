// L2-01: Frobenius series solution - Series and recurrence logic.
#include "ode_solver_frobenius_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/bigint.hpp"
#include "cas/calculus.hpp"
#include "cas/error.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <vector>
#include <string>

namespace cas::calculus {

// Recurrence:
//   c_0 = 1
//   c_n = -[ sum_{k=1..n} ( (n-k+r) * p_k + q_k ) * c_{n-k} ] / I(n+r)
Result<std::vector<ExprPtr>> compute_recurrence(
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

// Logarithmic Frobenius branch — Coddington-Levinson §4.8.
Result<ExprPtr> build_log_branch(
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
Result<ExprPtr> build_series(
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

}  // namespace cas::calculus
