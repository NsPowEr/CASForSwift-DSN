// ode_kovacic_case3_helpers.cpp — Step 3 (polynomial P search via the Kovacic
// §5 recurrence) and Step 4 (minimal polynomial of ω) for the Case 3 algorithm.
//
// Reference: Kovacic J.J. (1986), §5, "The Algorithm for Case 3".
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Kovacic_Case3_SL2C.md

#include "ode_kovacic_case3_helpers.hpp"
#include "../algebra/polynomial_internal.hpp"
#include <chrono>
#include <vector>
#include <optional>
#include <string>

namespace cas::calculus::kovacic_impl {

namespace {

using algebra::PolyExpr;

// Factorial helper (Kovacic §5 uses (n−i)! with n ≤ 12, so values fit in BigInt
// trivially).
[[nodiscard]] BigInt factorial(unsigned n) {
    BigInt r(1);
    for (unsigned k = 2U; k <= n; ++k) r = r * BigInt(static_cast<long long>(k));
    return r;
}

// d/dx of a PolyExpr with ExprPtr coefficients that are CONSTANT in x (the
// ansatz coefficients a_i are fresh symbols standing for unknowns to solve
// for, not functions of x) — plain term-by-term power-rule shift, no generic
// tree diff() needed.
[[nodiscard]] PolyExpr poly_derivative_const_coeffs(const PolyExpr& p, AstArena& a) {
    if (p.size() <= 1U) return PolyExpr{};
    std::vector<ExprPtr> out;
    out.reserve(p.size() - 1U);
    for (std::size_t k = 1U; k < p.size(); ++k) {
        ExprPtr ck = p[k];
        ExprPtr term = (k == 1U) ? ck : kv_mul(a, kv_int(a, static_cast<long long>(k)), ck);
        out.push_back(term);
    }
    return PolyExpr(std::move(out));
}

// Reduce a fixed (ansatz-independent) rational-function expression to a
// genuine polynomial in x.  Kovacic §5 guarantees S·θ and S²·r ARE
// polynomials under the Case 3 pole-order precondition already checked in
// case3_omega (S clears exactly the poles introduced by θ/r); a failure
// here signals that precondition was violated upstream, so it is reported
// as an explicit Unimplemented rather than silently mishandling a residual
// denominator.
[[nodiscard]] Result<PolyExpr> reduce_to_polynomial(
    ExprPtr expr, const Symbol& x, const char* what, symbolic::CASContext& ctx) {
    auto tog = algebra::together(expr, ctx);
    ExprPtr e = tog.is_ok() ? tog.value() : expr;
    auto s = ctx.simplify(e);
    if (s.is_ok()) e = s.value();
    auto poly_res = algebra::parse_polynomial(e, x, ctx);
    if (poly_res.is_error()) {
        return fail<PolyExpr>(kv_unimpl(
            std::string("Kovacic Case 3 recurrence: ") + what
            + " did not reduce to a polynomial in x (Case 3 pole-order "
              "precondition violated): " + poly_res.error().message));
    }
    return poly_res;
}

// Monic ansatz  P = x^d + a_{d−1}·x^{d−1} + ... + a_0  with fresh symbols a_i.
struct PAnsatz {
    ExprPtr             poly;
    std::vector<Symbol> coeffs;
};

[[nodiscard]] PAnsatz build_P_ansatz(
    long long d, const Symbol& x, symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();
    PAnsatz pa;
    if (d == 0) {
        pa.poly = kv_int(a, 1);
        return pa;
    }
    pa.coeffs.reserve(static_cast<std::size_t>(d));
    ExprPtr poly = (d == 1)
        ? static_cast<ExprPtr>(a.make<Symbol>(x.name))
        : a.make<Binary>(BinaryOp::Pow, a.make<Symbol>(x.name), kv_int(a, d));
    for (long long i = d - 1; i >= 0; --i) {
        Symbol s(ctx.make_fresh_symbol("p"));
        pa.coeffs.push_back(s);
        ExprPtr coeff_expr = a.make<Symbol>(s.name);
        ExprPtr term = (i == 0)
            ? coeff_expr
            : (i == 1)
                ? kv_mul(a, coeff_expr, a.make<Symbol>(x.name))
                : kv_mul(a, coeff_expr,
                    a.make<Binary>(BinaryOp::Pow,
                        a.make<Symbol>(x.name), kv_int(a, i)));
        poly = kv_add(a, poly, term);
    }
    pa.poly = poly;
    return pa;
}

// Compute  P_n, P_{n−1}, ..., P_{−1}  given (P, θ, S, r, n).  Stores into out
// indexed by  out[i]  for i ∈ {−1, 0, ..., n}; we map index  i → i + 1  to
// stay non-negative.
//
// Recurrence:  P_n = −P;
//              P_{i−1} = −S·P'_i + ((n−i)·S' − S·θ)·P_i − (n−i)·(i+1)·S²·r·P_{i+1}.
//
// HC-KV-06 fix: every term of this recurrence is a FIXED polynomial in x
// (S, S', S·θ, S²·r — all ansatz-independent, and genuine polynomials under
// the Case 3 pole-order precondition) times P_i/P_i'/P_{i+1}, which are
// themselves polynomials in x with coefficients that are CONSTANT in x (the
// ansatz symbols a_i).  Representing P_i as `algebra::PolyExpr` and doing the
// step in polynomial (coefficient-array) arithmetic — instead of building a
// generic rational-function AST and calling diff()+together()+simplify() on
// the whole tree at every step — avoids re-deriving/re-normalising an
// increasingly large symbolic tree n times.  S, S', S·θ, S²·r are reduced to
// PolyExpr ONCE up front (their cost is no longer multiplied by n); each
// recurrence step is then plain polynomial multiply/add/subtract, whose
// per-coefficient simplify calls stay cheap because the coefficients are
// short linear combinations of the a_i, not deep rational-function trees.
//
// Returns vector v of length n + 2 with v[i + 1] = P_i.
[[nodiscard]] Result<std::vector<ExprPtr>> compute_P_sequence(
    ExprPtr P, ExprPtr theta, ExprPtr S, ExprPtr r, unsigned n,
    const Symbol& x, symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();

    auto S_poly_res = algebra::parse_polynomial(S, x, ctx);
    if (S_poly_res.is_error()) {
        return fail<std::vector<ExprPtr>>(kv_unimpl(
            "Kovacic Case 3 recurrence: S is not a polynomial in x: "
            + S_poly_res.error().message));
    }
    const PolyExpr& S_poly = S_poly_res.value();
    const PolyExpr S_prime_poly = poly_derivative_const_coeffs(S_poly, a);

    auto S_theta_poly_res = reduce_to_polynomial(
        kv_mul(a, S, theta), x, "S*theta", ctx);
    if (S_theta_poly_res.is_error()) return fail<std::vector<ExprPtr>>(S_theta_poly_res.error());
    const PolyExpr& S_theta_poly = S_theta_poly_res.value();

    auto S2_r_poly_res = reduce_to_polynomial(
        kv_mul(a, kv_mul(a, S, S), r), x, "S^2*r", ctx);
    if (S2_r_poly_res.is_error()) return fail<std::vector<ExprPtr>>(S2_r_poly_res.error());
    const PolyExpr& S2_r_poly = S2_r_poly_res.value();

    auto P_poly_res = algebra::parse_polynomial(P, x, ctx);
    if (P_poly_res.is_error()) {
        return fail<std::vector<ExprPtr>>(kv_unimpl(
            "Kovacic Case 3 recurrence: ansatz P is not a polynomial in x: "
            + P_poly_res.error().message));
    }
    auto Pn_poly_res = algebra::poly_negate(P_poly_res.value(), ctx);
    if (Pn_poly_res.is_error()) return fail<std::vector<ExprPtr>>(Pn_poly_res.error());

    std::vector<PolyExpr> vp(n + 2U);
    // vp[i + 1] = P_i  →  vp[n + 1] = P_n = −P.
    vp[n + 1U] = Pn_poly_res.value();

    // Per-call wall-clock budget (ctx.kovacic_case3_budget_ms(), 0 = off):
    // retained as an outer safety net (the polynomial-mode arithmetic makes
    // blow-up far less likely, but the P_i degree still grows with each step
    // for multi-pole inputs).  Soft-fail return is std::vector<ExprPtr>{}
    // (empty), read by the dispatcher as "ansatz unverifiable; advance family".
    const auto start = std::chrono::steady_clock::now();
    const auto budget = std::chrono::milliseconds(ctx.kovacic_case3_budget_ms());
    const bool budget_enabled = ctx.kovacic_case3_budget_ms() > 0U;
    // i runs from n down to 0, computing P_{i−1} = vp[i].
    for (long long i = static_cast<long long>(n); i >= 0; --i) {
        if (budget_enabled && std::chrono::steady_clock::now() - start > budget)
            return ok(std::vector<ExprPtr>{});
        const PolyExpr& Pi = vp[static_cast<std::size_t>(i + 1)];
        const PolyExpr Pip1 = (static_cast<std::size_t>(i + 2) < vp.size())
            ? vp[static_cast<std::size_t>(i + 2)]
            : PolyExpr{};

        const PolyExpr Pi_prime = poly_derivative_const_coeffs(Pi, a);
        const long long n_minus_i = static_cast<long long>(n) - i;

        // term1 = -S * Pi'
        auto S_Pip_res = algebra::poly_multiply(S_poly, Pi_prime, ctx);
        if (S_Pip_res.is_error()) return fail<std::vector<ExprPtr>>(S_Pip_res.error());
        auto term1_res = algebra::poly_negate(S_Pip_res.value(), ctx);
        if (term1_res.is_error()) return fail<std::vector<ExprPtr>>(term1_res.error());

        // coeffA = (n-i)*S' - S*theta ; term2 = coeffA * Pi
        auto nmi_Sprime_res = algebra::poly_multiply(
            PolyExpr{{kv_int(a, n_minus_i)}}, S_prime_poly, ctx);
        if (nmi_Sprime_res.is_error()) return fail<std::vector<ExprPtr>>(nmi_Sprime_res.error());
        auto coeffA_res = algebra::poly_subtract(nmi_Sprime_res.value(), S_theta_poly, ctx);
        if (coeffA_res.is_error()) return fail<std::vector<ExprPtr>>(coeffA_res.error());
        auto term2_res = algebra::poly_multiply(coeffA_res.value(), Pi, ctx);
        if (term2_res.is_error()) return fail<std::vector<ExprPtr>>(term2_res.error());

        // term3 = -(n-i)(i+1) * S^2*r * Pi+1
        const long long coeff_B = n_minus_i * (i + 1);
        auto S2r_Pip1_res = algebra::poly_multiply(S2_r_poly, Pip1, ctx);
        if (S2r_Pip1_res.is_error()) return fail<std::vector<ExprPtr>>(S2r_Pip1_res.error());
        auto scaled_res = algebra::poly_multiply(
            PolyExpr{{kv_int(a, coeff_B)}}, S2r_Pip1_res.value(), ctx);
        if (scaled_res.is_error()) return fail<std::vector<ExprPtr>>(scaled_res.error());
        auto term3_res = algebra::poly_negate(scaled_res.value(), ctx);
        if (term3_res.is_error()) return fail<std::vector<ExprPtr>>(term3_res.error());

        auto sum12_res = algebra::poly_add(term1_res.value(), term2_res.value(), ctx);
        if (sum12_res.is_error()) return fail<std::vector<ExprPtr>>(sum12_res.error());
        auto Pim1_res = algebra::poly_add(sum12_res.value(), term3_res.value(), ctx);
        if (Pim1_res.is_error()) return fail<std::vector<ExprPtr>>(Pim1_res.error());

        vp[static_cast<std::size_t>(i)] = Pim1_res.value();
    }

    // Convert every P_i back to ExprPtr for the callers (build_omega_minpoly
    // needs the full sequence; search_polynomial_P_case3 needs only P_{-1}).
    std::vector<ExprPtr> v;
    v.reserve(vp.size());
    for (const auto& poly_i : vp) {
        auto e_res = algebra::polynomial_to_expr(poly_i, x, ctx);
        if (e_res.is_error()) return fail<std::vector<ExprPtr>>(e_res.error());
        v.push_back(e_res.value());
    }
    return ok(std::move(v));
}

} // anonymous namespace

Result<std::optional<ExprPtr>> search_polynomial_P_case3(
    ExprPtr theta, ExprPtr S, ExprPtr r,
    long long d, unsigned n,
    const Symbol& x, symbolic::CASContext& ctx) {

    if (d < 0) return ok(std::optional<ExprPtr>{});
    if (static_cast<std::size_t>(d) > ctx.kovacic_case3_max_poly_degree()) {
        return fail<std::optional<ExprPtr>>(kv_unimpl(
            "Kovacic Case 3 Step 3: required degree d = " + std::to_string(d)
            + " exceeds ctx.kovacic_case3_max_poly_degree() = "
            + std::to_string(ctx.kovacic_case3_max_poly_degree())));
    }

    AstArena& a = ctx.arena();
    PAnsatz pa = build_P_ansatz(d, x, ctx);
    ExprPtr P = pa.poly;

    auto seq_res = compute_P_sequence(P, theta, S, r, n, x, ctx);
    if (seq_res.is_error()) {
        return fail<std::optional<ExprPtr>>(seq_res.error());
    }
    if (seq_res.value().empty()) {
        // Sequence aborted (timeout) → ansatz cannot be checked at this
        // family; advance the family loop.
        return ok(std::optional<ExprPtr>{});
    }
    ExprPtr P_minus_one = seq_res.value()[0];

    // Clear denominators and extract numerator polynomial in x.  For large
    // n (typically 12) and non-trivial pole sets, together() and the
    // subsequent normalisation can blow up; treat such overflows as "ansatz
    // cannot be checked" rather than hard-failing the whole Case 3 attempt.
    auto together_res = algebra::together(P_minus_one, ctx);
    if (together_res.is_error()) return ok(std::optional<ExprPtr>{});
    auto split_res = algebra::apart_num_den(together_res.value(), ctx);
    if (split_res.is_error()) return ok(std::optional<ExprPtr>{});
    ExprPtr numerator = split_res.value().numerator;

    auto coeffs_res = algebra::univariate_coefficients(numerator, x, ctx);
    if (coeffs_res.is_error()) {
        // Numerator not polynomial in x → ansatz cannot satisfy ODE.
        return ok(std::optional<ExprPtr>{});
    }
    const auto& coeffs = coeffs_res.value();

    std::vector<ExprPtr> equations;
    equations.reserve(coeffs.size());
    for (const auto& c : coeffs) {
        auto cs = ctx.simplify(c);
        if (cs.is_error()) {
            return fail<std::optional<ExprPtr>>(kv_unimpl(
                "Kovacic Case 3 Step 3: coefficient simplify failed."));
        }
        if (!kv_is_zero(cs.value(), ctx)) {
            equations.push_back(
                a.make<Binary>(BinaryOp::Equal, cs.value(), kv_int(a, 0)));
        }
    }

    if (pa.coeffs.empty()) {
        if (equations.empty()) return ok(std::optional<ExprPtr>{P});
        return ok(std::optional<ExprPtr>{});
    }
    if (equations.empty()) {
        ExprPtr P_zero = P;
        for (const auto& s : pa.coeffs) {
            auto sub_res = ctx.substitute(P_zero, s, kv_int(a, 0));
            if (sub_res.is_error()) {
                return fail<std::optional<ExprPtr>>(kv_unimpl(
                    "Kovacic Case 3 Step 3: substitute(a_i=0) failed."));
            }
            P_zero = sub_res.value();
        }
        auto Pz_simp = ctx.simplify(P_zero);
        return ok(std::optional<ExprPtr>{
            Pz_simp.is_ok() ? Pz_simp.value() : P_zero});
    }

    std::vector<ExprPtr> var_exprs;
    var_exprs.reserve(pa.coeffs.size());
    for (const auto& s : pa.coeffs) var_exprs.push_back(a.make<Symbol>(s.name));
    ExprPtr eqs_mat = a.make<Matrix>(
        static_cast<int>(equations.size()), 1, std::move(equations));
    ExprPtr vars_mat = a.make<Matrix>(
        static_cast<int>(var_exprs.size()), 1, std::move(var_exprs));

    auto sol_res = algebra::csolve(eqs_mat, vars_mat, ctx);
    if (sol_res.is_error()) {
        return ok(std::optional<ExprPtr>{});
    }
    const auto* sol_mat = expr_cast<Matrix>(sol_res.value());
    if (!sol_mat || sol_mat->rows == 0 || sol_mat->cols == 0
        || sol_mat->elements.empty()) {
        return ok(std::optional<ExprPtr>{});
    }
    const std::size_t ncols = static_cast<std::size_t>(sol_mat->cols);
    if (sol_mat->elements.size() < ncols) return ok(std::optional<ExprPtr>{});

    ExprPtr P_sub = P;
    for (std::size_t i = 0U; i < pa.coeffs.size() && i < ncols; ++i) {
        ExprPtr value = sol_mat->elements[i];
        auto sub_res = ctx.substitute(P_sub, pa.coeffs[i], value);
        if (sub_res.is_error()) {
            return fail<std::optional<ExprPtr>>(kv_unimpl(
                "Kovacic Case 3 Step 3: substitute(a_i=sol) failed."));
        }
        P_sub = sub_res.value();
    }
    auto Pf = ctx.simplify(P_sub);
    return ok(std::optional<ExprPtr>{Pf.is_ok() ? Pf.value() : P_sub});
}

Result<ExprPtr> build_omega_minpoly_case3(
    ExprPtr theta, ExprPtr S, ExprPtr P, ExprPtr r, unsigned n,
    const Symbol& x, const Symbol& omega_var, symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();

    auto seq_res = compute_P_sequence(P, theta, S, r, n, x, ctx);
    if (seq_res.is_error()) return fail<ExprPtr>(seq_res.error());
    if (seq_res.value().empty()) {
        return fail<ExprPtr>(kv_unimpl(
            "Kovacic Case 3 Step 4: P-sequence aborted (timeout); cannot "
            "build minimal polynomial."));
    }
    const std::vector<ExprPtr>& v = seq_res.value();

    // Build  Σ_{i=0}^n  (S^i · P_i(x) / (n−i)!) · ω^i.
    ExprPtr omega = a.make<Symbol>(omega_var.name);
    ExprPtr S_pow_acc = kv_int(a, 1);   // S^0 = 1
    ExprPtr omega_pow_acc = kv_int(a, 1); // ω^0 = 1
    std::vector<ExprPtr> terms;
    terms.reserve(n + 1U);
    for (unsigned i = 0U; i <= n; ++i) {
        ExprPtr Pi = v[i + 1U];
        BigInt fac = factorial(n - i);
        ExprPtr fac_expr = a.make<IntegerLit>(fac);
        ExprPtr term = kv_div(a,
            kv_mul(a, S_pow_acc, kv_mul(a, Pi, omega_pow_acc)),
            fac_expr);
        terms.push_back(term);
        if (i < n) {
            S_pow_acc = kv_mul(a, S_pow_acc, S);
            omega_pow_acc = kv_mul(a, omega_pow_acc, omega);
        }
    }
    ExprPtr sum = terms.empty() ? kv_int(a, 0) : terms[0];
    for (std::size_t i = 1U; i < terms.size(); ++i)
        sum = kv_add(a, sum, terms[i]);
    auto s = ctx.simplify(sum);
    return ok(s.is_ok() ? s.value() : sum);
}

} // namespace cas::calculus::kovacic_impl
