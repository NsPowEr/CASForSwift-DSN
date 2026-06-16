// ode_kovacic_case3_helpers.cpp — Step 3 (polynomial P search via the Kovacic
// §5 recurrence) and Step 4 (minimal polynomial of ω) for the Case 3 algorithm.
//
// Reference: Kovacic J.J. (1986), §5, "The Algorithm for Case 3".
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Kovacic_Case3_SL2C.md

#include "ode_kovacic_case3_helpers.hpp"
#include <chrono>
#include <vector>
#include <optional>
#include <string>

namespace cas::calculus::kovacic_impl {

namespace {

// Factorial helper (Kovacic §5 uses (n−i)! with n ≤ 12, so values fit in BigInt
// trivially).
[[nodiscard]] BigInt factorial(unsigned n) {
    BigInt r(1);
    for (unsigned k = 2U; k <= n; ++k) r = r * BigInt(static_cast<long long>(k));
    return r;
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
// Returns vector v of length n + 2 with v[i + 1] = P_i.
[[nodiscard]] Result<std::vector<ExprPtr>> compute_P_sequence(
    ExprPtr P, ExprPtr theta, ExprPtr S, ExprPtr r, unsigned n,
    const Symbol& x, symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();

    auto S_prime_res = diff(S, x, 1U, ctx);
    if (S_prime_res.is_error()) {
        return fail<std::vector<ExprPtr>>(kv_unimpl(
            "Kovacic Case 3 recurrence: diff(S) failed: "
            + S_prime_res.error().message));
    }
    ExprPtr S_prime = S_prime_res.value();
    ExprPtr S_squared = kv_mul(a, S, S);

    std::vector<ExprPtr> v(n + 2U, nullptr);
    // v[i + 1] = P_i  →  v[n + 1] = P_n = −P;  v[n + 2 missing] never used.
    v[n + 1U] = kv_neg(a, P);

    // Per-call wall-clock budget: aborts recurrence promptly when the
    // intermediate symbolic forms explode (typical for n=12 with multi-pole
    // inputs).  Soft-fail return type is std::vector<ExprPtr>{} (empty),
    // which the dispatcher reads as "ansatz unverifiable; advance family".
    const auto start = std::chrono::steady_clock::now();
    const auto budget = std::chrono::milliseconds(500);
    // i runs from n down to 0, computing P_{i−1} = v[i].
    for (long long i = static_cast<long long>(n); i >= 0; --i) {
        if (std::chrono::steady_clock::now() - start > budget)
            return ok(std::vector<ExprPtr>{});
        ExprPtr Pi   = v[static_cast<std::size_t>(i + 1)];
        ExprPtr Pip1 = (static_cast<std::size_t>(i + 2) < v.size())
            ? v[static_cast<std::size_t>(i + 2)]
            : kv_int(a, 0);

        auto Pi_p_res = diff(Pi, x, 1U, ctx);
        if (Pi_p_res.is_error()) {
            // diff() can hit symbolic-op timeout on the deepest P_i for
            // large recurrence chains.  Treat as "ansatz unverifiable" and
            // bail back to the family loop rather than abort Case 3.
            return ok(std::vector<ExprPtr>{});
        }

        const long long n_minus_i = static_cast<long long>(n) - i;
        ExprPtr coeff_A = kv_sub(a,
            kv_mul(a, kv_int(a, n_minus_i), S_prime),
            kv_mul(a, S, theta));
        ExprPtr coeff_B = kv_int(a, n_minus_i * (i + 1));

        ExprPtr term1 = kv_neg(a, kv_mul(a, S, Pi_p_res.value()));
        ExprPtr term2 = kv_mul(a, coeff_A, Pi);
        ExprPtr term3 = kv_neg(a,
            kv_mul(a, coeff_B,
                kv_mul(a, S_squared,
                    kv_mul(a, r, Pip1))));
        ExprPtr Pim1 = kv_add(a, term1, kv_add(a, term2, term3));
        // Per-step canonicalisation: together() collapses sums of rational
        // functions into a single N/D, simplify trims content.  Sufficient
        // for n ∈ {4, 6}.  n = 12 hits the wall-clock budget in case3.cpp;
        // tracked HC-KV-06 (polynomial-mode rep is the true fix).
        auto tog = algebra::together(Pim1, ctx);
        if (tog.is_ok()) Pim1 = tog.value();
        auto s = ctx.simplify(Pim1);
        v[static_cast<std::size_t>(i)] = s.is_ok() ? s.value() : Pim1;
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
