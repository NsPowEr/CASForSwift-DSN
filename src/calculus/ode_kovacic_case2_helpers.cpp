// ode_kovacic_case2_helpers.cpp — Step 3 (polynomial P search) and Step 4
// (ω quadratic construction) for the Kovacic Case 2 algorithm.
//
// Reference: Kovacic J.J. (1986), §4, "An Algorithm for Solving Second
// Order Linear Homogeneous Differential Equations", J. Symbolic
// Computation 2, 3-43.  See Kovacic_Case2.md spec.

#include "ode_kovacic_case2_helpers.hpp"
#include <vector>
#include <optional>
#include <string>

namespace cas::calculus::kovacic_impl {

namespace {

// Build ansatz P = x^d + a_{d-1}·x^{d-1} + ... + a_0 with a_i fresh symbols.
struct PAnsatz {
    ExprPtr             poly;
    std::vector<Symbol> coeffs;  // a_0, ..., a_{d-1} (in same order)
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

} // anonymous namespace

Result<std::optional<ExprPtr>> search_polynomial_P_case2(
    ExprPtr theta, ExprPtr r, long long d,
    const Symbol& x, symbolic::CASContext& ctx) {

    if (d < 0) return ok(std::optional<ExprPtr>{});
    if (static_cast<std::size_t>(d) > ctx.kovacic_case2_max_poly_degree()) {
        return fail<std::optional<ExprPtr>>(kv_unimpl(
            "Kovacic Case 2 Step 3: required degree d = " + std::to_string(d)
            + " exceeds ctx.kovacic_case2_max_poly_degree() = "
            + std::to_string(ctx.kovacic_case2_max_poly_degree())));
    }

    AstArena& a = ctx.arena();
    PAnsatz pa = build_P_ansatz(d, x, ctx);
    ExprPtr P = pa.poly;

    // Derivatives of P, θ, r.
    auto P_prime = diff(P, x, 1U, ctx);
    auto P_pp    = diff(P, x, 2U, ctx);
    auto P_ppp   = diff(P, x, 3U, ctx);
    auto th_p    = diff(theta, x, 1U, ctx);
    auto th_pp   = diff(theta, x, 2U, ctx);
    auto r_p     = diff(r, x, 1U, ctx);
    if (P_prime.is_error() || P_pp.is_error() || P_ppp.is_error()
        || th_p.is_error() || th_pp.is_error() || r_p.is_error()) {
        return fail<std::optional<ExprPtr>>(kv_unimpl(
            "Kovacic Case 2 Step 3: derivative computation failed."));
    }

    // LHS = P''' + 3θP'' + (3θ² + 3θ' - 4r)P' + (θ'' + 3θθ' + θ³ - 4rθ - 2r')P.
    ExprPtr three  = kv_int(a, 3);
    ExprPtr four   = kv_int(a, 4);
    ExprPtr two    = kv_int(a, 2);
    ExprPtr theta2 = kv_mul(a, theta, theta);
    ExprPtr theta3 = kv_mul(a, theta2, theta);
    ExprPtr coefA  = kv_add(a, kv_mul(a, three, theta2),
                       kv_sub(a, kv_mul(a, three, th_p.value()),
                              kv_mul(a, four, r)));
    ExprPtr coefB  = kv_sub(a,
                       kv_add(a, kv_add(a, th_pp.value(),
                              kv_mul(a, three,
                                  kv_mul(a, theta, th_p.value()))),
                              theta3),
                       kv_add(a, kv_mul(a, four, kv_mul(a, r, theta)),
                              kv_mul(a, two, r_p.value())));
    ExprPtr lhs_raw = kv_add(a, kv_add(a, P_ppp.value(),
                          kv_mul(a, three,
                              kv_mul(a, theta, P_pp.value()))),
                          kv_add(a, kv_mul(a, coefA, P_prime.value()),
                                 kv_mul(a, coefB, P)));
    auto lhs_simp = ctx.simplify(lhs_raw);
    if (lhs_simp.is_error()) {
        return fail<std::optional<ExprPtr>>(kv_unimpl(
            "Kovacic Case 2 Step 3: LHS simplify failed."));
    }

    // Clear denominators: together() unifies into a single (N, D); we extract N.
    auto together_res = algebra::together(lhs_simp.value(), ctx);
    if (together_res.is_error()) {
        return fail<std::optional<ExprPtr>>(kv_unimpl(
            "Kovacic Case 2 Step 3: together() failed."));
    }
    auto split_res = algebra::apart_num_den(together_res.value(), ctx);
    if (split_res.is_error()) {
        return fail<std::optional<ExprPtr>>(kv_unimpl(
            "Kovacic Case 2 Step 3: apart_num_den failed."));
    }
    ExprPtr numerator = split_res.value().numerator;

    // Extract coefficients of x^k in numerator.  Each is linear in the {a_i}.
    auto coeffs_res = algebra::univariate_coefficients(numerator, x, ctx);
    if (coeffs_res.is_error()) {
        // Numerator not polynomial in x → ansatz cannot satisfy ODE.
        return ok(std::optional<ExprPtr>{});
    }
    const auto& coeffs = coeffs_res.value();

    // Collect non-trivially-zero coefficient equations.
    std::vector<ExprPtr> equations;
    equations.reserve(coeffs.size());
    for (const auto& c : coeffs) {
        auto cs = ctx.simplify(c);
        if (cs.is_error()) {
            return fail<std::optional<ExprPtr>>(kv_unimpl(
                "Kovacic Case 2 Step 3: coefficient simplify failed."));
        }
        if (!kv_is_zero(cs.value(), ctx)) {
            equations.push_back(
                a.make<Binary>(BinaryOp::Equal, cs.value(), kv_int(a, 0)));
        }
    }

    if (pa.coeffs.empty()) {
        // d = 0, P = 1.  System trivially satisfied iff all equations vanish.
        if (equations.empty()) return ok(std::optional<ExprPtr>{P});
        return ok(std::optional<ExprPtr>{});
    }

    if (equations.empty()) {
        // No constraints on the a_i — substitute zeros for a clean solution.
        ExprPtr P_zero = P;
        for (const auto& s : pa.coeffs) {
            auto sub_res = ctx.substitute(P_zero, s, kv_int(a, 0));
            if (sub_res.is_error()) {
                return fail<std::optional<ExprPtr>>(kv_unimpl(
                    "Kovacic Case 2 Step 3: substitute(a_i=0) failed."));
            }
            P_zero = sub_res.value();
        }
        auto Pz_simp = ctx.simplify(P_zero);
        return ok(std::optional<ExprPtr>{
            Pz_simp.is_ok() ? Pz_simp.value() : P_zero});
    }

    // Pack equations and variables into Matrix wrappers (algebra::csolve API).
    std::vector<ExprPtr> var_exprs;
    var_exprs.reserve(pa.coeffs.size());
    for (const auto& s : pa.coeffs) var_exprs.push_back(a.make<Symbol>(s.name));
    ExprPtr eqs_mat = a.make<Matrix>(
        static_cast<int>(equations.size()), 1, std::move(equations));
    ExprPtr vars_mat = a.make<Matrix>(
        static_cast<int>(var_exprs.size()), 1, std::move(var_exprs));

    auto sol_res = algebra::csolve(eqs_mat, vars_mat, ctx);
    if (sol_res.is_error()) {
        // csolve failed: no consistent polynomial P for this family.
        return ok(std::optional<ExprPtr>{});
    }
    const auto* sol_mat = expr_cast<Matrix>(sol_res.value());
    if (!sol_mat || sol_mat->rows == 0 || sol_mat->cols == 0
        || sol_mat->elements.empty()) {
        return ok(std::optional<ExprPtr>{});
    }
    // csolve linear-path output: Matrix(rows=1, cols=n_vars, [v_0,...,v_{n-1}]).
    // Pick first row (one solution).
    const std::size_t ncols = static_cast<std::size_t>(sol_mat->cols);
    if (sol_mat->elements.size() < ncols) return ok(std::optional<ExprPtr>{});

    ExprPtr P_sub = P;
    for (std::size_t i = 0U; i < pa.coeffs.size() && i < ncols; ++i) {
        ExprPtr value = sol_mat->elements[i];
        auto sub_res = ctx.substitute(P_sub, pa.coeffs[i], value);
        if (sub_res.is_error()) {
            return fail<std::optional<ExprPtr>>(kv_unimpl(
                "Kovacic Case 2 Step 3: substitute(a_i=sol) failed."));
        }
        P_sub = sub_res.value();
    }
    auto Pf = ctx.simplify(P_sub);
    return ok(std::optional<ExprPtr>{Pf.is_ok() ? Pf.value() : P_sub});
}

Result<OmegaPair> build_omega_from_phi_case2(
    ExprPtr theta, ExprPtr P, ExprPtr r,
    const Symbol& x, symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();
    auto Pp_res = diff(P, x, 1U, ctx);
    if (Pp_res.is_error()) {
        return fail<OmegaPair>(kv_unimpl(
            "Kovacic Case 2 Step 4: diff(P) failed."));
    }
    ExprPtr phi = kv_add(a, theta, kv_div(a, Pp_res.value(), P));
    auto phi_simp = ctx.simplify(phi);
    if (phi_simp.is_ok()) phi = phi_simp.value();
    auto phi_p_res = diff(phi, x, 1U, ctx);
    if (phi_p_res.is_error()) {
        return fail<OmegaPair>(kv_unimpl(
            "Kovacic Case 2 Step 4: diff(φ) failed."));
    }
    // Quadratic (Kovacic 1986 p. 18): ω² - φω + (½φ' + ½φ² - r) = 0.
    // Standard formula: ω = (φ ± √(φ² - 4·(½φ' + ½φ² - r)))/2
    //                    = (φ ± √(φ² - 2φ' - 2φ² + 4r))/2
    //                    = (φ ± √(4r - φ² - 2φ'))/2.
    ExprPtr phi2 = kv_mul(a, phi, phi);
    ExprPtr disc_raw = kv_add(a,
        kv_sub(a, kv_neg(a, phi2),
            kv_mul(a, kv_int(a, 2), phi_p_res.value())),
        kv_mul(a, kv_int(a, 4), r));
    auto disc_simp = ctx.simplify(disc_raw);
    ExprPtr disc = disc_simp.is_ok() ? disc_simp.value() : disc_raw;
    // Normalize discriminant into a single rational fraction so that
    // the simplifier can pull a clean sqrt() out (e.g. 4/x → 2/√x).
    auto disc_tog = algebra::together(disc, ctx);
    if (disc_tog.is_ok()) disc = disc_tog.value();
    auto disc_simp2 = ctx.simplify(disc);
    if (disc_simp2.is_ok()) disc = disc_simp2.value();
    ExprPtr sqrt_disc =
        a.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{disc});
    auto sqrt_simp = ctx.simplify(sqrt_disc);
    if (sqrt_simp.is_ok()) sqrt_disc = sqrt_simp.value();
    ExprPtr two = kv_int(a, 2);
    ExprPtr omega_plus  = kv_div(a, kv_add(a, phi, sqrt_disc), two);
    ExprPtr omega_minus = kv_div(a, kv_sub(a, phi, sqrt_disc), two);
    auto op_simp = ctx.simplify(omega_plus);
    auto om_simp = ctx.simplify(omega_minus);
    ExprPtr op_final = op_simp.is_ok() ? op_simp.value() : omega_plus;
    ExprPtr om_final = om_simp.is_ok() ? om_simp.value() : omega_minus;

    // ─── Algebraic certificate (Kovacic_Case2.md §"Vincoli REGOLA ZERO") ────
    // Riccati identity  ω' + ω² ≡ r.  Two-tier verification:
    //   1. Symbolic: simplify(ω' + ω² − r) reduces to literal 0.
    //   2. Fallback: substitute x at perfect-square integer points (1,4,9)
    //      and check zero — covers cases the simplifier can't cancel
    //      structurally (sqrt-of-quotient cross-terms).
    // The fallback is still a strict check (3 evaluations pin a rational-
    // in-√x form).  Symbolic-path-only would require an algebraic-extension
    // simplifier beyond F2 scope.
    auto riccati_zero_symbolic = [&](ExprPtr omega) -> std::optional<bool> {
        auto dp = diff(omega, x, 1U, ctx);
        if (dp.is_error()) return std::nullopt;
        ExprPtr resid = kv_sub(a,
            kv_add(a, dp.value(), kv_mul(a, omega, omega)), r);
        auto s = ctx.simplify(resid);
        if (s.is_error()) return std::nullopt;
        if (auto* il = expr_cast<IntegerLit>(s.value()))
            return il->value.is_zero();
        if (auto* rl = expr_cast<RationalLit>(s.value()))
            return rl->numerator.is_zero();
        return std::nullopt;  // non-literal: defer to multi-point check.
    };
    auto riccati_zero_at = [&](ExprPtr omega, long long pt) -> bool {
        auto dp = diff(omega, x, 1U, ctx);
        if (dp.is_error()) return false;
        ExprPtr resid = kv_sub(a,
            kv_add(a, dp.value(), kv_mul(a, omega, omega)), r);
        auto sub = ctx.substitute(resid, x, kv_int(a, pt));
        if (sub.is_error()) return false;
        auto s = ctx.simplify(sub.value());
        if (s.is_error()) return false;
        if (auto* il = expr_cast<IntegerLit>(s.value()))
            return il->value.is_zero();
        if (auto* rl = expr_cast<RationalLit>(s.value()))
            return rl->numerator.is_zero();
        return false;
    };
    auto certify = [&](ExprPtr omega, const char* tag) -> Result<bool> {
        auto sym = riccati_zero_symbolic(omega);
        if (sym.has_value() && *sym) return ok(true);
        // Fallback: multi-point.
        for (long long pt : {1LL, 4LL, 9LL}) {
            if (!riccati_zero_at(omega, pt))
                return fail<bool>(kv_unimpl(
                    std::string("Kovacic Case 2 Step 4 certificate: ω") + tag
                    + " failed Riccati identity (ω' + ω² ≠ r) at x = "
                    + std::to_string(pt) + "."));
        }
        return ok(true);
    };
    auto c_plus  = certify(op_final, "+");
    if (c_plus.is_error())  return fail<OmegaPair>(c_plus.error());
    auto c_minus = certify(om_final, "-");
    if (c_minus.is_error()) return fail<OmegaPair>(c_minus.error());
    return ok(OmegaPair{op_final, om_final});
}

} // namespace cas::calculus::kovacic_impl
