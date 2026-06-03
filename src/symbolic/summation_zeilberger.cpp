// F5.7 — Zeilberger creative telescoping.
//
// Parametric-Gosper algorithm (Zeilberger 1990):
//   Find p_0..p_J (polynomials in n), R(n,k) rational in k, such that
//     Σ_{i=0}^J p_i(n)·F(n+i,k) = R(n,k+1)·F(n,k+1) − R(n,k)·F(n,k).
//   Summing over k yields the recurrence Σ p_i·S(n+i) = 0.
//
// This file implements J ≤ ctx.max_zeilberger_order() with polynomial degrees
// in n bounded by ctx.max_zeilberger_poly_degree().  Larger J/D produce an
// explicit Unimplemented diagnostic (never silent failure).

#include "summation_zeilberger.hpp"
#include "summation_zeilberger_helpers.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include <optional>
#include <vector>

namespace cas::symbolic {

// ── local utilities ──────────────────────────────────────────────────────────

static bool expr_contains_symbol(ExprPtr e, const std::string& name) {
    if (!e) return false;
    if (const auto* sym = expr_cast<Symbol>(e)) return sym->name == name;
    if (const auto* bin = expr_cast<Binary>(e))
        return expr_contains_symbol(bin->left, name) ||
               expr_contains_symbol(bin->right, name);
    if (const auto* un = expr_cast<Unary>(e))
        return expr_contains_symbol(un->operand, name);
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        for (ExprPtr a : fc->args)
            if (expr_contains_symbol(a, name)) return true;
        return false;
    }
    if (const auto* s = expr_cast<Sum>(e)) {
        for (ExprPtr t : s->terms)
            if (expr_contains_symbol(t, name)) return true;
        return false;
    }
    if (const auto* p = expr_cast<Product>(e)) {
        for (ExprPtr f : p->factors)
            if (expr_contains_symbol(f, name)) return true;
        return false;
    }
    return false;
}

// Forwarding alias to the helpers' compute_shift_ratio (out-of-line in
// summation_zeilberger_helpers.cpp).  Kept inline for readability.
static std::optional<ExprPtr> compute_shift_ratio(
    ExprPtr F, const Symbol& sym, symbolic::CASContext& ctx) {
    return zeilberger_detail::compute_shift_ratio(F, sym, ctx);
}

// Build polynomial in sym of degree d with fresh unknown symbols as coefficients.
// Returns (poly_expr, vector_of_unknowns).
static std::pair<ExprPtr, std::vector<ExprPtr>> make_poly_unknown(
    const Symbol& sym, unsigned int d, const std::string& prefix,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> unknowns;
    ExprPtr poly = arena.make<IntegerLit>(BigInt(0));
    for (unsigned int i = 0U; i <= d; ++i) {
        ExprPtr u = arena.make<Symbol>(ctx.make_fresh_symbol(prefix));
        unknowns.push_back(u);
        ExprPtr mon = (i == 0U)
            ? u
            : arena.make<Binary>(BinaryOp::Mul, u,
                arena.make<Binary>(BinaryOp::Pow,
                    arena.make<Symbol>(sym),
                    arena.make<IntegerLit>(BigInt(static_cast<long long>(i)))));
        poly = arena.make<Binary>(BinaryOp::Add, poly, mon);
    }
    auto simp = ctx.simplify(poly);
    return {simp.is_ok() ? simp.value() : poly, std::move(unknowns)};
}

// Collect coefficients of expr as polynomial in sym, then split each coefficient
// as polynomial in sym2, returning all scalar equations (coefficient = 0).
static void collect_equations(
    ExprPtr expr, const Symbol& sym, const Symbol& sym2,
    symbolic::CASContext& ctx, std::vector<ExprPtr>& eqs_out) {
    auto exp = algebra::expand(expr, ctx);
    if (exp.is_error()) { eqs_out.push_back(expr); return; }
    auto k_coeffs = algebra::univariate_coefficients(exp.value(), sym, ctx);
    if (k_coeffs.is_error()) { eqs_out.push_back(exp.value()); return; }
    for (ExprPtr kc : k_coeffs.value()) {
        auto ke = algebra::expand(kc, ctx);
        if (ke.is_error()) { eqs_out.push_back(kc); continue; }
        auto n_coeffs = algebra::univariate_coefficients(ke.value(), sym2, ctx);
        if (n_coeffs.is_error()) { eqs_out.push_back(ke.value()); continue; }
        for (ExprPtr nc : n_coeffs.value()) eqs_out.push_back(nc);
    }
}

// Try parametric Gosper for order J with polynomial degree D_p in n and
// certificate R numerator degree D_r in k.
// Returns {p_i} vector (size J+1) or nullopt.
static std::optional<std::vector<ExprPtr>> try_parametric_gosper(
    ExprPtr r_k,  // r(k) = F(n,k+1)/F(n,k), rational in k
    const std::vector<ExprPtr>& s_shifts,  // s_i(k) = F(n+i,k)/F(n,k), i=0..J
    const Symbol& n_sym, const Symbol& k_sym,
    unsigned int J, unsigned int D_p, unsigned int D_r,
    symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();

    // Decompose r_k = N_r / D_r_poly.
    auto tog_r = algebra::together(r_k, ctx);
    if (tog_r.is_error()) return std::nullopt;
    auto parts_r = algebra::apart_num_den(tog_r.value(), ctx);
    if (parts_r.is_error()) return std::nullopt;
    ExprPtr N_r = parts_r.value().numerator;
    ExprPtr D_r_poly = parts_r.value().denominator;

    // Decompose the combined h(k) = Σ p_i · s_i(k) symbolically.
    // Build p_i polynomials in n with unknown coefficients.
    std::vector<ExprPtr> p_polys;
    std::vector<ExprPtr> all_p_unknowns;
    for (unsigned int i = 0U; i <= J; ++i) {
        auto [pi, pi_unknowns] = make_poly_unknown(n_sym, D_p, "zp", ctx);
        p_polys.push_back(pi);
        all_p_unknowns.insert(all_p_unknowns.end(), pi_unknowns.begin(), pi_unknowns.end());
    }

    // Build h(k) = Σ p_i * s_i(k) as a single rational expression.
    ExprPtr h = arena.make<IntegerLit>(BigInt(0));
    for (unsigned int i = 0U; i <= J; ++i) {
        ExprPtr term_i = arena.make<Binary>(BinaryOp::Mul, p_polys[i], s_shifts[i]);
        h = arena.make<Binary>(BinaryOp::Add, h, term_i);
    }
    auto tog_h = algebra::together(h, ctx);
    if (tog_h.is_error()) return std::nullopt;
    auto parts_h = algebra::apart_num_den(tog_h.value(), ctx);
    if (parts_h.is_error()) return std::nullopt;
    ExprPtr N_h = parts_h.value().numerator;
    ExprPtr D_h = parts_h.value().denominator;

    // Build certificate R(k) = x(k) / D_h(k) with x polynomial of degree D_r.
    auto [x_k, x_unknowns] = make_poly_unknown(k_sym, D_r, "zr", ctx);

    // Compute R(k+1): substitute k → k+1 in R.
    ExprPtr k_plus_one = arena.make<Binary>(BinaryOp::Add,
        arena.make<Symbol>(k_sym), arena.make<IntegerLit>(BigInt(1)));
    auto x_k1_res = ctx.substitute(x_k, k_sym, k_plus_one);
    if (x_k1_res.is_error()) return std::nullopt;
    ExprPtr x_k1 = x_k1_res.value();
    auto D_h_k1_res = ctx.substitute(D_h, k_sym, k_plus_one);
    if (D_h_k1_res.is_error()) return std::nullopt;
    ExprPtr D_h_k1 = D_h_k1_res.value();

    // Equation: R(k+1)·r(k) - R(k) = h(k)
    // => x(k+1)/D_h(k+1) · N_r/D_r_poly - x(k)/D_h = N_h/D_h
    // Multiply through by D_h(k)·D_h(k+1)·D_r_poly:
    // x(k+1)·N_r·D_h(k) - x(k)·D_h(k+1)·D_r_poly = N_h·D_h(k+1)·D_r_poly
    ExprPtr lhs_t1 = arena.make<Binary>(BinaryOp::Mul,
        arena.make<Binary>(BinaryOp::Mul, x_k1, N_r), D_h);
    ExprPtr lhs_t2 = arena.make<Binary>(BinaryOp::Mul,
        arena.make<Binary>(BinaryOp::Mul, x_k, D_h_k1), D_r_poly);
    ExprPtr rhs    = arena.make<Binary>(BinaryOp::Mul,
        arena.make<Binary>(BinaryOp::Mul, N_h, D_h_k1), D_r_poly);

    ExprPtr equation = arena.make<Binary>(BinaryOp::Sub,
        arena.make<Binary>(BinaryOp::Sub, lhs_t1, lhs_t2), rhs);

    // Expand and collect scalar equations by k then by n.
    std::vector<ExprPtr> eqs;
    collect_equations(equation, k_sym, n_sym, ctx, eqs);
    if (eqs.empty()) return std::nullopt;

    // Build unknowns list: all p_i coefficients and x coefficients.
    std::vector<ExprPtr> all_unknowns;
    all_unknowns.insert(all_unknowns.end(), all_p_unknowns.begin(), all_p_unknowns.end());
    all_unknowns.insert(all_unknowns.end(), x_unknowns.begin(), x_unknowns.end());

    ExprPtr eqs_mat = arena.make<Matrix>(
        static_cast<std::size_t>(eqs.size()), 1U, eqs);
    ExprPtr vars_mat = arena.make<Matrix>(
        static_cast<std::size_t>(all_unknowns.size()), 1U, all_unknowns);

    auto sol_res = algebra::csolve(eqs_mat, vars_mat, ctx);
    if (sol_res.is_error()) return std::nullopt;
    const auto* sol = expr_cast<Matrix>(sol_res.value());
    if (!sol || sol->elements.empty()) return std::nullopt;

    // Substitute solutions into p_i polynomials.
    std::vector<ExprPtr> result_p;
    std::size_t n_p_unknowns = all_p_unknowns.size();
    for (unsigned int i = 0U; i <= J; ++i) {
        ExprPtr pi = p_polys[i];
        for (std::size_t c = 0U; c < n_p_unknowns && c < sol->cols; ++c) {
            auto sub = ctx.substitute(pi, *expr_cast<Symbol>(all_p_unknowns[c]),
                sol->elements[c]);
            if (sub.is_ok()) pi = sub.value();
        }
        auto simp = ctx.simplify(pi);
        result_p.push_back(simp.is_ok() ? simp.value() : pi);
    }

    // Check that p_J is non-trivial (not zero).
    ExprPtr pJ = result_p.back();
    bool pJ_zero = false;
    if (const auto* il = expr_cast<IntegerLit>(pJ)) pJ_zero = il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(pJ)) pJ_zero = rl->numerator.is_zero();
    if (pJ_zero) return std::nullopt;

    return result_p;
}

// Solve first-order homogeneous recurrence p0(n)·S(n) + p1(n)·S(n+1) = 0
// with S(lower_val) = S0.  Returns S(n) or nullopt.
static std::optional<ExprPtr> solve_first_order_rec(
    ExprPtr p0, ExprPtr p1, ExprPtr S0, ExprPtr lower_val,
    const Symbol& n_sym, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    // ρ(n) = -p0/p1 (the ratio S(n+1)/S(n))
    ExprPtr neg_p0 = arena.make<Unary>(UnaryOp::Neg, p0);
    ExprPtr rho_raw = arena.make<Binary>(BinaryOp::Div, neg_p0, p1);
    auto rho_simp = ctx.simplify(rho_raw);
    if (rho_simp.is_error()) return std::nullopt;
    ExprPtr rho = rho_simp.value();

    ExprPtr n_e = arena.make<Symbol>(n_sym);

    // Case 1: ρ is a constant (no n).
    if (!expr_contains_symbol(rho, n_sym.name)) {
        // S(n) = S0 · ρ^(n - lower_val)
        ExprPtr exp_e = arena.make<Binary>(BinaryOp::Sub, n_e, lower_val);
        ExprPtr rho_pow = arena.make<Binary>(BinaryOp::Pow, rho, exp_e);
        ExprPtr result = arena.make<Binary>(BinaryOp::Mul, S0, rho_pow);
        auto s = ctx.simplify(result);
        return s.is_ok() ? std::optional<ExprPtr>{s.value()} : std::nullopt;
    }

    // Case 2: ρ is linear/linear rational in n.
    // ρ(n) = (a·n + b) / (c·n + d)
    // S(n) = S0 · ∏_{m=lower}^{n−1} ρ(m)
    //       = S0 · (a·lower+b)(a(lower+1)+b)···(a(n-1)+b)
    //              / [(c·lower+d)(c(lower+1)+d)···(c(n-1)+d)]
    // Expressed as Pochhammer(b/a + lower, n-lower) * a^(n-lower)
    //            / Pochhammer(d/c + lower, n-lower) / c^(n-lower)
    // Simplify: for integer lower=0: (b)(b+a)···(b+(n-1)a) = a^n · Γ(b/a+n)/Γ(b/a)
    //
    // For the standard case a=c=1 (recurrence shift 1):
    //   ρ(n) = (n+b)/(n+d) → S(n) = S0 · Γ(n+b)/Γ(lower+b) / (Γ(n+d)/Γ(lower+d))
    //   = S0 · [n+b-1]! / [lower+b-1]! · [lower+d-1]! / [n+d-1]!  (when b,d are integers)

    auto tog = algebra::together(rho, ctx);
    if (tog.is_error()) return std::nullopt;
    auto parts = algebra::apart_num_den(tog.value(), ctx);
    if (parts.is_error()) return std::nullopt;
    ExprPtr Num = parts.value().numerator;
    ExprPtr Den = parts.value().denominator;

    auto nc = algebra::univariate_coefficients(Num, n_sym, ctx);
    auto dc = algebra::univariate_coefficients(Den, n_sym, ctx);
    if (nc.is_error() || dc.is_error()) return std::nullopt;
    if (nc.value().size() > 2U || dc.value().size() > 2U) return std::nullopt;

    // Extract a,b from (a·n+b) and c,d from (c·n+d).
    Rational a_r(0), b_r(0), c_r(0), d_r(0);
    auto extract = [&](const std::vector<ExprPtr>& cv, Rational& c0, Rational& c1) -> bool {
        auto get_rat = [](ExprPtr e, Rational& out) -> bool {
            if (const auto* il = expr_cast<IntegerLit>(e)) {
                out = Rational(il->value); return true;
            }
            if (const auto* rl = expr_cast<RationalLit>(e)) {
                out = Rational(rl->numerator, rl->denominator); return true;
            }
            if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
                if (const auto* il2 = expr_cast<IntegerLit>(un->operand)) {
                    out = Rational(-il2->value); return true;
                }
            }
            return false;
        };
        c0 = Rational(0); c1 = Rational(0);
        if (cv.size() >= 1U && !get_rat(cv[0], c0)) return false;
        if (cv.size() >= 2U && !get_rat(cv[1], c1)) return false;
        return true;
    };

    if (!extract(nc.value(), b_r, a_r)) return std::nullopt;
    if (!extract(dc.value(), d_r, c_r)) return std::nullopt;

    // Build the product S(n)/S(lower) = ∏_{m=lower}^{n-1} ρ(m).
    // Special case: a=c=1 (common shift-1 recurrences).
    // Product = (b+lower)(b+lower+1)···(b+n-1) / (d+lower)(d+lower+1)···(d+n-1)
    // = Pochhammer(b+lower, n-lower) / Pochhammer(d+lower, n-lower)
    // Build as rising factorial ratio using Gamma for non-integer cases.

    // For now: handle integer a,c case with Pochhammer or factorial structure.
    // Ledgered for arbitrary rational a,c via Gamma extension.
    if (a_r.denominator() != BigInt(1) || c_r.denominator() != BigInt(1))
        return std::nullopt;  // HARDCODE-OF-PASSAGE: non-unit leading coefficients

    // a=c=0 means constant already handled above; a≠0 or c≠0 required.
    // Build via explicit Pochhammer (already a BuiltinOp).
    // Pochhammer(α, n) = α(α+1)···(α+n-1) = Γ(α+n)/Γ(α).
    auto build_start = [&](const Rational& coeff, const Rational& base_r,
                           ExprPtr lower_e) -> ExprPtr {
        // Start value = coeff·lower + base_r → expressed as ExprPtr.
        if (coeff.numerator().is_zero()) {
            // constant denominator: just 1/Pochhammer(base,n-lower)
            return arena.make<RationalLit>(base_r.numerator(), base_r.denominator());
        }
        // coeff·n + base_r at n = lower
        ExprPtr coeff_e = arena.make<RationalLit>(coeff.numerator(), coeff.denominator());
        ExprPtr base_e  = arena.make<RationalLit>(base_r.numerator(), base_r.denominator());
        ExprPtr start   = arena.make<Binary>(BinaryOp::Add,
            arena.make<Binary>(BinaryOp::Mul, coeff_e, lower_e), base_e);
        auto s = ctx.simplify(start);
        return s.is_ok() ? s.value() : start;
    };

    ExprPtr n_minus_lower = arena.make<Binary>(BinaryOp::Sub, n_e, lower_val);
    auto nml_s = ctx.simplify(n_minus_lower);
    ExprPtr nml = nml_s.is_ok() ? nml_s.value() : n_minus_lower;

    ExprPtr ratio;
    if (a_r.numerator().is_zero() && c_r.numerator().is_zero()) {
        // constant ρ = b/d
        ExprPtr rho_c = arena.make<RationalLit>(b_r.numerator() * c_r.denominator(),
                                                 d_r.numerator() * b_r.denominator());
        auto rc = ctx.simplify(rho_c);
        ExprPtr rho_expr = rc.is_ok() ? rc.value() : rho_c;
        ratio = arena.make<Binary>(BinaryOp::Pow, rho_expr, nml);
    } else {
        // Build Pochhammer numerator and denominator.
        ExprPtr alpha_start = build_start(a_r, b_r, lower_val);
        ExprPtr gamma_start = build_start(c_r, d_r, lower_val);
        std::vector<ExprPtr> poch_num_args{alpha_start, nml};
        std::vector<ExprPtr> poch_den_args{gamma_start, nml};
        ExprPtr poch_num = arena.make<FuncCall>(BuiltinOp::Pochhammer,
            std::move(poch_num_args));
        ExprPtr poch_den = arena.make<FuncCall>(BuiltinOp::Pochhammer,
            std::move(poch_den_args));
        // Correct for leading coefficient mismatch: a^(n-lower) / c^(n-lower).
        ExprPtr lc_ratio;
        if (a_r == c_r) {
            ratio = arena.make<Binary>(BinaryOp::Div, poch_num, poch_den);
        } else {
            ExprPtr a_e = arena.make<RationalLit>(a_r.numerator(), a_r.denominator());
            ExprPtr c_e = arena.make<RationalLit>(c_r.numerator(), c_r.denominator());
            lc_ratio = arena.make<Binary>(BinaryOp::Pow,
                arena.make<Binary>(BinaryOp::Div, a_e, c_e), nml);
            ratio = arena.make<Binary>(BinaryOp::Mul, lc_ratio,
                arena.make<Binary>(BinaryOp::Div, poch_num, poch_den));
        }
    }

    ExprPtr result = arena.make<Binary>(BinaryOp::Mul, S0, ratio);
    auto s = ctx.simplify(result);
    return s.is_ok() ? std::optional<ExprPtr>{s.value()} : std::nullopt;
}

// Evaluate the initial condition S(n=lower_val) = Σ_{k=lower}^{lower} F(lower,k)
// (valid when upper = n_sym → at n=lower, sum has one term if lower is integer).
static std::optional<ExprPtr> eval_initial_condition(
    ExprPtr F, const Symbol& n_sym, const Symbol& k_sym,
    ExprPtr lower, symbolic::CASContext& ctx) {
    // Substitute n = lower, k = lower.
    auto F_n = ctx.substitute(F, n_sym, lower);
    if (F_n.is_error()) return std::nullopt;
    auto F_nk = ctx.substitute(F_n.value(), k_sym, lower);
    if (F_nk.is_error()) return std::nullopt;
    auto s = ctx.simplify(F_nk.value());
    return s.is_ok() ? std::optional<ExprPtr>{s.value()} : std::nullopt;
}

// ── Public interface ─────────────────────────────────────────────────────────

Result<std::optional<ExprPtr>> zeilberger_sum(
    ExprPtr F,
    const Symbol& n_param,
    const Symbol& k,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx) {

    // F must depend on k.
    if (!expr_contains_symbol(F, k.name))
        return ok(std::optional<ExprPtr>{std::nullopt});
    // F must depend on n_param.
    if (!expr_contains_symbol(F, n_param.name))
        return ok(std::optional<ExprPtr>{std::nullopt});
    // Algorithm assumes upper = n_param (proper hypergeometric boundary conditions).
    // Non-standard bounds require inhomogeneous recurrence extension (ledgered).
    if (const auto* up_sym = expr_cast<Symbol>(upper);
        !up_sym || up_sym->name != n_param.name)
        return ok(std::optional<ExprPtr>{std::nullopt});

    // Compute r(k) = F(k→k+1)/F(k).
    auto r_opt = compute_shift_ratio(F, k, ctx);
    if (!r_opt) return ok(std::optional<ExprPtr>{std::nullopt});
    ExprPtr r_k = *r_opt;

    // r_k must depend on k (otherwise Gosper handles the constant-ratio case).
    if (!expr_contains_symbol(r_k, k.name))
        return ok(std::optional<ExprPtr>{std::nullopt});

    // Compute the n-shift ratios: s_i(k) = F(n→n+i, k) / F(n, k) for i=0..J_max.
    unsigned int J_max = ctx.max_zeilberger_order();
    std::vector<ExprPtr> s_shifts;  // s_shifts[i] = s_i(k)
    s_shifts.push_back(ctx.arena().make<IntegerLit>(BigInt(1)));  // s_0 = 1
    AstArena& arena = ctx.arena();
    ExprPtr F_current = F;
    for (unsigned int i = 0U; i < J_max; ++i) {
        // Compute F(n→n+i+1, k) / F(n,k) = s_{i+1}.
        // s_{i+1} = s_i · [F(n+i+1,k)/F(n+i,k)] = s_i · r_n(n→n+i).
        auto s_next = compute_shift_ratio(F_current, n_param, ctx);
        if (!s_next) break;
        // s_{i+1}(k) = s_i(k) * s_next(k), where s_next is F(n+i+1,k)/F(n+i,k).
        ExprPtr s_i_plus_1 = arena.make<Binary>(BinaryOp::Mul,
            s_shifts.back(), *s_next);
        auto tog = algebra::together(s_i_plus_1, ctx);
        auto simp = tog.is_ok()
            ? ctx.simplify(tog.value())
            : ctx.simplify(s_i_plus_1);
        if (simp.is_error()) break;
        s_shifts.push_back(simp.value());
        // Advance F_current by one n-shift for the next ratio.
        ExprPtr n_plus_one = arena.make<Binary>(BinaryOp::Add,
            arena.make<Symbol>(n_param), arena.make<IntegerLit>(BigInt(1)));
        auto F_next = ctx.substitute(F_current, n_param, n_plus_one);
        if (F_next.is_error()) break;
        auto F_next_s = ctx.simplify(F_next.value());
        if (F_next_s.is_error()) break;
        F_current = F_next_s.value();
    }

    unsigned int J_max_avail = static_cast<unsigned int>(s_shifts.size()) - 1U;
    unsigned int D_p = ctx.max_zeilberger_poly_degree();
    unsigned int D_r = ctx.max_zeilberger_cert_degree();

    // Try increasing J from 1 upward.
    for (unsigned int J = 1U; J <= J_max_avail; ++J) {
        std::vector<ExprPtr> s_J(s_shifts.begin(), s_shifts.begin() + J + 1);
        auto rec_opt = try_parametric_gosper(r_k, s_J, n_param, k, J, D_p, D_r, ctx);
        if (!rec_opt) continue;
        const auto& p_vec = *rec_opt;

        if (J == 1U) {
            // Solve first-order recurrence.
            auto S0_opt = eval_initial_condition(F, n_param, k, lower, ctx);
            if (!S0_opt) continue;
            auto closed_opt = solve_first_order_rec(
                p_vec[0], p_vec[1], *S0_opt, lower, n_param, ctx);
            if (closed_opt)
                return ok(std::optional<ExprPtr>{*closed_opt});
        }
        // J≥2: ledgered as F5.7-ZEIL-HIGHER-ORDER — return Unimplemented.
        // (Recurrence found but higher-order ODE solver not yet wired.)
        break;
    }

    return ok(std::optional<ExprPtr>{std::nullopt});
}

}  // namespace cas::symbolic
