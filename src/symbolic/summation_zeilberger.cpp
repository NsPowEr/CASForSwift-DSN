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
#include "summation_hyper.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include <cstdio>
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
static std::optional<std::pair<ExprPtr, ExprPtr>> compute_shift_ratio(
    ExprPtr F, const Symbol& sym, symbolic::CASContext& ctx, long long delta = 1) {
    return zeilberger_detail::compute_shift_ratio(F, sym, ctx, delta);
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
    auto simp_expr = ctx.simplify(expr);
    ExprPtr e = simp_expr.is_ok() ? simp_expr.value() : expr;
    auto exp = algebra::expand(e, ctx);
    if (exp.is_error()) { 
        eqs_out.push_back(e); 
        return; 
    }
    auto k_coeffs = algebra::univariate_coefficients(exp.value(), sym, ctx);
    if (k_coeffs.is_error()) { 
        eqs_out.push_back(exp.value()); 
        return; 
    }
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
    const algebra::RationalParts& r_parts,           // r(k) = N_r / D_r_poly
    const std::vector<algebra::RationalParts>& s_parts, // s_i(k) = N_s_i / D_s_i
    const Symbol& n_sym, const Symbol& k_sym,
    unsigned int J, unsigned int D_p, unsigned int D_r,
    symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();

    ExprPtr N_r = r_parts.numerator;
    ExprPtr D_r_poly = r_parts.denominator;

    // Decompose the combined h(k) = Σ p_i · s_i(k) symbolically.
    // Build p_i polynomials in n with unknown coefficients.
    std::vector<ExprPtr> p_polys;
    std::vector<ExprPtr> all_p_unknowns;
    for (unsigned int i = 0U; i <= J; ++i) {
        auto [pi, pi_unknowns] = make_poly_unknown(n_sym, D_p, "zp", ctx);
        p_polys.push_back(pi);
        all_p_unknowns.insert(all_p_unknowns.end(), pi_unknowns.begin(), pi_unknowns.end());
    }

    // Build h(k) = Σ p_i * (N_s_i / D_s_i) as a single rational expression.
    ExprPtr common_den_h = arena.make<IntegerLit>(BigInt(1));
    for (unsigned int i = 0U; i <= J; ++i) {
        common_den_h = arena.make<Binary>(BinaryOp::Mul, common_den_h, s_parts[i].denominator);
    }
    auto cd_simp = ctx.simplify(common_den_h);
    if (cd_simp.is_ok()) common_den_h = cd_simp.value();

    ExprPtr N_h = arena.make<IntegerLit>(BigInt(0));
    for (unsigned int i = 0U; i <= J; ++i) {
        ExprPtr ratio = arena.make<Binary>(BinaryOp::Div, common_den_h, s_parts[i].denominator);
        auto rs = ctx.simplify(ratio);
        ExprPtr term_i = arena.make<Binary>(BinaryOp::Mul, p_polys[i],
            arena.make<Binary>(BinaryOp::Mul, s_parts[i].numerator,
                rs.is_ok() ? rs.value() : ratio));
        N_h = arena.make<Binary>(BinaryOp::Add, N_h, term_i);
    }
    ExprPtr D_h = common_den_h;
    
    auto [x_k, x_unknowns] = make_poly_unknown(k_sym, D_r, "zr", ctx);
    
    ExprPtr k_plus_one = arena.make<Binary>(BinaryOp::Add,
        arena.make<Symbol>(k_sym), arena.make<IntegerLit>(BigInt(1)));
    auto x_k1_res = ctx.substitute(x_k, k_sym, k_plus_one);
    if (x_k1_res.is_error()) return std::nullopt;
    ExprPtr x_k1 = x_k1_res.value();
    auto D_h_k1_res = ctx.substitute(D_h, k_sym, k_plus_one);
    if (D_h_k1_res.is_error()) return std::nullopt;
    ExprPtr D_h_k1 = D_h_k1_res.value();

    ExprPtr lhs_t1 = arena.make<Binary>(BinaryOp::Mul,
        arena.make<Binary>(BinaryOp::Mul, x_k1, N_r), D_h);
    ExprPtr lhs_t2 = arena.make<Binary>(BinaryOp::Mul,
        arena.make<Binary>(BinaryOp::Mul, x_k, D_h_k1), D_r_poly);
    ExprPtr rhs    = arena.make<Binary>(BinaryOp::Mul,
        arena.make<Binary>(BinaryOp::Mul, N_h, D_h_k1), D_r_poly);

    ExprPtr equation = arena.make<Binary>(BinaryOp::Sub,
        arena.make<Binary>(BinaryOp::Sub, lhs_t1, lhs_t2), rhs);
    
    std::vector<ExprPtr> eqs;
    collect_equations(equation, k_sym, n_sym, ctx, eqs);
    if (eqs.empty()) return std::nullopt;

    std::vector<ExprPtr> all_unknowns;
    all_unknowns.insert(all_unknowns.end(), all_p_unknowns.begin(), all_p_unknowns.end());
    all_unknowns.insert(all_unknowns.end(), x_unknowns.begin(), x_unknowns.end());

    // The system is homogeneous. csolve will just return 0.
    // We want a non-trivial solution where p_J != 0.
    // We probe by setting one coefficient of p_J to 1.
    std::size_t pJ_start_idx = J * (D_p + 1);
    std::size_t pJ_end_idx   = pJ_start_idx + D_p + 1;
    
    std::vector<ExprPtr> sol_vals;
    bool found_sol = false;

    for (std::size_t probe_idx = pJ_start_idx; probe_idx < pJ_end_idx && probe_idx < all_unknowns.size(); ++probe_idx) {
        ExprPtr probe_var = all_unknowns[probe_idx];
        ExprPtr one = arena.make<IntegerLit>(BigInt(1));
        
        std::vector<ExprPtr> sub_eqs;
        for (ExprPtr eq : eqs) {
            auto sub_eq = ctx.substitute(eq, *expr_cast<Symbol>(probe_var), one);
            sub_eqs.push_back(sub_eq.is_ok() ? sub_eq.value() : eq);
        }
        
        std::vector<ExprPtr> sub_vars;
        for (std::size_t i = 0; i < all_unknowns.size(); ++i) {
            if (i != probe_idx) sub_vars.push_back(all_unknowns[i]);
        }
        
        ExprPtr eqs_mat = arena.make<Matrix>(sub_eqs.size(), 1U, sub_eqs);
        ExprPtr vars_mat = arena.make<Matrix>(sub_vars.size(), 1U, sub_vars);
        
        auto sol_res = algebra::csolve(eqs_mat, vars_mat, ctx);
        if (sol_res.is_ok()) {
            if (const auto* sol_mat = expr_cast<Matrix>(sol_res.value())) {
                if (!sol_mat->elements.empty()) {
                    // Reconstruct full solution vector.
                    sol_vals.resize(all_unknowns.size());
                    std::size_t sol_idx = 0;
                    for (std::size_t i = 0; i < all_unknowns.size(); ++i) {
                        if (i == probe_idx) {
                            sol_vals[i] = one;
                        } else {
                            sol_vals[i] = sol_mat->elements[sol_idx++];
                        }
                    }
                    found_sol = true;
                    break;
                }
            }
        }
    }
    
    if (!found_sol) return std::nullopt;

    std::vector<ExprPtr> result_p;
    std::size_t n_p_unknowns = all_p_unknowns.size();
    for (unsigned int i = 0U; i <= J; ++i) {
        ExprPtr pi = p_polys[i];
        for (std::size_t c = 0U; c < n_p_unknowns && c < sol_vals.size(); ++c) {
            auto sub = ctx.substitute(pi, *expr_cast<Symbol>(all_p_unknowns[c]),
                sol_vals[c]);
            if (sub.is_ok()) pi = sub.value();
        }
        auto simp = ctx.simplify(pi);
        result_p.push_back(simp.is_ok() ? simp.value() : pi);
    }

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

    if (a_r.denominator() != BigInt(1) || c_r.denominator() != BigInt(1))
        return std::nullopt;

    auto build_start = [&](const Rational& coeff, const Rational& base_r,
                           ExprPtr lower_e) -> ExprPtr {
        if (coeff.numerator().is_zero()) {
            return arena.make<RationalLit>(base_r.numerator(), base_r.denominator());
        }
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
        ExprPtr rho_c = arena.make<RationalLit>(b_r.numerator() * c_r.denominator(),
                                                 d_r.numerator() * b_r.denominator());
        auto rc = ctx.simplify(rho_c);
        ExprPtr rho_expr = rc.is_ok() ? rc.value() : rho_c;
        ratio = arena.make<Binary>(BinaryOp::Pow, rho_expr, nml);
    } else {
        ExprPtr alpha_start = build_start(a_r, b_r, lower_val);
        ExprPtr gamma_start = build_start(c_r, d_r, lower_val);
        std::vector<ExprPtr> poch_num_args{alpha_start, nml};
        std::vector<ExprPtr> poch_den_args{gamma_start, nml};
        ExprPtr poch_num = arena.make<FuncCall>(BuiltinOp::Pochhammer,
            std::move(poch_num_args));
        ExprPtr poch_den = arena.make<FuncCall>(BuiltinOp::Pochhammer,
            std::move(poch_den_args));
        if (a_r == c_r) {
            ratio = arena.make<Binary>(BinaryOp::Div, poch_num, poch_den);
        } else {
            ExprPtr a_e = arena.make<RationalLit>(a_r.numerator(), a_r.denominator());
            ExprPtr c_e = arena.make<RationalLit>(c_r.numerator(), c_r.denominator());
            ExprPtr lc_ratio = arena.make<Binary>(BinaryOp::Pow,
                arena.make<Binary>(BinaryOp::Div, a_e, c_e), nml);
            ratio = arena.make<Binary>(BinaryOp::Mul, lc_ratio,
                arena.make<Binary>(BinaryOp::Div, poch_num, poch_den));
        }
    }

    ExprPtr result = arena.make<Binary>(BinaryOp::Mul, S0, ratio);
    auto s = ctx.simplify(result);
    return s.is_ok() ? std::optional<ExprPtr>{s.value()} : std::nullopt;
}

static std::optional<ExprPtr> eval_initial_condition(
    ExprPtr F, const Symbol& n_sym, const Symbol& k_sym,
    ExprPtr lower, symbolic::CASContext& ctx) {
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

    if (!expr_contains_symbol(F, k.name)) return ok(std::optional<ExprPtr>{std::nullopt});
    if (!expr_contains_symbol(F, n_param.name)) return ok(std::optional<ExprPtr>{std::nullopt});
    if (const auto* up_sym = expr_cast<Symbol>(upper);
        !up_sym || up_sym->name != n_param.name) return ok(std::optional<ExprPtr>{std::nullopt});

    auto r_pair_opt = compute_shift_ratio(F, k, ctx);
    if (!r_pair_opt) return ok(std::optional<ExprPtr>{std::nullopt});
    algebra::RationalParts r_parts = { .numerator = r_pair_opt->first, .denominator = r_pair_opt->second };

    unsigned int J_max = ctx.max_zeilberger_order();
    std::vector<algebra::RationalParts> s_parts;
    s_parts.push_back({
        .numerator = ctx.arena().make<IntegerLit>(BigInt(1)),
        .denominator = ctx.arena().make<IntegerLit>(BigInt(1))
    });

    for (unsigned int i = 1U; i <= J_max; ++i) {
        auto si_pair_opt = compute_shift_ratio(F, n_param, ctx, static_cast<long long>(i));
        if (!si_pair_opt) break;
        s_parts.push_back({ .numerator = si_pair_opt->first, .denominator = si_pair_opt->second });
    }

    unsigned int J_max_avail = static_cast<unsigned int>(s_parts.size()) - 1U;
    unsigned int D_p = ctx.max_zeilberger_poly_degree();
    unsigned int D_r = ctx.max_zeilberger_cert_degree();

    for (unsigned int J = 1U; J <= J_max_avail; ++J) {
        std::vector<algebra::RationalParts> s_J(s_parts.begin(), s_parts.begin() + J + 1);
        auto rec_opt = try_parametric_gosper(r_parts, s_J, n_param, k, J, D_p, D_r, ctx);
        if (!rec_opt) continue;
        const auto& p_vec = *rec_opt;

        if (J == 1U) {
            auto S0_opt = eval_initial_condition(F, n_param, k, lower, ctx);
            if (!S0_opt) continue;
            auto closed_opt = solve_first_order_rec(
                p_vec[0], p_vec[1], *S0_opt, lower, n_param, ctx);
            if (closed_opt) return ok(std::optional<ExprPtr>{*closed_opt});
            break;
        }

        // J ≥ 2: solve the recurrence into a verified Petkovšek closed form
        // (a ℚ-linear combination of hypergeometric terms), cross-checked
        // against directly-computed sums.  ok(nullopt) ⇒ no closed form.
        auto closed = sum_closed_form_from_recurrence(p_vec, F, n_param, k, lower, ctx);
        if (closed.is_ok() && closed.value().has_value())
            return ok(std::optional<ExprPtr>{*closed.value()});
        break;
    }

    return ok(std::optional<ExprPtr>{std::nullopt});
}

}  // namespace cas::symbolic
