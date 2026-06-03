#include "summation_gosper.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include <algorithm>

namespace cas::symbolic {

Result<std::optional<ExprPtr>> gosper_sum(
    const ExprPtr& term, 
    const Symbol& k, 
    symbolic::CASContext& ctx) 
{
    AstArena& arena = ctx.arena();
    auto k_plus_1 = linalg::add_expr(ctx, arena.make<Symbol>(k.name), linalg::integer(ctx, 1)).value();
    
    auto t_next_res = substitute(term, k, k_plus_1, ctx);
    if (t_next_res.is_error()) return fail<std::optional<ExprPtr>>(t_next_res.error());
    auto t_next = t_next_res.value();

    auto ratio = linalg::div_expr(ctx, t_next, term).value();
    
    auto r_k_res = simplify(ratio, ctx);
    if (r_k_res.is_error()) return fail<std::optional<ExprPtr>>(r_k_res.error());
    auto r_k = r_k_res.value();

    // A/B
    auto parts_res = algebra::apart_num_den(r_k, ctx);
    if (parts_res.is_error()) return fail<std::optional<ExprPtr>>(parts_res.error());
    auto parts = parts_res.value();
    
    ExprPtr A = parts.numerator;
    ExprPtr B = parts.denominator;

    auto j_sym = arena.make<Symbol>(ctx.make_fresh_symbol("j"));
    auto k_plus_j = linalg::add_expr(ctx, arena.make<Symbol>(k.name), j_sym).value();
    
    auto B_k_plus_j_res = substitute(B, k, k_plus_j, ctx);
    if (B_k_plus_j_res.is_error()) return fail<std::optional<ExprPtr>>(B_k_plus_j_res.error());
    auto B_k_plus_j = B_k_plus_j_res.value();
    
    auto res_j_res = algebra::polynomial_resultant(A, B_k_plus_j, k, ctx);
    if (res_j_res.is_error()) return fail<std::optional<ExprPtr>>(res_j_res.error());
    auto res_j = res_j_res.value();
    
    auto roots_res = algebra::solve_polynomial(res_j, *expr_cast<Symbol>(j_sym), ctx);
    if (roots_res.is_error()) return fail<std::optional<ExprPtr>>(roots_res.error());
    auto roots = roots_res.value();
    
    std::vector<int> J;
    for (const auto& root : roots) {
        if (const auto* i = expr_cast<IntegerLit>(root)) {
            if (!i->value.is_negative()) {
                J.push_back(static_cast<int>(i->value.to_u64()));
            }
        }
    }
    std::sort(J.rbegin(), J.rend());

    ExprPtr p = linalg::integer(ctx, 1);
    ExprPtr q = A;
    ExprPtr r = B;

    for (int j : J) {
        auto k_plus_j_val = linalg::add_expr(ctx, arena.make<Symbol>(k.name), linalg::integer(ctx, j)).value();
        
        auto r_shift_res = substitute(r, k, k_plus_j_val, ctx);
        if (r_shift_res.is_error()) return fail<std::optional<ExprPtr>>(r_shift_res.error());
        auto r_shift = r_shift_res.value();
        
        auto g_res = algebra::polynomial_gcd(q, r_shift, k, ctx);
        if (g_res.is_error()) return fail<std::optional<ExprPtr>>(g_res.error());
        auto g = g_res.value();
        
        auto deg_g_res = algebra::polynomial_degree(g, k, ctx);
        if (deg_g_res.is_error()) return fail<std::optional<ExprPtr>>(deg_g_res.error());
        auto deg_g = deg_g_res.value();
        
        if (deg_g == 0) continue;

        auto q_div_res = algebra::polynomial_exact_divide(q, g, k, ctx);
        if (q_div_res.is_error()) return fail<std::optional<ExprPtr>>(q_div_res.error());
        q = q_div_res.value();

        auto k_minus_j = linalg::sub_expr(ctx, arena.make<Symbol>(k.name), linalg::integer(ctx, j)).value();
        
        auto g_shift_back_res = substitute(g, k, k_minus_j, ctx);
        if (g_shift_back_res.is_error()) return fail<std::optional<ExprPtr>>(g_shift_back_res.error());
        auto g_shift_back = g_shift_back_res.value();
        
        auto r_div_res = algebra::polynomial_exact_divide(r, g_shift_back, k, ctx);
        if (r_div_res.is_error()) return fail<std::optional<ExprPtr>>(r_div_res.error());
        r = r_div_res.value();

        // Petkovšek polynomial-part update: p(k) *= ∏_{i=1..j} g(k-i).
        // (Bug fix: loop was 0..j inclusive, producing one extra factor and
        // mis-decomposing t_{k+1}/t_k. Verified against t_k=k → s=k(k-1)/2.)
        for (int i = 1; i <= j; ++i) {
            auto k_minus_i = linalg::sub_expr(ctx, arena.make<Symbol>(k.name), linalg::integer(ctx, i)).value();

            auto g_i_res = substitute(g, k, k_minus_i, ctx);
            if (g_i_res.is_error()) return fail<std::optional<ExprPtr>>(g_i_res.error());
            auto g_i = g_i_res.value();

            p = linalg::mul_expr(ctx, p, g_i).value();
        }
        
        auto p_simp_res = simplify(p, ctx);
        if (p_simp_res.is_error()) return fail<std::optional<ExprPtr>>(p_simp_res.error());
        p = p_simp_res.value();
    }

    auto q_coeffs_res = algebra::univariate_coefficients(q, k, ctx);
    if (q_coeffs_res.is_error()) return fail<std::optional<ExprPtr>>(q_coeffs_res.error());
    auto q_coeffs = q_coeffs_res.value();
    
    auto r_coeffs_res = algebra::univariate_coefficients(r, k, ctx);
    if (r_coeffs_res.is_error()) return fail<std::optional<ExprPtr>>(r_coeffs_res.error());
    auto r_coeffs = r_coeffs_res.value();
    
    auto p_coeffs_res = algebra::univariate_coefficients(p, k, ctx);
    if (p_coeffs_res.is_error()) return fail<std::optional<ExprPtr>>(p_coeffs_res.error());
    auto p_coeffs = p_coeffs_res.value();

    int D_q = std::max(0, static_cast<int>(q_coeffs.size()) - 1);
    int D_r = std::max(0, static_cast<int>(r_coeffs.size()) - 1);
    int D_p = std::max(0, static_cast<int>(p_coeffs.size()) - 1);

    long long d = -1;
    if (D_q != D_r) {
        d = D_p - std::max(D_q, D_r);
    } else {
        int D = D_q;
        ExprPtr q_D = D < static_cast<int>(q_coeffs.size()) ? q_coeffs[D] : linalg::integer(ctx, 0);
        ExprPtr r_D = D < static_cast<int>(r_coeffs.size()) ? r_coeffs[D] : linalg::integer(ctx, 0);
        
        auto q_D_eq_r_D_res = mathematically_equal(q_D, r_D, ctx);
        if (q_D_eq_r_D_res.is_error()) return fail<std::optional<ExprPtr>>(q_D_eq_r_D_res.error());
        bool q_D_eq_r_D = q_D_eq_r_D_res.value();
        
        long long d1 = D_p - D + 1;
        
        if (!q_D_eq_r_D) {
            d = d1;
        } else {
            ExprPtr q_Dm1 = (D >= 1 && D - 1 < static_cast<int>(q_coeffs.size())) ? q_coeffs[D-1] : linalg::integer(ctx, 0);
            ExprPtr r_Dm1 = (D >= 1 && D - 1 < static_cast<int>(r_coeffs.size())) ? r_coeffs[D-1] : linalg::integer(ctx, 0);
            
            auto diff = linalg::sub_expr(ctx, r_Dm1, q_Dm1).value();
            auto frac = linalg::div_expr(ctx, diff, q_D).value();
            
            auto frac_simp_res = simplify(frac, ctx);
            if (frac_simp_res.is_error()) return fail<std::optional<ExprPtr>>(frac_simp_res.error());
            auto frac_simp = frac_simp_res.value();
            
            long long d2 = -1;
            if (const auto* frac_int = expr_cast<IntegerLit>(frac_simp)) {
                d2 = static_cast<long long>(frac_int->value.to_u64()) - D;
            } else if (const auto* frac_rat = expr_cast<RationalLit>(frac_simp)) {
                if (frac_rat->denominator == BigInt(1)) {
                    d2 = static_cast<long long>(frac_rat->numerator.to_u64()) - D;
                }
            }
            d = std::max(d1, d2);
        }
    }

    if (d < 0) {
        return ok(std::optional<ExprPtr>{std::nullopt});
    }

    std::vector<ExprPtr> u_vars;
    ExprPtr x_k = linalg::integer(ctx, 0);
    for (long long i = 0; i <= d; ++i) {
        auto u_i = arena.make<Symbol>(ctx.make_fresh_symbol("u"));
        u_vars.push_back(u_i);
        auto pow_k = arena.make<Binary>(BinaryOp::Pow, arena.make<Symbol>(k.name), linalg::integer(ctx, i));
        auto term_i = linalg::mul_expr(ctx, u_i, pow_k).value();
        x_k = linalg::add_expr(ctx, x_k, term_i).value();
    }
    
    auto x_k_simp_res = simplify(x_k, ctx);
    if (x_k_simp_res.is_error()) return fail<std::optional<ExprPtr>>(x_k_simp_res.error());
    x_k = x_k_simp_res.value();
    
    auto x_k_plus_1_sub_res = substitute(x_k, k, k_plus_1, ctx);
    if (x_k_plus_1_sub_res.is_error()) return fail<std::optional<ExprPtr>>(x_k_plus_1_sub_res.error());
    auto x_k_plus_1 = x_k_plus_1_sub_res.value();
    
    auto x_k_plus_1_simp_res = simplify(x_k_plus_1, ctx);
    if (x_k_plus_1_simp_res.is_error()) return fail<std::optional<ExprPtr>>(x_k_plus_1_simp_res.error());
    x_k_plus_1 = x_k_plus_1_simp_res.value();

    auto k_minus_1 = linalg::sub_expr(ctx, arena.make<Symbol>(k.name), linalg::integer(ctx, 1)).value();
    
    auto r_k_minus_1_res = substitute(r, k, k_minus_1, ctx);
    if (r_k_minus_1_res.is_error()) return fail<std::optional<ExprPtr>>(r_k_minus_1_res.error());
    auto r_k_minus_1 = r_k_minus_1_res.value();

    auto t1 = linalg::mul_expr(ctx, q, x_k_plus_1).value();
    auto t2 = linalg::mul_expr(ctx, r_k_minus_1, x_k).value();
    auto L_term = linalg::sub_expr(ctx, t1, t2).value();
    auto L = linalg::sub_expr(ctx, L_term, p).value();
    
    auto L_simp_res = simplify(L, ctx);
    if (L_simp_res.is_error()) return fail<std::optional<ExprPtr>>(L_simp_res.error());
    auto L_simp = L_simp_res.value();
    
    auto L_expand_res = algebra::expand(L_simp, ctx);
    if (L_expand_res.is_error()) return fail<std::optional<ExprPtr>>(L_expand_res.error());
    auto L_expand = L_expand_res.value();
    
    auto L_coeffs_res = algebra::univariate_coefficients(L_expand, k, ctx);
    if (L_coeffs_res.is_error()) return fail<std::optional<ExprPtr>>(L_coeffs_res.error());
    auto L_coeffs = L_coeffs_res.value();

    std::vector<ExprPtr> eqs;
    for (auto c : L_coeffs) {
        eqs.push_back(c);
    }

    auto eqs_list = arena.make<Matrix>(eqs.size(), 1, eqs);
    auto vars_list = arena.make<Matrix>(u_vars.size(), 1, u_vars);

    auto sol_res = algebra::csolve(eqs_list, vars_list, ctx);
    if (sol_res.is_error()) return fail<std::optional<ExprPtr>>(sol_res.error());
    auto sol = sol_res.value();
    
    const auto* sol_matrix = expr_cast<Matrix>(sol);
    if (!sol_matrix || sol_matrix->elements.empty()) {
        return ok(std::optional<ExprPtr>{std::nullopt});
    }

    ExprPtr x_sol = x_k;
    for (size_t c = 0; c < u_vars.size() && c < sol_matrix->cols; ++c) {
        x_sol = substitute(x_sol, *expr_cast<Symbol>(u_vars[c]), sol_matrix->elements[c], ctx).value();
    }
    
    auto x_sol_simp_res = simplify(x_sol, ctx);
    if (x_sol_simp_res.is_error()) return fail<std::optional<ExprPtr>>(x_sol_simp_res.error());
    x_sol = x_sol_simp_res.value();

    // Petkovšek closing formula: s(k) = r(k-1)·x(k)·t(k)/p(k).
    // Divide x(k) by p(k) FIRST as exact polynomial division (Gosper's theorem
    // guarantees p | x): the Binary(Div) simplifier currently drops rational
    // coefficients on Div of rational-coeff polynomial by polynomial. Routing
    // through polynomial_exact_divide preserves Q[k] coefficient field.
    auto x_div_p_res = algebra::polynomial_exact_divide(x_sol, p, k, ctx);
    ExprPtr x_div_p;
    if (x_div_p_res.is_ok()) {
        x_div_p = x_div_p_res.value();
    } else {
        // Fallback: residual rational expression (still subject to upstream
        // Binary(Div) limitation). Use linalg::div_expr as last resort.
        auto fallback = linalg::div_expr(ctx, x_sol, p);
        if (fallback.is_error()) return fail<std::optional<ExprPtr>>(fallback.error());
        x_div_p = fallback.value();
    }
    auto s_term1 = linalg::mul_expr(ctx, r_k_minus_1, x_div_p).value();
    auto s = linalg::mul_expr(ctx, s_term1, term).value();
    auto s_simp_res = simplify(s, ctx);
    if (s_simp_res.is_error()) return fail<std::optional<ExprPtr>>(s_simp_res.error());
    ExprPtr s_final = s_simp_res.value();
    // Canonicalize as single rational fraction: cancels structural sharing
    // duplicates and reduces common factors that Binary(Div) simplifier
    // misses (notably `(k+1)^2 / (k+1)` and Sum-of-rationals forms).
    auto s_tog = algebra::together(s_final, ctx);
    if (s_tog.is_ok()) {
        auto s_tog_simp = simplify(s_tog.value(), ctx);
        if (s_tog_simp.is_ok()) s_final = s_tog_simp.value();
    }

    // F5.7-GOSPER-K2-NORMALIZATION fix — `csolve` on the over-determined
    // rational system at lines 218-226 occasionally returns the u_i vector
    // pre-scaled by the LCM of denominators (verified for term = k²:
    // expected u = (1/6, −1/2, 1/3), returned u = (1, −3, 2)).
    //
    // Run the rescaling check ONLY when `term` is a polynomial in `k` (the
    // case where the csolve over-scaling has been observed).  For
    // non-polynomial Gosper inputs the heavy substitute + simplify pass
    // would add measurable cost across long-running suites without ever
    // triggering the rescale, so we cheaply gate on
    // `algebra::univariate_coefficients(term, k)` succeeding first.
    {
        auto term_coeffs = algebra::univariate_coefficients(term, k, ctx);
        bool term_is_poly = term_coeffs.is_ok() && !term_coeffs.value().empty();
        if (term_is_poly) {
            auto s_next = substitute(s_final, k, k_plus_1, ctx);
            if (s_next.is_ok()) {
                auto delta_raw = linalg::sub_expr(ctx, s_next.value(), s_final);
                if (delta_raw.is_ok()) {
                    auto delta_simp = simplify(delta_raw.value(), ctx);
                    if (delta_simp.is_ok()) {
                        auto delta_exp = algebra::expand(delta_simp.value(), ctx);
                        auto term_exp  = algebra::expand(term, ctx);
                        if (delta_exp.is_ok() && term_exp.is_ok()) {
                            auto ratio = linalg::div_expr(ctx, delta_exp.value(), term_exp.value());
                            if (ratio.is_ok()) {
                                auto ratio_simp = simplify(ratio.value(), ctx);
                                if (ratio_simp.is_ok()) {
                                    ExprPtr r_val = ratio_simp.value();
                                    bool is_unit = false;
                                    if (const auto* il = expr_cast<IntegerLit>(r_val)) {
                                        is_unit = (il->value == BigInt(1));
                                    } else if (const auto* rl = expr_cast<RationalLit>(r_val)) {
                                        is_unit = (rl->numerator == BigInt(1) &&
                                                   rl->denominator == BigInt(1));
                                    }
                                    if (!is_unit) {
                                        bool is_rational_const =
                                            expr_cast<IntegerLit>(r_val) != nullptr ||
                                            expr_cast<RationalLit>(r_val) != nullptr;
                                        if (is_rational_const) {
                                            auto rescaled = linalg::div_expr(ctx, s_final, r_val);
                                            if (rescaled.is_ok()) {
                                                auto rescaled_simp = simplify(rescaled.value(), ctx);
                                                if (rescaled_simp.is_ok()) {
                                                    s_final = rescaled_simp.value();
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return ok(std::optional<ExprPtr>{s_final});
}

} // namespace cas::symbolic