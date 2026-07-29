#include "calculus_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"
#include "../algebra/polynomial_internal.hpp"
#include <vector>
#include <algorithm>

namespace cas::calculus {

namespace {

[[nodiscard]] Result<ExprPtr> solve_risch_de_field(
    ExprPtr f,
    ExprPtr g,
    std::size_t ext_idx,
    const DifferentialField& field,
    symbolic::CASContext& ctx) {
    
    AstArena& arena = ctx.arena();
    
    // Base case: Q(x)
    if (ext_idx == 0U) {
        return solve_risch_de_q(f, g, field.base_var(), ctx);
    }
    
    // Topmost extension at this level
    const auto& ext = field.extensions()[ext_idx - 1U];
    const Symbol& t = ext.t_var;
    
    // Compute Denominator bound: LCM of denominators of f and g in t
    auto f_parts_res = algebra::apart_num_den(f, ctx);
    if (f_parts_res.is_error()) return fail<ExprPtr>(f_parts_res.error());
    ExprPtr B = f_parts_res.value().denominator;
    
    auto g_parts_res = algebra::apart_num_den(g, ctx);
    if (g_parts_res.is_error()) return fail<ExprPtr>(g_parts_res.error());
    ExprPtr Q = g_parts_res.value().denominator;
    
    // D = LCM(Q, B) = (Q * B) / GCD(Q, B)
    auto gcd_res = algebra::polynomial_gcd(Q, B, t, ctx);
    if (gcd_res.is_error()) return fail<ExprPtr>(gcd_res.error());
    ExprPtr QB = arena.make<Binary>(BinaryOp::Mul, Q, B);
    ExprPtr D = arena.make<Binary>(BinaryOp::Div, QB, gcd_res.value());
    if (auto s = ctx.simplify(D); s.is_ok()) D = s.value();
    
    // f_new = f - D'/D,  g_new = g * D
    auto D_prime_res = field.derive(D, ctx);
    if (D_prime_res.is_error()) return fail<ExprPtr>(D_prime_res.error());
    ExprPtr D_prime = D_prime_res.value();
    
    ExprPtr D_prime_over_D = arena.make<Binary>(BinaryOp::Div, D_prime, D);
    ExprPtr f_new = arena.make<Binary>(BinaryOp::Sub, f, D_prime_over_D);
    ExprPtr g_new = arena.make<Binary>(BinaryOp::Mul, g, D);
    
    if (auto s = ctx.simplify(f_new); s.is_ok()) f_new = s.value();
    if (auto s = ctx.simplify(g_new); s.is_ok()) g_new = s.value();
    
    // Now solve polynomial Risch DE: p' + f_new * p = g_new
    auto f_poly_res = algebra::parse_polynomial(f_new, t, ctx);
    auto g_poly_res = algebra::parse_polynomial(g_new, t, ctx);
    if (f_poly_res.is_error() || g_poly_res.is_error()) {
        return make_unimplemented<ExprPtr>(
            "calculus", "solve_risch_de_field",
            "f or g is not polynomial after denominator bound",
            cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
            "Risch DE solver: make sure inputs are polynomials after bounding");
    }
    
    const auto& f_poly = f_poly_res.value();
    const auto& g_poly = g_poly_res.value();
    
    int df = static_cast<int>(poly_degree(f_poly));
    int dg = static_cast<int>(poly_degree(g_poly));
    
    // Handle leading-coefficient zero cases
    if (is_zero_poly(f_poly)) df = -1;
    if (is_zero_poly(g_poly)) dg = -1;
    
    // Degree bound N
    int N = 0;
    if (ext.type == ExtensionType::Logarithmic) {
        if (df > 0) {
            N = std::max(0, dg - df);
        } else if (df == 0) {
            N = dg;
        } else { // df == -1 (f_new = 0)
            N = dg + 1;
        }
    } else { // Exponential
        if (df > 0) {
            N = std::max(0, dg - df);
        } else if (df == 0) {
            N = std::max(dg, 0);
        } else { // df == -1 (f_new = 0)
            N = dg + 1;
        }
    }
    
    if (N < 0) N = 0;
    
    std::vector<ExprPtr> p_coeffs(static_cast<std::size_t>(N) + 1U, ExprPtr{});
    
    // If df > 0: solve algebraically top-down
    if (df > 0) {
        ExprPtr f_ld = leading_coefficient(f_poly);
        for (int i = N; i >= 0; --i) {
            std::vector<ExprPtr> p_known_coeffs(p_coeffs.begin() + i + 1, p_coeffs.end());
            p_known_coeffs.insert(p_known_coeffs.begin(), i + 1, arena.make<IntegerLit>(BigInt(0)));
            algebra::PolyExpr p_known(p_known_coeffs);
            
            auto p_known_expr_res = algebra::polynomial_to_expr(p_known, t, ctx);
            if (p_known_expr_res.is_error()) return fail<ExprPtr>(p_known_expr_res.error());
            ExprPtr pk = p_known_expr_res.value();
            
            auto Dpk_res = field.derive(pk, ctx);
            if (Dpk_res.is_error()) return fail<ExprPtr>(Dpk_res.error());
            ExprPtr Dpk = Dpk_res.value();
            
            ExprPtr term = arena.make<Binary>(BinaryOp::Add, Dpk, arena.make<Binary>(BinaryOp::Mul, f_new, pk));
            ExprPtr R = arena.make<Binary>(BinaryOp::Sub, g_new, term);
            if (auto s = ctx.simplify(R); s.is_ok()) R = s.value();
            
            auto R_poly_res = algebra::parse_polynomial(R, t, ctx);
            if (R_poly_res.is_error()) return fail<ExprPtr>(R_poly_res.error());
            
            ExprPtr R_coeff = R_poly_res.value()[static_cast<std::size_t>(i + df)];
            if (!R_coeff) R_coeff = arena.make<IntegerLit>(BigInt(0));
            
            ExprPtr p_i = arena.make<Binary>(BinaryOp::Div, R_coeff, f_ld);
            if (auto s = ctx.simplify(p_i); s.is_ok()) p_i = s.value();
            p_coeffs[static_cast<std::size_t>(i)] = p_i;
        }
        
        // Verify solution
        algebra::PolyExpr p_sol(p_coeffs);
        auto p_sol_expr_res = algebra::polynomial_to_expr(p_sol, t, ctx);
        if (p_sol_expr_res.is_error()) return fail<ExprPtr>(p_sol_expr_res.error());
        ExprPtr pk = p_sol_expr_res.value();
        
        auto Dpk_res = field.derive(pk, ctx);
        if (Dpk_res.is_error()) return fail<ExprPtr>(Dpk_res.error());
        ExprPtr Dpk = Dpk_res.value();
        
        ExprPtr LHS = arena.make<Binary>(BinaryOp::Add, Dpk, arena.make<Binary>(BinaryOp::Mul, f_new, pk));
        ExprPtr diff = arena.make<Binary>(BinaryOp::Sub, LHS, g_new);
        auto diff_tog = algebra::together(diff, ctx);
        ExprPtr diff_simp = diff;
        if (diff_tog.is_ok()) {
            if (auto s = ctx.simplify(diff_tog.value()); s.is_ok()) diff_simp = s.value();
        }
        
        bool verified = false;
        if (const auto* il = expr_cast<IntegerLit>(diff_simp)) verified = il->value.is_zero();
        
        if (!verified) {
            return make_unimplemented<ExprPtr>(
                "calculus", "solve_risch_de_field",
                "LHS != RHS in algebraic Risch DE solver",
                cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
                "Risch DE solver: check algebraic coefficients verification");
        }
        
        ExprPtr p_expr = p_sol_expr_res.value();
        ExprPtr y = arena.make<Binary>(BinaryOp::Div, p_expr, D);
        auto y_tog = algebra::together(y, ctx);
        if (y_tog.is_ok()) {
            if (auto s = ctx.simplify(y_tog.value()); s.is_ok()) return ok(s.value());
        }
        return ok(y);
    } else {
        // df <= 0. f_new is a constant in t (f_new = f_0 in the lower field).
        ExprPtr f_0 = leading_coefficient(f_poly);
        if (!f_0) f_0 = arena.make<IntegerLit>(BigInt(0));
        
        if (ext.type == ExtensionType::Logarithmic) {
            // Logarithmic case: recursive solving top-down
            auto theta_prime_res = field.derive(arena.make<Symbol>(t.name), ctx);
            if (theta_prime_res.is_error()) return fail<ExprPtr>(theta_prime_res.error());
            ExprPtr theta_prime = theta_prime_res.value();
            
            for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(N); i >= 0; --i) {
                const std::size_t iz = static_cast<std::size_t>(i);
                ExprPtr g_i = g_poly[iz];
                if (!g_i) g_i = arena.make<IntegerLit>(BigInt(0));
                
                ExprPtr rhs = g_i;
                if (iz + 1U <= static_cast<std::size_t>(N)) {
                    ExprPtr coef = arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(iz + 1U)));
                    ExprPtr term = arena.make<Product>(std::vector<ExprPtr>{coef, theta_prime, p_coeffs[iz + 1U]});
                    rhs = arena.make<Binary>(BinaryOp::Sub, rhs, term);
                }
                if (auto s = ctx.simplify(rhs); s.is_ok()) rhs = s.value();
                
                auto p_i_res = solve_risch_de_field(f_0, rhs, ext_idx - 1U, field, ctx);
                if (p_i_res.is_error()) return p_i_res;
                p_coeffs[iz] = p_i_res.value();
            }
        } else {
            // Exponential case: decoupled solving
            auto u_prime_res = field.derive(ext.argument, ctx);
            if (u_prime_res.is_error()) return fail<ExprPtr>(u_prime_res.error());
            ExprPtr u_prime = u_prime_res.value();
            
            for (std::size_t i = 0; i <= static_cast<std::size_t>(N); ++i) {
                ExprPtr g_i = g_poly[i];
                if (!g_i) g_i = arena.make<IntegerLit>(BigInt(0));
                
                ExprPtr f_eff = f_0;
                if (i > 0U) {
                    ExprPtr i_coef = arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(i)));
                    ExprPtr term = arena.make<Binary>(BinaryOp::Mul, i_coef, u_prime);
                    f_eff = arena.make<Binary>(BinaryOp::Add, term, f_0);
                    if (auto s = ctx.simplify(f_eff); s.is_ok()) f_eff = s.value();
                }
                
                auto p_i_res = solve_risch_de_field(f_eff, g_i, ext_idx - 1U, field, ctx);
                if (p_i_res.is_error()) return p_i_res;
                p_coeffs[i] = p_i_res.value();
            }
        }
        
        algebra::PolyExpr p_sol(p_coeffs);
        auto p_sol_expr_res = algebra::polynomial_to_expr(p_sol, t, ctx);
        if (p_sol_expr_res.is_error()) return fail<ExprPtr>(p_sol_expr_res.error());
        
        ExprPtr p_expr = p_sol_expr_res.value();
        ExprPtr y = arena.make<Binary>(BinaryOp::Div, p_expr, D);
        auto y_tog = algebra::together(y, ctx);
        if (y_tog.is_ok()) {
            if (auto s = ctx.simplify(y_tog.value()); s.is_ok()) return ok(s.value());
        }
        return ok(y);
    }
}

} // namespace

Result<ExprPtr> solve_risch_de_general(
    ExprPtr f,
    ExprPtr g,
    const Symbol& var,
    const DifferentialField& field,
    symbolic::CASContext& ctx) {
    (void)var;
    return solve_risch_de_field(f, g, field.extensions().size(), field, ctx);
}

// solve_risch_de_parametric_field lives in risch_rde_parametric_field.cpp
// (A29 anti-monolith split).

} // namespace cas::calculus
