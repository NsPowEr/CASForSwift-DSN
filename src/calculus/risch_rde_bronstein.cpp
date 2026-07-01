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

Result<std::vector<ParametricRischDeQSolution>> solve_risch_de_parametric_field(
    ExprPtr f,
    const std::vector<ExprPtr>& g_vec,
    std::size_t ext_idx,
    const DifferentialField& field,
    symbolic::CASContext& ctx) {
    
    AstArena& arena = ctx.arena();
    
    // Base case: Q(x)
    if (ext_idx == 0U) {
        auto poly_res = solve_risch_de_parametric_q(f, g_vec, field.base_var(), ctx);
        if (poly_res.is_ok()) return poly_res;
        // The Q[x] solver bailed (typically a rational f or g, which the
        // primitive tower descent produces).  Split by f:
        //   f == 0 → rational limited integration (integrate each g_i, split
        //            rational + log/arctan atoms);
        //   f != 0 → full parametric Risch DE over Q(x) via the P/D ansatz.
        // Both are sound by construction (back-substitution verified).
        // (A26 / HC-A26-PRIMITIVE-PARAMQ-RATIONAL.)
        ExprPtr f_s = f;
        if (auto s = ctx.simplify(f); s.is_ok()) f_s = s.value();
        bool f_is_zero = false;
        if (const auto* il = expr_cast<IntegerLit>(f_s)) f_is_zero = il->value.is_zero();
        if (const auto* rl = expr_cast<RationalLit>(f_s)) f_is_zero = rl->numerator.is_zero();
        if (f_is_zero)
            return solve_param_limited_integration_rational_q(g_vec, field.base_var(), ctx);
        auto rat = solve_param_risch_de_rational_q(f, g_vec, field.base_var(), ctx);
        if (rat.is_ok()) return rat;
        return poly_res;
    }
    
    // Topmost extension at this level
    const auto& ext = field.extensions()[ext_idx - 1U];
    const Symbol& t = ext.t_var;
    const std::size_t m = g_vec.size();
    
    // Denominator bound: LCM of denominators of f and all g_i
    auto f_parts_res = algebra::apart_num_den(f, ctx);
    if (f_parts_res.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(f_parts_res.error());
    ExprPtr B = f_parts_res.value().denominator;
    
    ExprPtr Q = arena.make<IntegerLit>(BigInt(1));
    for (ExprPtr g_expr : g_vec) {
        auto g_parts_res = algebra::apart_num_den(g_expr, ctx);
        if (g_parts_res.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(g_parts_res.error());
        ExprPtr den = g_parts_res.value().denominator;
        auto gcd_den = algebra::polynomial_gcd(Q, den, t, ctx);
        if (gcd_den.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(gcd_den.error());
        ExprPtr Q_den = arena.make<Binary>(BinaryOp::Mul, Q, den);
        Q = arena.make<Binary>(BinaryOp::Div, Q_den, gcd_den.value());
        if (auto s = ctx.simplify(Q); s.is_ok()) Q = s.value();
    }
    
    auto gcd_qb = algebra::polynomial_gcd(Q, B, t, ctx);
    if (gcd_qb.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(gcd_qb.error());
    ExprPtr QB = arena.make<Binary>(BinaryOp::Mul, Q, B);
    ExprPtr D = arena.make<Binary>(BinaryOp::Div, QB, gcd_qb.value());
    if (auto s = ctx.simplify(D); s.is_ok()) D = s.value();
    
    // f_new = f - D'/D,  g_new_i = g_i * D
    auto D_prime_res = field.derive(D, ctx);
    if (D_prime_res.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(D_prime_res.error());
    ExprPtr D_prime = D_prime_res.value();
    
    ExprPtr D_prime_over_D = arena.make<Binary>(BinaryOp::Div, D_prime, D);
    ExprPtr f_new = arena.make<Binary>(BinaryOp::Sub, f, D_prime_over_D);
    if (auto s = ctx.simplify(f_new); s.is_ok()) f_new = s.value();
    
    std::vector<ExprPtr> g_new_vec;
    g_new_vec.reserve(m);
    for (ExprPtr g_expr : g_vec) {
        ExprPtr g_new = arena.make<Binary>(BinaryOp::Mul, g_expr, D);
        if (auto s = ctx.simplify(g_new); s.is_ok()) g_new = s.value();
        g_new_vec.push_back(g_new);
    }
    
    // Parse f_new and all g_new_i as polynomials in t
    auto f_poly_res = algebra::parse_polynomial(f_new, t, ctx);
    if (f_poly_res.is_error()) {
        return make_unimplemented<std::vector<ParametricRischDeQSolution>>(
            "calculus", "solve_risch_de_parametric_field",
            "f_new is not polynomial",
            cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
            "Risch DE solver: make sure f_new is polynomial");
    }
    const auto& f_poly = f_poly_res.value();
    
    std::vector<algebra::PolyExpr> g_polys;
    g_polys.reserve(m);
    int dg_max = -1;
    for (ExprPtr g_new : g_new_vec) {
        auto g_poly_res = algebra::parse_polynomial(g_new, t, ctx);
        if (g_poly_res.is_error()) {
            return make_unimplemented<std::vector<ParametricRischDeQSolution>>(
                "calculus", "solve_risch_de_parametric_field",
                "g_new_i is not polynomial",
                cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
                "Risch DE solver: make sure g_new_i is polynomial");
        }
        g_polys.push_back(g_poly_res.value());
        int dg = static_cast<int>(poly_degree(g_poly_res.value()));
        if (is_zero_poly(g_poly_res.value())) dg = -1;
        dg_max = std::max(dg_max, dg);
    }
    
    int df = static_cast<int>(poly_degree(f_poly));
    if (is_zero_poly(f_poly)) df = -1;
    
    // Degree bound N
    int N = 0;
    if (ext.type == ExtensionType::Logarithmic) {
        if (df > 0) {
            N = std::max(0, dg_max - df);
        } else if (df == 0) {
            N = dg_max;
        } else {
            N = dg_max + 1;
        }
    } else { // Exponential
        if (df > 0) {
            N = std::max(0, dg_max - df);
        } else if (df == 0) {
            N = std::max(dg_max, 0);
        } else {
            N = dg_max + 1;
        }
    }
    if (N < 0) N = 0;
    
    // df > 0 — non-cancellation parametric PolyRischDE (A1).  For the log/exp
    // monomials here, deg_t(f_new) > 0 ⇒ deg(b) > max(0, δ(t)−1), i.e. the
    // "deg(b) is Large Enough" case: Bronstein Symbolic Integration I §7.1,
    // ParamPolyRischDENoCancel1.  Solve  D(q) + f_new·q = Σ c_i·g_new_i  in
    // K[t], then divide the polynomial solutions q by D (as the df≤0 tail does).
    // Sound: solve_param_poly_risch_de_nocancel1 verifies each candidate by
    // field back-substitution.  The residual constant system is solved for any
    // tower K ⊇ Q(x) via ConstantSystem (Bronstein §7.1, Lemma 7.1.2).
    if (df > 0) {
        auto nc = solve_param_poly_risch_de_nocancel1(f_new, g_new_vec, N, t, field, ctx);
        if (nc.is_error()) return nc;
        std::vector<ParametricRischDeQSolution> out_sols;
        out_sols.reserve(nc.value().size());
        for (auto& sol : nc.value()) {
            ExprPtr y = arena.make<Binary>(BinaryOp::Div, sol.y, D);
            auto y_tog = algebra::together(y, ctx);
            if (y_tog.is_ok()) {
                if (auto s = ctx.simplify(y_tog.value()); s.is_ok()) y = s.value();
            }
            out_sols.push_back({y, std::move(sol.c)});
        }
        return ok(std::move(out_sols));
    }

    // df <= 0. f_new is f_0 in lower field.
    ExprPtr f_0 = leading_coefficient(f_poly);
    if (!f_0) f_0 = arena.make<IntegerLit>(BigInt(0));
    
    // We start at i = N.  Bind a const& so the bounds-checked const operator[]
    // is selected (the non-const overload is unchecked, std-vector style):
    // N = dg_max+1 routinely exceeds deg(g_s), and OOB must read as the zero
    // coefficient (nullptr → 0), not a heap-buffer-overflow.
    std::vector<ExprPtr> H_vec(m);
    for (std::size_t s = 0; s < m; ++s) {
        const algebra::PolyExpr& gp = g_polys[s];
        ExprPtr coeff = gp[static_cast<std::size_t>(N)];
        H_vec[s] = coeff ? coeff : arena.make<IntegerLit>(BigInt(0));
    }
    
    auto solve_recursive = [&](auto& self, int i, const std::vector<ExprPtr>& H_current) -> Result<std::vector<ParametricRischDeQSolution>> {
        ExprPtr F_eff = f_0;
        if (ext.type == ExtensionType::Exponential && i > 0) {
            auto u_prime_res = field.derive(ext.argument, ctx);
            if (u_prime_res.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(u_prime_res.error());
            ExprPtr i_coef = arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(i)));
            ExprPtr term = arena.make<Binary>(BinaryOp::Mul, i_coef, u_prime_res.value());
            F_eff = arena.make<Binary>(BinaryOp::Add, term, f_0);
            if (auto s = ctx.simplify(F_eff); s.is_ok()) F_eff = s.value();
        }
        
        auto sols_res = solve_risch_de_parametric_field(F_eff, H_current, ext_idx - 1U, field, ctx);
        if (sols_res.is_error()) return sols_res;
        
        const auto& sols = sols_res.value();
        if (i == 0) {
            return ok(sols);
        }
        
        const std::size_t num_sols = sols.size();
        std::vector<ExprPtr> H_next(num_sols);
        
        ExprPtr theta_prime = nullptr;
        if (ext.type == ExtensionType::Logarithmic) {
            auto theta_prime_res = field.derive(arena.make<Symbol>(t.name), ctx);
            if (theta_prime_res.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(theta_prime_res.error());
            theta_prime = theta_prime_res.value();
        }
        
        for (std::size_t r = 0; r < num_sols; ++r) {
            const auto& sol = sols[r];
            std::vector<ExprPtr> terms;
            for (std::size_t s = 0; s < H_current.size(); ++s) {
                const algebra::PolyExpr& gp = g_polys[s];  // const& → bounds-checked operator[]
                ExprPtr g_val = gp[static_cast<std::size_t>(i - 1)];
                if (!g_val) g_val = arena.make<IntegerLit>(BigInt(0));
                
                if (sol.c[s].numerator().is_zero()) continue;
                ExprPtr c_val = arena.make<RationalLit>(sol.c[s].numerator(), sol.c[s].denominator());
                terms.push_back(arena.make<Binary>(BinaryOp::Mul, c_val, g_val));
            }
            ExprPtr sum_g = terms.empty() ? arena.make<IntegerLit>(BigInt(0)) :
                            (terms.size() == 1U ? terms[0] : arena.make<Sum>(std::move(terms)));
            
            if (ext.type == ExtensionType::Logarithmic) {
                // PRIMITIVE-CASE GAP (HC-A26-PRIMITIVE-PARAMQ-RATIONAL): the
                // correction i·y·θ' with θ' = D(t) = u'/u re-introduces a
                // denominator, so H_next can be rational in the lower field.
                // The base case solve_risch_de_parametric_q is polynomial-only
                // (Q[x]) and will return a diagnostic Unimplemented for such a
                // forcing.  Completing this needs ParamRischDE over Q(x)
                // (weak-normalizer + denominator bound, Bronstein §5.12/§6.5).
                ExprPtr i_coef = arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(i)));
                ExprPtr corr = arena.make<Product>(std::vector<ExprPtr>{i_coef, sol.y, theta_prime});
                ExprPtr val = arena.make<Binary>(BinaryOp::Sub, sum_g, corr);
                if (auto s = ctx.simplify(val); s.is_ok()) val = s.value();
                H_next[r] = val;
            } else {
                if (auto s = ctx.simplify(sum_g); s.is_ok()) sum_g = s.value();
                H_next[r] = sum_g;
            }
        }
        
        auto next_sols_res = self(self, i - 1, H_next);
        if (next_sols_res.is_error()) return next_sols_res;
        
        const auto& next_sols = next_sols_res.value();
        std::vector<ParametricRischDeQSolution> combined_sols;
        combined_sols.reserve(next_sols.size());
        
        for (const auto& next_sol : next_sols) {
            std::vector<ExprPtr> y_i_terms;
            for (std::size_t r = 0; r < num_sols; ++r) {
                const Rational& cr = next_sol.c[r];
                if (cr.numerator().is_zero()) continue;
                ExprPtr cr_e = arena.make<RationalLit>(cr.numerator(), cr.denominator());
                y_i_terms.push_back(arena.make<Binary>(BinaryOp::Mul, cr_e, sols[r].y));
            }
            ExprPtr y_i = y_i_terms.empty() ? arena.make<IntegerLit>(BigInt(0)) :
                          (y_i_terms.size() == 1U ? y_i_terms[0] : arena.make<Sum>(std::move(y_i_terms)));
            if (auto s = ctx.simplify(y_i); s.is_ok()) y_i = s.value();
            
            ExprPtr ti = (i == 1) ? static_cast<ExprPtr>(arena.make<Symbol>(t.name)) :
                         static_cast<ExprPtr>(arena.make<Binary>(BinaryOp::Pow, arena.make<Symbol>(t.name),
                             arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(i)))));
            ExprPtr term_ti = arena.make<Binary>(BinaryOp::Mul, y_i, ti);
            ExprPtr y_comb = arena.make<Binary>(BinaryOp::Add, next_sol.y, term_ti);
            if (auto s = ctx.simplify(y_comb); s.is_ok()) y_comb = s.value();
            
            std::vector<Rational> c_comb(m, Rational(BigInt(0)));
            for (std::size_t r = 0; r < num_sols; ++r) {
                const Rational& cr = next_sol.c[r];
                for (std::size_t s = 0; s < m; ++s) {
                    c_comb[s] = c_comb[s] + cr * sols[r].c[s];
                }
            }
            
            combined_sols.push_back({y_comb, std::move(c_comb)});
        }
        
        return ok(std::move(combined_sols));
    };
    
    auto sols_res = solve_recursive(solve_recursive, N, H_vec);
    if (sols_res.is_error()) return sols_res;
    
    std::vector<ParametricRischDeQSolution> out_sols;
    out_sols.reserve(sols_res.value().size());
    for (auto& sol : sols_res.value()) {
        ExprPtr y = arena.make<Binary>(BinaryOp::Div, sol.y, D);
        auto y_tog = algebra::together(y, ctx);
        if (y_tog.is_ok()) {
            if (auto s = ctx.simplify(y_tog.value()); s.is_ok()) y = s.value();
        }
        out_sols.push_back({y, std::move(sol.c)});
    }
    
    return ok(std::move(out_sols));
}

} // namespace cas::calculus
