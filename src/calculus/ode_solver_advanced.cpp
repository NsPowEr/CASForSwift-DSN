#include "cas/ode.hpp"
#include "calculus_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include <iostream>

namespace cas::calculus {

[[nodiscard]] static CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] static bool is_zero_expr(ExprPtr expr, symbolic::CASContext& ctx) {
    auto res = ctx.simplify(expr);
    if (res.is_error()) return false;
    if (const auto* il = expr_cast<IntegerLit>(res.value())) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(res.value())) return rl->numerator.is_zero();
    return false;
}

[[nodiscard]] static ExprPtr make_exp_rx(ExprPtr r, const Symbol& x, AstArena& arena) {
    return arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{
        arena.make<Binary>(BinaryOp::Mul, r, arena.make<Symbol>(x.name))
    });
}

[[nodiscard]] static Result<std::vector<ExprPtr>> solve_vop_system(const std::vector<ExprPtr>& y_funcs, const Symbol& x, symbolic::CASContext& ctx) {
    size_t n = y_funcs.size();
    AstArena& arena = ctx.arena();
    
    // Matrice aumentata [Wronskiano | B]
    std::vector<std::vector<ExprPtr>> aug(n, std::vector<ExprPtr>(n + 1));
    for (size_t j = 0; j < n; ++j) {
        aug[0][j] = y_funcs[j];
        for (size_t i = 1; i < n; ++i) {
            auto d_res = diff(aug[i-1][j], x, 1, ctx);
            if (d_res.is_error()) return fail<std::vector<ExprPtr>>(d_res.error());
            aug[i][j] = d_res.value();
        }
    }
    
    for (size_t i = 0; i < n; ++i) {
        aug[i][n] = (i == n - 1) ? arena.make<IntegerLit>(1) : arena.make<IntegerLit>(0);
    }
    
    // Eliminazione di Gauss per il Wronskiano (simbolica)
    for (size_t i = 0; i < n; ++i) {
        size_t pivot = i;
        for (size_t j = i + 1; j < n; ++j) {
            if (!is_zero_expr(aug[j][i], ctx)) {
                pivot = j;
                break;
            }
        }
        std::swap(aug[i], aug[pivot]);
        
        ExprPtr pivot_val = aug[i][i];
        if (is_zero_expr(pivot_val, ctx)) continue;
        
        for (size_t j = i + 1; j <= n; ++j) {
            auto div_res = ctx.simplify(arena.make<Binary>(BinaryOp::Div, aug[i][j], pivot_val));
            if (div_res.is_error()) return fail<std::vector<ExprPtr>>(div_res.error());
            aug[i][j] = div_res.value();
        }
        aug[i][i] = arena.make<IntegerLit>(1);
        
        for (size_t k = 0; k < n; ++k) {
            if (k != i && !is_zero_expr(aug[k][i], ctx)) {
                ExprPtr factor = aug[k][i];
                for (size_t j = i + 1; j <= n; ++j) {
                    auto sub_res = ctx.simplify(arena.make<Binary>(BinaryOp::Sub, aug[k][j],
                        arena.make<Binary>(BinaryOp::Mul, factor, aug[i][j])));
                    if (sub_res.is_error()) return fail<std::vector<ExprPtr>>(sub_res.error());
                    aug[k][j] = sub_res.value();
                }
                aug[k][i] = arena.make<IntegerLit>(0);
            }
        }
    }
    
    std::vector<ExprPtr> u_prime(n);
    for (size_t i = 0; i < n; ++i) {
        u_prime[i] = aug[i][n];
    }
    return ok(u_prime);
}

[[nodiscard]] static Result<ExprPtr> solve_linear_nth_order_constant_coeffs(
    const OdeClassification& ode,
    symbolic::CASContext& ctx) {
    
    AstArena& arena = ctx.arena();
    if (ode.components.size() < 2) return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Invalid ODE components"));
    size_t n = ode.components.size() - 2;
    ExprPtr an = ode.components[0];
    ExprPtr f = ode.components.back();

    // Equazione caratteristica: sum a_i r^i = 0
    Symbol r_sym("r");
    ExprPtr r_ptr = arena.make<Symbol>("r");
    std::vector<ExprPtr> char_terms;
    for (size_t i = 0; i <= n; ++i) {
        size_t power = n - i;
        ExprPtr term = ode.components[i];
        if (power > 0) {
            ExprPtr r_pow = (power == 1) ? r_ptr : arena.make<Binary>(BinaryOp::Pow, r_ptr, arena.make<IntegerLit>(power));
            term = arena.make<Binary>(BinaryOp::Mul, term, r_pow);
        }
        char_terms.push_back(term);
    }
    ExprPtr char_poly = arena.make<Sum>(std::move(char_terms));

    auto roots_res = algebra::solve_polynomial(char_poly, r_sym, ctx);
    if (roots_res.is_error()) return fail<ExprPtr>(roots_res.error());
    auto roots = roots_res.value();

    if (roots.size() != n) return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Expected " + std::to_string(n) + " roots for characteristic equation"));

    // Raggruppamento radici per gestire le molteplicità
    std::vector<std::pair<ExprPtr, size_t>> grouped_roots;
    for (auto r : roots) {
        bool found = false;
        for (auto& gr : grouped_roots) {
            auto eq_res = symbolic::mathematically_equal(r, gr.first, ctx);
            if (eq_res.is_ok() && eq_res.value()) {
                gr.second++;
                found = true;
                break;
            }
        }
        if (!found) {
            grouped_roots.push_back({r, 1});
        }
    }

    // Costruzione delle soluzioni fondamentali omogenee
    std::vector<ExprPtr> y_funcs;
    for (const auto& gr : grouped_roots) {
        ExprPtr root = gr.first;
        size_t mult = gr.second;
        ExprPtr exp_part = make_exp_rx(root, ode.x, arena);
        for (size_t i = 0; i < mult; ++i) {
            if (i == 0) {
                y_funcs.push_back(exp_part);
            } else {
                ExprPtr x_pow = (i == 1) ? arena.make<Symbol>(ode.x.name) : arena.make<Binary>(BinaryOp::Pow, arena.make<Symbol>(ode.x.name), arena.make<IntegerLit>(i));
                y_funcs.push_back(arena.make<Binary>(BinaryOp::Mul, x_pow, exp_part));
            }
        }
    }

    // Soluzione particolare f != 0 tramite Variazione delle Costanti N-esima (Wronskiano)
    ExprPtr yp = nullptr;
    if (!is_zero_expr(f, ctx)) {
        if (n == 2) {
            // Fast-path per ODE di 2° ordine con calcolo chiuso del Wronskiano
            ExprPtr y1 = y_funcs[0];
            ExprPtr y2 = y_funcs[1];
            
            ExprPtr r1 = grouped_roots[0].first;
            ExprPtr W = nullptr;
            
            if (grouped_roots.size() == 2) {
                // Radici distinte
                ExprPtr r2 = grouped_roots[1].first;
                ExprPtr r2_minus_r1 = ctx.simplify(arena.make<Binary>(BinaryOp::Sub, r2, r1)).value();
                ExprPtr r1_plus_r2 = ctx.simplify(arena.make<Binary>(BinaryOp::Add, r1, r2)).value();
                W = arena.make<Binary>(BinaryOp::Mul, r2_minus_r1, make_exp_rx(r1_plus_r2, ode.x, arena));
            } else {
                // Radice ripetuta
                ExprPtr two_r1 = ctx.simplify(arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(2), r1)).value();
                W = make_exp_rx(two_r1, ode.x, arena);
            }
                
            ExprPtr a2W = ctx.simplify(arena.make<Binary>(BinaryOp::Mul, an, W)).value();
            
            auto unsimplified1 = arena.make<Unary>(UnaryOp::Neg, arena.make<Binary>(BinaryOp::Div, arena.make<Binary>(BinaryOp::Mul, y2, f), a2W));
            auto exp1 = algebra::expand(unsimplified1, ctx);
            ExprPtr integrand1 = ctx.simplify(exp1.is_ok() ? exp1.value() : unsimplified1).value();
            
            auto unsimplified2 = arena.make<Binary>(BinaryOp::Div, arena.make<Binary>(BinaryOp::Mul, y1, f), a2W);
            auto exp2 = algebra::expand(unsimplified2, ctx);
            ExprPtr integrand2 = ctx.simplify(exp2.is_ok() ? exp2.value() : unsimplified2).value();
            
            auto u1_res = integrate(integrand1, ode.x, ctx);
            auto u2_res = integrate(integrand2, ode.x, ctx);
            
            ExprPtr u1 = u1_res.is_ok() ? u1_res.value() : arena.make<Integral>(integrand1, ode.x, std::nullopt, std::nullopt);
            ExprPtr u2 = u2_res.is_ok() ? u2_res.value() : arena.make<Integral>(integrand2, ode.x, std::nullopt, std::nullopt);
            
            yp = arena.make<Sum>(std::vector<ExprPtr>{
                arena.make<Binary>(BinaryOp::Mul, u1, y1),
                arena.make<Binary>(BinaryOp::Mul, u2, y2)
            });
        } else {
            auto u_prime_res = solve_vop_system(y_funcs, ode.x, ctx);
            if (u_prime_res.is_error()) return fail<ExprPtr>(u_prime_res.error());
            
            auto u_prime = u_prime_res.value();
            std::vector<ExprPtr> yp_terms;
            for (size_t i = 0; i < n; ++i) {
                // u_i' = u_prime[i] * f / an
                auto integrand_unsimplified = arena.make<Binary>(BinaryOp::Mul,
                                    u_prime[i],
                                    arena.make<Binary>(BinaryOp::Div, f, an));
                auto integrand_res = ctx.simplify(integrand_unsimplified);
                if (integrand_res.is_error()) return fail<ExprPtr>(integrand_res.error());
                ExprPtr integrand = integrand_res.value();
                
                auto u_i_res = integrate(integrand, ode.x, ctx);

                ExprPtr u_i;
                if (u_i_res.is_ok()) {
                    u_i = u_i_res.value();
                } else {
                    u_i = arena.make<Integral>(integrand, ode.x, std::nullopt, std::nullopt);
                }
                yp_terms.push_back(arena.make<Binary>(BinaryOp::Mul, u_i, y_funcs[i]));
            }
            yp = arena.make<Sum>(std::move(yp_terms));
        }
    }

    std::vector<ExprPtr> total_terms;
    for (size_t i = 0; i < n; ++i) {
        ExprPtr Ci = arena.make<Symbol>("C" + std::to_string(i + 1));
        total_terms.push_back(arena.make<Binary>(BinaryOp::Mul, Ci, y_funcs[i]));
    }
    if (yp) total_terms.push_back(yp);
    
    auto sol = ctx.simplify(arena.make<Sum>(std::move(total_terms)));
    if (sol.is_error()) return fail<ExprPtr>(sol.error());
    
    return ok(arena.make<Binary>(BinaryOp::Equal, arena.make<Symbol>(ode.y.name), sol.value()));
}

[[nodiscard]] Result<ExprPtr> solve_ode_advanced(const OdeClassification& classification, symbolic::CASContext& ctx) {
    if (classification.type == OdeType::Linear2ndOrderConstantCoeff ||
        classification.type == OdeType::LinearNthOrderConstantCoeff) {
        return solve_linear_nth_order_constant_coeffs(classification, ctx);
    }
    
    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Risolutore avanzato non ancora implementato per questo tipo di ODE."));
}

} // namespace cas::calculus
