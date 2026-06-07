#include "cas/ode.hpp"
#include "calculus_internal.hpp"
#include "cas/algebra.hpp"

namespace cas::calculus {

[[nodiscard]] static CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] Result<ExprPtr> solve_ode_1st_order(const OdeClassification& classification, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    if (classification.type == OdeType::Linear1stOrder) {
        // y' + P(x)y = Q(x)
        if (classification.components.size() < 2) {
            return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Componenti P(x) o Q(x) mancanti per ODE lineare."));
        }
        ExprPtr P = classification.components[0];
        ExprPtr Q = classification.components[1];
        
        // mu(x) = exp(integral(P(x) dx))
        auto int_P_res = integrate_risch(P, classification.x, ctx);
        if (int_P_res.is_error()) return fail<ExprPtr>(int_P_res.error());
        
        ExprPtr mu = arena.make<FuncCall>("exp", std::vector<ExprPtr>{int_P_res.value()});
        auto mu_simp = ctx.simplify(mu);
        if (mu_simp.is_error()) return fail<ExprPtr>(mu_simp.error());
        mu = mu_simp.value();
        
        // integral(mu(x) * Q(x) dx)
        ExprPtr muQ = arena.make<Binary>(BinaryOp::Mul, mu, Q);
        auto int_muQ_res = integrate_risch(muQ, classification.x, ctx);
        if (int_muQ_res.is_error()) return fail<ExprPtr>(int_muQ_res.error());
        
        // y = (integral(muQ) + C) / mu — C generata fresh per evitare collisione utente.
        Symbol C_fresh = ctx.make_fresh_symbol("C");
        ExprPtr C1 = arena.make<Symbol>(C_fresh.name);
        ExprPtr num = arena.make<Binary>(BinaryOp::Add, int_muQ_res.value(), C1);
        ExprPtr sol = arena.make<Binary>(BinaryOp::Div, num, mu);
        
        return ctx.simplify(sol);
    }
    
    if (classification.type == OdeType::Separable) {
        // N(y) dy = M(x) dx
        if (classification.components.size() < 2) {
            return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Componenti N(y) o M(x) mancanti per ODE separabile."));
        }
        ExprPtr Ny = classification.components[0];
        ExprPtr Mx = classification.components[1];
        
        auto int_Ny_res = integrate_risch(Ny, classification.y, ctx);
        auto int_Mx_res = integrate_risch(Mx, classification.x, ctx);
        
        if (int_Ny_res.is_error()) return fail<ExprPtr>(int_Ny_res.error());
        if (int_Mx_res.is_error()) return fail<ExprPtr>(int_Mx_res.error());

        auto int_Ny_simp = ctx.simplify(int_Ny_res.value());
        if (int_Ny_simp.is_error()) return fail<ExprPtr>(int_Ny_simp.error());
        
        auto int_Mx_simp = ctx.simplify(int_Mx_res.value());
        if (int_Mx_simp.is_error()) return fail<ExprPtr>(int_Mx_simp.error());
        
        Symbol C_fresh = ctx.make_fresh_symbol("C");
        ExprPtr C1 = arena.make<Symbol>(C_fresh.name);
        ExprPtr rhs = arena.make<Binary>(BinaryOp::Add, int_Mx_simp.value(), C1);
        
        return ok(arena.make<FuncCall>("equal", std::vector<ExprPtr>{int_Ny_simp.value(), rhs}));
    }

    if (classification.type == OdeType::Exact) {
        // M(x,y) dx + N(x,y) dy = 0
        if (classification.components.size() < 2) {
            return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Componenti M o N mancanti per ODE esatta."));
        }
        ExprPtr M = classification.components[0];
        ExprPtr N = classification.components[1];
        
        auto int_M_res = integrate_risch(M, classification.x, ctx);
        if (int_M_res.is_error()) return fail<ExprPtr>(int_M_res.error());
        auto int_M_simp = ctx.simplify(int_M_res.value());
        if (int_M_simp.is_error()) return int_M_simp;
        ExprPtr U1 = int_M_simp.value();
        
        auto dU1_dy_res = partial_diff(U1, classification.y, ctx);
        if (dU1_dy_res.is_error()) return dU1_dy_res;
        auto dU1_dy_simp = ctx.simplify(dU1_dy_res.value());
        if (dU1_dy_simp.is_error()) return dU1_dy_simp;
        ExprPtr dU1_dy = dU1_dy_simp.value();
        
        ExprPtr diff_N = arena.make<Binary>(BinaryOp::Sub, N, dU1_dy);
        auto diff_simp_res = ctx.simplify(diff_N);
        if (diff_simp_res.is_error()) return diff_simp_res;
        diff_N = diff_simp_res.value();
        
        auto int_N_res = integrate_risch(diff_N, classification.y, ctx);
        if (int_N_res.is_error()) return fail<ExprPtr>(int_N_res.error());
        auto int_N_simp = ctx.simplify(int_N_res.value());
        if (int_N_simp.is_error()) return int_N_simp;
        ExprPtr U2 = int_N_simp.value();
        
        ExprPtr U = arena.make<Binary>(BinaryOp::Add, U1, U2);
        auto U_simp = ctx.simplify(U);
        if (U_simp.is_error()) return U_simp;
        U = U_simp.value();
        
        Symbol C_fresh = ctx.make_fresh_symbol("C");
        ExprPtr C1 = arena.make<Symbol>(C_fresh.name);
        
        return ok(arena.make<FuncCall>("equal", std::vector<ExprPtr>{U, C1}));
    }

    if (classification.type == OdeType::Bernoulli) {
        // y' + P(x)y = Q(x)y^n
        if (classification.components.size() < 3) {
            return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Componenti mancanti per ODE di Bernoulli."));
        }
        ExprPtr P = classification.components[0];
        ExprPtr Q = classification.components[1];
        ExprPtr n_expr = classification.components[2];

        ExprPtr one = arena.make<IntegerLit>(BigInt(1));
        ExprPtr one_minus_n = arena.make<Binary>(BinaryOp::Sub, one, n_expr);
        auto one_minus_n_simp_res = ctx.simplify(one_minus_n);
        if (one_minus_n_simp_res.is_error()) return fail<ExprPtr>(one_minus_n_simp_res.error());
        ExprPtr one_minus_n_simp = one_minus_n_simp_res.value();

        ExprPtr P_new = arena.make<Binary>(BinaryOp::Mul, one_minus_n_simp, P);
        ExprPtr Q_new = arena.make<Binary>(BinaryOp::Mul, one_minus_n_simp, Q);

        auto int_P_res = integrate_risch(P_new, classification.x, ctx);
        if (int_P_res.is_error()) return fail<ExprPtr>(int_P_res.error());
        
        ExprPtr mu = arena.make<FuncCall>("exp", std::vector<ExprPtr>{int_P_res.value()});
        auto mu_simp = ctx.simplify(mu);
        if (mu_simp.is_error()) return fail<ExprPtr>(mu_simp.error());
        mu = mu_simp.value();
        
        ExprPtr muQ = arena.make<Binary>(BinaryOp::Mul, mu, Q_new);
        auto int_muQ_res = integrate_risch(muQ, classification.x, ctx);
        if (int_muQ_res.is_error()) return fail<ExprPtr>(int_muQ_res.error());
        
        Symbol C_fresh = ctx.make_fresh_symbol("C");
        ExprPtr C1 = arena.make<Symbol>(C_fresh.name);
        ExprPtr num = arena.make<Binary>(BinaryOp::Add, int_muQ_res.value(), C1);
        ExprPtr v_sol = arena.make<Binary>(BinaryOp::Div, num, mu);
        
        auto v_sol_simp = ctx.simplify(v_sol);
        if (v_sol_simp.is_error()) return fail<ExprPtr>(v_sol_simp.error());
        
        ExprPtr y_pow = arena.make<Binary>(BinaryOp::Pow, arena.make<Symbol>(classification.y.name), one_minus_n_simp);
        
        return ok(arena.make<FuncCall>("equal", std::vector<ExprPtr>{y_pow, v_sol_simp.value()}));
    }

    if (classification.type == OdeType::Homogeneous) {
        // y' = F(v), v = y/x
        if (classification.components.size() < 1) {
            return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Componente F(v) mancante per ODE omogenea."));
        }
        ExprPtr F_v = classification.components[0];
        
        if (!classification.parameter.has_value()) {
            return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Parametro v mancante per ODE omogenea."));
        }
        Symbol v_sym = classification.parameter.value();
        ExprPtr v = arena.make<Symbol>(v_sym.name);
        
        ExprPtr denom = arena.make<Binary>(BinaryOp::Sub, F_v, v);
        
        ExprPtr one = arena.make<IntegerLit>(BigInt(1));
        ExprPtr LHS_integrand = arena.make<Binary>(BinaryOp::Div, one, denom);
        
        auto LHS_integrand_simp = ctx.simplify(LHS_integrand);
        if (LHS_integrand_simp.is_error()) return fail<ExprPtr>(LHS_integrand_simp.error());
        
        auto int_LHS_res = integrate_risch(LHS_integrand_simp.value(), v_sym, ctx);
        if (int_LHS_res.is_error()) return fail<ExprPtr>(int_LHS_res.error());
        
        auto LHS_simp = ctx.simplify(int_LHS_res.value());
        if (LHS_simp.is_error()) return fail<ExprPtr>(LHS_simp.error());
        ExprPtr LHS = LHS_simp.value();
        
        ExprPtr RHS_integrand = arena.make<Binary>(BinaryOp::Div, one, arena.make<Symbol>(classification.x.name));
        auto int_RHS_res = integrate_risch(RHS_integrand, classification.x, ctx);
        if (int_RHS_res.is_error()) return fail<ExprPtr>(int_RHS_res.error());
        
        auto RHS_simp = ctx.simplify(int_RHS_res.value());
        if (RHS_simp.is_error()) return fail<ExprPtr>(RHS_simp.error());
        ExprPtr RHS_int = RHS_simp.value();
        
        Symbol C_fresh = ctx.make_fresh_symbol("C");
        ExprPtr C1 = arena.make<Symbol>(C_fresh.name);
        ExprPtr RHS = arena.make<Binary>(BinaryOp::Add, RHS_int, C1);
        
        ExprPtr y_over_x = arena.make<Binary>(BinaryOp::Div, arena.make<Symbol>(classification.y.name), arena.make<Symbol>(classification.x.name));
        auto LHS_subst_res = ctx.substitute(LHS, v_sym, y_over_x);
        if (LHS_subst_res.is_error()) return fail<ExprPtr>(LHS_subst_res.error());
        
        return ok(arena.make<FuncCall>("equal", std::vector<ExprPtr>{LHS_subst_res.value(), RHS}));
    }
    
    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Metodo non implementato per questo tipo di ODE del primo ordine."));
}

} // namespace cas::calculus
