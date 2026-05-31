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
        
        Symbol C_fresh = ctx.make_fresh_symbol("C");
        ExprPtr C1 = arena.make<Symbol>(C_fresh.name);
        ExprPtr rhs = arena.make<Binary>(BinaryOp::Add, int_Mx_res.value(), C1);
        
        // Restituiamo una relazione implicita tramite FuncCall "equal"
        return ok(arena.make<FuncCall>("equal", std::vector<ExprPtr>{int_Ny_res.value(), rhs}));
    }
    
    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Metodo non implementato per questo tipo di ODE del primo ordine."));
}

} // namespace cas::calculus
