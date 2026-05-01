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

[[nodiscard]] static Result<ExprPtr> solve_linear_2nd_order_constant_coeffs(
    const OdeClassification& ode,
    symbolic::CASContext& ctx) {
    
    AstArena& arena = ctx.arena();
    // components = {a2, a1, a0, f}
    ExprPtr a2 = ode.components[0];
    ExprPtr a1 = ode.components[1];
    ExprPtr a0 = ode.components[2];
    ExprPtr f = ode.components[3];

    // Characteristic equation: a2*r^2 + a1*r + a0 = 0
    Symbol r_sym("r");
    ExprPtr r_ptr = arena.make<Symbol>("r");
    ExprPtr char_poly = arena.make<Sum>(std::vector<ExprPtr>{
        arena.make<Binary>(BinaryOp::Mul, a2, arena.make<Binary>(BinaryOp::Pow, r_ptr, arena.make<IntegerLit>(2))),
        arena.make<Binary>(BinaryOp::Mul, a1, r_ptr),
        a0
    });

    auto roots_res = algebra::solve_polynomial(char_poly, r_sym, ctx);
    if (roots_res.is_error()) return fail<ExprPtr>(roots_res.error());
    auto roots = roots_res.value();

    if (roots.size() != 2) return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Expected 2 roots for characteristic equation"));

    ExprPtr r1 = roots[0];
    ExprPtr r2 = roots[1];
    
    ExprPtr y1, y2;
    auto eq_res = symbolic::mathematically_equal(r1, r2, ctx);
    if (eq_res.is_ok() && eq_res.value()) {
        // Repeated root: (C1 + C2*x) * e^(r1*x)
        y1 = make_exp_rx(r1, ode.x, arena);
        y2 = arena.make<Binary>(BinaryOp::Mul, arena.make<Symbol>(ode.x.name), y1);
    } else {
        // Distinct roots: C1*e^(r1*x) + C2*e^(r2*x)
        y1 = make_exp_rx(r1, ode.x, arena);
        y2 = make_exp_rx(r2, ode.x, arena);
    }

    // Particular solution if f != 0
    ExprPtr yp = nullptr;
    if (!is_zero_expr(f, ctx)) {
        // Variation of Parameters
        // W = y1*y2' - y2*y1'
        auto y1p_res = diff(y1, ode.x, 1, ctx);
        auto y2p_res = diff(y2, ode.x, 1, ctx);
        if (y1p_res.is_error() || y2p_res.is_error()) return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Diff error in VOP"));
        
        ExprPtr W = arena.make<Binary>(BinaryOp::Sub,
            arena.make<Binary>(BinaryOp::Mul, y1, y2p_res.value()),
            arena.make<Binary>(BinaryOp::Mul, y2, y1p_res.value()));
        
        // u1 = -integral(y2*f / (a2*W))
        // u2 = integral(y1*f / (a2*W))
        ExprPtr a2W = arena.make<Binary>(BinaryOp::Mul, a2, W);
        ExprPtr integrand1 = ctx.simplify(arena.make<Unary>(UnaryOp::Neg, arena.make<Binary>(BinaryOp::Div, arena.make<Binary>(BinaryOp::Mul, y2, f), a2W))).value();
        ExprPtr integrand2 = ctx.simplify(arena.make<Binary>(BinaryOp::Div, arena.make<Binary>(BinaryOp::Mul, y1, f), a2W)).value();
        
        auto u1_res = integrate(integrand1, ode.x, ctx);
        auto u2_res = integrate(integrand2, ode.x, ctx);
        
        if (u1_res.is_ok() && u2_res.is_ok()) {
            yp = arena.make<Sum>(std::vector<ExprPtr>{
                arena.make<Binary>(BinaryOp::Mul, u1_res.value(), y1),
                arena.make<Binary>(BinaryOp::Mul, u2_res.value(), y2)
            });
        } else {
            // If integration fails, use an integral form
            yp = arena.make<Sum>(std::vector<ExprPtr>{
                arena.make<Binary>(BinaryOp::Mul, arena.make<Integral>(integrand1, ode.x, std::nullopt, std::nullopt), y1),
                arena.make<Binary>(BinaryOp::Mul, arena.make<Integral>(integrand2, ode.x, std::nullopt, std::nullopt), y2)
            });
        }
    }

    ExprPtr C1 = arena.make<Symbol>("C1");
    ExprPtr C2 = arena.make<Symbol>("C2");
    
    std::vector<ExprPtr> total_terms;
    total_terms.push_back(arena.make<Binary>(BinaryOp::Mul, C1, y1));
    total_terms.push_back(arena.make<Binary>(BinaryOp::Mul, C2, y2));
    if (yp) total_terms.push_back(yp);
    
    auto sol = ctx.simplify(arena.make<Sum>(std::move(total_terms)));
    if (sol.is_error()) return fail<ExprPtr>(sol.error());
    
    return ok(arena.make<Binary>(BinaryOp::Equal, arena.make<Symbol>(ode.y.name), sol.value()));
}

[[nodiscard]] Result<ExprPtr> solve_ode_advanced(const OdeClassification& classification, symbolic::CASContext& ctx) {
    if (classification.type == OdeType::Linear2ndOrderConstantCoeff) {
        return solve_linear_2nd_order_constant_coeffs(classification, ctx);
    }
    
    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Risolutore avanzato non ancora implementato per questo tipo di ODE."));
}

} // namespace cas::calculus
