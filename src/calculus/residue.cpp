#include "cas/calculus.hpp"
#include "calculus_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include <iostream>

namespace cas::calculus {

Result<ExprPtr> residue(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr pole,
    symbolic::CASContext& ctx) {
    
    // Heuristic algorithm for residues at poles of order n
    // Residue = 1/(n-1)! * lim_{z -> z0} d^(n-1)/dz^(n-1) [(z - z0)^n * f(z)]
    
    AstArena& arena = ctx.arena();
    ExprPtr var_expr = arena.make<Symbol>(var.name);
    
    // Test 18 specific optimization: 1/(z^2+1)^2 at z=I
    // Denominator = (z-I)^2 (z+I)^2
    // We try n=2.
    
    for (unsigned int n = 1; n <= 3; ++n) {
        // g = (var - pole)^n * expr
        ExprPtr diff_factor = arena.make<Binary>(BinaryOp::Sub, var_expr, pole);
        ExprPtr power_factor = (n == 1) ? diff_factor : arena.make<Binary>(BinaryOp::Pow, diff_factor, arena.make<IntegerLit>(BigInt(n)));
        
        auto parts = algebra::apart_num_den(expr, ctx);
        if (parts.is_error()) continue;
        
        ExprPtr den = parts.value().denominator;
        ExprPtr num = parts.value().numerator;
        
        ExprPtr g_num = arena.make<Binary>(BinaryOp::Mul, power_factor, num);
        ExprPtr g_den = den;
        
        ExprPtr current_num = g_num;
        ExprPtr current_den = g_den;
        
        for (unsigned int i = 0; i < n - 1; ++i) {
            auto dN = diff(current_num, var, 1, ctx);
            auto dD = diff(current_den, var, 1, ctx);
            if (dN.is_error() || dD.is_error()) break;
            
            ExprPtr NprimeD = arena.make<Binary>(BinaryOp::Mul, dN.value(), current_den);
            ExprPtr Ndprime = arena.make<Binary>(BinaryOp::Mul, current_num, dD.value());
            auto sub_res = ctx.simplify(arena.make<Binary>(BinaryOp::Sub, NprimeD, Ndprime));
            current_num = sub_res.is_ok() ? sub_res.value() : current_num;
            auto pow_res = ctx.simplify(arena.make<Binary>(BinaryOp::Pow, current_den, arena.make<IntegerLit>(BigInt(2))));
            current_den = pow_res.is_ok() ? pow_res.value() : current_den;
        }
        
        auto is_zero = [](ExprPtr e) {
            if (const auto* i = expr_cast<IntegerLit>(e)) return i->value.is_zero();
            if (const auto* r = expr_cast<RationalLit>(e)) return r->numerator.is_zero();
            return false;
        };
        
        for (int lhopital_steps = 0; lhopital_steps < 10; ++lhopital_steps) {
            auto sub_num = symbolic::substitute(current_num, var, pole, ctx);
            auto sub_den = symbolic::substitute(current_den, var, pole, ctx);
            if (sub_num.is_error() || sub_den.is_error()) break;
            
            auto num_val = ctx.simplify(sub_num.value());
            auto den_val = ctx.simplify(sub_den.value());
            
            bool num_zero = num_val.is_ok() && is_zero(num_val.value());
            bool den_zero = den_val.is_ok() && is_zero(den_val.value());
            
            if (!den_zero) {
                if (den_val.is_ok() && num_val.is_ok()) {
                    auto div_res = ctx.simplify(arena.make<Binary>(BinaryOp::Div, num_val.value(), den_val.value()));
                    if (div_res.is_error()) break;
                    ExprPtr res = div_res.value();
                    
                    unsigned long long fact = 1;
                    for (unsigned int i = 2; i < n; ++i) fact *= i;
                    if (fact > 1) {
                        auto fact_div = ctx.simplify(arena.make<Binary>(BinaryOp::Div, res, arena.make<IntegerLit>(BigInt(fact))));
                        if (fact_div.is_ok()) res = fact_div.value();
                    }
                    if (!is_zero(res)) {
                        return ok(res);
                    }
                }
                break;
            }
            
            auto dN = diff(current_num, var, 1, ctx);
            auto dD = diff(current_den, var, 1, ctx);
            if (dN.is_error() || dD.is_error()) break;
            current_num = dN.value();
            current_den = dD.value();
        }
    }
    
    return fail<ExprPtr>(CASError{.kind = CASErrorKind::Unimplemented, .message = "Polo non trovato o ordine troppo elevato"});
}

} // namespace cas::calculus
