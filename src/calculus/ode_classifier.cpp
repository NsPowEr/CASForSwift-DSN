#include "cas/ode.hpp"
#include "calculus_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include "../algebra/polynomial_internal.hpp"
#include <functional>

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

[[nodiscard]] static Result<ExprPtr> substitute_any(ExprPtr expr, ExprPtr target, ExprPtr replacement, symbolic::CASContext& ctx) {
    if (structural_equal(expr, target)) return ok(replacement);
    
    return visit_expr(expr, [&](const auto& node) -> Result<ExprPtr> {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, Unary>) {
            auto op = substitute_any(node.operand, target, replacement, ctx);
            if (op.is_error()) return op;
            return ok(ctx.arena().make<Unary>(node.op, op.value()));
        } else if constexpr (std::is_same_v<T, Binary>) {
            auto l = substitute_any(node.left, target, replacement, ctx);
            if (l.is_error()) return l;
            auto r = substitute_any(node.right, target, replacement, ctx);
            if (r.is_error()) return r;
            return ok(ctx.arena().make<Binary>(node.op, l.value(), r.value()));
        } else if constexpr (std::is_same_v<T, Sum>) {
            std::vector<ExprPtr> terms;
            for (auto t : node.terms) {
                auto r = substitute_any(t, target, replacement, ctx);
                if (r.is_error()) return r;
                terms.push_back(r.value());
            }
            return ok(ctx.arena().make<Sum>(std::move(terms)));
        } else if constexpr (std::is_same_v<T, Product>) {
            std::vector<ExprPtr> factors;
            for (auto f : node.factors) {
                auto r = substitute_any(f, target, replacement, ctx);
                if (r.is_error()) return r;
                factors.push_back(r.value());
            }
            return ok(ctx.arena().make<Product>(std::move(factors)));
        } else if constexpr (std::is_same_v<T, FuncCall>) {
            std::vector<ExprPtr> args;
            for (auto a : node.args) {
                auto r = substitute_any(a, target, replacement, ctx);
                if (r.is_error()) return r;
                args.push_back(r.value());
            }
            return ok(ctx.arena().make<FuncCall>(node.name, std::move(args)));
        } else if constexpr (std::is_same_v<T, Derivative>) {
            auto e = substitute_any(node.expression, target, replacement, ctx);
            if (e.is_error()) return e;
            return ok(ctx.arena().make<Derivative>(e.value(), Symbol(node.variable.name), node.order));
        } else {
            return ok(expr);
        }
    });
}

[[nodiscard]] static std::optional<uint32_t> find_max_order(ExprPtr expr, const Symbol& y, const Symbol& x) {
    uint32_t max_order = 0;
    bool found = false;

    std::function<void(ExprPtr)> visit = [&](ExprPtr e) {
        if (const auto* deriv = expr_cast<Derivative>(e)) {
            if (const auto* sym = expr_cast<Symbol>(deriv->expression)) {
                if (sym->name == y.name && deriv->variable.name == x.name) {
                    if (deriv->order > max_order) max_order = deriv->order;
                    found = true;
                }
            }
        } else if (const auto* sym = expr_cast<Symbol>(e)) {
            if (sym->name == y.name) {
                found = true;
            }
        }
        
        if (const auto* un = expr_cast<Unary>(e)) visit(un->operand);
        else if (const auto* bin = expr_cast<Binary>(e)) { visit(bin->left); visit(bin->right); }
        else if (const auto* sum = expr_cast<Sum>(e)) { for (auto t : sum->terms) visit(t); }
        else if (const auto* prod = expr_cast<Product>(e)) { for (auto f : prod->factors) visit(f); }
        else if (const auto* func = expr_cast<FuncCall>(e)) { for (auto a : func->args) visit(a); }
        else if (const auto* intg = expr_cast<Integral>(e)) { visit(intg->integrand); }
        else if (const auto* der = expr_cast<Derivative>(e)) { visit(der->expression); }
    };
    visit(expr);
    
    if (found) return max_order;
    return std::nullopt;
}

[[nodiscard]] static Result<ExprPtr> substitute_y_derivatives(ExprPtr E, const std::vector<ExprPtr>& y_ders, const std::vector<ExprPtr>& vals, symbolic::CASContext& ctx) {
    ExprPtr res = E;
    // Sostituzione in ordine inverso: partiamo dalle derivate più alte
    for (int i = static_cast<int>(y_ders.size()) - 1; i >= 0; --i) {
        auto sub_res = substitute_any(res, y_ders[i], vals[i], ctx);
        if (sub_res.is_error()) return sub_res;
        res = sub_res.value();
    }
    return ctx.simplify(res);
}

[[nodiscard]] Result<OdeClassification> classify_ode(
    ExprPtr equation,
    const Symbol& y,
    const Symbol& x,
    symbolic::CASContext& ctx) {
    
    AstArena& arena = ctx.arena();
    ExprPtr eq_lhs = equation;
    if (const auto* bin = expr_cast<Binary>(equation)) {
        if (bin->op == BinaryOp::Equal) {
            auto sub = ctx.simplify(arena.make<Binary>(BinaryOp::Sub, bin->left, bin->right));
            if (sub.is_error()) return fail<OdeClassification>(sub.error());
            eq_lhs = sub.value();
        }
    }

    auto E_res = algebra::expand(eq_lhs, ctx);
    if (E_res.is_error()) return fail<OdeClassification>(E_res.error());
    ExprPtr E = E_res.value();

    auto max_order_opt = find_max_order(E, y, x);
    if (!max_order_opt) return ok(OdeClassification(OdeType::Unknown, equation, y, x));
    
    uint32_t n = *max_order_opt;
    
    std::vector<ExprPtr> y_ders;
    y_ders.push_back(arena.make<Symbol>(y.name));
    for (uint32_t i = 1; i <= n; ++i) {
        y_ders.push_back(arena.make<Derivative>(y_ders[0], Symbol(x.name), i));
    }
    
    auto get_val = [&](const std::vector<long long>& vals) -> Result<ExprPtr> {
        std::vector<ExprPtr> sub_vals;
        for (long long v : vals) {
            sub_vals.push_back(algebra::poly_make_integer(arena, v));
        }
        return substitute_y_derivatives(E, y_ders, sub_vals, ctx);
    };
    
    std::vector<long long> zero_vals(n + 1, 0);
    auto v0_res = get_val(zero_vals);
    if (v0_res.is_error()) return fail<OdeClassification>(v0_res.error());
    ExprPtr v0 = v0_res.value();
    
    std::vector<ExprPtr> coeffs(n + 1);
    for (uint32_t i = 0; i <= n; ++i) {
        std::vector<long long> one_vals(n + 1, 0);
        one_vals[i] = 1;
        auto vi_res = get_val(one_vals);
        if (vi_res.is_error()) return fail<OdeClassification>(vi_res.error());
        
        auto a_i = ctx.simplify(arena.make<Binary>(BinaryOp::Sub, vi_res.value(), v0));
        if (a_i.is_error()) return fail<OdeClassification>(a_i.error());
        coeffs[i] = a_i.value();
    }
    
    std::vector<ExprPtr> L_terms;
    for (uint32_t i = 0; i <= n; ++i) {
        L_terms.push_back(arena.make<Binary>(BinaryOp::Mul, coeffs[i], y_ders[i]));
    }
    L_terms.push_back(v0);
    ExprPtr L = arena.make<Sum>(std::move(L_terms));
    
    auto diff_res = algebra::expand(arena.make<Binary>(BinaryOp::Sub, E, L), ctx);
    if (diff_res.is_error() || !is_zero_expr(diff_res.value(), ctx)) {
        return ok(OdeClassification(OdeType::Unknown, equation, y, x));
    }
    
    auto deps_y = [&](ExprPtr e) { return depends_on(e, y); };
    if (deps_y(v0)) return ok(OdeClassification(OdeType::Unknown, equation, y, x));
    for (const auto& a_i : coeffs) {
        if (deps_y(a_i)) return ok(OdeClassification(OdeType::Unknown, equation, y, x));
    }
    
    bool constant_coeffs = true;
    for (const auto& a_i : coeffs) {
        if (depends_on(a_i, x)) {
            constant_coeffs = false;
            break;
        }
    }
    
    if (n >= 1 && !is_zero_expr(coeffs[n], ctx)) {
        if (constant_coeffs) {
            OdeType t = OdeType::LinearNthOrderConstantCoeff;
            if (n == 1) t = OdeType::Linear1stOrder;
            else if (n == 2) t = OdeType::Linear2ndOrderConstantCoeff;
            
            OdeClassification res(t, equation, y, x);
            // Salva come a_n, a_{n-1}, ..., a_0, e infine f(x) = -v0
            for (int i = static_cast<int>(n); i >= 0; --i) {
                res.components.push_back(coeffs[i]);
            }
            res.components.push_back(ctx.simplify(arena.make<Unary>(UnaryOp::Neg, v0)).value());
            return ok(res);
        } else {
             if (n == 2) {
                 OdeClassification res(OdeType::Linear2ndOrderRationalCoeff, equation, y, x);
                 for (int i = 2; i >= 0; --i) {
                     res.components.push_back(coeffs[i]);
                 }
                 res.components.push_back(ctx.simplify(arena.make<Unary>(UnaryOp::Neg, v0)).value());
                 return ok(res);
             }
        }
    }
    
    return ok(OdeClassification(OdeType::Unknown, equation, y, x));
}

[[nodiscard]] Result<ExprPtr> solve_ode(ExprPtr equation, const Symbol& y, const Symbol& x, symbolic::CASContext& ctx) {
    auto class_res = classify_ode(equation, y, x, ctx);
    if (class_res.is_error()) return fail<ExprPtr>(class_res.error());
    
    const auto& classification = class_res.value();
    
    switch (classification.type) {
        case OdeType::Separable:
        case OdeType::Linear1stOrder:
        case OdeType::Bernoulli:
        case OdeType::Exact:
            return solve_ode_1st_order(classification, ctx);
            
        case OdeType::Linear2ndOrderConstantCoeff:
        case OdeType::Linear2ndOrderRationalCoeff:
        case OdeType::LinearNthOrderConstantCoeff:
            return solve_ode_advanced(classification, ctx);
            
        default:
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Tipo di ODE non riconosciuto o non supportato analiticamente."));
    }
}

} // namespace cas::calculus
