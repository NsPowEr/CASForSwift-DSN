#include "cas/ode.hpp"
#include "calculus_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include "../algebra/polynomial_internal.hpp"

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

    ExprPtr y_ptr = arena.make<Symbol>(y.name);
    ExprPtr y_p = arena.make<Derivative>(y_ptr, Symbol(x.name), 1);
    ExprPtr y_pp = arena.make<Derivative>(y_ptr, Symbol(x.name), 2);

    auto get_val = [&](long long v2, long long v1, long long v0) -> Result<ExprPtr> {
        auto res = substitute_any(E, y_pp, algebra::poly_make_integer(arena, v2), ctx);
        if (res.is_error()) return res;
        res = substitute_any(res.value(), y_p, algebra::poly_make_integer(arena, v1), ctx);
        if (res.is_error()) return res;
        res = substitute_any(res.value(), y_ptr, algebra::poly_make_integer(arena, v0), ctx);
        if (res.is_error()) return res;
        return ctx.simplify(res.value());
    };

    auto v000 = get_val(0, 0, 0); if (v000.is_error()) return fail<OdeClassification>(v000.error());
    auto v100 = get_val(1, 0, 0); if (v100.is_error()) return fail<OdeClassification>(v100.error());
    auto v010 = get_val(0, 1, 0); if (v010.is_error()) return fail<OdeClassification>(v010.error());
    auto v001 = get_val(0, 0, 1); if (v001.is_error()) return fail<OdeClassification>(v001.error());

    auto a2 = ctx.simplify(arena.make<Binary>(BinaryOp::Sub, v100.value(), v000.value()));
    auto a1 = ctx.simplify(arena.make<Binary>(BinaryOp::Sub, v010.value(), v000.value()));
    auto a0 = ctx.simplify(arena.make<Binary>(BinaryOp::Sub, v001.value(), v000.value()));
    
    if (a2.is_error() || a1.is_error() || a0.is_error()) return fail<OdeClassification>(make_error(CASErrorKind::InternalError, "Extraction error"));

    ExprPtr L = arena.make<Sum>(std::vector<ExprPtr>{
        arena.make<Binary>(BinaryOp::Mul, a2.value(), y_pp),
        arena.make<Binary>(BinaryOp::Mul, a1.value(), y_p),
        arena.make<Binary>(BinaryOp::Mul, a0.value(), y_ptr),
        v000.value()
    });

    auto diff_res = algebra::expand(arena.make<Binary>(BinaryOp::Sub, E, L), ctx);
    if (diff_res.is_error() || !is_zero_expr(diff_res.value(), ctx)) {
        return ok(OdeClassification(OdeType::Unknown, equation, y, x));
    }

    auto deps_y = [&](ExprPtr e) { return depends_on(e, y); };
    if (deps_y(a2.value()) || deps_y(a1.value()) || deps_y(a0.value()) || deps_y(v000.value())) {
        return ok(OdeClassification(OdeType::Unknown, equation, y, x));
    }

    if (depends_on(a2.value(), x) || depends_on(a1.value(), x) || depends_on(a0.value(), x)) {
        OdeClassification res(OdeType::Linear2ndOrderRationalCoeff, equation, y, x);
        res.components = {a2.value(), a1.value(), a0.value(), ctx.simplify(arena.make<Unary>(UnaryOp::Neg, v000.value())).value()};
        return ok(res);
    }

    if (!is_zero_expr(a2.value(), ctx)) {
        OdeClassification res(OdeType::Linear2ndOrderConstantCoeff, equation, y, x);
        res.components = {a2.value(), a1.value(), a0.value(), ctx.simplify(arena.make<Unary>(UnaryOp::Neg, v000.value())).value()};
        return ok(res);
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
            return solve_ode_advanced(classification, ctx);
            
        default:
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Tipo di ODE non riconosciuto o non supportato analiticamente."));
    }
}

} // namespace cas::calculus
