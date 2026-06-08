#include "ode_classifier_internal.hpp"

// Public ODE classification entry point + solve_ode dispatcher.
// Detector implementations live in:
//   ode_classifier_first_order.cpp   (separable, exact, Bernoulli, homogeneous, Riccati)
//   ode_classifier_higher_order.cpp  (Clairaut, d'Alembert / Lagrange)

namespace cas::calculus {

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
        // The linear template did not match.  For order-1 ODEs we still have
        // a chance to recognise nonlinear families (F5.3): Riccati first,
        // then Lagrange (Clairaut / d'Alembert).
        if (n == 1U) {
            auto separable = try_separable(E, y, x, y_ders[1], ctx);
            if (separable.has_value()) return ok(std::move(*separable));

            auto exact = try_exact(E, y, x, y_ders[1], ctx);
            if (exact.has_value()) return ok(std::move(*exact));

            auto bernoulli = try_bernoulli(E, y, x, y_ders[1], ctx);
            if (bernoulli.has_value()) return ok(std::move(*bernoulli));

            auto homogeneous = try_homogeneous(E, y, x, y_ders[1], ctx);
            if (homogeneous.has_value()) return ok(std::move(*homogeneous));

            auto riccati = try_riccati(E, y, x, y_ders[0], y_ders[1], ctx);
            if (riccati.has_value()) return ok(std::move(*riccati));
            auto lagrange = try_lagrange(E, y, x, y_ders[0], y_ders[1], ctx);
            if (lagrange.has_value()) return ok(std::move(*lagrange));
        }
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

[[nodiscard]] Result<ExprPtr> solve_ode(
    ExprPtr equation, const Symbol& y, const Symbol& x, symbolic::CASContext& ctx) {
    auto class_res = classify_ode(equation, y, x, ctx);
    if (class_res.is_error()) return fail<ExprPtr>(class_res.error());

    const auto& classification = class_res.value();

    switch (classification.type) {
        case OdeType::Separable:
        case OdeType::Linear1stOrder:
        case OdeType::Bernoulli:
        case OdeType::Exact:
            return solve_ode_1st_order(classification, ctx);

        case OdeType::Riccati:
        case OdeType::Clairaut:
        case OdeType::DAlembert:
            return solve_ode_nonlinear(classification, ctx);

        case OdeType::Linear2ndOrderConstantCoeff:
        case OdeType::Linear2ndOrderRationalCoeff:
        case OdeType::LinearNthOrderConstantCoeff:
            return solve_ode_advanced(classification, ctx);

        default:
            return fail<ExprPtr>(ode_make_error(
                CASErrorKind::Unimplemented,
                "Tipo di ODE non riconosciuto o non supportato analiticamente."));
    }
}

} // namespace cas::calculus
