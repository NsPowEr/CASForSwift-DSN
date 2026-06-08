#include "ode_classifier_internal.hpp"

// Higher-order / special-form ODE type detectors:
//   Clairaut, d'Alembert (Lagrange family)

namespace cas::calculus {

// Detect Lagrange family:  E ≡ 0  ⇔  y = x·F(p) + G(p),  p ≡ y'.
// Sub-cases:
//   F(p) ≡ p  →  Clairaut, components = [G(p)], parameter = p.
//   F(p) ≢ p  →  d'Alembert, components = [F(p), G(p)], parameter = p.
//
// Algorithm (no hardcode, no shortcut):
//   1. View E as polynomial in y of degree ≤ 1: c_1·y + c_0(x, y').
//      c_1 = E|_{y=1} - E|_{y=0}.  Degree witness E|_{y=2} = 2c_1 + c_0.
//   2. Reject if c_1 depends on x, y, or y' (Lagrange has constant y-coeff).
//   3. Reject if c_1 ≡ 0 (no y term).
//   4. View c_0 as affine in x: c_0 = α(y')·x + β(y').
//      α = c_0|_{x=1} - c_0|_{x=0}.  Witness c_0|_{x=2} = 2α + β.
//   5. Reject if α or β depends on x or y.
//   6. F(p) = -α/c_1, G(p) = -β/c_1, with y' substituted by a fresh p.
//   7. Reject if F or G depends on x or y, or if F ≡ 0 (degenerate).
//   8. Classify: F ≡ p → Clairaut, else → d'Alembert.
[[nodiscard]] std::optional<OdeClassification> try_lagrange(
    ExprPtr E,
    const Symbol& y,
    const Symbol& x,
    ExprPtr y_symbol,
    ExprPtr y_prime,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    auto eval_y = [&](long long k) -> Result<ExprPtr> {
        ExprPtr val = arena.make<IntegerLit>(BigInt(k));
        auto s = substitute_y_shielded(E, y_symbol, val, arena);
        if (s.is_error()) return s;
        return ctx.simplify(s.value());
    };
    auto e0 = eval_y(0); if (e0.is_error()) return std::nullopt;
    auto e1 = eval_y(1); if (e1.is_error()) return std::nullopt;
    auto e2 = eval_y(2); if (e2.is_error()) return std::nullopt;

    ExprPtr c0 = e0.value();
    auto c1_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Sub,
        e1.value(), e0.value()));
    if (c1_raw.is_error()) return std::nullopt;
    ExprPtr c1 = c1_raw.value();

    auto two = arena.make<IntegerLit>(BigInt(2));
    auto expected_e2 = ctx.simplify(arena.make<Binary>(BinaryOp::Add,
        arena.make<Binary>(BinaryOp::Mul, two, c1), c0));
    if (expected_e2.is_error()) return std::nullopt;
    auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub,
        e2.value(), expected_e2.value()), ctx);
    if (diff.is_error() || !is_zero_expr(diff.value(), ctx)) return std::nullopt;

    if (is_zero_expr(c1, ctx)) return std::nullopt;

    auto depends_on_yprime = [&](ExprPtr e) {
        auto mo = find_max_order(e, y, x);
        return mo.has_value() && *mo >= 1U;
    };
    if (depends_on(c1, x) || depends_on(c1, y) || depends_on_yprime(c1))
        return std::nullopt;

    // c_0 = α(y')·x + β(y') affine-in-x extraction.
    auto c0_at_x = [&](long long k) -> Result<ExprPtr> {
        ExprPtr val = arena.make<IntegerLit>(BigInt(k));
        auto s = substitute_any(c0, arena.make<Symbol>(x.name), val, ctx);
        if (s.is_error()) return s;
        return ctx.simplify(s.value());
    };
    auto beta_raw = c0_at_x(0); if (beta_raw.is_error()) return std::nullopt;
    auto c0_at_1  = c0_at_x(1); if (c0_at_1.is_error())  return std::nullopt;
    auto c0_at_2  = c0_at_x(2); if (c0_at_2.is_error())  return std::nullopt;
    ExprPtr beta = beta_raw.value();

    auto alpha_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Sub,
        c0_at_1.value(), beta));
    if (alpha_raw.is_error()) return std::nullopt;
    ExprPtr alpha = alpha_raw.value();

    auto expected_c0_at_2 = ctx.simplify(arena.make<Binary>(BinaryOp::Add,
        arena.make<Binary>(BinaryOp::Mul, two, alpha), beta));
    if (expected_c0_at_2.is_error()) return std::nullopt;
    auto c0_diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub,
        c0_at_2.value(), expected_c0_at_2.value()), ctx);
    if (c0_diff.is_error() || !is_zero_expr(c0_diff.value(), ctx))
        return std::nullopt;

    if (depends_on(alpha, x) || depends_on(alpha, y)) return std::nullopt;
    if (depends_on(beta,  x) || depends_on(beta,  y)) return std::nullopt;

    // F(p) = -α/c_1, G(p) = -β/c_1 with y' replaced by fresh p.
    Symbol p_sym = ctx.make_fresh_symbol("p");
    ExprPtr p_expr = arena.make<Symbol>(p_sym.name);

    auto neg_div = [&](ExprPtr num) -> Result<ExprPtr> {
        ExprPtr q = arena.make<Binary>(BinaryOp::Div,
            arena.make<Unary>(UnaryOp::Neg, num), c1);
        return ctx.simplify(q);
    };
    auto F_in_yp = neg_div(alpha); if (F_in_yp.is_error()) return std::nullopt;
    auto G_in_yp = neg_div(beta);  if (G_in_yp.is_error()) return std::nullopt;

    auto F_in_p = substitute_any(F_in_yp.value(), y_prime, p_expr, ctx);
    if (F_in_p.is_error()) return std::nullopt;
    auto G_in_p = substitute_any(G_in_yp.value(), y_prime, p_expr, ctx);
    if (G_in_p.is_error()) return std::nullopt;
    auto F_p = ctx.simplify(F_in_p.value());
    auto G_p = ctx.simplify(G_in_p.value());
    if (F_p.is_error() || G_p.is_error()) return std::nullopt;

    if (depends_on(F_p.value(), x) || depends_on(F_p.value(), y))
        return std::nullopt;
    if (depends_on(G_p.value(), x) || depends_on(G_p.value(), y))
        return std::nullopt;
    if (is_zero_expr(F_p.value(), ctx)) return std::nullopt;  // degenerate

    auto f_minus_p = algebra::expand(arena.make<Binary>(BinaryOp::Sub,
        F_p.value(), p_expr), ctx);
    bool is_clairaut = f_minus_p.is_ok() && is_zero_expr(f_minus_p.value(), ctx);

    if (is_clairaut) {
        OdeClassification res(OdeType::Clairaut, /*equation*/E, y, x);
        res.parameter = p_sym;
        res.components.push_back(G_p.value());
        return res;
    }
    OdeClassification res(OdeType::DAlembert, /*equation*/E, y, x);
    res.parameter = p_sym;
    res.components.push_back(F_p.value());
    res.components.push_back(G_p.value());
    return res;
}

} // namespace cas::calculus
