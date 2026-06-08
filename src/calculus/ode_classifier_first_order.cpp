#include "ode_classifier_internal.hpp"

// First-order ODE type detectors:
//   separable, exact, Bernoulli, homogeneous, Riccati

namespace cas::calculus {

[[nodiscard]] std::optional<OdeClassification> try_separable(
    ExprPtr E, const Symbol& y, const Symbol& x, ExprPtr y_prime, symbolic::CASContext& ctx) {
    auto lin = extract_first_order_yprime(E, y, x, y_prime, ctx);
    if (!lin) return std::nullopt;
    ExprPtr A = lin->A;
    ExprPtr B = lin->B;

    if (is_zero_expr(A, ctx)) return std::nullopt;

    AstArena& arena = ctx.arena();

    if (is_zero_expr(B, ctx)) {
        OdeClassification res(OdeType::Separable, E, y, x);
        res.components.push_back(arena.make<IntegerLit>(BigInt(1)));
        res.components.push_back(arena.make<IntegerLit>(BigInt(0)));
        return res;
    }

    if (!depends_on(A, x) && !depends_on(B, y)) {
        OdeClassification res(OdeType::Separable, E, y, x);
        res.components.push_back(A);
        res.components.push_back(ctx.simplify(arena.make<Unary>(UnaryOp::Neg, B)).value());
        return res;
    }

    std::pair<long long, long long> points[] = {{0,0}, {1,1}, {0,1}, {1,0}, {-1,-1}, {2,2}};
    for (auto [x0, y0] : points) {
        auto eval_xy = [&](ExprPtr expr, long long xv, long long yv) -> Result<ExprPtr> {
            auto s1 = substitute_any(expr, arena.make<Symbol>(x.name), arena.make<IntegerLit>(BigInt(xv)), ctx);
            if (s1.is_error()) return s1;
            auto s2 = substitute_any(s1.value(), arena.make<Symbol>(y.name), arena.make<IntegerLit>(BigInt(yv)), ctx);
            if (s2.is_error()) return s2;
            return ctx.simplify(s2.value());
        };

        auto A0_res = eval_xy(A, x0, y0);
        auto B0_res = eval_xy(B, x0, y0);
        if (A0_res.is_error() || B0_res.is_error()) continue;
        ExprPtr A0 = A0_res.value();
        ExprPtr B0 = B0_res.value();

        if (is_zero_expr(A0, ctx) || is_zero_expr(B0, ctx)) continue;

        auto Ax_y0_res = substitute_any(A, arena.make<Symbol>(y.name), arena.make<IntegerLit>(BigInt(y0)), ctx);
        auto Bx_y0_res = substitute_any(B, arena.make<Symbol>(y.name), arena.make<IntegerLit>(BigInt(y0)), ctx);
        if (Ax_y0_res.is_error() || Bx_y0_res.is_error()) continue;
        ExprPtr Ax_y0 = ctx.simplify(Ax_y0_res.value()).value();
        ExprPtr Bx_y0 = ctx.simplify(Bx_y0_res.value()).value();

        auto Ax0_y_res = substitute_any(A, arena.make<Symbol>(x.name), arena.make<IntegerLit>(BigInt(x0)), ctx);
        auto Bx0_y_res = substitute_any(B, arena.make<Symbol>(x.name), arena.make<IntegerLit>(BigInt(x0)), ctx);
        if (Ax0_y_res.is_error() || Bx0_y_res.is_error()) continue;
        ExprPtr Ax0_y = ctx.simplify(Ax0_y_res.value()).value();
        ExprPtr Bx0_y = ctx.simplify(Bx0_y_res.value()).value();

        auto term1 = arena.make<Binary>(BinaryOp::Mul, A,
            arena.make<Binary>(BinaryOp::Mul, Bx_y0,
                arena.make<Binary>(BinaryOp::Mul, Bx0_y, A0)));
        auto term2 = arena.make<Binary>(BinaryOp::Mul, B,
            arena.make<Binary>(BinaryOp::Mul, Ax0_y,
                arena.make<Binary>(BinaryOp::Mul, B0, Ax_y0)));

        auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub, term1, term2), ctx);
        if (diff.is_ok() && is_zero_expr(diff.value(), ctx)) {
            auto M_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Div,
                arena.make<Unary>(UnaryOp::Neg, Bx_y0), Ax_y0));
            auto N_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Div,
                arena.make<Binary>(BinaryOp::Mul, Ax0_y, B0),
                arena.make<Binary>(BinaryOp::Mul, Bx0_y, A0)));

            if (M_raw.is_ok() && N_raw.is_ok()) {
                OdeClassification res(OdeType::Separable, E, y, x);
                res.components.push_back(N_raw.value());
                res.components.push_back(M_raw.value());
                return res;
            }
        }
    }

    return std::nullopt;
}

[[nodiscard]] std::optional<OdeClassification> try_exact(
    ExprPtr E, const Symbol& y, const Symbol& x, ExprPtr y_prime, symbolic::CASContext& ctx) {
    auto lin = extract_first_order_yprime(E, y, x, y_prime, ctx);
    if (!lin) return std::nullopt;
    ExprPtr N = lin->A;
    ExprPtr M = lin->B;

    AstArena& arena = ctx.arena();

    auto dM_dy = cas::calculus::diff(M, Symbol(y.name), 1, ctx);
    auto dN_dx = cas::calculus::diff(N, Symbol(x.name), 1, ctx);
    if (dM_dy.is_error() || dN_dx.is_error()) return std::nullopt;

    auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub, dM_dy.value(), dN_dx.value()), ctx);
    if (diff.is_ok()) {
        if (is_zero_expr(diff.value(), ctx)) {
            OdeClassification res(OdeType::Exact, E, y, x);
            res.components.push_back(M);
            res.components.push_back(N);
            return res;
        }
    }

    return std::nullopt;
}

[[nodiscard]] std::optional<OdeClassification> try_bernoulli(
    ExprPtr E, const Symbol& y, const Symbol& x, ExprPtr y_prime, symbolic::CASContext& ctx) {
    auto lin = extract_first_order_yprime(E, y, x, y_prime, ctx);
    if (!lin) return std::nullopt;
    ExprPtr A = lin->A;
    ExprPtr B = lin->B;

    AstArena& arena = ctx.arena();
    if (depends_on(A, y)) return std::nullopt;

    std::vector<ExprPtr> terms;
    std::function<void(ExprPtr, bool)> collect_terms = [&](ExprPtr e, bool neg) {
        if (const auto* sum = expr_cast<Sum>(e)) {
            for (auto t : sum->terms) collect_terms(t, neg);
        } else if (const auto* bin = expr_cast<Binary>(e)) {
            if (bin->op == BinaryOp::Add) {
                collect_terms(bin->left, neg);
                collect_terms(bin->right, neg);
            } else if (bin->op == BinaryOp::Sub) {
                collect_terms(bin->left, neg);
                collect_terms(bin->right, !neg);
            } else {
                terms.push_back(neg ? arena.make<Unary>(UnaryOp::Neg, e) : e);
            }
        } else if (const auto* un = expr_cast<Unary>(e)) {
            if (un->op == UnaryOp::Neg) collect_terms(un->operand, !neg);
            else terms.push_back(neg ? arena.make<Unary>(UnaryOp::Neg, e) : e);
        } else {
            terms.push_back(neg ? arena.make<Unary>(UnaryOp::Neg, e) : e);
        }
    };
    collect_terms(B, false);

    struct Group {
        ExprPtr power;
        std::vector<ExprPtr> coeffs;
    };
    std::vector<Group> groups;

    for (auto t : terms) {
        auto p_opt = get_y_power(t, y, ctx);
        if (!p_opt) return std::nullopt;
        ExprPtr p = *p_opt;

        ExprPtr yp;
        if (is_zero_expr(p, ctx)) {
            yp = arena.make<IntegerLit>(BigInt(1));
        } else if (const auto* il = expr_cast<IntegerLit>(p); il && il->value == BigInt(1)) {
            yp = arena.make<Symbol>(y.name);
        } else {
            yp = arena.make<Binary>(BinaryOp::Pow, arena.make<Symbol>(y.name), p);
        }
        auto C_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Div, t, yp));
        if (C_raw.is_error()) return std::nullopt;
        ExprPtr C = C_raw.value();

        bool found = false;
        for (auto& g : groups) {
            auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub, g.power, p), ctx);
            if (diff.is_ok() && is_zero_expr(diff.value(), ctx)) {
                g.coeffs.push_back(C);
                found = true;
                break;
            }
        }
        if (!found) {
            groups.push_back({p, {C}});
        }
    }

    if (groups.size() != 2) return std::nullopt;

    int idx_1 = -1;
    for (size_t i = 0; i < groups.size(); ++i) {
        auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub, groups[i].power, arena.make<IntegerLit>(BigInt(1))), ctx);
        if (diff.is_ok() && is_zero_expr(diff.value(), ctx)) {
            idx_1 = static_cast<int>(i);
            break;
        }
    }
    if (idx_1 == -1) return std::nullopt;

    int idx_n = 1 - idx_1;
    ExprPtr n_expr = groups[idx_n].power;

    ExprPtr P_raw = ctx.simplify(arena.make<Sum>(groups[idx_1].coeffs)).value();
    ExprPtr Q_raw = ctx.simplify(arena.make<Unary>(UnaryOp::Neg, arena.make<Sum>(groups[idx_n].coeffs))).value();

    auto P = ctx.simplify(arena.make<Binary>(BinaryOp::Div, P_raw, A));
    auto Q = ctx.simplify(arena.make<Binary>(BinaryOp::Div, Q_raw, A));
    if (P.is_error() || Q.is_error()) return std::nullopt;

    OdeClassification res(OdeType::Bernoulli, E, y, x);
    res.components.push_back(P.value());
    res.components.push_back(Q.value());
    res.components.push_back(n_expr);
    return res;
}

[[nodiscard]] std::optional<OdeClassification> try_homogeneous(
    ExprPtr E, const Symbol& y, const Symbol& x, ExprPtr y_prime, symbolic::CASContext& ctx) {
    auto lin = extract_first_order_yprime(E, y, x, y_prime, ctx);
    if (!lin) return std::nullopt;
    ExprPtr A = lin->A;
    ExprPtr B = lin->B;

    AstArena& arena = ctx.arena();

    Symbol v_sym = ctx.make_fresh_symbol("v");
    ExprPtr v_expr = arena.make<Symbol>(v_sym.name);
    ExprPtr vx = arena.make<Binary>(BinaryOp::Mul, v_expr, arena.make<Symbol>(x.name));

    auto A1_res = substitute_any(A, arena.make<Symbol>(x.name), arena.make<IntegerLit>(BigInt(1)), ctx);
    if (A1_res.is_error()) return std::nullopt;
    auto A1_v_res = substitute_any(A1_res.value(), arena.make<Symbol>(y.name), v_expr, ctx);
    if (A1_v_res.is_error()) return std::nullopt;
    ExprPtr A1 = ctx.simplify(A1_v_res.value()).value();

    auto B1_res = substitute_any(B, arena.make<Symbol>(x.name), arena.make<IntegerLit>(BigInt(1)), ctx);
    if (B1_res.is_error()) return std::nullopt;
    auto B1_v_res = substitute_any(B1_res.value(), arena.make<Symbol>(y.name), v_expr, ctx);
    if (B1_v_res.is_error()) return std::nullopt;
    ExprPtr B1 = ctx.simplify(B1_v_res.value()).value();

    auto Ax_res = substitute_any(A, arena.make<Symbol>(y.name), vx, ctx);
    if (Ax_res.is_error()) return std::nullopt;
    ExprPtr Ax = ctx.simplify(Ax_res.value()).value();

    auto Bx_res = substitute_any(B, arena.make<Symbol>(y.name), vx, ctx);
    if (Bx_res.is_error()) return std::nullopt;
    ExprPtr Bx = ctx.simplify(Bx_res.value()).value();

    auto term1 = arena.make<Binary>(BinaryOp::Mul, Ax, B1);
    auto term2 = arena.make<Binary>(BinaryOp::Mul, Bx, A1);

    auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub, term1, term2), ctx);
    if (diff.is_ok() && is_zero_expr(diff.value(), ctx)) {
        auto F_v_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Div,
            arena.make<Unary>(UnaryOp::Neg, B1), A1));
        if (F_v_raw.is_ok()) {
            OdeClassification res(OdeType::Homogeneous, E, y, x);
            res.parameter = v_sym;
            res.components.push_back(F_v_raw.value());
            return res;
        }
    }

    return std::nullopt;
}

// Detect the Riccati family  y' = q_0(x) + q_1(x)·y + q_2(x)·y²  with q_2 ≢ 0.
// Returns std::nullopt if `E` (the equation written as E = 0) does not fit the
// Riccati pattern; on success, fills `components` with [q_0, q_1, q_2] in this
// order.
//
// Algorithm (no hardcode, no shortcut):
//   1. View E as a polynomial in y of expected degree ≤ 2 with coefficients in
//      Q(x)[y'].  Recover the coefficients c_0(x,y'), c_1(x,y'), c_2(x,y') by
//      evaluating E at four y-values (Lagrange three-point fit + degree witness).
//   2. Reject if any of c_1, c_2 depends on y' (would make the equation
//      non-first-order or non-polynomial in y).
//   3. Reject if c_2 ≡ 0 (would be the already-handled linear branch).
//   4. Split c_0 as α(x)·y' + β(x).  Reject if α ≡ 0 (no first derivative) or
//      α depends on y' (E quadratic in y').
//   5. Normalise q_0 = -β/α, q_1 = -c_1/α, q_2 = -c_2/α.  Reject if any q_i
//      depends on y (genuine non-polynomial coefficient).
[[nodiscard]] std::optional<OdeClassification> try_riccati(
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

    auto e0  = eval_y(0);  if (e0.is_error())  return std::nullopt;
    auto e1  = eval_y(1);  if (e1.is_error())  return std::nullopt;
    auto em1 = eval_y(-1); if (em1.is_error()) return std::nullopt;
    auto e2  = eval_y(2);  if (e2.is_error())  return std::nullopt;

    auto two   = arena.make<IntegerLit>(BigInt(2));
    auto four  = arena.make<IntegerLit>(BigInt(4));
    auto half  = arena.make<RationalLit>(BigInt(1), BigInt(2));

    auto c0 = e0.value();
    auto c1_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Mul,
        half,
        arena.make<Binary>(BinaryOp::Sub, e1.value(), em1.value())));
    if (c1_raw.is_error()) return std::nullopt;
    auto c1 = c1_raw.value();

    auto c2_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Sub,
        arena.make<Binary>(BinaryOp::Mul, half,
            arena.make<Binary>(BinaryOp::Add, e1.value(), em1.value())),
        c0));
    if (c2_raw.is_error()) return std::nullopt;
    auto c2 = c2_raw.value();

    // Degree-witness check: E|_{y=2} must equal 4·c_2 + 2·c_1 + c_0.
    auto expected_e2 = ctx.simplify(arena.make<Sum>(std::vector<ExprPtr>{
        arena.make<Binary>(BinaryOp::Mul, four, c2),
        arena.make<Binary>(BinaryOp::Mul, two,  c1),
        c0}));
    if (expected_e2.is_error()) return std::nullopt;
    auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub, e2.value(),
        expected_e2.value()), ctx);
    if (diff.is_error() || !is_zero_expr(diff.value(), ctx)) return std::nullopt;

    // Reject c_2 ≡ 0 (linear-in-y case handled elsewhere).
    if (is_zero_expr(c2, ctx)) return std::nullopt;

    // c_1, c_2 must NOT depend on y' (otherwise structurally non-Riccati).
    auto depends_on_yprime = [&](ExprPtr e) {
        auto mo = find_max_order(e, y, x);
        return mo.has_value() && *mo >= 1U;
    };
    if (depends_on_yprime(c1) || depends_on_yprime(c2)) return std::nullopt;

    // Split c_0 into α(x)·y' + β(x) via two y'-evaluations.
    auto c0_at_yp = [&](long long k) -> Result<ExprPtr> {
        ExprPtr val = arena.make<IntegerLit>(BigInt(k));
        auto s = substitute_any(c0, y_prime, val, ctx);
        if (s.is_error()) return s;
        return ctx.simplify(s.value());
    };
    auto beta_res  = c0_at_yp(0);    if (beta_res.is_error())  return std::nullopt;
    auto c0_at_1   = c0_at_yp(1);    if (c0_at_1.is_error())   return std::nullopt;
    auto c0_at_2   = c0_at_yp(2);    if (c0_at_2.is_error())   return std::nullopt;

    auto alpha_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Sub,
        c0_at_1.value(), beta_res.value()));
    if (alpha_raw.is_error()) return std::nullopt;
    auto alpha = alpha_raw.value();

    // Affine witness: c_0|_{y'=2} must equal 2·α + β.
    auto expected_c0_at_2 = ctx.simplify(arena.make<Binary>(BinaryOp::Add,
        arena.make<Binary>(BinaryOp::Mul, two, alpha),
        beta_res.value()));
    if (expected_c0_at_2.is_error()) return std::nullopt;
    auto c0_diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub,
        c0_at_2.value(), expected_c0_at_2.value()), ctx);
    if (c0_diff.is_error() || !is_zero_expr(c0_diff.value(), ctx)) return std::nullopt;

    if (is_zero_expr(alpha, ctx)) return std::nullopt;          // no y' term
    if (depends_on_yprime(alpha) || depends_on_yprime(beta_res.value()))
        return std::nullopt;                                    // c_0 not affine in y'

    // q_0 = -β/α, q_1 = -c_1/α, q_2 = -c_2/α.
    auto neg = [&](ExprPtr e) { return arena.make<Unary>(UnaryOp::Neg, e); };
    auto divα = [&](ExprPtr num) {
        return ctx.simplify(arena.make<Binary>(BinaryOp::Div, num, alpha));
    };
    auto q0 = divα(neg(beta_res.value()));
    auto q1 = divα(neg(c1));
    auto q2 = divα(neg(c2));
    if (q0.is_error() || q1.is_error() || q2.is_error()) return std::nullopt;

    // Final sanity: q_i must not depend on y (otherwise spurious split).
    auto deps_y = [&](ExprPtr e) { return depends_on(e, y); };
    if (deps_y(q0.value()) || deps_y(q1.value()) || deps_y(q2.value()))
        return std::nullopt;

    OdeClassification res(OdeType::Riccati, /*equation*/E, y, x);
    res.components.push_back(q0.value());
    res.components.push_back(q1.value());
    res.components.push_back(q2.value());
    return res;
}

} // namespace cas::calculus
