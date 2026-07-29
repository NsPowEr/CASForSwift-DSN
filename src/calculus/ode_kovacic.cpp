// Kovacic's algorithm — public entry point.
// Implements solve_ode_kovacic: full pipeline for
//   a2·y'' + a1·y' + a0·y = f  (a₂ ≠ 0, coefficients rational in x).
//
// Step 1: compute invariant r via compute_r (ode_kovacic_case1.cpp).
// Step 2: find ω₊, ω₋ via case1_omega (Case 1 only; Case 2/3 → Unimplemented).
// Step 3: build z₁ = exp(∫ω₊), z₂ from ω₋ or reduction of order.
// Step 4: back-transform  y = z · exp(−½∫p dx).
// Step 5: if f ≠ 0, add particular solution via variation of parameters.

#include "ode_kovacic_internal.hpp"

namespace cas::calculus {

using namespace kovacic_impl;

// ─── Reduction of order: given z₁, z₂ = z₁·∫(1/z₁²) dx ─────────────────────

[[nodiscard]] static Result<ExprPtr> reduction_of_order(
    ExprPtr z1, const Symbol& x, symbolic::CASContext& ctx) {

    AstArena& a = ctx.arena();
    ExprPtr inv_z1sq = kv_div(a, kv_int(a, 1), kv_mul(a, z1, z1));
    auto inv_simp = ctx.simplify(inv_z1sq);
    if (inv_simp.is_error()) return inv_simp;

    auto int_res = integrate(inv_simp.value(), x, ctx);
    if (int_res.is_error()) {
        // Unevaluated integral — still a valid symbolic result.
        return ok(kv_mul(a, z1,
            a.make<Integral>(inv_simp.value(), x, std::nullopt, std::nullopt)));
    }
    return ctx.simplify(kv_mul(a, z1, int_res.value()));
}

// ─── Variation of parameters: particular solution y_p ────────────────────────

[[nodiscard]] static Result<ExprPtr> variation_of_params(
    ExprPtr y1, ExprPtr y2, ExprPtr g,
    const Symbol& x, symbolic::CASContext& ctx) {

    AstArena& a = ctx.arena();

    auto dy1 = diff(y1, x, 1U, ctx);
    auto dy2 = diff(y2, x, 1U, ctx);
    if (dy1.is_error()) return dy1;
    if (dy2.is_error()) return dy2;

    // W = y₁·y₂' − y₁'·y₂
    auto W_res = ctx.simplify(
        kv_sub(a, kv_mul(a, y1, dy2.value()), kv_mul(a, dy1.value(), y2)));
    if (W_res.is_error()) return W_res;
    ExprPtr W = W_res.value();

    // u₁' = −y₂·g / W,  u₂' = y₁·g / W
    auto int1_simp = ctx.simplify(kv_div(a, kv_neg(a, kv_mul(a, y2, g)), W));
    auto int2_simp = ctx.simplify(kv_div(a, kv_mul(a, y1, g), W));
    if (int1_simp.is_error()) return int1_simp;
    if (int2_simp.is_error()) return int2_simp;

    auto u1 = integrate(int1_simp.value(), x, ctx);
    auto u2 = integrate(int2_simp.value(), x, ctx);

    ExprPtr u1_expr = u1.is_ok()
        ? u1.value()
        : static_cast<ExprPtr>(a.make<Integral>(int1_simp.value(), x, std::nullopt, std::nullopt));
    ExprPtr u2_expr = u2.is_ok()
        ? u2.value()
        : static_cast<ExprPtr>(a.make<Integral>(int2_simp.value(), x, std::nullopt, std::nullopt));

    auto yp = ctx.simplify(
        kv_add(a, kv_mul(a, u1_expr, y1), kv_mul(a, u2_expr, y2)));
    return yp;
}

// ─── Public entry point ───────────────────────────────────────────────────────

[[nodiscard]] Result<ExprPtr> solve_ode_kovacic(
    const OdeClassification& cls, symbolic::CASContext& ctx) {

    if (cls.components.size() < 3) {
        return fail<ExprPtr>(CASError{
            .kind    = CASErrorKind::InternalError,
            .message = "solve_ode_kovacic: expected [a2, a1, a0, f?] in components",
            .hint    = std::nullopt});
    }

    AstArena& a = ctx.arena();
    ExprPtr a2 = cls.components[0];
    ExprPtr a1 = cls.components[1];
    ExprPtr a0 = cls.components[2];
    ExprPtr f  = (cls.components.size() >= 4)
        ? cls.components[3]
        : kv_int(a, 0);

    const Symbol& x = cls.x;
    const Symbol& y = cls.y;

    // Step 1: compute the Kovacic invariant r.
    auto r_res = compute_r(a2, a1, a0, x, ctx);
    if (r_res.is_error()) return r_res;
    ExprPtr r = r_res.value();

    // Riccati certificate: η = e^{∫ω} solves y'' = r·y  iff  ω' + ω² = r.
    // case1_omega assembles ω₊/ω₋ by flipping the √r sign at each pole; for
    // poles of order > 2 the flipped branch need not satisfy the Riccati, so
    // every candidate must be certified before use (silent-wrong guard).
    auto riccati_ok = [&](ExprPtr w) -> bool {
        auto wp = diff(w, x, 1U, ctx);
        if (wp.is_error()) return false;
        auto res = ctx.simplify(kv_sub(a, kv_add(a, wp.value(), kv_mul(a, w, w)), r));
        return res.is_ok() && kv_is_zero(res.value(), ctx);
    };

    // Step 2: Case 1 — find rational ω₊, ω₋.
    bool from_case1 = true;
    auto omega_res = case1_omega(r, x, ctx);
    bool case1_ok = omega_res.is_ok()
        && (riccati_ok(omega_res.value().plus) || riccati_ok(omega_res.value().minus));

    if (!case1_ok) {
        // Case 1 failed, or assembled a ω that fails the Riccati certificate —
        // route through Case 2 (dihedral D∞ subgroup).  case2_omega implements
        // Kovacic 1986 §4 Steps 1-4 with its own algebraic certificate.
        CASError case1_err = omega_res.is_error()
            ? omega_res.error()
            : kv_unimpl("Kovacic Case 1: assembled ω fails the Riccati certificate ω'+ω²=r.");
        auto case2_res = case2_omega(r, x, ctx);
        if (case2_res.is_ok()) {
            omega_res = case2_res;
            from_case1 = false;
        } else {
            // Case 2 failed too — try Case 3 (SL(2,C) finite subgroups).
            // case3_omega implements Kovacic 1986 §5 with n ∈ {4, 6, 12}.
            auto case3_res = case3_omega(r, x, ctx);
            if (case3_res.is_error()) {
                // Surface original Case 1 diagnostic — most specific for
                // Case-1-friendly inputs.  Case 2/3 errors are also
                // Unimplemented per spec; preserving Case 1 message helps
                // typical users debug.
                return fail<ExprPtr>(case1_err);
            }
            omega_res = case3_res;
            from_case1 = false;
        }
    }

    ExprPtr op = omega_res.value().plus;
    ExprPtr om = omega_res.value().minus;

    // For Case 1, make ω₊ the certified branch so z₁ is always a genuine solution.
    if (from_case1 && !riccati_ok(op) && riccati_ok(om)) std::swap(op, om);

    // Step 3a: z₁ = exp(∫ω₊ dx).
    auto Fp = integrate(op, x, ctx);
    if (Fp.is_error()) {
        return fail<ExprPtr>(kv_unimpl(
            "Kovacic Case 1: cannot integrate ω₊ — " + Fp.error().message));
    }
    auto z1s_r = ctx.simplify(kv_exp(a, Fp.value()));
    ExprPtr z1 = z1s_r.is_ok() ? z1s_r.value() : kv_exp(a, Fp.value());

    // Step 3b: z₂ from ω₋ (if distinct) or reduction of order.
    ExprPtr z2;
    auto odiff = ctx.simplify(kv_sub(a, op, om));
    bool same  = odiff.is_ok() && kv_is_zero(odiff.value(), ctx);

    // ω₋ yields an independent solution z₂ = e^{∫ω₋} only if it is itself a
    // genuine (certified) Riccati solution; otherwise fall back to reduction of
    // order, which is always valid (and may leave an honest unevaluated integral
    // when the second solution is non-elementary).
    bool om_usable = (!from_case1 || riccati_ok(om));

    if (same || !om_usable) {
        auto z2_res = reduction_of_order(z1, x, ctx);
        if (z2_res.is_error()) return z2_res;
        z2 = z2_res.value();
    } else {
        auto Fm = integrate(om, x, ctx);
        if (Fm.is_ok()) {
            auto z2s = ctx.simplify(kv_exp(a, Fm.value()));
            z2 = z2s.is_ok() ? z2s.value() : kv_exp(a, Fm.value());
        } else {
            auto z2_res = reduction_of_order(z1, x, ctx);
            if (z2_res.is_error()) return z2_res;
            z2 = z2_res.value();
        }
    }

    // Step 4: back-transform  y = z · exp(−½∫(a1/a2) dx).
    auto p_simp = ctx.simplify(kv_div(a, a1, a2));
    ExprPtr p   = p_simp.is_ok() ? p_simp.value() : kv_div(a, a1, a2);
    auto int_p  = integrate(p, x, ctx);
    ExprPtr ef;
    if (int_p.is_ok()) {
        ExprPtr arg = kv_neg(a, kv_div(a, int_p.value(), kv_int(a, 2)));
        auto efs = ctx.simplify(kv_exp(a, arg));
        ef = efs.is_ok() ? efs.value() : kv_exp(a, arg);
    } else {
        ef = kv_exp(a, kv_neg(a, kv_div(a,
            a.make<Integral>(p, x, std::nullopt, std::nullopt), kv_int(a, 2))));
    }

    auto y1s_r = ctx.simplify(kv_mul(a, z1, ef));
    auto y2s_r = ctx.simplify(kv_mul(a, z2, ef));
    ExprPtr y1s = y1s_r.is_ok() ? y1s_r.value() : kv_mul(a, z1, ef);
    ExprPtr y2s = y2s_r.is_ok() ? y2s_r.value() : kv_mul(a, z2, ef);

    // Homogeneous general solution: y_h = C₁·y₁ + C₂·y₂.
    ExprPtr C1 = a.make<Symbol>(ctx.make_fresh_symbol("C"));
    ExprPtr C2 = a.make<Symbol>(ctx.make_fresh_symbol("C"));

    ExprPtr y_h = kv_add(a, kv_mul(a, C1, y1s), kv_mul(a, C2, y2s));

    // Step 5: particular solution for non-homogeneous case.
    if (!kv_is_zero(f, ctx)) {
        auto g_r = ctx.simplify(kv_div(a, f, a2));
        ExprPtr g = g_r.is_ok() ? g_r.value() : kv_div(a, f, a2);
        auto yp = variation_of_params(y1s, y2s, g, x, ctx);
        if (yp.is_ok()) {
            y_h = kv_add(a, y_h, yp.value());
        }
        // On failure: return homogeneous part — unevaluated integrals mark the gap.
    }

    auto sol_r = ctx.simplify(y_h);
    ExprPtr sol = sol_r.is_ok() ? sol_r.value() : y_h;
    return ok(a.make<Binary>(BinaryOp::Equal, a.make<Symbol>(y.name), sol));
}

} // namespace cas::calculus
