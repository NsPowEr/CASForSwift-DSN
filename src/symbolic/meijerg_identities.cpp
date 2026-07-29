// A7 step 5 — Meijer G rewrite identities (Meijer_G_Slater.md §6):
//   §6.2 power shift        z^mu G(z|a;b) = G(z|a+mu;b+mu)     (DLMF 16.19.2)
//   §6.3 parameter cancel   G(z|a0,a;b,a0) = G(z|a;b)          (DLMF 16.19.3)
//   §6.6 antiderivative     int G(t)dt = z G^{m,n+1}(z|0,a;b,-1) (h=-1 case)
//   §6.1 argument inversion G(1/z|a;b) = G^{n,m}_{q,p}(z|1-b;1-a) (DLMF 16.19.1)
// All four were numerically verified in the spec (Appendice C) — every
// rebuild goes back through make_meijerg so the §2.2 guards re-run.

#include "cas/error.hpp"
#include "cas/meijerg.hpp"
#include "cas/symbolic.hpp"

#include <utility>
#include <vector>

namespace cas::symbolic {

namespace {

[[nodiscard]] ExprPtr int_lit(AstArena& arena, long long v) {
    return arena.make<IntegerLit>(BigInt(v));
}

// Each parameter shifted by mu, engine-simplified.
[[nodiscard]] Result<std::vector<ExprPtr>> shift_all(
    CASContext& ctx, const std::vector<ExprPtr>& params, ExprPtr mu) {
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> out;
    out.reserve(params.size());
    for (ExprPtr p : params) {
        auto shifted = ctx.simplify(arena.make<Sum>(std::vector<ExprPtr>{p, mu}));
        if (shifted.is_error()) return fail<std::vector<ExprPtr>>(shifted.error());
        out.push_back(shifted.value());
    }
    return ok(std::move(out));
}

// Decidable structural equality after simplification.
[[nodiscard]] Result<bool> decidably_equal(CASContext& ctx, ExprPtr x, ExprPtr y) {
    AstArena& arena = ctx.arena();
    auto diff = ctx.simplify(arena.make<Sum>(std::vector<ExprPtr>{
        x, arena.make<Unary>(UnaryOp::Neg, y)}));
    if (diff.is_error()) return fail<bool>(diff.error());
    const auto* lit = expr_cast<IntegerLit>(diff.value());
    return ok(lit != nullptr && lit->value.is_zero());
}

}  // namespace

Result<ExprPtr> meijerg_power_shift(CASContext& ctx, const FuncCall& g, ExprPtr mu) {
    auto view_res = view_meijerg(g);
    if (view_res.is_error()) return fail<ExprPtr>(view_res.error());
    const MeijerGView& v = view_res.value();
    auto a_shifted = shift_all(ctx, v.a, mu);
    if (a_shifted.is_error()) return fail<ExprPtr>(a_shifted.error());
    auto b_shifted = shift_all(ctx, v.b, mu);
    if (b_shifted.is_error()) return fail<ExprPtr>(b_shifted.error());
    return make_meijerg(ctx, v.m, v.n,
                        std::move(a_shifted.value()),
                        std::move(b_shifted.value()), v.z);
}

Result<ExprPtr> meijerg_cancel_common_param(CASContext& ctx, const FuncCall& g) {
    auto view_res = view_meijerg(g);
    if (view_res.is_error()) return fail<ExprPtr>(view_res.error());
    MeijerGView v = view_res.value();

    // §6.3: a parameter that is FIRST in the upper n-group and LAST in the
    // lower (m+1..q)-group cancels. Applied repeatedly until fixpoint. The
    // within-group canonical order (§7.2) makes "first"/"last" positional
    // choices immaterial: we scan every n-group upper against every
    // outside-m lower (the identity holds for any such pair after the
    // group-invariant permutations).
    bool changed = true;
    bool any = false;
    while (changed) {
        changed = false;
        for (std::size_t i = 0; i < v.n && !changed; ++i) {
            for (std::size_t j = v.m; j < v.q && !changed; ++j) {
                auto eq = decidably_equal(ctx, v.a[i], v.b[j]);
                if (eq.is_error()) return fail<ExprPtr>(eq.error());
                if (!eq.value()) continue;
                v.a.erase(v.a.begin() + static_cast<std::ptrdiff_t>(i));
                v.b.erase(v.b.begin() + static_cast<std::ptrdiff_t>(j));
                v.n -= 1U;
                v.p -= 1U;
                v.q -= 1U;
                changed = true;
                any = true;
            }
        }
    }
    (void)any;  // no cancellation -> rebuild is structurally identical
    return make_meijerg(ctx, v.m, v.n, std::move(v.a), std::move(v.b), v.z);
}

Result<ExprPtr> meijerg_invert_argument(CASContext& ctx, const FuncCall& g) {
    auto view_res = view_meijerg(g);
    if (view_res.is_error()) return fail<ExprPtr>(view_res.error());
    const MeijerGView& v = view_res.value();
    AstArena& arena = ctx.arena();

    // §6.1: G^{m,n}_{p,q}(1/z | a; b) = G^{n,m}_{q,p}(z | 1-b ; 1-a) — the
    // caller passes the node whose argument it wants rewritten from z to
    // 1/z; here we build the swapped node with argument inverted.
    auto one_minus = [&](const std::vector<ExprPtr>& params)
        -> Result<std::vector<ExprPtr>> {
        std::vector<ExprPtr> out;
        out.reserve(params.size());
        for (ExprPtr p : params) {
            auto val = ctx.simplify(arena.make<Sum>(std::vector<ExprPtr>{
                int_lit(arena, 1), arena.make<Unary>(UnaryOp::Neg, p)}));
            if (val.is_error()) return fail<std::vector<ExprPtr>>(val.error());
            out.push_back(val.value());
        }
        return ok(std::move(out));
    };
    auto new_a = one_minus(v.b);  // 1-b becomes the upper list (size q)
    if (new_a.is_error()) return fail<ExprPtr>(new_a.error());
    auto new_b = one_minus(v.a);  // 1-a becomes the lower list (size p)
    if (new_b.is_error()) return fail<ExprPtr>(new_b.error());
    auto inv_z = ctx.simplify(arena.make<Binary>(BinaryOp::Pow, v.z,
                                                 int_lit(arena, -1)));
    if (inv_z.is_error()) return inv_z;
    return make_meijerg(ctx, v.n, v.m,
                        std::move(new_a.value()), std::move(new_b.value()),
                        inv_z.value());
}

namespace {

// Shared builder for the theta-shift family (§6.5/§6.6):
//   z^h (d/dz)^h G^{m,n}_{p,q}(z|a;b) = G^{m,n+1}_{p+1,q+1}(z | 0,a ; b, h).
// The prepended 0 joins the upper n-group; the appended h joins the lower
// (m+1..q)-group. make_meijerg re-runs the §2.2 pole checks on the new
// lists — an overlap there is a structured refusal, never a silently wrong
// rewrite.
[[nodiscard]] Result<ExprPtr> make_theta_shifted(
    CASContext& ctx, const MeijerGView& v, long long h) {
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> new_a;
    new_a.reserve(v.p + 1U);
    new_a.push_back(int_lit(arena, 0));
    new_a.insert(new_a.end(), v.a.begin(), v.a.end());
    std::vector<ExprPtr> new_b;
    new_b.reserve(v.q + 1U);
    new_b.insert(new_b.end(), v.b.begin(),
                 v.b.begin() + static_cast<std::ptrdiff_t>(v.m));
    std::vector<ExprPtr> b_tail(v.b.begin() + static_cast<std::ptrdiff_t>(v.m),
                                v.b.end());
    b_tail.push_back(int_lit(arena, h));
    new_b.insert(new_b.end(), b_tail.begin(), b_tail.end());
    return make_meijerg(ctx, v.m, v.n + 1U,
                        std::move(new_a), std::move(new_b), v.z);
}

}  // namespace

Result<ExprPtr> meijerg_antiderivative(CASContext& ctx, const FuncCall& g) {
    auto view_res = view_meijerg(g);
    if (view_res.is_error()) return fail<ExprPtr>(view_res.error());
    const MeijerGView& v = view_res.value();

    // §6.6 (h = -1): int G(t|a;b) dt = z * G^{m,n+1}_{p+1,q+1}(z|0,a;b,-1) + C.
    auto shifted = make_theta_shifted(ctx, v, -1);
    if (shifted.is_error()) return shifted;
    return ok(ctx.arena().make<Product>(
        std::vector<ExprPtr>{v.z, shifted.value()}));
}

Result<ExprPtr> meijerg_derivative_shift(CASContext& ctx, const FuncCall& g) {
    auto view_res = view_meijerg(g);
    if (view_res.is_error()) return fail<ExprPtr>(view_res.error());
    // §6.5 (h = +1, certified): z * G'(z) = G^{m,n+1}_{p+1,q+1}(z|0,a;b,+1).
    return make_theta_shifted(ctx, view_res.value(), 1);
}

}  // namespace cas::symbolic
