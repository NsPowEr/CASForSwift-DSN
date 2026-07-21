// A7 step 4 — from_meijerg (table §5 inverse fast path, shape-gated) and
// expand_meijerg_nodes (tree map). General fallback = slater_expand
// (meijerg_slater.cpp) — the table is ONLY a fast path (CLAUDE.md cat. 8).
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Meijer_G_Slater.md §5, §9.4.
//
// Sign convention (Mellin half-line): the §5 pairs are exact for a POSITIVE
// elementary argument u (G's argument z = c*u^2 is even in u, so the sign of
// u is not encoded in the node — the standard Mellin-transform convention on
// (0, infinity), the domain the Adamchik-Marichev pipeline §9 operates on).
// The same convention already applies in the forward direction (§5 was
// verified numerically on positive sample points, Appendice C).

#include "cas/error.hpp"
#include "cas/meijerg.hpp"
#include "cas/symbolic.hpp"
#include "simplify_impl.hpp"

#include <optional>
#include <utility>
#include <vector>

namespace cas::symbolic {

using detail::is_zero_expr;

namespace {

[[nodiscard]] ExprPtr int_lit(AstArena& arena, long long v) {
    return arena.make<IntegerLit>(BigInt(v));
}

[[nodiscard]] ExprPtr sqrt_pi(AstArena& arena) {
    return arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{arena.make<Constant>(MathConstant::Pi)});
}

[[nodiscard]] ExprPtr inv_sqrt_pi(AstArena& arena) {
    return arena.make<Binary>(BinaryOp::Pow, sqrt_pi(arena), int_lit(arena, -1));
}

[[nodiscard]] bool is_rat(ExprPtr e, long long num, long long den) {
    if (den == 1) {
        const auto* il = expr_cast<IntegerLit>(e);
        if (il != nullptr) return il->value == BigInt(num);
    }
    const auto* rl = expr_cast<RationalLit>(e);
    return rl != nullptr && rl->numerator == BigInt(num)
        && rl->denominator == BigInt(den);
}

// Matches z == c * u^2 with c a literal rational, returning {c, u}.
// Recognizes the canonical shapes the engine produces for u^2, -u^2,
// u^2/4, -u^2/4 (bare Pow, Product with rational-literal coefficients,
// Neg wrapper). Anything else: no match (the caller falls back to Slater —
// never a wrong extraction).
struct ScaledSquare {
    Rational coeff;
    ExprPtr base;
};
[[nodiscard]] std::optional<ScaledSquare> match_scaled_square(ExprPtr z) {
    if (z == nullptr) return std::nullopt;
    if (const auto* un = expr_cast<Unary>(z); un != nullptr && un->op == UnaryOp::Neg) {
        auto inner = match_scaled_square(un->operand);
        if (!inner.has_value()) return std::nullopt;
        inner->coeff = Rational(BigInt(0)) - inner->coeff;
        return inner;
    }
    if (const auto* pw = expr_cast<Binary>(z);
        pw != nullptr && pw->op == BinaryOp::Pow) {
        const auto* e2 = expr_cast<IntegerLit>(pw->right);
        if (e2 != nullptr && e2->value == BigInt(2)) {
            return ScaledSquare{Rational(BigInt(1)), pw->left};
        }
        return std::nullopt;
    }
    if (const auto* prod = expr_cast<Product>(z)) {
        Rational coeff(BigInt(1));
        ExprPtr base = nullptr;
        for (ExprPtr f : prod->factors) {
            if (const auto* il = expr_cast<IntegerLit>(f)) {
                coeff = coeff * Rational(il->value);
                continue;
            }
            if (const auto* rl = expr_cast<RationalLit>(f)) {
                coeff = coeff * Rational(rl->numerator, rl->denominator);
                continue;
            }
            const auto* pw = expr_cast<Binary>(f);
            if (pw != nullptr && pw->op == BinaryOp::Pow && base == nullptr) {
                const auto* e2 = expr_cast<IntegerLit>(pw->right);
                if (e2 != nullptr && e2->value == BigInt(2)) {
                    base = pw->left;
                    continue;
                }
            }
            return std::nullopt;  // non-square, non-literal factor
        }
        if (base == nullptr) return std::nullopt;
        return ScaledSquare{coeff, base};
    }
    return std::nullopt;
}

[[nodiscard]] bool coeff_is(const ScaledSquare& s, long long num, long long den) {
    return s.coeff == Rational(BigInt(num), BigInt(den));
}

// Inverse table entries for the trig/Bessel family
// G^{1,0}_{0,2}(c*u^2 | ; b1 ; b2) (§5.2/5.3/5.8). Returns the value of the
// BARE G node (the caller's own prefactor completes the elementary form).
[[nodiscard]] Result<ExprPtr> inverse_0_2(
    CASContext& ctx, const MeijerGView& v, const FuncCall& g) {
    AstArena& arena = ctx.arena();
    auto shape = match_scaled_square(v.z);
    if (!shape.has_value()) return slater_expand(ctx, g);
    ExprPtr u = shape->base;

    // sin u = sqrt(pi) G(u^2/4 | ; 1/2 ; 0)  =>  G = sin(u)/sqrt(pi)  (§5.2)
    if (is_rat(v.b[0], 1, 2) && is_zero_expr(v.b[1]) && coeff_is(*shape, 1, 4)) {
        return ok(arena.make<Product>(std::vector<ExprPtr>{
            inv_sqrt_pi(arena),
            arena.make<FuncCall>(BuiltinOp::Sin, std::vector<ExprPtr>{u})}));
    }
    // cos u = sqrt(pi) G(u^2/4 | ; 0 ; 1/2)  =>  G = cos(u)/sqrt(pi)  (§5.2)
    if (is_zero_expr(v.b[0]) && is_rat(v.b[1], 1, 2) && coeff_is(*shape, 1, 4)) {
        return ok(arena.make<Product>(std::vector<ExprPtr>{
            inv_sqrt_pi(arena),
            arena.make<FuncCall>(BuiltinOp::Cos, std::vector<ExprPtr>{u})}));
    }
    // cosh u = sqrt(pi) G(-u^2/4 | ; 0 ; 1/2)  =>  G = cosh(u)/sqrt(pi)  (§5.3)
    if (is_zero_expr(v.b[0]) && is_rat(v.b[1], 1, 2) && coeff_is(*shape, -1, 4)) {
        return ok(arena.make<Product>(std::vector<ExprPtr>{
            inv_sqrt_pi(arena),
            arena.make<FuncCall>(BuiltinOp::Cosh, std::vector<ExprPtr>{u})}));
    }
    // sinh u = (sqrt(pi)/2) u G(-u^2/4 | ; 0 ; -1/2)
    //   =>  G = 2 sinh(u) / (sqrt(pi) u)  (§5.3 real primary form)
    if (is_zero_expr(v.b[0]) && is_rat(v.b[1], -1, 2) && coeff_is(*shape, -1, 4)) {
        return ok(arena.make<Product>(std::vector<ExprPtr>{
            int_lit(arena, 2), inv_sqrt_pi(arena),
            arena.make<Binary>(BinaryOp::Pow, u, int_lit(arena, -1)),
            arena.make<FuncCall>(BuiltinOp::Sinh, std::vector<ExprPtr>{u})}));
    }
    // J_nu(u) = G(u^2/4 | ; nu/2 ; -nu/2)  (§5.8): b2 == -b1 decidably.
    if (coeff_is(*shape, 1, 4)) {
        auto sum = ctx.simplify(arena.make<Sum>(std::vector<ExprPtr>{
            v.b[0], v.b[1]}));
        if (sum.is_error()) return sum;
        if (is_zero_expr(sum.value())) {
            auto nu = ctx.simplify(arena.make<Product>(std::vector<ExprPtr>{
                int_lit(arena, 2), v.b[0]}));
            if (nu.is_error()) return nu;
            return ok(arena.make<FuncCall>(BuiltinOp::BesselJ,
                std::vector<ExprPtr>{nu.value(), u}));
        }
    }
    return slater_expand(ctx, g);
}

// Inverse entries for G^{1,2}_{2,2} (§5.4/5.6): ln (table-ONLY — Slater's
// series would keep a 2F1 the engine does not fold back), arctan, arcsin.
[[nodiscard]] Result<ExprPtr> inverse_2_2(
    CASContext& ctx, const MeijerGView& v, const FuncCall& g) {
    AstArena& arena = ctx.arena();
    // ln(1+z) = G(z | 1,1 ; 1 ; 0)  =>  G = ln(1+z), z as-is (§5.4).
    if (is_rat(v.a[0], 1, 1) && is_rat(v.a[1], 1, 1)
        && is_rat(v.b[0], 1, 1) && is_zero_expr(v.b[1])) {
        ExprPtr one_plus_z = arena.make<Sum>(std::vector<ExprPtr>{
            int_lit(arena, 1), v.z});
        return ok(arena.make<FuncCall>(BuiltinOp::Ln,
            std::vector<ExprPtr>{one_plus_z}));
    }
    auto shape = match_scaled_square(v.z);
    if (shape.has_value()) {
        ExprPtr u = shape->base;
        // arctan u = (1/2) G(u^2 | 1/2,1 ; 1/2 ; 0)  =>  G = 2 arctan(u) (§5.6)
        // (factory canonicalization sorts the a-group: order-insensitive
        // membership check via the two literal tests below).
        const bool a_half_one =
            (is_rat(v.a[0], 1, 2) && is_rat(v.a[1], 1, 1))
            || (is_rat(v.a[0], 1, 1) && is_rat(v.a[1], 1, 2));
        if (a_half_one && is_rat(v.b[0], 1, 2) && is_zero_expr(v.b[1])
            && coeff_is(*shape, 1, 1)) {
            return ok(arena.make<Product>(std::vector<ExprPtr>{
                int_lit(arena, 2),
                arena.make<FuncCall>(BuiltinOp::Atan, std::vector<ExprPtr>{u})}));
        }
        // arcsin u = u/(2 sqrt(pi)) G(-u^2 | 1/2,1/2 ; 0 ; -1/2)
        //   =>  G = 2 sqrt(pi) arcsin(u) / u  (§5.6)
        if (is_rat(v.a[0], 1, 2) && is_rat(v.a[1], 1, 2)
            && is_zero_expr(v.b[0]) && is_rat(v.b[1], -1, 2)
            && coeff_is(*shape, -1, 1)) {
            return ok(arena.make<Product>(std::vector<ExprPtr>{
                int_lit(arena, 2), sqrt_pi(arena),
                arena.make<Binary>(BinaryOp::Pow, u, int_lit(arena, -1)),
                arena.make<FuncCall>(BuiltinOp::Asin, std::vector<ExprPtr>{u})}));
        }
    }
    return slater_expand(ctx, g);
}

}  // namespace

Result<ExprPtr> from_meijerg(CASContext& ctx, const FuncCall& g) {
    auto view_res = view_meijerg(g);
    if (view_res.is_error()) return fail<ExprPtr>(view_res.error());
    const MeijerGView& v = view_res.value();
    AstArena& arena = ctx.arena();

    Result<ExprPtr> raw = fail<ExprPtr>(CASError{
        .kind = CASErrorKind::InternalError, .message = "unset"});
    if (v.m == 1U && v.n == 0U && v.p == 0U && v.q == 2U) {
        raw = inverse_0_2(ctx, v, g);
    } else if (v.m == 1U && v.n == 2U && v.p == 2U && v.q == 2U) {
        raw = inverse_2_2(ctx, v, g);
    } else if (v.m == 2U && v.n == 0U && v.p == 1U && v.q == 2U
               && is_rat(v.a[0], 1, 1)
               && (is_zero_expr(v.b[0]) || is_zero_expr(v.b[1]))) {
        // Γ(a,z) = G^{2,0}_{1,2}(z | 1 ; 0,a)  (§5.9, verified mpmath).
        // The b-group {0, a} is canonicalized by the factory ⇒ order-insensitive:
        // the literal 0 is the "0" of the formula, the other member is the order a
        // (a = 0 ⇒ both members zero ⇒ Γ(0,z), still exact).
        ExprPtr order = is_zero_expr(v.b[0]) ? v.b[1] : v.b[0];
        raw = ok(arena.make<FuncCall>(BuiltinOp::GammaIncomplete,
            std::vector<ExprPtr>{order, v.z}));
    } else if (v.m == 1U && v.n == 1U && v.p == 1U && v.q == 2U
               && is_rat(v.a[0], 1, 1) && is_zero_expr(v.b[1])) {
        // Shape G^{1,1}_{1,2}(· | 1 ; b0 ; 0). Two §5 entries share it:
        //   erf u = (1/sqrt(pi)) G(u^2 | 1 ; 1/2 ; 0) ⇒ sqrt(pi) erf(u)  (§5.7),
        //   γ(a,z) =              G(z   | 1 ;  a  ; 0)                     (§5.9).
        // erf is the elementary special case (b0 = 1/2 AND a square argument);
        // otherwise the exact inverse is the lower incomplete gamma.
        auto shape = match_scaled_square(v.z);
        if (is_rat(v.b[0], 1, 2) && shape.has_value() && coeff_is(*shape, 1, 1)) {
            raw = ok(arena.make<Product>(std::vector<ExprPtr>{
                sqrt_pi(arena),
                arena.make<FuncCall>(BuiltinOp::Erf,
                    std::vector<ExprPtr>{shape->base})}));
        } else {
            raw = ok(arena.make<FuncCall>(BuiltinOp::GammaIncompleteLower,
                std::vector<ExprPtr>{v.b[0], v.z}));
        }
    } else {
        // No table signature: the GENERAL path (§3.2).
        raw = slater_expand(ctx, g);
    }
    if (raw.is_error()) return raw;
    return ctx.simplify(raw.value());
}

Result<ExprPtr> expand_meijerg_nodes(CASContext& ctx, ExprPtr expr) {
    if (expr == nullptr) return ok(expr);
    AstArena& arena = ctx.arena();

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if (call->func_id == BuiltinOp::MeijerG) {
            auto expanded = from_meijerg(ctx, *call);
            if (expanded.is_ok()) return expanded;
            // Unimplemented expansion: the G node legitimately stays (§9.4).
            if (expanded.error().kind == CASErrorKind::Unimplemented) {
                return ok(expr);
            }
            return expanded;
        }
        std::vector<ExprPtr> new_args;
        new_args.reserve(call->args.size());
        bool changed = false;
        for (ExprPtr arg : call->args) {
            auto mapped = expand_meijerg_nodes(ctx, arg);
            if (mapped.is_error()) return mapped;
            changed = changed || mapped.value() != arg;
            new_args.push_back(mapped.value());
        }
        if (!changed) return ok(expr);
        return ok(arena.make<FuncCall>(call->func_id, std::move(new_args)));
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        std::vector<ExprPtr> new_factors;
        new_factors.reserve(prod->factors.size());
        bool changed = false;
        for (ExprPtr f : prod->factors) {
            auto mapped = expand_meijerg_nodes(ctx, f);
            if (mapped.is_error()) return mapped;
            changed = changed || mapped.value() != f;
            new_factors.push_back(mapped.value());
        }
        if (!changed) return ok(expr);
        return ok(arena.make<Product>(std::move(new_factors)));
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> new_terms;
        new_terms.reserve(sum->terms.size());
        bool changed = false;
        for (ExprPtr t : sum->terms) {
            auto mapped = expand_meijerg_nodes(ctx, t);
            if (mapped.is_error()) return mapped;
            changed = changed || mapped.value() != t;
            new_terms.push_back(mapped.value());
        }
        if (!changed) return ok(expr);
        return ok(arena.make<Sum>(std::move(new_terms)));
    }
    if (const auto* bin = expr_cast<Binary>(expr)) {
        auto left = expand_meijerg_nodes(ctx, bin->left);
        if (left.is_error()) return left;
        auto right = expand_meijerg_nodes(ctx, bin->right);
        if (right.is_error()) return right;
        if (left.value() == bin->left && right.value() == bin->right)
            return ok(expr);
        return ok(arena.make<Binary>(bin->op, left.value(), right.value()));
    }
    if (const auto* un = expr_cast<Unary>(expr)) {
        auto inner = expand_meijerg_nodes(ctx, un->operand);
        if (inner.is_error()) return inner;
        if (inner.value() == un->operand) return ok(expr);
        return ok(arena.make<Unary>(un->op, inner.value()));
    }
    return ok(expr);  // leaves (literals, symbols, constants)
}

}  // namespace cas::symbolic
