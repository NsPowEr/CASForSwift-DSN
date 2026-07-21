// A7 step 3 — elementary/hypergeometric → Meijer G conversion.
// to_meijerg(): table §5 fast path + general pFq bridge §3.1 (DLMF 16.18.1).
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Meijer_G_Slater.md.
// Every table entry cites its §5.x line; the formulas were numerically
// verified against mpmath 1.3.0 in the spec (Appendice C) — do not "adjust"
// a prefactor here without re-running that verification.

#include "cas/error.hpp"
#include "cas/meijerg.hpp"
#include "cas/symbolic.hpp"
#include "simplify_impl.hpp"

#include <string>
#include <utility>
#include <vector>

namespace cas::symbolic {

using detail::is_constant_expr;
using detail::is_zero_expr;

namespace {

[[nodiscard]] ExprPtr int_lit(AstArena& arena, long long v) {
    return arena.make<IntegerLit>(BigInt(v));
}

[[nodiscard]] ExprPtr rat_lit(AstArena& arena, long long num, long long den) {
    return arena.make<RationalLit>(BigInt(num), BigInt(den));
}

[[nodiscard]] ExprPtr negate(AstArena& arena, ExprPtr e) {
    return arena.make<Unary>(UnaryOp::Neg, e);
}

[[nodiscard]] ExprPtr sqrt_pi(AstArena& arena) {
    return arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{arena.make<Constant>(MathConstant::Pi)});
}

// Canonical form of an argument expression built for a G node (z, z^2/4,
// w-1, nu/2, ...): run it through the engine so equivalent constructions
// share structure (Regola 2) and downstream pattern matching sees one shape.
[[nodiscard]] Result<ExprPtr> simp(CASContext& ctx, ExprPtr e) {
    return ctx.simplify(e);
}

// z^2 / 4 (argument of the sin/cos/sinh/cosh/BesselJ entries §5.2/5.3/5.8).
[[nodiscard]] Result<ExprPtr> square_over_four(CASContext& ctx, ExprPtr z) {
    AstArena& arena = ctx.arena();
    ExprPtr sq = arena.make<Binary>(BinaryOp::Pow, z, int_lit(arena, 2));
    return simp(ctx, arena.make<Product>(std::vector<ExprPtr>{
        rat_lit(arena, 1, 4), sq}));
}

// Decidably a nonpositive integer (0, -1, -2, ...) after simplification —
// the exact Gamma-pole guard of §3.1. Symbolic values that do not reduce to
// an integer literal are NOT flagged (no false rejection; matches the
// exact-when-decidable stance of make_meijerg's pole check).
[[nodiscard]] Result<bool> is_decidably_nonpositive_integer(
    CASContext& ctx, ExprPtr e) {
    auto simplified = ctx.simplify(e);
    if (simplified.is_error()) return fail<bool>(simplified.error());
    const auto* lit = expr_cast<IntegerLit>(simplified.value());
    return ok(lit != nullptr && !lit->value.is_positive());
}

[[nodiscard]] CASError unsupported(const std::string& shape) {
    return make_unimplemented_error(
        UnimplementedInfo{
            .module = "symbolic",
            .function = "to_meijerg",
            .input_shape = shape,
            .reason = error::reason_codes::GENERIC,
            .suggestion = "No Meijer G representation implemented for this "
                          "shape (spec Meijer_G_Slater.md §5/§3.1)",
            .ticket = "A7"},
        "to_meijerg: unsupported expression shape");
}

// prefactor * G, with a plain G when there is no prefactor.
[[nodiscard]] ExprPtr with_prefactor(
    AstArena& arena, std::vector<ExprPtr> prefactor, ExprPtr g) {
    if (prefactor.empty()) return g;
    prefactor.push_back(g);
    return arena.make<Product>(std::move(prefactor));
}

// §5.1: e^w = G^{1,0}_{0,1}(-w | ; 0).
[[nodiscard]] Result<ExprPtr> exp_to_g(CASContext& ctx, ExprPtr w) {
    AstArena& arena = ctx.arena();
    auto z = simp(ctx, negate(arena, w));
    if (z.is_error()) return z;
    return make_meijerg(ctx, 1U, 0U, {}, {int_lit(arena, 0)}, z.value());
}

// §5.2: sin w = sqrt(pi) * G^{1,0}_{0,2}(w^2/4 | ; 1/2 ; 0)
//       cos w = sqrt(pi) * G^{1,0}_{0,2}(w^2/4 | ; 0 ; 1/2)
// §5.3: cosh w = sqrt(pi) * G^{1,0}_{0,2}(-w^2/4 | ; 0 ; 1/2)
//       sinh w = (sqrt(pi)/2) * w * G^{1,0}_{0,2}(-w^2/4 | ; 0 ; -1/2)
//       (real primary forms; the -i*sqrt(pi) sinh variant is FORBIDDEN, §11)
[[nodiscard]] Result<ExprPtr> trig_pair_to_g(
    CASContext& ctx, ExprPtr w, bool negate_arg,
    ExprPtr b1, ExprPtr b2, std::vector<ExprPtr> prefactor) {
    AstArena& arena = ctx.arena();
    auto arg = square_over_four(ctx, w);
    if (arg.is_error()) return arg;
    ExprPtr z = arg.value();
    if (negate_arg) {
        auto neg = simp(ctx, negate(arena, z));
        if (neg.is_error()) return neg;
        z = neg.value();
    }
    auto g = make_meijerg(ctx, 1U, 0U, {}, {b1, b2}, z);
    if (g.is_error()) return g;
    return ok(with_prefactor(arena, std::move(prefactor), g.value()));
}

}  // namespace

Result<ExprPtr> pfq_to_meijerg(
    CASContext& ctx,
    std::vector<ExprPtr> alpha, std::vector<ExprPtr> beta, ExprPtr z) {
    AstArena& arena = ctx.arena();
    const std::size_t p = alpha.size();
    const std::size_t q = beta.size();

    // §3.1 precondition: the series/bridge is for p <= q+1.
    if (p > q + 1U) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "pfq_to_meijerg: requires p <= q+1 (got p="
                + std::to_string(p) + ", q=" + std::to_string(q) + ")"});
    }
    // §3.1 Gamma-pole guards: alpha_k in Z<=0 -> pFq is a polynomial and
    // Gamma(alpha_k) in the prefactor is a pole (formula inapplicable);
    // beta_k in Z<=0 -> the pFq itself is undefined.
    for (ExprPtr a_k : alpha) {
        auto bad = is_decidably_nonpositive_integer(ctx, a_k);
        if (bad.is_error()) return fail<ExprPtr>(bad.error());
        if (bad.value()) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::InvalidArgument,
                .message = "pfq_to_meijerg: upper parameter is a nonpositive "
                           "integer (polynomial pFq; DLMF 16.18.1 does not "
                           "apply)"});
        }
    }
    for (ExprPtr b_k : beta) {
        auto bad = is_decidably_nonpositive_integer(ctx, b_k);
        if (bad.is_error()) return fail<ExprPtr>(bad.error());
        if (bad.value()) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::InvalidArgument,
                .message = "pfq_to_meijerg: lower parameter is a nonpositive "
                           "integer (pFq undefined)"});
        }
    }

    // a_G = (1-alpha_1 .. 1-alpha_p), all in the FIRST upper group (n = p);
    // b_G = (0 | 1-beta_1 .. 1-beta_q), first lower group = (0) (m = 1).
    std::vector<ExprPtr> a_g;
    a_g.reserve(p);
    for (ExprPtr a_k : alpha) {
        auto shifted = simp(ctx, arena.make<Sum>(std::vector<ExprPtr>{
            int_lit(arena, 1), negate(arena, a_k)}));
        if (shifted.is_error()) return shifted;
        a_g.push_back(shifted.value());
    }
    std::vector<ExprPtr> b_g;
    b_g.reserve(q + 1U);
    b_g.push_back(int_lit(arena, 0));
    for (ExprPtr b_k : beta) {
        auto shifted = simp(ctx, arena.make<Sum>(std::vector<ExprPtr>{
            int_lit(arena, 1), negate(arena, b_k)}));
        if (shifted.is_error()) return shifted;
        b_g.push_back(shifted.value());
    }

    auto z_g = simp(ctx, negate(arena, z));
    if (z_g.is_error()) return z_g;

    auto g = make_meijerg(ctx, 1U, p, std::move(a_g), std::move(b_g),
                          z_g.value());
    if (g.is_error()) return g;

    // Prefactor prod Gamma(beta_k) / prod Gamma(alpha_k).
    std::vector<ExprPtr> prefactor;
    prefactor.reserve(p + q);
    for (ExprPtr b_k : beta) {
        prefactor.push_back(arena.make<FuncCall>(BuiltinOp::Gamma,
            std::vector<ExprPtr>{b_k}));
    }
    for (ExprPtr a_k : alpha) {
        prefactor.push_back(arena.make<Binary>(BinaryOp::Pow,
            arena.make<FuncCall>(BuiltinOp::Gamma, std::vector<ExprPtr>{a_k}),
            int_lit(arena, -1)));
    }
    return ok(with_prefactor(arena, std::move(prefactor), g.value()));
}

Result<ExprPtr> to_meijerg(CASContext& ctx, ExprPtr expr) {
    if (expr == nullptr) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "to_meijerg: null expression"});
    }
    AstArena& arena = ctx.arena();

    // Pow shapes: E^w (routes to §5.1) and the binomial (1+z)^{-a} (§5.5).
    if (const auto* pow = expr_cast<Binary>(expr);
        pow != nullptr && pow->op == BinaryOp::Pow) {
        if (is_constant_expr(pow->left, MathConstant::E)) {
            return exp_to_g(ctx, pow->right);
        }
        // §5.5: (1+z)^{-a} = (1/Gamma(a)) * G^{1,1}_{1,1}(z | 1-a ; 0),
        // with a = -exponent, z = base-1. A nonnegative-integer exponent is
        // a plain polynomial (and a = -e <= 0 hits the Gamma pole): refuse.
        auto a_val = simp(ctx, negate(arena, pow->right));
        if (a_val.is_error()) return a_val;
        auto a_pole = is_decidably_nonpositive_integer(ctx, a_val.value());
        if (a_pole.is_error()) return fail<ExprPtr>(a_pole.error());
        if (a_pole.value()) {
            return fail<ExprPtr>(unsupported(
                "Pow with nonnegative-integer exponent (polynomial)"));
        }
        auto z_g = simp(ctx, arena.make<Sum>(std::vector<ExprPtr>{
            pow->left, int_lit(arena, -1)}));
        if (z_g.is_error()) return z_g;
        auto one_minus_a = simp(ctx, arena.make<Sum>(std::vector<ExprPtr>{
            int_lit(arena, 1), negate(arena, a_val.value())}));
        if (one_minus_a.is_error()) return one_minus_a;
        auto g = make_meijerg(ctx, 1U, 1U, {one_minus_a.value()},
                              {int_lit(arena, 0)}, z_g.value());
        if (g.is_error()) return g;
        ExprPtr inv_gamma = arena.make<Binary>(BinaryOp::Pow,
            arena.make<FuncCall>(BuiltinOp::Gamma,
                std::vector<ExprPtr>{a_val.value()}),
            int_lit(arena, -1));
        return ok(with_prefactor(arena, {inv_gamma}, g.value()));
    }

    // §5.5 first form: w^a * e^{-w} = G^{1,0}_{0,1}(w | ; a) — strict
    // two-factor fast path (the general route is per-factor conversion plus
    // the §6.2 power-shift identity, a later step).
    if (const auto* prod = expr_cast<Product>(expr);
        prod != nullptr && prod->factors.size() == 2U) {
        for (std::size_t i = 0; i < 2U; ++i) {
            const auto* pw = expr_cast<Binary>(prod->factors[i]);
            ExprPtr other = prod->factors[1U - i];
            if (pw == nullptr || pw->op != BinaryOp::Pow) continue;
            ExprPtr exp_arg = nullptr;
            if (const auto* call = expr_cast<FuncCall>(other);
                call != nullptr && call->func_id == BuiltinOp::Exp
                && call->args.size() == 1U) {
                exp_arg = call->args.front();
            } else if (const auto* epow = expr_cast<Binary>(other);
                       epow != nullptr && epow->op == BinaryOp::Pow
                       && is_constant_expr(epow->left, MathConstant::E)) {
                exp_arg = epow->right;
            }
            if (exp_arg == nullptr) continue;
            // Require exp argument == -(pow base) exactly (engine-decided).
            auto residue = simp(ctx, arena.make<Sum>(std::vector<ExprPtr>{
                exp_arg, pw->left}));
            if (residue.is_error()) return residue;
            if (!is_zero_expr(residue.value())) continue;
            return make_meijerg(ctx, 1U, 0U, {}, {pw->right}, pw->left);
        }
    }

    const auto* call = expr_cast<FuncCall>(expr);
    if (call == nullptr) return fail<ExprPtr>(unsupported("non-FuncCall"));

    switch (call->func_id) {
        case BuiltinOp::Exp:
            if (call->args.size() != 1U) break;
            return exp_to_g(ctx, call->args.front());

        case BuiltinOp::Sin:  // §5.2
            if (call->args.size() != 1U) break;
            return trig_pair_to_g(ctx, call->args.front(), false,
                rat_lit(arena, 1, 2), int_lit(arena, 0), {sqrt_pi(arena)});

        case BuiltinOp::Cos:  // §5.2
            if (call->args.size() != 1U) break;
            return trig_pair_to_g(ctx, call->args.front(), false,
                int_lit(arena, 0), rat_lit(arena, 1, 2), {sqrt_pi(arena)});

        case BuiltinOp::Cosh:  // §5.3 (b1 = 0: branch-safe on negative arg)
            if (call->args.size() != 1U) break;
            return trig_pair_to_g(ctx, call->args.front(), true,
                int_lit(arena, 0), rat_lit(arena, 1, 2), {sqrt_pi(arena)});

        case BuiltinOp::Sinh:  // §5.3 real primary form (never the -i form)
            if (call->args.size() != 1U) break;
            return trig_pair_to_g(ctx, call->args.front(), true,
                int_lit(arena, 0), rat_lit(arena, -1, 2),
                {rat_lit(arena, 1, 2), sqrt_pi(arena), call->args.front()});

        case BuiltinOp::Ln:
        case BuiltinOp::Log: {  // §5.4: ln(1+z) = G^{1,2}_{2,2}(z | 1,1 ; 1 ; 0)
            if (call->args.size() != 1U) break;
            auto z_g = simp(ctx, arena.make<Sum>(std::vector<ExprPtr>{
                call->args.front(), int_lit(arena, -1)}));
            if (z_g.is_error()) return z_g;
            return make_meijerg(ctx, 1U, 2U,
                {int_lit(arena, 1), int_lit(arena, 1)},
                {int_lit(arena, 1), int_lit(arena, 0)}, z_g.value());
        }

        case BuiltinOp::Atan: {  // §5.6: (1/2) G^{1,2}_{2,2}(w^2 | 1/2,1 ; 1/2 ; 0)
            if (call->args.size() != 1U) break;
            auto z_g = simp(ctx, arena.make<Binary>(BinaryOp::Pow,
                call->args.front(), int_lit(arena, 2)));
            if (z_g.is_error()) return z_g;
            auto g = make_meijerg(ctx, 1U, 2U,
                {rat_lit(arena, 1, 2), int_lit(arena, 1)},
                {rat_lit(arena, 1, 2), int_lit(arena, 0)}, z_g.value());
            if (g.is_error()) return g;
            return ok(with_prefactor(arena, {rat_lit(arena, 1, 2)}, g.value()));
        }

        case BuiltinOp::Asin: {  // §5.6: w/(2 sqrt(pi)) G^{1,2}_{2,2}(-w^2 | ...)
            if (call->args.size() != 1U) break;
            auto z_g = simp(ctx, negate(arena,
                arena.make<Binary>(BinaryOp::Pow, call->args.front(),
                                   int_lit(arena, 2))));
            if (z_g.is_error()) return z_g;
            auto g = make_meijerg(ctx, 1U, 2U,
                {rat_lit(arena, 1, 2), rat_lit(arena, 1, 2)},
                {int_lit(arena, 0), rat_lit(arena, -1, 2)}, z_g.value());
            if (g.is_error()) return g;
            ExprPtr inv_two_sqrt_pi = arena.make<Binary>(BinaryOp::Pow,
                arena.make<Product>(std::vector<ExprPtr>{
                    int_lit(arena, 2), sqrt_pi(arena)}),
                int_lit(arena, -1));
            return ok(with_prefactor(arena,
                {call->args.front(), inv_two_sqrt_pi}, g.value()));
        }

        case BuiltinOp::Erf: {  // §5.7: (1/sqrt(pi)) G^{1,1}_{1,2}(w^2 | 1 ; 1/2 ; 0)
            if (call->args.size() != 1U) break;
            auto z_g = simp(ctx, arena.make<Binary>(BinaryOp::Pow,
                call->args.front(), int_lit(arena, 2)));
            if (z_g.is_error()) return z_g;
            auto g = make_meijerg(ctx, 1U, 1U, {int_lit(arena, 1)},
                {rat_lit(arena, 1, 2), int_lit(arena, 0)}, z_g.value());
            if (g.is_error()) return g;
            ExprPtr inv_sqrt_pi = arena.make<Binary>(BinaryOp::Pow,
                sqrt_pi(arena), int_lit(arena, -1));
            return ok(with_prefactor(arena, {inv_sqrt_pi}, g.value()));
        }

        case BuiltinOp::GammaIncomplete: {  // §5.9: Γ(a,z) = G^{2,0}_{1,2}(z | 1 ; 0,a)
            if (call->args.size() != 2U) break;
            return make_meijerg(ctx, 2U, 0U, {int_lit(arena, 1)},
                {int_lit(arena, 0), call->args[0]}, call->args[1]);
        }

        case BuiltinOp::GammaIncompleteLower: {  // §5.9: γ(a,z) = G^{1,1}_{1,2}(z | 1 ; a,0)
            if (call->args.size() != 2U) break;
            return make_meijerg(ctx, 1U, 1U, {int_lit(arena, 1)},
                {call->args[0], int_lit(arena, 0)}, call->args[1]);
        }

        case BuiltinOp::BesselJ: {  // §5.8: G^{1,0}_{0,2}(w^2/4 | ; nu/2 ; -nu/2)
            if (call->args.size() != 2U) break;
            ExprPtr nu = call->args[0];
            auto half_nu = simp(ctx, arena.make<Product>(std::vector<ExprPtr>{
                rat_lit(arena, 1, 2), nu}));
            if (half_nu.is_error()) return half_nu;
            auto neg_half_nu = simp(ctx, negate(arena, half_nu.value()));
            if (neg_half_nu.is_error()) return neg_half_nu;
            auto arg = square_over_four(ctx, call->args[1]);
            if (arg.is_error()) return arg;
            return make_meijerg(ctx, 1U, 0U, {},
                {half_nu.value(), neg_half_nu.value()}, arg.value());
        }

        case BuiltinOp::Hypergeometric0F1:  // §3.1 general bridge
            if (call->args.size() != 2U) break;
            return pfq_to_meijerg(ctx, {}, {call->args[0]}, call->args[1]);

        case BuiltinOp::Hypergeometric1F1:
            if (call->args.size() != 3U) break;
            return pfq_to_meijerg(ctx, {call->args[0]}, {call->args[1]},
                                  call->args[2]);

        case BuiltinOp::Hypergeometric2F1:
            if (call->args.size() != 4U) break;
            return pfq_to_meijerg(ctx, {call->args[0], call->args[1]},
                                  {call->args[2]}, call->args[3]);

        default:
            break;
    }
    return fail<ExprPtr>(unsupported(
        "FuncCall " + std::string(builtin_op_name(call->func_id))));
}

}  // namespace cas::symbolic
