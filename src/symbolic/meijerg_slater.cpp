// A7 step 4 — Slater expansion (DLMF 16.17.2), the GENERAL inverse path
// from a Meijer G node to a sum of hypergeometric terms.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Meijer_G_Slater.md §3.2
// (formula numerically verified against mpmath, Appendice C).

#include "cas/error.hpp"
#include "cas/meijerg.hpp"
#include "cas/symbolic.hpp"
#include "simplify_impl.hpp"

#include <string>
#include <utility>
#include <vector>

namespace cas::symbolic {

using detail::is_zero_expr;

namespace {

[[nodiscard]] ExprPtr int_lit(AstArena& arena, long long v) {
    return arena.make<IntegerLit>(BigInt(v));
}

[[nodiscard]] ExprPtr gamma_of(AstArena& arena, ExprPtr arg) {
    return arena.make<FuncCall>(BuiltinOp::Gamma, std::vector<ExprPtr>{arg});
}

[[nodiscard]] ExprPtr inverse_of(AstArena& arena, ExprPtr e) {
    return arena.make<Binary>(BinaryOp::Pow, e, int_lit(arena, -1));
}

// 1 + x - y, engine-simplified.
[[nodiscard]] Result<ExprPtr> one_plus_diff(
    CASContext& ctx, ExprPtr x, ExprPtr y) {
    AstArena& arena = ctx.arena();
    return ctx.simplify(arena.make<Sum>(std::vector<ExprPtr>{
        int_lit(arena, 1), x, arena.make<Unary>(UnaryOp::Neg, y)}));
}

// b_j - b_k decidably an integer (confluent poles, §3.2 precondition).
[[nodiscard]] Result<bool> is_decidably_integer_diff(
    CASContext& ctx, ExprPtr x, ExprPtr y) {
    AstArena& arena = ctx.arena();
    auto diff = ctx.simplify(arena.make<Sum>(std::vector<ExprPtr>{
        x, arena.make<Unary>(UnaryOp::Neg, y)}));
    if (diff.is_error()) return fail<bool>(diff.error());
    return ok(expr_cast<IntegerLit>(diff.value()) != nullptr);
}

[[nodiscard]] CASError slater_unimplemented(const std::string& shape) {
    return make_unimplemented_error(
        UnimplementedInfo{
            .module = "symbolic",
            .function = "slater_expand",
            .input_shape = shape,
            .reason = error::reason_codes::GENERIC,
            .suggestion = "Slater expansion (DLMF 16.17.2) not applicable; "
                          "the G node stays as-is (Meijer_G_Slater.md §9.4)",
            .ticket = "A7"},
        "slater_expand: expansion not applicable");
}

// The k-th summand's pF(q-1) with upper (1+b_k-a_l, l=1..p), lower
// (1+b_k-b_l, l=1..q, l != k), argument w. Closed forms for the arities
// with no engine node: 0F0(;;w) = e^w, 1F0(a;;w) = (1-w)^{-a} (both exact;
// DLMF 16.2.1 degenerate cases). Otherwise a 0F1/1F1/2F1 FuncCall.
[[nodiscard]] Result<ExprPtr> build_pfq_term(
    CASContext& ctx, const std::vector<ExprPtr>& upper,
    const std::vector<ExprPtr>& lower, ExprPtr w) {
    AstArena& arena = ctx.arena();
    const std::size_t up = upper.size();
    const std::size_t low = lower.size();
    if (up == 0U && low == 0U) {
        return ok(arena.make<FuncCall>(BuiltinOp::Exp,
            std::vector<ExprPtr>{w}));
    }
    if (up == 1U && low == 0U) {
        ExprPtr one_minus_w = arena.make<Sum>(std::vector<ExprPtr>{
            int_lit(arena, 1), arena.make<Unary>(UnaryOp::Neg, w)});
        return ok(arena.make<Binary>(BinaryOp::Pow, one_minus_w,
            arena.make<Unary>(UnaryOp::Neg, upper[0])));
    }
    if (up == 0U && low == 1U) {
        return ok(arena.make<FuncCall>(BuiltinOp::Hypergeometric0F1,
            std::vector<ExprPtr>{lower[0], w}));
    }
    if (up == 1U && low == 1U) {
        return ok(arena.make<FuncCall>(BuiltinOp::Hypergeometric1F1,
            std::vector<ExprPtr>{upper[0], lower[0], w}));
    }
    if (up == 2U && low == 1U) {
        return ok(arena.make<FuncCall>(BuiltinOp::Hypergeometric2F1,
            std::vector<ExprPtr>{upper[0], upper[1], lower[0], w}));
    }
    return fail<ExprPtr>(slater_unimplemented(
        "pF(q-1) arity (" + std::to_string(up) + "," + std::to_string(low)
        + ") has no engine node"));
}

}  // namespace

Result<ExprPtr> slater_expand(CASContext& ctx, const FuncCall& g) {
    auto view_res = view_meijerg(g);
    if (view_res.is_error()) return fail<ExprPtr>(view_res.error());
    const MeijerGView& v = view_res.value();
    AstArena& arena = ctx.arena();

    // §3.2 preconditions.
    if (v.p > v.q) {
        return fail<ExprPtr>(slater_unimplemented(
            "p > q (use the argument-inversion identity 16.19.1 first)"));
    }
    if (v.m == 0U) {
        return fail<ExprPtr>(slater_unimplemented(
            "m = 0 (no residue series over the b-group)"));
    }
    for (std::size_t j = 0; j < v.m; ++j) {
        for (std::size_t k = j + 1U; k < v.m; ++k) {
            auto confluent = is_decidably_integer_diff(ctx, v.b[j], v.b[k]);
            if (confluent.is_error()) return fail<ExprPtr>(confluent.error());
            if (confluent.value()) {
                return fail<ExprPtr>(slater_unimplemented(
                    "b_j - b_k integer within the m-group (confluent poles "
                    "need logarithmic terms)"));
            }
        }
    }

    // Argument (-1)^{p-m-n} z: parity of p-m-n == parity of p+m+n.
    const bool negate_arg = ((v.p + v.m + v.n) % 2U) == 1U;
    ExprPtr w = v.z;
    if (negate_arg) w = arena.make<Unary>(UnaryOp::Neg, v.z);
    auto w_simplified = ctx.simplify(w);
    if (w_simplified.is_error()) return w_simplified;
    w = w_simplified.value();

    std::vector<ExprPtr> summands;
    summands.reserve(v.m);
    for (std::size_t k = 0; k < v.m; ++k) {
        std::vector<ExprPtr> factors;

        // A_k numerator: prod_{l != k, l <= m} Gamma(b_l - b_k)
        //                * prod_{l <= n} Gamma(1 + b_k - a_l).
        for (std::size_t l = 0; l < v.m; ++l) {
            if (l == k) continue;
            auto diff = ctx.simplify(arena.make<Sum>(std::vector<ExprPtr>{
                v.b[l], arena.make<Unary>(UnaryOp::Neg, v.b[k])}));
            if (diff.is_error()) return diff;
            factors.push_back(gamma_of(arena, diff.value()));
        }
        for (std::size_t l = 0; l < v.n; ++l) {
            auto arg = one_plus_diff(ctx, v.b[k], v.a[l]);
            if (arg.is_error()) return arg;
            factors.push_back(gamma_of(arena, arg.value()));
        }
        // A_k denominator: prod_{l = m+1..q} Gamma(1 + b_k - b_l)
        //                  * prod_{l = n+1..p} Gamma(a_l + 1 - b_k)  — as
        // Pow(Gamma, -1) factors (no Div in canonical products).
        for (std::size_t l = v.m; l < v.q; ++l) {
            auto arg = one_plus_diff(ctx, v.b[k], v.b[l]);
            if (arg.is_error()) return arg;
            factors.push_back(inverse_of(arena, gamma_of(arena, arg.value())));
        }
        for (std::size_t l = v.n; l < v.p; ++l) {
            auto arg = one_plus_diff(ctx, v.a[l], v.b[k]);
            if (arg.is_error()) return arg;
            factors.push_back(inverse_of(arena, gamma_of(arena, arg.value())));
        }
        // z^{b_k} prefactor (omitted for b_k = 0: z^0 = 1 and adding the
        // factor would only inject a spurious NonZero(z) side condition).
        if (!is_zero_expr(v.b[k])) {
            factors.push_back(arena.make<Binary>(BinaryOp::Pow, v.z, v.b[k]));
        }

        // pF(q-1) parameter lists for this k.
        std::vector<ExprPtr> upper;
        upper.reserve(v.p);
        for (std::size_t l = 0; l < v.p; ++l) {
            auto arg = one_plus_diff(ctx, v.b[k], v.a[l]);
            if (arg.is_error()) return arg;
            upper.push_back(arg.value());
        }
        std::vector<ExprPtr> lower;
        lower.reserve(v.q > 0U ? v.q - 1U : 0U);
        for (std::size_t l = 0; l < v.q; ++l) {
            if (l == k) continue;  // the omitted 1 + b_k - b_k factor (§3.2)
            auto arg = one_plus_diff(ctx, v.b[k], v.b[l]);
            if (arg.is_error()) return arg;
            lower.push_back(arg.value());
        }
        auto pfq = build_pfq_term(ctx, upper, lower, w);
        if (pfq.is_error()) return pfq;
        factors.push_back(pfq.value());

        summands.push_back(factors.size() == 1U
            ? factors.front()
            : arena.make<Product>(std::move(factors)));
    }

    ExprPtr result = summands.size() == 1U
        ? summands.front()
        : arena.make<Sum>(std::move(summands));
    return ctx.simplify(result);
}

}  // namespace cas::symbolic
