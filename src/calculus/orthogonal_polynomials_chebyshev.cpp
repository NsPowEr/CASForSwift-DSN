#include "orthogonal_polynomials_internal.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/error.hpp"
#include "cas/extended_real.hpp"

#include <optional>
#include <vector>

namespace cas::calculus {
namespace {

// Does `expr` simplify to (1 − x²) where x is `var`?
[[nodiscard]] bool match_one_minus_xsq(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    if (!expr) return false;
    AstArena& arena = ctx.arena();
    ExprPtr x = arena.make<Symbol>(var);
    ExprPtr x_sq = arena.make<Binary>(BinaryOp::Pow, x, arena.make<IntegerLit>(BigInt(2)));
    ExprPtr expected = arena.make<Binary>(BinaryOp::Sub, arena.make<IntegerLit>(BigInt(1)), x_sq);
    ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, expr, expected);
    auto s = ctx.simplify(delta);
    if (s.is_error()) return false;
    const auto* lit = expr_cast<IntegerLit>(s.value());
    return lit != nullptr && lit->value.is_zero();
}

// Match Sqrt(1 − x²)  ↔  the U weight.
[[nodiscard]] bool match_sqrt_one_minus_xsq(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    const auto* call = expr_cast<FuncCall>(expr);
    if (!call || call->func_id != BuiltinOp::Sqrt || call->args.size() != 1U) return false;
    return match_one_minus_xsq(call->args.front(), var, ctx);
}

// Match 1 / Sqrt(1 − x²)  ↔  the T weight, in any of:
//   Binary(Div, 1, Sqrt(…))
//   Binary(Pow, Sqrt(…), -1)
//   Binary(Pow, (1 − x²), -1/2)
[[nodiscard]] bool match_inv_sqrt_one_minus_xsq(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    if (const auto* div = expr_cast<Binary>(expr); div && div->op == BinaryOp::Div) {
        const auto* one = expr_cast<IntegerLit>(div->left);
        if (one && one->value == BigInt(1)) {
            return match_sqrt_one_minus_xsq(div->right, var, ctx);
        }
    }
    if (const auto* pw = expr_cast<Binary>(expr); pw && pw->op == BinaryOp::Pow) {
        // Sqrt(…)^(-1)
        if (match_sqrt_one_minus_xsq(pw->left, var, ctx)) {
            const auto* il = expr_cast<IntegerLit>(pw->right);
            if (il && il->value == BigInt(-1)) return true;
            if (const auto* un = expr_cast<Unary>(pw->right); un && un->op == UnaryOp::Neg) {
                const auto* il2 = expr_cast<IntegerLit>(un->operand);
                if (il2 && il2->value == BigInt(1)) return true;
            }
        }
        // (1 − x²)^(-1/2)
        if (match_one_minus_xsq(pw->left, var, ctx)) {
            if (const auto* rl = expr_cast<RationalLit>(pw->right);
                rl && rl->numerator == BigInt(-1) && rl->denominator == BigInt(2)) return true;
        }
    }
    return false;
}

}  // namespace

Result<std::optional<ExprPtr>> pattern_chebyshev_t_orthogonality(const DefiniteContext& dc) {
    // ∫_{-1}^{1} T_m(x)·T_n(x) / √(1−x²) dx
    //   = π        if m = n = 0
    //   = π/2      if m = n > 0
    //   = 0        otherwise
    if (!is_literal_neg_one(dc.lower)) return ok(std::optional<ExprPtr>{});
    if (!is_literal_rational(dc.upper, 1, 1)) return ok(std::optional<ExprPtr>{});

    // Two structural shapes for the weight 1/√(1−x²) :
    //   (A) integrand = Binary(Div, numerator, √(1−x²))
    //   (B) explicit reciprocal factor inside a Mul/Product tree
    ExprPtr numerator_form = dc.integrand;
    bool weight_found = false;
    if (const auto* div = expr_cast<Binary>(dc.integrand); div && div->op == BinaryOp::Div) {
        if (match_sqrt_one_minus_xsq(div->right, dc.var, dc.ctx)) {
            numerator_form = div->left;
            weight_found = true;
        }
    }
    std::vector<ExprPtr> factors;
    flatten_mul_factors(numerator_form, factors);
    std::vector<ExprPtr> non_weight;
    for (ExprPtr f : factors) {
        if (!weight_found && match_inv_sqrt_one_minus_xsq(f, dc.var, dc.ctx)) {
            weight_found = true;
            continue;
        }
        non_weight.push_back(f);
    }
    if (!weight_found) return ok(std::optional<ExprPtr>{});

    AstArena& arena = dc.ctx.arena();
    ExprPtr residual = arena.make<IntegerLit>(BigInt(1));
    if (!non_weight.empty()) {
        if (non_weight.size() == 1U) residual = non_weight.front();
        else residual = arena.make<Product>(std::move(non_weight));
    }
    auto match = match_two_poly_product(residual, BuiltinOp::ChebyshevT, dc.var);
    if (!match.has_value()) return ok(std::optional<ExprPtr>{});

    ExprPtr base;
    if (match->m == match->n) {
        ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
        if (match->n.is_zero()) {
            base = pi;
        } else {
            base = arena.make<Binary>(BinaryOp::Div, pi, arena.make<IntegerLit>(BigInt(2)));
        }
    } else {
        base = arena.make<IntegerLit>(BigInt(0));
    }
    return apply_constant_factors(base, match->other_factors, dc.var, dc.ctx);
}

Result<std::optional<ExprPtr>> pattern_chebyshev_u_orthogonality(const DefiniteContext& dc) {
    // ∫_{-1}^{1} U_m(x)·U_n(x) · √(1−x²) dx = (π/2) · δ_{mn}.
    if (!is_literal_neg_one(dc.lower)) return ok(std::optional<ExprPtr>{});
    if (!is_literal_rational(dc.upper, 1, 1)) return ok(std::optional<ExprPtr>{});

    std::vector<ExprPtr> factors;
    flatten_mul_factors(dc.integrand, factors);
    std::vector<ExprPtr> non_weight;
    bool weight_found = false;
    for (ExprPtr f : factors) {
        if (!weight_found && match_sqrt_one_minus_xsq(f, dc.var, dc.ctx)) {
            weight_found = true;
            continue;
        }
        non_weight.push_back(f);
    }
    if (!weight_found) return ok(std::optional<ExprPtr>{});

    AstArena& arena = dc.ctx.arena();
    ExprPtr residual = arena.make<IntegerLit>(BigInt(1));
    if (!non_weight.empty()) {
        if (non_weight.size() == 1U) residual = non_weight.front();
        else residual = arena.make<Product>(std::move(non_weight));
    }
    auto match = match_two_poly_product(residual, BuiltinOp::ChebyshevU, dc.var);
    if (!match.has_value()) return ok(std::optional<ExprPtr>{});

    ExprPtr base;
    if (match->m == match->n) {
        ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
        base = arena.make<Binary>(BinaryOp::Div, pi, arena.make<IntegerLit>(BigInt(2)));
    } else {
        base = arena.make<IntegerLit>(BigInt(0));
    }
    return apply_constant_factors(base, match->other_factors, dc.var, dc.ctx);
}

}  // namespace cas::calculus
