#include "cas/algebraic_number_bridge.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"
#include "cas/formatter.hpp"

#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <cstddef>
#include <cstdio>
#include <optional>
#include <utility>
#include <vector>

namespace cas {
namespace algebra {

namespace {

// Walk an expression tree, collecting structurally distinct RootOf nodes.
// Stops accumulating once two distinct RootOfs are encountered (caller only
// needs to know whether there is exactly one).
void collect_distinct_rootofs(ExprPtr expr, std::vector<ExprPtr>& out) {
    if (!expr) return;
    if (out.size() >= 2U) return;
    if (expr_is<RootOf>(expr)) {
        for (ExprPtr existing : out) {
            if (structural_equal(existing, expr)) return;
        }
        out.push_back(expr);
        // Continue: the RootOf's polynomial subexpression itself could contain
        // other RootOfs in pathological inputs.  Walk it anyway for correctness.
        const auto& root = expr_ref<RootOf>(expr);
        collect_distinct_rootofs(root.polynomial, out);
        return;
    }

    visit_expr(expr, [&](const auto& node) {
        using Node = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<Node, Unary>) {
            collect_distinct_rootofs(node.operand, out);
        } else if constexpr (std::is_same_v<Node, Binary>) {
            collect_distinct_rootofs(node.left, out);
            collect_distinct_rootofs(node.right, out);
        } else if constexpr (std::is_same_v<Node, FuncCall>) {
            for (ExprPtr arg : node.args) collect_distinct_rootofs(arg, out);
        } else if constexpr (std::is_same_v<Node, Sum>) {
            for (ExprPtr t : node.terms) collect_distinct_rootofs(t, out);
        } else if constexpr (std::is_same_v<Node, Product>) {
            for (ExprPtr f : node.factors) collect_distinct_rootofs(f, out);
        } else if constexpr (std::is_same_v<Node, Integral>) {
            collect_distinct_rootofs(node.integrand, out);
            if (node.lower.has_value()) collect_distinct_rootofs(*node.lower, out);
            if (node.upper.has_value()) collect_distinct_rootofs(*node.upper, out);
        } else if constexpr (std::is_same_v<Node, Derivative>) {
            collect_distinct_rootofs(node.expression, out);
        } else if constexpr (std::is_same_v<Node, Limit>) {
            collect_distinct_rootofs(node.expression, out);
            collect_distinct_rootofs(node.point, out);
        } else if constexpr (std::is_same_v<Node, Matrix>) {
            for (ExprPtr e : node.elements) collect_distinct_rootofs(e, out);
        }
        // Leaf kinds and unhandled aggregates: no descent.
    });
}

}  // namespace

[[nodiscard]] Result<ExprPtr> try_reduce_in_q_alpha(
    ExprPtr expr,
    symbolic::CASContext& ctx) {
    if (!expr) return ok(expr);

    std::vector<ExprPtr> roots;
    collect_distinct_rootofs(expr, roots);
    if (roots.size() != 1U) return ok(expr);

    ExprPtr alpha_expr = roots.front();
    const auto& root_node = expr_ref<RootOf>(alpha_expr);

    auto mp_res = rootof_min_poly(root_node, ctx);
    if (mp_res.is_error()) return ok(expr);  // non-rational min_poly: no-op

    auto an_res = try_express_in_q_alpha(expr, alpha_expr, mp_res.value(), ctx);
    if (an_res.is_error()) return ok(expr);
    if (!an_res.value().has_value()) return ok(expr);

    ExprPtr rendered = algebraic_number_to_expr_raw(an_res.value().value(), alpha_expr, ctx.arena());
    if (structural_equal(rendered, expr)) return ok(expr);
    return ok(rendered);
}

[[nodiscard]] static bool is_algebraic_generator(ExprPtr expr) {
    if (expr_is<RootOf>(expr)) return true;
    if (const auto* call = expr_cast<FuncCall>(expr); call && call->func_id == BuiltinOp::Sqrt && call->args.size() == 1U) {
        if (const auto* il = expr_cast<IntegerLit>(call->args[0])) return il->value.is_positive();
        if (const auto* rl = expr_cast<RationalLit>(call->args[0])) return rl->numerator.is_positive();
    }
    if (const auto* bin = expr_cast<Binary>(expr); bin && bin->op == BinaryOp::Pow) {
        bool base_ok = false;
        if (const auto* il = expr_cast<IntegerLit>(bin->left)) base_ok = il->value.is_positive();
        else if (const auto* rl = expr_cast<RationalLit>(bin->left)) base_ok = rl->numerator.is_positive();
        if (base_ok) {
            if (const auto* rl = expr_cast<RationalLit>(bin->right)) {
                return rl->numerator == BigInt(1) && rl->denominator > BigInt(1);
            }
        }
    }
    return false;
}

static void collect_algebraic_generators(ExprPtr expr, std::vector<ExprPtr>& out, symbolic::CASContext& ctx) {
    if (!expr) return;
    if (out.size() >= 2U) return;
    if (is_algebraic_generator(expr)) {
        for (ExprPtr existing : out) {
            if (same_generator_expr(existing, expr, ctx)) return;
        }
        out.push_back(expr);
        if (expr_is<RootOf>(expr)) {
            const auto& root = expr_ref<RootOf>(expr);
            collect_algebraic_generators(root.polynomial, out, ctx);
        }
        return;
    }

    visit_expr(expr, [&](const auto& node) {
        using Node = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<Node, Unary>) {
            collect_algebraic_generators(node.operand, out, ctx);
        } else if constexpr (std::is_same_v<Node, Binary>) {
            collect_algebraic_generators(node.left, out, ctx);
            collect_algebraic_generators(node.right, out, ctx);
        } else if constexpr (std::is_same_v<Node, FuncCall>) {
            for (ExprPtr arg : node.args) collect_algebraic_generators(arg, out, ctx);
        } else if constexpr (std::is_same_v<Node, Sum>) {
            for (ExprPtr t : node.terms) collect_algebraic_generators(t, out, ctx);
        } else if constexpr (std::is_same_v<Node, Product>) {
            for (ExprPtr f : node.factors) collect_algebraic_generators(f, out, ctx);
        } else if constexpr (std::is_same_v<Node, Integral>) {
            collect_algebraic_generators(node.integrand, out, ctx);
            if (node.lower.has_value()) collect_algebraic_generators(*node.lower, out, ctx);
            if (node.upper.has_value()) collect_algebraic_generators(*node.upper, out, ctx);
        } else if constexpr (std::is_same_v<Node, Derivative>) {
            collect_algebraic_generators(node.expression, out, ctx);
        } else if constexpr (std::is_same_v<Node, Limit>) {
            collect_algebraic_generators(node.expression, out, ctx);
            collect_algebraic_generators(node.point, out, ctx);
        } else if constexpr (std::is_same_v<Node, Matrix>) {
            for (ExprPtr e : node.elements) collect_algebraic_generators(e, out, ctx);
        }
    });
}

[[nodiscard]] static Result<std::vector<Rational>> generator_min_poly(ExprPtr expr, symbolic::CASContext& ctx) {
    if (const auto* root = expr_cast<RootOf>(expr)) {
        return rootof_min_poly(*root, ctx);
    }
    Rational c;
    BigInt n;
    if (const auto* call = expr_cast<FuncCall>(expr); call && call->func_id == BuiltinOp::Sqrt) {
        if (const auto* il = expr_cast<IntegerLit>(call->args[0])) c = Rational(il->value);
        else if (const auto* rl = expr_cast<RationalLit>(call->args[0])) c = Rational(rl->numerator, rl->denominator);
        n = BigInt(2);
    } else if (const auto* bin = expr_cast<Binary>(expr); bin && bin->op == BinaryOp::Pow) {
        if (const auto* il = expr_cast<IntegerLit>(bin->left)) c = Rational(il->value);
        else if (const auto* rl = expr_cast<RationalLit>(bin->left)) c = Rational(rl->numerator, rl->denominator);
        n = expr_cast<RationalLit>(bin->right)->denominator;
    }
    // minimal polynomial x^n - c -> coefficients are [-c, 0, ..., 0, 1]
    std::vector<Rational> coeffs(n.to_u64() + 1, Rational(0));
    coeffs[0] = -c;
    coeffs[n.to_u64()] = Rational(1);
    return ok(std::move(coeffs));
}

[[nodiscard]] Result<ExprPtr> simplify_in_q_alpha(
    ExprPtr expr,
    symbolic::CASContext& ctx) {
    auto first = ctx.simplify(expr);
    if (first.is_error()) return first;
    auto reduced = try_reduce_in_q_alpha(first.value(), ctx);
    if (reduced.is_error()) return reduced;
    if (reduced.value() == first.value()) return first;
    return ctx.simplify(reduced.value());
}

void register_algebraic_simplify_hook(symbolic::CASContext& ctx) {
    ctx.set_post_simplify_hook([](ExprPtr e, symbolic::CASContext& c) -> Result<ExprPtr> {
        std::vector<ExprPtr> gens;
        collect_algebraic_generators(e, gens, c);
        if (gens.size() == 1U) {
            ExprPtr alpha_expr = gens.front();
            if (!expr_is<RootOf>(alpha_expr)) {
                auto mp_res = generator_min_poly(alpha_expr, c);
                if (mp_res.is_ok()) {
                    auto an_res = try_express_in_q_alpha(e, alpha_expr, mp_res.value(), c);
                    if (an_res.is_ok() && an_res.value().has_value()) {
                        ExprPtr rendered = algebraic_number_to_expr_raw(an_res.value().value(), alpha_expr, c.arena());
                        if (!structural_equal(rendered, e)) {
                            return ok(rendered);
                        }
                    }
                }
            }
        }
        return try_reduce_in_q_alpha(e, c);
    });
}

[[nodiscard]] Result<ExprPtr> simplify_polynomial_in_x_over_q_alpha(
    ExprPtr expr,
    const Symbol& poly_var,
    symbolic::CASContext& ctx) {
    if (!expr) return ok(expr);

    // 1. Run global simplify + expand so the expression is a flat polynomial in poly_var.
    auto simplified = ctx.simplify(expr);
    if (simplified.is_error()) return simplified;
    auto expanded = expand(simplified.value(), ctx);
    if (expanded.is_error()) return expanded;

    // 2. Parse as polynomial in poly_var.  If parse fails, fall back to direct
    //    try_reduce_in_q_alpha on the whole expression.
    auto poly_res = parse_polynomial(expanded.value(), poly_var, ctx);
    if (poly_res.is_error()) {
        return try_reduce_in_q_alpha(expanded.value(), ctx);
    }
    const auto& poly = poly_res.value();
    if (poly.empty()) return ok(expanded.value());

    // 3. For each coefficient, reduce in Q(alpha).  If a coefficient is
    //    already a pure rational, try_reduce_in_q_alpha is a no-op.
    PolyExpr reduced_poly;
    reduced_poly.reserve(poly.size());
    for (std::size_t k = 0; k < poly.size(); ++k) {
        ExprPtr coeff = poly[k];
        if (!coeff) {
            reduced_poly.push_back(coeff);
            continue;
        }
        auto reduced_coeff = try_reduce_in_q_alpha(coeff, ctx);
        if (reduced_coeff.is_error()) {
            reduced_poly.push_back(coeff);
        } else {
            auto canon = ctx.simplify(reduced_coeff.value());
            reduced_poly.push_back(canon.is_ok() ? canon.value() : reduced_coeff.value());
        }
    }

    // 4. Rebuild polynomial expression.
    auto rebuilt = polynomial_to_expr(reduced_poly, poly_var, ctx);
    if (rebuilt.is_error()) return rebuilt;
    return ctx.simplify(rebuilt.value());
}

}  // namespace algebra
}  // namespace cas
