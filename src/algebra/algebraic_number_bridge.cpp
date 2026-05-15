#include "cas/algebraic_number_bridge.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"

#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace cas {
namespace algebra {

namespace {

[[nodiscard]] Result<std::size_t> bigint_to_size_t_nonneg(const BigInt& value) {
    if (value.is_negative()) {
        return fail<std::size_t>(make_error(
            CASErrorKind::InvalidArgument,
            "AlgebraicNumber bridge: negative integer exponent not representable as size_t"));
    }
    if (value.bit_length() > 63U) {
        return fail<std::size_t>(make_error(
            CASErrorKind::Unimplemented,
            "AlgebraicNumber bridge: integer exponent exceeds 64 bits"));
    }
    return ok(static_cast<std::size_t>(value.to_u64()));
}

[[nodiscard]] AlgebraicNumber make_zero(const AlgebraicNumber::CoeffVec& min_poly) {
    return AlgebraicNumber(AlgebraicNumber::CoeffVec{Rational(BigInt(0))}, min_poly);
}

[[nodiscard]] AlgebraicNumber make_one(const AlgebraicNumber::CoeffVec& min_poly) {
    return AlgebraicNumber(AlgebraicNumber::CoeffVec{Rational(BigInt(1))}, min_poly);
}

[[nodiscard]] AlgebraicNumber make_alpha(const AlgebraicNumber::CoeffVec& min_poly) {
    // value(x) = x  →  coefficients [0, 1] (ascending degree).
    return AlgebraicNumber(
        AlgebraicNumber::CoeffVec{Rational(BigInt(0)), Rational(BigInt(1))},
        min_poly);
}

[[nodiscard]] AlgebraicNumber from_rational(
    const Rational& r,
    const AlgebraicNumber::CoeffVec& min_poly) {
    return AlgebraicNumber(AlgebraicNumber::CoeffVec{r}, min_poly);
}

// Forward declaration: recursive worker.
[[nodiscard]] Result<std::optional<AlgebraicNumber>> express_recursive(
    ExprPtr e,
    ExprPtr alpha_expr,
    const AlgebraicNumber::CoeffVec& min_poly,
    symbolic::CASContext& ctx,
    unsigned int depth);

[[nodiscard]] Result<std::optional<AlgebraicNumber>> express_integer_power(
    const AlgebraicNumber& base,
    const BigInt& exponent,
    const AlgebraicNumber::CoeffVec& min_poly) {
    if (exponent.is_zero()) {
        return ok(std::optional<AlgebraicNumber>(make_one(min_poly)));
    }
    const bool negative = exponent.is_negative();
    BigInt abs_exp = exponent.abs();
    auto n_res = bigint_to_size_t_nonneg(abs_exp);
    if (n_res.is_error()) return fail<std::optional<AlgebraicNumber>>(n_res.error());
    const std::size_t n = n_res.value();

    auto pow_res = base.pow(n);
    if (pow_res.is_error()) return fail<std::optional<AlgebraicNumber>>(pow_res.error());
    AlgebraicNumber result = pow_res.value();
    if (negative) {
        if (result.is_zero()) {
            return ok(std::optional<AlgebraicNumber>{});  // 0^(-k): not in Q(alpha); fail gracefully
        }
        auto inv = result.inverse();
        if (inv.is_error()) return fail<std::optional<AlgebraicNumber>>(inv.error());
        result = inv.value();
    }
    return ok(std::optional<AlgebraicNumber>(std::move(result)));
}

// Recursion depth limit prevents accidental runaway on adversarial trees.
// Real Q(alpha) expressions are shallow; 256 is more than ample.
constexpr unsigned int kMaxBridgeDepth = 256U;

[[nodiscard]] Result<std::optional<AlgebraicNumber>> express_recursive(
    ExprPtr e,
    ExprPtr alpha_expr,
    const AlgebraicNumber::CoeffVec& min_poly,
    symbolic::CASContext& ctx,
    unsigned int depth) {
    if (depth > kMaxBridgeDepth) {
        return ok(std::optional<AlgebraicNumber>{});
    }
    if (!e) {
        return ok(std::optional<AlgebraicNumber>{});
    }

    // Structural match against the supplied generator expression.
    if (alpha_expr && structural_equal(e, alpha_expr)) {
        return ok(std::optional<AlgebraicNumber>(make_alpha(min_poly)));
    }

    if (const auto* lit = expr_cast<IntegerLit>(e)) {
        return ok(std::optional<AlgebraicNumber>(from_rational(Rational(lit->value), min_poly)));
    }
    if (const auto* rlit = expr_cast<RationalLit>(e)) {
        return ok(std::optional<AlgebraicNumber>(
            from_rational(Rational(rlit->numerator, rlit->denominator), min_poly)));
    }
    if (const auto* un = expr_cast<Unary>(e)) {
        if (un->op == UnaryOp::Neg) {
            auto inner = express_recursive(un->operand, alpha_expr, min_poly, ctx, depth + 1U);
            if (inner.is_error()) return inner;
            if (!inner.value().has_value()) return ok(std::optional<AlgebraicNumber>{});
            return ok(std::optional<AlgebraicNumber>(-inner.value().value()));
        }
        return ok(std::optional<AlgebraicNumber>{});  // Factorial etc. not in Q(alpha).
    }
    if (const auto* bin = expr_cast<Binary>(e)) {
        if (bin->op == BinaryOp::Pow) {
            auto base = express_recursive(bin->left, alpha_expr, min_poly, ctx, depth + 1U);
            if (base.is_error()) return base;
            if (!base.value().has_value()) return ok(std::optional<AlgebraicNumber>{});
            const auto* exp_lit = expr_cast<IntegerLit>(bin->right);
            if (!exp_lit) {
                return ok(std::optional<AlgebraicNumber>{});
            }
            return express_integer_power(base.value().value(), exp_lit->value, min_poly);
        }
        auto lhs = express_recursive(bin->left, alpha_expr, min_poly, ctx, depth + 1U);
        if (lhs.is_error()) return lhs;
        if (!lhs.value().has_value()) return ok(std::optional<AlgebraicNumber>{});
        auto rhs = express_recursive(bin->right, alpha_expr, min_poly, ctx, depth + 1U);
        if (rhs.is_error()) return rhs;
        if (!rhs.value().has_value()) return ok(std::optional<AlgebraicNumber>{});

        const AlgebraicNumber& l = lhs.value().value();
        const AlgebraicNumber& r = rhs.value().value();
        switch (bin->op) {
            case BinaryOp::Add:
                return ok(std::optional<AlgebraicNumber>(l + r));
            case BinaryOp::Sub:
                return ok(std::optional<AlgebraicNumber>(l - r));
            case BinaryOp::Mul:
                return ok(std::optional<AlgebraicNumber>(l * r));
            case BinaryOp::Div: {
                if (r.is_zero()) {
                    return ok(std::optional<AlgebraicNumber>{});  // 1/0 outside Q(alpha).
                }
                auto div_res = l.div(r);
                if (div_res.is_error()) {
                    // Non-invertible: minimal polynomial may be reducible; treat as "not in Q(alpha)".
                    return ok(std::optional<AlgebraicNumber>{});
                }
                return ok(std::optional<AlgebraicNumber>(div_res.value()));
            }
            default:
                return ok(std::optional<AlgebraicNumber>{});  // Mod / comparison: not field ops.
        }
    }
    if (const auto* sum = expr_cast<Sum>(e)) {
        AlgebraicNumber acc = make_zero(min_poly);
        for (ExprPtr term : sum->terms) {
            auto t = express_recursive(term, alpha_expr, min_poly, ctx, depth + 1U);
            if (t.is_error()) return t;
            if (!t.value().has_value()) return ok(std::optional<AlgebraicNumber>{});
            acc = acc + t.value().value();
        }
        return ok(std::optional<AlgebraicNumber>(std::move(acc)));
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        AlgebraicNumber acc = make_one(min_poly);
        for (ExprPtr factor : prod->factors) {
            auto f = express_recursive(factor, alpha_expr, min_poly, ctx, depth + 1U);
            if (f.is_error()) return f;
            if (!f.value().has_value()) return ok(std::optional<AlgebraicNumber>{});
            acc = acc * f.value().value();
        }
        return ok(std::optional<AlgebraicNumber>(std::move(acc)));
    }

    // Anything else (Symbol, Constant, FuncCall, Matrix, Integral, ..., other RootOf):
    // not expressible in Q(alpha) via this bridge.
    return ok(std::optional<AlgebraicNumber>{});
}

[[nodiscard]] ExprPtr rational_to_expr(AstArena& arena, const Rational& r) {
    if (r.denominator() == BigInt(1)) {
        return arena.make<IntegerLit>(r.numerator());
    }
    return arena.make<RationalLit>(r.numerator(), r.denominator());
}

}  // namespace

[[nodiscard]] Result<AlgebraicNumber::CoeffVec> rootof_min_poly(
    const RootOf& root,
    symbolic::CASContext& ctx) {
    auto parsed = parse_polynomial(root.polynomial, root.variable, ctx);
    if (parsed.is_error()) return fail<AlgebraicNumber::CoeffVec>(parsed.error());

    auto rational_poly = poly_to_rational_poly(parsed.value());
    if (rational_poly.is_error()) {
        return fail<AlgebraicNumber::CoeffVec>(make_error(
            CASErrorKind::Unimplemented,
            "RootOf bridge: minimal polynomial requires rational coefficients"));
    }

    RatPoly min_poly = rational_poly.value();
    normalize_rational_coefficients(min_poly);
    if (min_poly.empty() || min_poly.is_zero() || min_poly.degree() == 0U) {
        return fail<AlgebraicNumber::CoeffVec>(make_error(
            CASErrorKind::InvalidArgument,
            "RootOf bridge: minimal polynomial must have positive degree"));
    }
    const Rational leading = min_poly.leading_coeff();
    if (leading.numerator().is_zero()) {
        return fail<AlgebraicNumber::CoeffVec>(make_error(
            CASErrorKind::InvalidArgument,
            "RootOf bridge: leading coefficient must be non-zero"));
    }
    if (!(leading == Rational(BigInt(1)))) {
        for (auto& coefficient : min_poly.coefficients()) {
            coefficient = coefficient / leading;
        }
        normalize_rational_coefficients(min_poly);
    }
    return ok(min_poly.coefficients());
}

[[nodiscard]] Result<AlgebraicNumber> alpha_from_rootof(
    const RootOf& root,
    symbolic::CASContext& ctx) {
    auto mp = rootof_min_poly(root, ctx);
    if (mp.is_error()) return fail<AlgebraicNumber>(mp.error());
    return ok(make_alpha(mp.value()));
}

[[nodiscard]] Result<std::optional<AlgebraicNumber>> try_express_in_q_alpha(
    ExprPtr e,
    ExprPtr alpha_expr,
    const AlgebraicNumber::CoeffVec& min_poly,
    symbolic::CASContext& ctx) {
    if (min_poly.empty() || min_poly.size() < 2U) {
        return fail<std::optional<AlgebraicNumber>>(make_error(
            CASErrorKind::InvalidArgument,
            "AlgebraicNumber bridge: minimal polynomial must have degree >= 1"));
    }
    return express_recursive(e, alpha_expr, min_poly, ctx, 0U);
}

[[nodiscard]] ExprPtr algebraic_number_to_expr_raw(
    const AlgebraicNumber& value,
    ExprPtr alpha_expr,
    AstArena& arena) {
    const auto& coeffs = value.value();
    if (coeffs.empty()) {
        return arena.make<IntegerLit>(BigInt(0));
    }
    std::vector<ExprPtr> terms;
    terms.reserve(coeffs.size());
    for (std::size_t k = 0U; k < coeffs.size(); ++k) {
        const Rational& c = coeffs[k];
        if (c.numerator().is_zero()) continue;
        ExprPtr coeff_expr = rational_to_expr(arena, c);
        ExprPtr alpha_power;
        if (k == 0U) {
            alpha_power = nullptr;
        } else if (k == 1U) {
            alpha_power = alpha_expr;
        } else if (alpha_expr) {
            alpha_power = arena.make<Binary>(
                BinaryOp::Pow,
                alpha_expr,
                arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(k))));
        } else {
            // No generator supplied but non-constant: render as integer 0
            // since we cannot represent the higher-degree term.  Caller is
            // expected to supply alpha_expr whenever value() has size > 1.
            continue;
        }
        ExprPtr term;
        if (!alpha_power) {
            term = coeff_expr;
        } else if (c == Rational(BigInt(1))) {
            term = alpha_power;
        } else if (c == Rational(BigInt(-1))) {
            term = arena.make<Unary>(UnaryOp::Neg, alpha_power);
        } else {
            term = arena.make<Binary>(BinaryOp::Mul, coeff_expr, alpha_power);
        }
        terms.push_back(term);
    }
    if (terms.empty()) {
        return arena.make<IntegerLit>(BigInt(0));
    }
    if (terms.size() == 1U) return terms.front();
    return arena.make<Sum>(std::move(terms));
}

[[nodiscard]] Result<ExprPtr> algebraic_number_to_expr(
    const AlgebraicNumber& value,
    ExprPtr alpha_expr,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    if (!alpha_expr && value.value().size() > 1U) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::InvalidArgument,
            "AlgebraicNumber bridge: non-constant element requires alpha_expr"));
    }
    ExprPtr raw = algebraic_number_to_expr_raw(value, alpha_expr, arena);
    return ctx.simplify(raw);
}

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
