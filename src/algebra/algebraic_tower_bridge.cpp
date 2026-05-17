#include "cas/algebraic_tower_bridge.hpp"

#include "cas/algebraic_number_bridge.hpp"
#include "cas/error.hpp"

#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace cas {
namespace algebra {

namespace {

class RootOfExplicitDegreeGuard {
public:
    RootOfExplicitDegreeGuard(symbolic::CASContext& ctx, std::size_t temp)
        : ctx_(ctx), saved_(ctx.max_rootof_explicit_degree()) {
        ctx_.set_max_rootof_explicit_degree(temp);
    }

    ~RootOfExplicitDegreeGuard() {
        ctx_.set_max_rootof_explicit_degree(saved_);
    }

private:
    symbolic::CASContext& ctx_;
    std::size_t saved_;
};

[[nodiscard]] Result<ExprPtr> clone_expr_raw(ExprPtr expr, symbolic::CASContext& ctx) {
    return symbolic::materialize_expr(expr, ctx.arena());
}

[[nodiscard]] Result<ExprPtr> canonicalize_root_expr(ExprPtr expr, symbolic::CASContext& ctx) {
    RootOfExplicitDegreeGuard rootof_guard(ctx, 1U);
    return ctx.simplify(expr);
}

[[nodiscard]] Result<PolyExpr> make_constant_poly_raw(ExprPtr coefficient, symbolic::CASContext& ctx) {
    auto simplified = poly_simplify_expr(coefficient, ctx);
    if (simplified.is_error()) return fail<PolyExpr>(simplified.error());
    PolyExpr poly;
    if (!poly_is_zero_expr(simplified.value())) poly.push_back(simplified.value());
    return ok(std::move(poly));
}

[[nodiscard]] Result<PolyExpr> poly_pow_raw(PolyExpr base, std::size_t exponent, symbolic::CASContext& ctx) {
    auto one = make_constant_poly_raw(poly_make_integer(ctx.arena(), 1), ctx);
    if (one.is_error()) return fail<PolyExpr>(one.error());
    PolyExpr result = one.value();
    while (exponent > 0U) {
        if ((exponent & 1U) != 0U) {
            auto multiplied = poly_multiply(result, base, ctx);
            if (multiplied.is_error()) return fail<PolyExpr>(multiplied.error());
            result = std::move(multiplied.value());
        }
        exponent >>= 1U;
        if (exponent == 0U) break;
        auto squared = poly_multiply(base, base, ctx);
        if (squared.is_error()) return fail<PolyExpr>(squared.error());
        base = std::move(squared.value());
    }
    return ok(std::move(result));
}

[[nodiscard]] Result<PolyExpr> parse_polynomial_raw(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<PolyExpr>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "Tower bridge: null polynomial expression",
        });
    }
    if (poly_contains_decimal_literal(expr)) {
        return fail<PolyExpr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Tower bridge: decimal literals are not supported in exact polynomial parsing",
        });
    }
    if (!poly_depends_on(expr, var.name)) {
        auto cloned = clone_expr_raw(expr, ctx);
        if (cloned.is_error()) return fail<PolyExpr>(cloned.error());
        return make_constant_poly_raw(cloned.value(), ctx);
    }
    if (const auto* symbol = expr_cast<Symbol>(expr)) {
        if (symbol->name == var.name) {
            return ok(poly_make_monomial(poly_make_integer(ctx.arena(), 1), 1U));
        }
    }
    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op != UnaryOp::Neg) {
            return fail<PolyExpr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Tower bridge: unsupported unary operator in polynomial parsing",
            });
        }
        auto operand = parse_polynomial_raw(unary->operand, var, ctx);
        if (operand.is_error()) return fail<PolyExpr>(operand.error());
        return poly_negate(operand.value(), ctx);
    }
    if (const auto* binary = expr_cast<Binary>(expr)) {
        switch (binary->op) {
            case BinaryOp::Add: {
                auto lhs = parse_polynomial_raw(binary->left, var, ctx);
                if (lhs.is_error()) return fail<PolyExpr>(lhs.error());
                auto rhs = parse_polynomial_raw(binary->right, var, ctx);
                if (rhs.is_error()) return fail<PolyExpr>(rhs.error());
                return poly_add(lhs.value(), rhs.value(), ctx);
            }
            case BinaryOp::Sub: {
                auto lhs = parse_polynomial_raw(binary->left, var, ctx);
                if (lhs.is_error()) return fail<PolyExpr>(lhs.error());
                auto rhs = parse_polynomial_raw(binary->right, var, ctx);
                if (rhs.is_error()) return fail<PolyExpr>(rhs.error());
                return poly_subtract(lhs.value(), rhs.value(), ctx);
            }
            case BinaryOp::Mul: {
                auto lhs = parse_polynomial_raw(binary->left, var, ctx);
                if (lhs.is_error()) return fail<PolyExpr>(lhs.error());
                auto rhs = parse_polynomial_raw(binary->right, var, ctx);
                if (rhs.is_error()) return fail<PolyExpr>(rhs.error());
                return poly_multiply(lhs.value(), rhs.value(), ctx);
            }
            case BinaryOp::Div: {
                if (poly_depends_on(binary->right, var.name)) {
                    return fail<PolyExpr>(CASError{
                        .kind = CASErrorKind::Unimplemented,
                        .message = "Tower bridge: denominator depending on tower variable is not a polynomial",
                    });
                }
                auto numerator = parse_polynomial_raw(binary->left, var, ctx);
                if (numerator.is_error()) return fail<PolyExpr>(numerator.error());
                auto denominator = clone_expr_raw(binary->right, ctx);
                if (denominator.is_error()) return fail<PolyExpr>(denominator.error());
                return poly_divide_by_scalar(numerator.value(), denominator.value(), ctx);
            }
            case BinaryOp::Pow: {
                auto base = parse_polynomial_raw(binary->left, var, ctx);
                if (base.is_error()) return fail<PolyExpr>(base.error());
                auto exponent = poly_parse_nonnegative_integer_exponent(binary->right);
                if (exponent.is_error()) return fail<PolyExpr>(exponent.error());
                return poly_pow_raw(base.value(), exponent.value(), ctx);
            }
            default:
                return fail<PolyExpr>(CASError{
                    .kind = CASErrorKind::Unimplemented,
                    .message = "Tower bridge: unsupported binary operator in polynomial parsing",
                });
        }
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        PolyExpr result;
        for (ExprPtr term : sum->terms) {
            auto parsed = parse_polynomial_raw(term, var, ctx);
            if (parsed.is_error()) return fail<PolyExpr>(parsed.error());
            auto next = poly_add(result, parsed.value(), ctx);
            if (next.is_error()) return fail<PolyExpr>(next.error());
            result = std::move(next.value());
        }
        return ok(std::move(result));
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        auto one = make_constant_poly_raw(poly_make_integer(ctx.arena(), 1), ctx);
        if (one.is_error()) return fail<PolyExpr>(one.error());
        PolyExpr result = one.value();
        for (ExprPtr factor : product->factors) {
            auto parsed = parse_polynomial_raw(factor, var, ctx);
            if (parsed.is_error()) return fail<PolyExpr>(parsed.error());
            auto next = poly_multiply(result, parsed.value(), ctx);
            if (next.is_error()) return fail<PolyExpr>(next.error());
            result = std::move(next.value());
        }
        return ok(std::move(result));
    }
    return fail<PolyExpr>(CASError{
        .kind = CASErrorKind::Unimplemented,
        .message = "Tower bridge: expression not supported by raw polynomial parser",
    });
}

[[nodiscard]] std::string rational_key(const Rational& value) {
    return value.numerator().decimal() + "/" + value.denominator().decimal();
}

[[nodiscard]] std::string min_poly_sort_key(const AlgebraicNumber::CoeffVec& min_poly) {
    std::string out;
    for (std::size_t i = min_poly.size(); i > 0U; --i) {
        if (!out.empty()) out += '|';
        out += rational_key(min_poly[i - 1U]);
    }
    return out;
}

[[nodiscard]] bool same_generator_expr(
    ExprPtr expr,
    ExprPtr generator_expr,
    symbolic::CASContext& ctx) {
    if (generator_expr && structural_equal(expr, generator_expr)) return true;

    const auto* lhs = expr_cast<RootOf>(expr);
    const auto* rhs = expr_cast<RootOf>(generator_expr);
    if (!lhs || !rhs) return false;
    if (lhs->root_index != rhs->root_index) return false;

    auto lhs_mp = rootof_min_poly(*lhs, ctx);
    if (lhs_mp.is_error()) return false;
    auto rhs_mp = rootof_min_poly(*rhs, ctx);
    if (rhs_mp.is_error()) return false;
    return lhs_mp.value() == rhs_mp.value();
}

[[nodiscard]] AlgebraicNumber inner_from_rational(
    const Rational& value,
    const AlgebraicNumber::CoeffVec& min_poly) {
    return AlgebraicNumber({value}, min_poly);
}

[[nodiscard]] AlgebraicNumber inner_zero(const AlgebraicNumber::CoeffVec& min_poly) {
    return inner_from_rational(Rational(BigInt(0)), min_poly);
}

[[nodiscard]] AlgebraicNumber inner_one(const AlgebraicNumber::CoeffVec& min_poly) {
    return inner_from_rational(Rational(BigInt(1)), min_poly);
}

[[nodiscard]] AlgebraicNumber inner_alpha(const AlgebraicNumber::CoeffVec& min_poly) {
    return AlgebraicNumber(
        AlgebraicNumber::CoeffVec{Rational(BigInt(0)), Rational(BigInt(1))},
        min_poly);
}

[[nodiscard]] AlgebraicTowerTwoLevel tower_constant(
    const AlgebraicNumber& coeff,
    const std::vector<AlgebraicNumber>& min_poly) {
    return AlgebraicTowerTwoLevel({coeff}, min_poly);
}

[[nodiscard]] AlgebraicTowerTwoLevel tower_one(
    const AlgebraicNumber::CoeffVec& inner_min_poly,
    const std::vector<AlgebraicNumber>& outer_min_poly) {
    return tower_constant(inner_one(inner_min_poly), outer_min_poly);
}

[[nodiscard]] Result<std::size_t> bigint_to_size_t_nonneg(const BigInt& value) {
    if (value.is_negative()) {
        return fail<std::size_t>(make_error(
            CASErrorKind::InvalidArgument,
            "Tower bridge: negative integer exponent not representable as size_t"));
    }
    if (value.bit_length() > 63U) {
        return fail<std::size_t>(make_error(
            CASErrorKind::Unimplemented,
            "Tower bridge: integer exponent exceeds 64 bits"));
    }
    return ok(static_cast<std::size_t>(value.to_u64()));
}

[[nodiscard]] Result<std::vector<AlgebraicNumber>> monic_outer_min_poly(
    std::vector<AlgebraicNumber> min_poly) {
    tower_detail::strip_trailing(min_poly);
    if (min_poly.empty() || min_poly.size() < 2U) {
        return fail<std::vector<AlgebraicNumber>>(make_error(
            CASErrorKind::InvalidArgument,
            "Tower bridge: outer minimal polynomial must have positive degree"));
    }
    if (min_poly.back().is_zero()) {
        return fail<std::vector<AlgebraicNumber>>(make_error(
            CASErrorKind::InvalidArgument,
            "Tower bridge: outer minimal polynomial leading coefficient is zero"));
    }
    auto lead_inv = min_poly.back().inverse();
    if (lead_inv.is_error()) {
        return fail<std::vector<AlgebraicNumber>>(lead_inv.error());
    }
    for (auto& coeff : min_poly) coeff = coeff * lead_inv.value();
    tower_detail::strip_trailing(min_poly);
    return ok(std::move(min_poly));
}

[[nodiscard]] Result<std::optional<std::vector<AlgebraicNumber>>> rootof_min_poly_over_q_alpha(
    const RootOf& root,
    ExprPtr alpha_1_expr,
    const AlgebraicNumber::CoeffVec& min_poly_1,
    symbolic::CASContext& ctx) {
    RootOfExplicitDegreeGuard rootof_guard(ctx, 1U);
    auto parsed = parse_polynomial_raw(root.polynomial, root.variable, ctx);
    if (parsed.is_error()) return ok(std::optional<std::vector<AlgebraicNumber>>{});

    std::vector<AlgebraicNumber> coeffs;
    coeffs.reserve(parsed.value().coefficients().size());
    for (ExprPtr coeff_expr : parsed.value().coefficients()) {
        if (!coeff_expr || poly_is_zero_expr(coeff_expr)) {
            coeffs.push_back(inner_zero(min_poly_1));
            continue;
        }
        auto reduced = try_express_in_q_alpha(coeff_expr, alpha_1_expr, min_poly_1, ctx);
        if (reduced.is_error()) return fail<std::optional<std::vector<AlgebraicNumber>>>(reduced.error());
        if (!reduced.value().has_value()) return ok(std::optional<std::vector<AlgebraicNumber>>{});
        coeffs.push_back(std::move(reduced.value().value()));
    }

    auto monic = monic_outer_min_poly(std::move(coeffs));
    if (monic.is_error()) return ok(std::optional<std::vector<AlgebraicNumber>>{});
    return ok(std::optional<std::vector<AlgebraicNumber>>(std::move(monic.value())));
}

[[nodiscard]] bool same_outer_generator_expr(
    ExprPtr expr,
    ExprPtr generator_expr,
    ExprPtr alpha_1_expr,
    const AlgebraicNumber::CoeffVec& min_poly_1,
    symbolic::CASContext& ctx) {
    if (same_generator_expr(expr, generator_expr, ctx)) return true;

    const auto* lhs = expr_cast<RootOf>(expr);
    const auto* rhs = expr_cast<RootOf>(generator_expr);
    if (!lhs || !rhs) return false;
    if (lhs->root_index != rhs->root_index) return false;

    auto lhs_mp = rootof_min_poly_over_q_alpha(*lhs, alpha_1_expr, min_poly_1, ctx);
    if (lhs_mp.is_error() || !lhs_mp.value().has_value()) return false;
    auto rhs_mp = rootof_min_poly_over_q_alpha(*rhs, alpha_1_expr, min_poly_1, ctx);
    if (rhs_mp.is_error() || !rhs_mp.value().has_value()) return false;
    return lhs_mp.value().value() == rhs_mp.value().value();
}

void collect_distinct_rootofs(ExprPtr expr, std::vector<ExprPtr>& out, symbolic::CASContext& ctx) {
    if (!expr || out.size() >= 3U) return;
    if (expr_is<RootOf>(expr)) {
        bool seen = false;
        for (ExprPtr existing : out) {
            if (same_generator_expr(existing, expr, ctx)) {
                seen = true;
                break;
            }
        }
        if (!seen) out.push_back(expr);
        const auto& root = expr_ref<RootOf>(expr);
        collect_distinct_rootofs(root.polynomial, out, ctx);
        return;
    }

    visit_expr(expr, [&](const auto& node) {
        using Node = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<Node, Unary>) {
            collect_distinct_rootofs(node.operand, out, ctx);
        } else if constexpr (std::is_same_v<Node, Binary>) {
            collect_distinct_rootofs(node.left, out, ctx);
            collect_distinct_rootofs(node.right, out, ctx);
        } else if constexpr (std::is_same_v<Node, FuncCall>) {
            for (ExprPtr arg : node.args) collect_distinct_rootofs(arg, out, ctx);
        } else if constexpr (std::is_same_v<Node, Sum>) {
            for (ExprPtr term : node.terms) collect_distinct_rootofs(term, out, ctx);
        } else if constexpr (std::is_same_v<Node, Product>) {
            for (ExprPtr factor : node.factors) collect_distinct_rootofs(factor, out, ctx);
        } else if constexpr (std::is_same_v<Node, Integral>) {
            collect_distinct_rootofs(node.integrand, out, ctx);
            if (node.lower.has_value()) collect_distinct_rootofs(*node.lower, out, ctx);
            if (node.upper.has_value()) collect_distinct_rootofs(*node.upper, out, ctx);
        } else if constexpr (std::is_same_v<Node, Derivative>) {
            collect_distinct_rootofs(node.expression, out, ctx);
        } else if constexpr (std::is_same_v<Node, Limit>) {
            collect_distinct_rootofs(node.expression, out, ctx);
            collect_distinct_rootofs(node.point, out, ctx);
        } else if constexpr (std::is_same_v<Node, Matrix>) {
            for (ExprPtr item : node.elements) collect_distinct_rootofs(item, out, ctx);
        }
    });
}

[[nodiscard]] Result<std::optional<TowerGenerators>> try_build_candidate(
    ExprPtr alpha_1,
    ExprPtr alpha_2,
    symbolic::CASContext& ctx) {
    auto alpha_1_norm = canonicalize_root_expr(alpha_1, ctx);
    if (alpha_1_norm.is_error()) return fail<std::optional<TowerGenerators>>(alpha_1_norm.error());
    auto alpha_2_norm = canonicalize_root_expr(alpha_2, ctx);
    if (alpha_2_norm.is_error()) return fail<std::optional<TowerGenerators>>(alpha_2_norm.error());

    const auto* root_1 = expr_cast<RootOf>(alpha_1_norm.value());
    const auto* root_2 = expr_cast<RootOf>(alpha_2_norm.value());
    if (!root_1 || !root_2) return ok(std::optional<TowerGenerators>{});

    auto min_poly_1 = rootof_min_poly(*root_1, ctx);
    if (min_poly_1.is_error()) return ok(std::optional<TowerGenerators>{});

    RootOfExplicitDegreeGuard rootof_guard(ctx, 1U);
    auto parsed_outer = parse_polynomial_raw(root_2->polynomial, root_2->variable, ctx);
    if (parsed_outer.is_error()) return ok(std::optional<TowerGenerators>{});

    std::vector<AlgebraicNumber> outer_coeffs;
    const auto& parsed_coeffs = parsed_outer.value().coefficients();
    outer_coeffs.reserve(parsed_coeffs.size());
    for (ExprPtr coeff_expr : parsed_coeffs) {
        if (!coeff_expr || poly_is_zero_expr(coeff_expr)) {
            outer_coeffs.push_back(inner_zero(min_poly_1.value()));
            continue;
        }
        auto coeff = try_express_in_q_alpha(coeff_expr, alpha_1_norm.value(), min_poly_1.value(), ctx);
        if (coeff.is_error()) return fail<std::optional<TowerGenerators>>(coeff.error());
        if (!coeff.value().has_value()) return ok(std::optional<TowerGenerators>{});
        outer_coeffs.push_back(std::move(coeff.value().value()));
    }

    auto min_poly_2 = monic_outer_min_poly(std::move(outer_coeffs));
    if (min_poly_2.is_error()) return ok(std::optional<TowerGenerators>{});

    TowerGenerators gens{
        .alpha_1 = alpha_1_norm.value(),
        .min_poly_1 = min_poly_1.value(),
        .alpha_2 = alpha_2_norm.value(),
        .min_poly_2 = min_poly_2.value(),
    };
    return ok(std::optional<TowerGenerators>(std::move(gens)));
}

[[nodiscard]] Result<std::optional<AlgebraicTowerTwoLevel>> express_integer_power(
    const AlgebraicTowerTwoLevel& base,
    const BigInt& exponent,
    const TowerGenerators& gens) {
    if (exponent.is_zero()) {
        return ok(std::optional<AlgebraicTowerTwoLevel>(
            tower_one(gens.min_poly_1, gens.min_poly_2)));
    }
    const bool negative = exponent.is_negative();
    auto power_res = bigint_to_size_t_nonneg(exponent.abs());
    if (power_res.is_error()) return fail<std::optional<AlgebraicTowerTwoLevel>>(power_res.error());

    AlgebraicTowerTwoLevel value = tower_one(gens.min_poly_1, gens.min_poly_2);
    AlgebraicTowerTwoLevel factor = base;
    std::size_t power = power_res.value();
    while (power > 0U) {
        if ((power & 1U) != 0U) value = value * factor;
        power >>= 1U;
        if (power > 0U) factor = factor * factor;
    }
    if (negative) {
        if (value.is_zero()) return ok(std::optional<AlgebraicTowerTwoLevel>{});
        auto inv_res = value.inverse();
        if (inv_res.is_error()) return fail<std::optional<AlgebraicTowerTwoLevel>>(inv_res.error());
        value = inv_res.value();
    }
    return ok(std::optional<AlgebraicTowerTwoLevel>(std::move(value)));
}

[[nodiscard]] Result<std::optional<AlgebraicTowerTwoLevel>> express_recursive(
    ExprPtr expr,
    const TowerGenerators& gens,
    symbolic::CASContext& ctx,
    unsigned int depth) {
    if (!expr) return ok(std::optional<AlgebraicTowerTwoLevel>{});
    if (static_cast<std::size_t>(depth) > ctx.max_q_alpha_bridge_depth()) {
        return ok(std::optional<AlgebraicTowerTwoLevel>{});
    }

    if (same_outer_generator_expr(expr, gens.alpha_2, gens.alpha_1, gens.min_poly_1, ctx)) {
        return ok(std::optional<AlgebraicTowerTwoLevel>(
            AlgebraicTowerTwoLevel(
                {inner_zero(gens.min_poly_1), inner_one(gens.min_poly_1)},
                gens.min_poly_2)));
    }
    if (same_generator_expr(expr, gens.alpha_1, ctx)) {
        return ok(std::optional<AlgebraicTowerTwoLevel>(
            tower_constant(inner_alpha(gens.min_poly_1), gens.min_poly_2)));
    }

    if (const auto* lit = expr_cast<IntegerLit>(expr)) {
        return ok(std::optional<AlgebraicTowerTwoLevel>(
            tower_constant(inner_from_rational(Rational(lit->value), gens.min_poly_1), gens.min_poly_2)));
    }
    if (const auto* lit = expr_cast<RationalLit>(expr)) {
        return ok(std::optional<AlgebraicTowerTwoLevel>(
            tower_constant(
                inner_from_rational(Rational(lit->numerator, lit->denominator), gens.min_poly_1),
                gens.min_poly_2)));
    }
    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op != UnaryOp::Neg) return ok(std::optional<AlgebraicTowerTwoLevel>{});
        auto inner = express_recursive(unary->operand, gens, ctx, depth + 1U);
        if (inner.is_error()) return inner;
        if (!inner.value().has_value()) return ok(std::optional<AlgebraicTowerTwoLevel>{});
        return ok(std::optional<AlgebraicTowerTwoLevel>(-inner.value().value()));
    }
    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Pow) {
            const auto* exponent = expr_cast<IntegerLit>(binary->right);
            if (!exponent) return ok(std::optional<AlgebraicTowerTwoLevel>{});
            auto base = express_recursive(binary->left, gens, ctx, depth + 1U);
            if (base.is_error()) return base;
            if (!base.value().has_value()) return ok(std::optional<AlgebraicTowerTwoLevel>{});
            return express_integer_power(base.value().value(), exponent->value, gens);
        }

        auto lhs = express_recursive(binary->left, gens, ctx, depth + 1U);
        if (lhs.is_error()) return lhs;
        if (!lhs.value().has_value()) return ok(std::optional<AlgebraicTowerTwoLevel>{});

        auto rhs = express_recursive(binary->right, gens, ctx, depth + 1U);
        if (rhs.is_error()) return rhs;
        if (!rhs.value().has_value()) return ok(std::optional<AlgebraicTowerTwoLevel>{});

        switch (binary->op) {
            case BinaryOp::Add:
                return ok(std::optional<AlgebraicTowerTwoLevel>(lhs.value().value() + rhs.value().value()));
            case BinaryOp::Sub:
                return ok(std::optional<AlgebraicTowerTwoLevel>(lhs.value().value() - rhs.value().value()));
            case BinaryOp::Mul:
                return ok(std::optional<AlgebraicTowerTwoLevel>(lhs.value().value() * rhs.value().value()));
            case BinaryOp::Div: {
                if (rhs.value().value().is_zero()) return ok(std::optional<AlgebraicTowerTwoLevel>{});
                auto div_res = lhs.value().value().div(rhs.value().value());
                if (div_res.is_error()) return ok(std::optional<AlgebraicTowerTwoLevel>{});
                return ok(std::optional<AlgebraicTowerTwoLevel>(div_res.value()));
            }
            default:
                return ok(std::optional<AlgebraicTowerTwoLevel>{});
        }
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        AlgebraicTowerTwoLevel acc({}, gens.min_poly_2);
        for (ExprPtr term : sum->terms) {
            auto item = express_recursive(term, gens, ctx, depth + 1U);
            if (item.is_error()) return item;
            if (!item.value().has_value()) return ok(std::optional<AlgebraicTowerTwoLevel>{});
            acc = acc + item.value().value();
        }
        return ok(std::optional<AlgebraicTowerTwoLevel>(std::move(acc)));
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        AlgebraicTowerTwoLevel acc = tower_one(gens.min_poly_1, gens.min_poly_2);
        for (ExprPtr factor : product->factors) {
            auto item = express_recursive(factor, gens, ctx, depth + 1U);
            if (item.is_error()) return item;
            if (!item.value().has_value()) return ok(std::optional<AlgebraicTowerTwoLevel>{});
            acc = acc * item.value().value();
        }
        return ok(std::optional<AlgebraicTowerTwoLevel>(std::move(acc)));
    }

    return ok(std::optional<AlgebraicTowerTwoLevel>{});
}

}  // namespace

[[nodiscard]] Result<std::optional<TowerGenerators>> detect_two_level_tower(
    ExprPtr expr,
    symbolic::CASContext& ctx) {
    std::vector<ExprPtr> roots;
    collect_distinct_rootofs(expr, roots, ctx);
    if (roots.size() != 2U) return ok(std::optional<TowerGenerators>{});

    auto first = try_build_candidate(roots[0], roots[1], ctx);
    if (first.is_error()) return first;
    auto second = try_build_candidate(roots[1], roots[0], ctx);
    if (second.is_error()) return second;

    if (!first.value().has_value()) return second;
    if (!second.value().has_value()) return first;

    const std::string first_key = min_poly_sort_key(first.value()->min_poly_1);
    const std::string second_key = min_poly_sort_key(second.value()->min_poly_1);
    if (second_key < first_key) return second;
    return first;
}

[[nodiscard]] Result<std::optional<AlgebraicTowerTwoLevel>> try_express_in_tower_two_level(
    ExprPtr expr,
    const TowerGenerators& gens,
    symbolic::CASContext& ctx) {
    if (gens.min_poly_1.size() < 2U || gens.min_poly_2.size() < 2U) {
        return fail<std::optional<AlgebraicTowerTwoLevel>>(make_error(
            CASErrorKind::InvalidArgument,
            "Tower bridge: both tower minimal polynomials must have positive degree"));
    }
    auto canonical_expr = canonicalize_root_expr(expr, ctx);
    if (canonical_expr.is_error()) return fail<std::optional<AlgebraicTowerTwoLevel>>(canonical_expr.error());
    return express_recursive(canonical_expr.value(), gens, ctx, 0U);
}

[[nodiscard]] ExprPtr tower_to_expr(
    const AlgebraicTowerTwoLevel& value,
    const TowerGenerators& gens,
    AstArena& arena) {
    const auto& coeffs = value.value();
    if (coeffs.empty()) return arena.make<IntegerLit>(BigInt(0));

    std::vector<ExprPtr> terms;
    terms.reserve(coeffs.size());
    for (std::size_t degree = 0; degree < coeffs.size(); ++degree) {
        if (coeffs[degree].is_zero()) continue;
        ExprPtr coeff_expr = algebraic_number_to_expr_raw(coeffs[degree], gens.alpha_1, arena);
        ExprPtr term = coeff_expr;
        if (degree > 0U) {
            ExprPtr alpha_power = (degree == 1U)
                ? gens.alpha_2
                : arena.make<Binary>(
                      BinaryOp::Pow,
                      gens.alpha_2,
                      arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(degree))));
            term = arena.make<Binary>(BinaryOp::Mul, coeff_expr, alpha_power);
        }
        terms.push_back(term);
    }
    if (terms.empty()) return arena.make<IntegerLit>(BigInt(0));
    if (terms.size() == 1U) return terms.front();
    return arena.make<Sum>(std::move(terms));
}

}  // namespace algebra
}  // namespace cas
