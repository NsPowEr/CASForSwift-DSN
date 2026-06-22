#include "cas/algebraic_number_bridge.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"
#include "cas/formatter.hpp"

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
// Real Q(alpha) expressions are shallow; configurable via
// ctx.max_q_alpha_bridge_depth() (default 256).
[[nodiscard]] Result<std::optional<AlgebraicNumber>> express_recursive(
    ExprPtr e,
    ExprPtr alpha_expr,
    const AlgebraicNumber::CoeffVec& min_poly,
    symbolic::CASContext& ctx,
    unsigned int depth) {
    if (static_cast<std::size_t>(depth) > ctx.max_q_alpha_bridge_depth()) {
        return ok(std::optional<AlgebraicNumber>{});
    }
    if (!e) {
        return ok(std::optional<AlgebraicNumber>{});
    }

    // Structural match against the supplied generator expression.
    if (alpha_expr && same_generator_expr(e, alpha_expr, ctx)) {
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

[[nodiscard]] bool same_generator_expr(
    ExprPtr expr,
    ExprPtr alpha_expr,
    symbolic::CASContext& ctx) {
    if (alpha_expr && structural_equal(expr, alpha_expr)) return true;

    const auto* lhs = expr_cast<RootOf>(expr);
    const auto* rhs = expr_cast<RootOf>(alpha_expr);
    if (!lhs || !rhs) return false;
    if (lhs->root_index != rhs->root_index) return false;

    auto lhs_mp = rootof_min_poly(*lhs, ctx);
    if (lhs_mp.is_error()) return false;
    auto rhs_mp = rootof_min_poly(*rhs, ctx);
    if (rhs_mp.is_error()) return false;
    return lhs_mp.value() == rhs_mp.value();
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

} // namespace algebra
} // namespace cas
