#include "integrate_definite_patterns.hpp"

#include "cas/extended_real.hpp"
#include "cas/rational.hpp"
#include "cas/residue_theorem.hpp"

#include <optional>

namespace cas::calculus {

namespace {

// Adopt the canonical extended-real predicates from cas::; legacy local
// copies missed Constant(NegInfinity).
using cas::is_pos_infinity;
using cas::is_neg_infinity;

// Match exp(-a*x^2) for positive rational a.  Returns a.
[[nodiscard]] std::optional<Rational> match_gaussian_exp(ExprPtr expr, const Symbol& var) {
    const auto* call = expr_cast<FuncCall>(expr);
    if (!call || call->func_id != BuiltinOp::Exp || call->args.size() != 1U)
        return std::nullopt;

    ExprPtr arg = call->args.front();

    if (const auto* neg = expr_cast<Unary>(arg)) {
        if (neg->op == UnaryOp::Neg) {
            if (const auto* pw = expr_cast<Binary>(neg->operand)) {
                if (pw->op == BinaryOp::Pow) {
                    if (const auto* sym = expr_cast<Symbol>(pw->left)) {
                        if (sym->name == var.name) {
                            if (const auto* e2 = expr_cast<IntegerLit>(pw->right)) {
                                if (e2->value == BigInt(2)) return Rational(BigInt(1));
                            }
                        }
                    }
                }
            }
        }
    }

    if (const auto* mul = expr_cast<Binary>(arg)) {
        if (mul->op == BinaryOp::Mul) {
            auto try_match = [&](ExprPtr coeff_expr, ExprPtr pow_expr) -> std::optional<Rational> {
                const auto* pw = expr_cast<Binary>(pow_expr);
                if (!pw || pw->op != BinaryOp::Pow) return std::nullopt;
                const auto* sym = expr_cast<Symbol>(pw->left);
                if (!sym || sym->name != var.name) return std::nullopt;
                const auto* e2 = expr_cast<IntegerLit>(pw->right);
                if (!e2 || e2->value != BigInt(2)) return std::nullopt;
                if (const auto* i = expr_cast<IntegerLit>(coeff_expr)) {
                    if (i->value.is_negative()) return Rational(-i->value);
                } else if (const auto* r = expr_cast<RationalLit>(coeff_expr)) {
                    Rational v(r->numerator, r->denominator);
                    if (v.numerator().is_negative()) return -v;
                }
                return std::nullopt;
            };
            auto m = try_match(mul->left, mul->right);
            if (!m) m = try_match(mul->right, mul->left);
            if (m) return m;
        }
    }

    return std::nullopt;
}

}  // namespace

[[nodiscard]] Result<std::optional<ExprPtr>> pattern_gaussian_full_line(const DefiniteContext& dc) {
    if (!is_neg_infinity(dc.lower) || !is_pos_infinity(dc.upper)) {
        return ok(std::optional<ExprPtr>{});
    }
    auto a = match_gaussian_exp(dc.integrand, dc.var);
    if (!a.has_value()) {
        // Also try the normalized form (post simplify).
        a = match_gaussian_exp(dc.integrand_normalized, dc.var);
        if (!a.has_value()) return ok(std::optional<ExprPtr>{});
    }

    AstArena& arena = dc.ctx.arena();
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    if (a->numerator() == BigInt(1) && a->denominator() == BigInt(1)) {
        ExprPtr sqrt_pi = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{pi});
        auto simplified = dc.ctx.simplify(sqrt_pi);
        if (simplified.is_error()) return fail<std::optional<ExprPtr>>(simplified.error());
        return ok(std::optional<ExprPtr>(simplified.value()));
    }
    ExprPtr a_expr = a->is_integer()
        ? static_cast<ExprPtr>(arena.make<IntegerLit>(a->numerator()))
        : static_cast<ExprPtr>(arena.make<RationalLit>(a->numerator(), a->denominator()));
    ExprPtr pi_over_a = arena.make<Binary>(BinaryOp::Div, pi, a_expr);
    ExprPtr sqrt_pi_a = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{pi_over_a});
    auto simplified = dc.ctx.simplify(sqrt_pi_a);
    if (simplified.is_error()) return fail<std::optional<ExprPtr>>(simplified.error());
    return ok(std::optional<ExprPtr>(simplified.value()));
}

[[nodiscard]] Result<std::optional<ExprPtr>> pattern_rational_full_real_line(const DefiniteContext& dc) {
    if (!is_neg_infinity(dc.lower) || !is_pos_infinity(dc.upper))
        return ok(std::optional<ExprPtr>{});
    auto result = integrate_rational_full_real_line(dc.integrand_normalized, dc.var, dc.ctx);
    if (result.is_error()) {
        if (result.error().kind == CASErrorKind::Unimplemented)
            return ok(std::optional<ExprPtr>{});
        return fail<std::optional<ExprPtr>>(result.error());
    }
    return ok(std::optional<ExprPtr>(result.value()));
}

const std::vector<DefinitePatternFn>& definite_patterns() {
    static const std::vector<DefinitePatternFn> registry = {
        &pattern_gaussian_full_line,
        &pattern_rational_full_real_line,
        &pattern_bessel_orthogonality,
        &pattern_legendre_orthogonality,
        &pattern_hermite_h_orthogonality,
        &pattern_hermite_he_orthogonality,
        &pattern_chebyshev_t_orthogonality,
        &pattern_chebyshev_u_orthogonality,
        &pattern_mellin_g_convolution,
    };
    return registry;
}

}  // namespace cas::calculus
