#include "cas/calculus.hpp"
#include "integrate_engine.hpp"

namespace cas::calculus {

namespace {

[[nodiscard]] bool is_pos_infinity(ExprPtr expr) {
    const auto* c = expr_cast<Constant>(expr);
    return c != nullptr && c->value == MathConstant::Infinity;
}

[[nodiscard]] bool is_neg_infinity(ExprPtr expr) {
    const auto* u = expr_cast<Unary>(expr);
    return u != nullptr && u->op == UnaryOp::Neg && is_pos_infinity(u->operand);
}

// Matches exp(-a*x^2) for positive rational a. Returns a on success.
[[nodiscard]] std::optional<Rational> match_gaussian_exp(ExprPtr expr, const Symbol& var) {
    const auto* call = expr_cast<FuncCall>(expr);
    if (!call || call->func_id != BuiltinOp::Exp || call->args.size() != 1U)
        return std::nullopt;

    ExprPtr arg = call->args.front();

    // Pattern: exp(-x^2) = exp(Neg(Pow(x, 2)))
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

    // Pattern: exp(-a*x^2) = exp(Mul(neg_rational, Pow(x,2)))
    if (const auto* mul = expr_cast<Binary>(arg)) {
        if (mul->op == BinaryOp::Mul) {
            // Try both orderings: -a * x^2 and x^2 * -a
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

} // anonymous namespace

Result<ExprPtr> integrate(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    auto primitive = integrate_detail::integrate_indefinite_impl(expr, var, ctx);
    if (primitive.is_error()) {
        return primitive;
    }
    return symbolic::materialize_expr(primitive.value(), ctx.arena());
}

Result<ExprPtr> definite_integral(ExprPtr expr, const Symbol& var, ExprPtr lower, ExprPtr upper, symbolic::CASContext& ctx) {
    // Gaussian integral: integral from -inf to +inf of exp(-a*x^2) dx = sqrt(pi/a)
    if (is_neg_infinity(lower) && is_pos_infinity(upper)) {
        if (auto a = match_gaussian_exp(expr, var)) {
            AstArena& arena = ctx.arena();
            ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
            if (a->numerator() == BigInt(1) && a->denominator() == BigInt(1)) {
                // a=1: sqrt(pi)
                ExprPtr sqrt_pi = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{pi});
                return ctx.simplify(sqrt_pi);
            }
            // General a: sqrt(pi/a) = sqrt(pi) / sqrt(a)
            ExprPtr a_expr = a->is_integer()
                ? static_cast<ExprPtr>(arena.make<IntegerLit>(a->numerator()))
                : static_cast<ExprPtr>(arena.make<RationalLit>(a->numerator(), a->denominator()));
            ExprPtr pi_over_a = arena.make<Binary>(BinaryOp::Div, pi, a_expr);
            ExprPtr sqrt_pi_a = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{pi_over_a});
            return ctx.simplify(sqrt_pi_a);
        }
        return fail<ExprPtr>(integrate_detail::make_error(CASErrorKind::Unimplemented,
            "Integrazione su dominio infinito: pattern non riconosciuto."));
    }

    auto primitive = integrate(expr, var, ctx);
    if (primitive.is_error()) {
        return primitive;
    }

    auto lower_value = ctx.substitute(primitive.value(), var, lower);
    if (lower_value.is_error()) {
        return lower_value;
    }

    auto upper_value = ctx.substitute(primitive.value(), var, upper);
    if (upper_value.is_error()) {
        return upper_value;
    }

    return ctx.simplify(integrate_detail::make_sum(ctx.arena(), {
        upper_value.value(),
        integrate_detail::make_unary(ctx.arena(), UnaryOp::Neg, lower_value.value()),
    }));
}

}  // namespace cas::calculus
