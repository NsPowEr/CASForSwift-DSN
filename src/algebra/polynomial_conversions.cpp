#include "polynomial_internal.hpp"
#include "cas/symbolic.hpp"
#include "cas/numtheory.hpp"
#include <algorithm>
#include <utility>

namespace cas {
namespace algebra {

[[nodiscard]] Result<Rational> expr_to_exact_rational_coefficient(ExprPtr expr) {
    if (!expr) {
        return ok(Rational(BigInt(0)));
    }
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return ok(Rational(integer->value));
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return ok(Rational(rational->numerator, rational->denominator));
    }
    return fail<Rational>(CASError{
        .kind = CASErrorKind::Unimplemented,
        .message = "partial_fractions richiede coefficienti razionali esatti",
    });
}

[[nodiscard]] Result<BigInt> expr_to_integer_coefficient(ExprPtr expr) {
    if (!expr) {
        return ok(BigInt(0));
    }
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return ok(integer->value);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        if (rational->denominator == BigInt(1)) {
            return ok(rational->numerator);
        }
    }
    return fail<BigInt>(CASError{
        .kind = CASErrorKind::Unimplemented,
        .message = "factor_over_integers supporta solo coefficienti interi esatti",
    });
}

Result<IntPoly> poly_to_integer_poly(const PolyExpr& poly) {
    if (poly.is_zero()) {
        return ok(IntPoly{});
    }

    // Pass 1: find all denominators and compute LCM
    BigInt common_lcm(1);
    for (ExprPtr coefficient : poly.coefficients()) {
        if (!coefficient || poly_is_zero_expr(coefficient)) {
            continue;
        }
        auto rat_res = expr_to_exact_rational_coefficient(coefficient);
        if (rat_res.is_error()) {
            return fail<IntPoly>(rat_res.error());
        }
        const Rational& rat = rat_res.value();
        if (!rat.denominator().is_zero()) {
            BigInt d = rat.denominator();
            common_lcm = (common_lcm * d) / gcd(common_lcm, d);
        }
    }

    // Pass 2: multiply all by LCM to get integers
    IntPoly coefficients;
    coefficients.reserve(poly.size());
    for (ExprPtr coefficient : poly.coefficients()) {
        if (!coefficient || poly_is_zero_expr(coefficient)) {
            coefficients.push_back(BigInt(0));
            continue;
        }
        auto rat_res = expr_to_exact_rational_coefficient(coefficient);
        // We already checked errors in pass 1, but for safety:
        if (rat_res.is_error()) return fail<IntPoly>(rat_res.error());
        
        Rational adjusted = rat_res.value() * Rational(common_lcm);
        if (adjusted.denominator() != BigInt(1)) {
            return fail<IntPoly>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "Errore interno nella rimozione dei denominatori",
            });
        }
        coefficients.push_back(adjusted.numerator());
    }
    
    normalize_integer_poly(coefficients);
    return ok(std::move(coefficients));
}

Result<std::vector<BigInt>> poly_to_integer_coefficients(const PolyExpr& poly) {
    auto coefficients = poly_to_integer_poly(poly);
    if (coefficients.is_error()) {
        return fail<std::vector<BigInt>>(coefficients.error());
    }
    return ok(coefficients.value().coefficients());
}

Result<RatPoly> poly_to_rational_poly(const PolyExpr& poly) {
    RatPoly coefficients;
    coefficients.reserve(poly.size());
    for (ExprPtr coefficient : poly.coefficients()) {
        auto rational = expr_to_exact_rational_coefficient(coefficient);
        if (rational.is_error()) {
            return fail<RatPoly>(rational.error());
        }
        coefficients.push_back(rational.value());
    }
    normalize_rational_coefficients(coefficients);
    return ok(std::move(coefficients));
}

Result<std::vector<Rational>> poly_to_rational_coefficients(const PolyExpr& poly) {
    auto coefficients = poly_to_rational_poly(poly);
    if (coefficients.is_error()) {
        return fail<std::vector<Rational>>(coefficients.error());
    }
    return ok(coefficients.value().coefficients());
}

ExprPtr poly_make_integer(AstArena& arena, long long value) {
    return arena.make<IntegerLit>(BigInt(value));
}

bool poly_is_zero_expr(ExprPtr expr) {
    if (!expr) return true;
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value.is_zero();
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator.is_zero();
    }
    return false;
}

bool poly_is_one_expr(ExprPtr expr) {
    if (!expr) return false;
    static const BigInt one(1);
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value == one;
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator == one && rational->denominator == one;
    }
    return false;
}

bool poly_is_minus_one_expr(ExprPtr expr) {
    if (!expr) return false;
    static const BigInt minus_one(-1);
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value == minus_one;
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator == minus_one && rational->denominator == BigInt(1);
    }
    return false;
}

bool poly_depends_on(ExprPtr expr, const std::string& variable_name) {
    if (!expr) {
        return false;
    }

    if (const auto* symbol = expr_cast<Symbol>(expr)) {
        return symbol->name == variable_name;
    }

    switch (expr->kind) {
    case ExprKind::Null:
    case ExprKind::IntegerLit:
    case ExprKind::RationalLit:
    case ExprKind::DecimalLit:
    case ExprKind::Symbol:
    case ExprKind::Constant:
    case ExprKind::ComplexLit:

        return false;
    case ExprKind::Unary:
        return poly_depends_on(expr_ref<Unary>(expr).operand, variable_name);
    case ExprKind::Binary: {
        const auto& binary = expr_ref<Binary>(expr);
        return poly_depends_on(binary.left, variable_name) || poly_depends_on(binary.right, variable_name);
    }
    case ExprKind::FuncCall: {
        for (ExprPtr arg : expr_ref<FuncCall>(expr).args) {
            if (poly_depends_on(arg, variable_name)) {
                return true;
            }
        }
        return false;
    }
    case ExprKind::Sum: {
        for (ExprPtr term : expr_ref<Sum>(expr).terms) {
            if (poly_depends_on(term, variable_name)) {
                return true;
            }
        }
        return false;
    }
    case ExprKind::Product: {
        for (ExprPtr factor : expr_ref<Product>(expr).factors) {
            if (poly_depends_on(factor, variable_name)) {
                return true;
            }
        }
        return false;
    }
    case ExprKind::Integral: {
        const auto& node = expr_ref<Integral>(expr);
        return poly_depends_on(node.integrand, variable_name) ||
               (node.lower.has_value() && poly_depends_on(*node.lower, variable_name)) ||
               (node.upper.has_value() && poly_depends_on(*node.upper, variable_name));
    }
    case ExprKind::Derivative:
        return poly_depends_on(expr_ref<Derivative>(expr).expression, variable_name);
    case ExprKind::Limit: {
        const auto& node = expr_ref<Limit>(expr);
        return poly_depends_on(node.expression, variable_name) || poly_depends_on(node.point, variable_name);
    }
    case ExprKind::RootOf:
        return poly_depends_on(expr_ref<RootOf>(expr).polynomial, variable_name);
    case ExprKind::Matrix: {
        for (ExprPtr element : expr_ref<Matrix>(expr).elements) {
            if (poly_depends_on(element, variable_name)) {
                return true;
            }
        }
        return false;
    }
    case ExprKind::SeriesExp:
    case ExprKind::Quantity:
        return false;
    }

    return false;
}

bool poly_contains_decimal_literal(ExprPtr expr) {
    if (!expr) {
        return false;
    }
    if (expr_is<DecimalLit>(expr)) {
        return true;
    }

    return visit_expr(expr, [](const auto& node) {
        using Node = std::decay_t<decltype(node)>;
        if constexpr (
            std::is_same_v<Node, IntegerLit> ||
            std::is_same_v<Node, RationalLit> ||
            std::is_same_v<Node, Symbol> ||
            std::is_same_v<Node, Constant>) {
            return false;
        } else if constexpr (std::is_same_v<Node, Unary>) {
            return poly_contains_decimal_literal(node.operand);
        } else if constexpr (std::is_same_v<Node, Binary>) {
            return poly_contains_decimal_literal(node.left) || poly_contains_decimal_literal(node.right);
        } else if constexpr (std::is_same_v<Node, FuncCall>) {
            return std::any_of(node.args.begin(), node.args.end(), [](ExprPtr arg) {
                return poly_contains_decimal_literal(arg);
            });
        } else if constexpr (std::is_same_v<Node, Sum>) {
            return std::any_of(node.terms.begin(), node.terms.end(), [](ExprPtr term) {
                return poly_contains_decimal_literal(term);
            });
        } else if constexpr (std::is_same_v<Node, Product>) {
            return std::any_of(node.factors.begin(), node.factors.end(), [](ExprPtr factor) {
                return poly_contains_decimal_literal(factor);
            });
        } else if constexpr (std::is_same_v<Node, Integral>) {
            return poly_contains_decimal_literal(node.integrand) ||
                   (node.lower.has_value() && poly_contains_decimal_literal(*node.lower)) ||
                   (node.upper.has_value() && poly_contains_decimal_literal(*node.upper));
        } else if constexpr (std::is_same_v<Node, Derivative>) {
            return poly_contains_decimal_literal(node.expression);
        } else if constexpr (std::is_same_v<Node, Limit>) {
            return poly_contains_decimal_literal(node.expression) || poly_contains_decimal_literal(node.point);
        } else if constexpr (std::is_same_v<Node, RootOf>) {
            return poly_contains_decimal_literal(node.polynomial);
        } else if constexpr (std::is_same_v<Node, Matrix>) {
            return std::any_of(node.elements.begin(), node.elements.end(), [](ExprPtr element) {
                return poly_contains_decimal_literal(element);
            });
        } else {
            return false;
        }
    });
}

Result<ExprPtr> integer_coefficients_to_expr(const IntPoly& poly, const Symbol& var, symbolic::CASContext& ctx) {
    if (poly.empty()) {
        return ok(ctx.arena().make<IntegerLit>(BigInt(0)));
    }
    
    std::vector<ExprPtr> terms;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        const BigInt& coeff = poly[i];
        if (coeff.is_zero()) continue;
        
        ExprPtr term;
        if (i == 0) {
            term = ctx.arena().make<IntegerLit>(coeff);
        } else {
            ExprPtr var_node = ctx.arena().make<Symbol>(var.name);
            ExprPtr pow_node;
            if (i == 1) {
                pow_node = var_node;
            } else {
                pow_node = ctx.arena().make<Binary>(BinaryOp::Pow, var_node, ctx.arena().make<IntegerLit>(BigInt(i)));
            }
            
            if (coeff == BigInt(1)) {
                term = pow_node;
            } else if (coeff == BigInt(-1)) {
                term = ctx.arena().make<Unary>(UnaryOp::Neg, pow_node);
            } else {
                term = ctx.arena().make<Product>(std::vector<ExprPtr>{
                    ctx.arena().make<IntegerLit>(coeff),
                    pow_node
                });
            }
        }
        terms.push_back(term);
    }
    
    if (terms.empty()) return ok(ctx.arena().make<IntegerLit>(BigInt(0)));
    if (terms.size() == 1) return ok(terms[0]);
    
    // Reverse to have descending order
    std::reverse(terms.begin(), terms.end());
    return ok(ctx.arena().make<Sum>(std::move(terms)));
}

Result<ExprPtr> poly_simplify_expr(ExprPtr expr, symbolic::CASContext& ctx) {
    auto simplified = ctx.simplify(expr);
    if (simplified.is_error()) {
        return fail<ExprPtr>(simplified.error());
    }
    return simplified;
}

Result<ExprPtr> poly_clone_into_context(ExprPtr expr, symbolic::CASContext& ctx) {
    auto cloned = symbolic::materialize_expr(expr, ctx.arena());
    if (cloned.is_error()) {
        return fail<ExprPtr>(cloned.error());
    }
    return poly_simplify_expr(cloned.value(), ctx);
}

}  // namespace algebra
}  // namespace cas
