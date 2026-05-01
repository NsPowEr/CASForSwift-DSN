#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <vector>

namespace cas::algebra {

[[nodiscard]] Result<IntegerExponent> parse_integer_exponent(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        const bool is_neg = integer->value.is_negative();
        BigInt magnitude_val = integer->value;
        if (is_neg) {
            magnitude_val = -magnitude_val;
        }
        
        const std::string decimal = magnitude_val.decimal();
        const auto parsed = parse_bounded_unsigned_decimal<std::size_t>(decimal);
        if (!parsed.has_value()) {
            return fail<IntegerExponent>(make_error(
                CASErrorKind::Unimplemented,
                "Esponente troppo grande"));
        }
        return ok(IntegerExponent{.magnitude = parsed.value(), .negative = is_neg});
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op == UnaryOp::Neg && expr_is<IntegerLit>(unary->operand)) {
            const auto& inner = expr_ref<IntegerLit>(unary->operand);
            const std::string decimal = inner.value.decimal();
            const auto parsed = parse_bounded_unsigned_decimal<std::size_t>(decimal);
            if (!parsed.has_value()) {
                return fail<IntegerExponent>(make_error(
                    CASErrorKind::Unimplemented,
                    "Esponente troppo grande"));
            }
            return ok(IntegerExponent{.magnitude = parsed.value(), .negative = true});
        }
    }

    return fail<IntegerExponent>(make_error(
        CASErrorKind::Unimplemented,
        "Le potenze razionali supportano solo esponenti interi espliciti"));
}

[[nodiscard]] bool is_explicit_integer_exponent_expr(ExprPtr expr) {
    if (expr_is<IntegerLit>(expr)) {
        return true;
    }

    const auto* unary = expr_cast<Unary>(expr);
    return unary != nullptr &&
        unary->op == UnaryOp::Neg &&
        expr_is<IntegerLit>(unary->operand);
}

[[nodiscard]] Result<ExprPtr> make_sum_expr(std::vector<ExprPtr> terms, symbolic::CASContext& ctx) {
    if (terms.empty()) {
        return ok(make_integer(ctx.arena(), 0));
    }
    if (terms.size() == 1U) {
        return simplify_expr(terms.front(), ctx);
    }
    return simplify_expr(ctx.arena().make<Sum>(std::move(terms)), ctx);
}

[[nodiscard]] Result<ExprPtr> make_product_expr(std::vector<ExprPtr> factors, symbolic::CASContext& ctx) {
    if (factors.empty()) {
        return ok(make_integer(ctx.arena(), 1));
    }
    if (factors.size() == 1U) {
        return simplify_expr(factors.front(), ctx);
    }
    return simplify_expr(ctx.arena().make<Product>(std::move(factors)), ctx);
}

void append_sum_terms(std::vector<ExprPtr>& out, ExprPtr expr) {
    if (const auto* sum = expr_cast<Sum>(expr)) {
        out.insert(out.end(), sum->terms.begin(), sum->terms.end());
        return;
    }
    out.push_back(expr);
}

void append_product_factors(std::vector<ExprPtr>& out, ExprPtr expr) {
    if (const auto* product = expr_cast<Product>(expr)) {
        out.insert(out.end(), product->factors.begin(), product->factors.end());
        return;
    }
    out.push_back(expr);
}

[[nodiscard]] Result<ExprPtr> distribute_product_terms(const std::vector<ExprPtr>& factors, symbolic::CASContext& ctx) {
    std::vector<ExprPtr> accumulated_terms{make_integer(ctx.arena(), 1)};

    for (ExprPtr factor : factors) {
        std::vector<ExprPtr> factor_terms;
        append_sum_terms(factor_terms, factor);

        std::vector<ExprPtr> next_terms;
        next_terms.reserve(accumulated_terms.size() * factor_terms.size());

        for (ExprPtr lhs_term : accumulated_terms) {
            for (ExprPtr rhs_term : factor_terms) {
                std::vector<ExprPtr> term_factors;
                if (!is_one_expr(lhs_term)) {
                    append_product_factors(term_factors, lhs_term);
                }
                if (!is_one_expr(rhs_term)) {
                    append_product_factors(term_factors, rhs_term);
                }

                auto product = make_product_expr(std::move(term_factors), ctx);
                if (product.is_error()) {
                    return fail<ExprPtr>(product.error());
                }
                next_terms.push_back(product.value());
            }
        }

        accumulated_terms = std::move(next_terms);
    }

    return make_sum_expr(std::move(accumulated_terms), ctx);
}

[[nodiscard]] Result<ExprPtr> expand_power(ExprPtr base, const IntegerExponent& exponent, symbolic::CASContext& ctx) {
    auto expanded_base = expand_expr_impl(base, ctx);
    if (expanded_base.is_error()) {
        return fail<ExprPtr>(expanded_base.error());
    }

    if (exponent.negative) {
        if (exponent.magnitude == 0U) {
            return ok(make_integer(ctx.arena(), 1));
        }

        auto positive_power = expand_power(base, IntegerExponent{.magnitude = exponent.magnitude, .negative = false}, ctx);
        if (positive_power.is_error()) {
            return fail<ExprPtr>(positive_power.error());
        }
        return divide_exprs(make_integer(ctx.arena(), 1), positive_power.value(), ctx);
    }

    if (exponent.magnitude == 0U) {
        return ok(make_integer(ctx.arena(), 1));
    }
    if (exponent.magnitude == 1U) {
        return expanded_base;
    }

    ExprPtr result = make_integer(ctx.arena(), 1);
    for (std::size_t index = 0; index < exponent.magnitude; ++index) {
        auto multiplied = distribute_product_terms({result, expanded_base.value()}, ctx);
        if (multiplied.is_error()) {
            return fail<ExprPtr>(multiplied.error());
        }
        result = multiplied.value();
    }
    return ok(result);
}

[[nodiscard]] Result<ExprPtr> expand_expr_impl(ExprPtr expr, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "expand richiede un'espressione non nulla"));
    }
    if (contains_decimal_literal(expr)) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::Unimplemented,
            "I literal decimali non sono supportati nell'espansione simbolica esatta"));
    }

    if (
        expr_is<IntegerLit>(expr) ||
        expr_is<RationalLit>(expr) ||
        expr_is<Symbol>(expr) ||
        expr_is<Constant>(expr)) {
        return clone_into_context(expr, ctx);
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        auto operand = expand_expr_impl(unary->operand, ctx);
        if (operand.is_error()) {
            return fail<ExprPtr>(operand.error());
        }

        if (unary->op == UnaryOp::Neg) {
            return distribute_product_terms({make_integer(ctx.arena(), -1), operand.value()}, ctx);
        }
        return simplify_expr(ctx.arena().make<Unary>(unary->op, operand.value()), ctx);
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        switch (binary->op) {
        case BinaryOp::Add: {
            auto lhs = expand_expr_impl(binary->left, ctx);
            if (lhs.is_error()) {
                return fail<ExprPtr>(lhs.error());
            }
            auto rhs = expand_expr_impl(binary->right, ctx);
            if (rhs.is_error()) {
                return fail<ExprPtr>(rhs.error());
            }
            return make_sum_expr({lhs.value(), rhs.value()}, ctx);
        }
        case BinaryOp::Sub: {
            auto lhs = expand_expr_impl(binary->left, ctx);
            if (lhs.is_error()) {
                return fail<ExprPtr>(lhs.error());
            }
            auto rhs = expand_expr_impl(binary->right, ctx);
            if (rhs.is_error()) {
                return fail<ExprPtr>(rhs.error());
            }
            auto negated = distribute_product_terms({make_integer(ctx.arena(), -1), rhs.value()}, ctx);
            if (negated.is_error()) {
                return fail<ExprPtr>(negated.error());
            }
            return make_sum_expr({lhs.value(), negated.value()}, ctx);
        }
        case BinaryOp::Mul: {
            auto lhs = expand_expr_impl(binary->left, ctx);
            if (lhs.is_error()) {
                return fail<ExprPtr>(lhs.error());
            }
            auto rhs = expand_expr_impl(binary->right, ctx);
            if (rhs.is_error()) {
                return fail<ExprPtr>(rhs.error());
            }
            return distribute_product_terms({lhs.value(), rhs.value()}, ctx);
        }
        case BinaryOp::Div: {
            auto lhs = expand_expr_impl(binary->left, ctx);
            if (lhs.is_error()) {
                return fail<ExprPtr>(lhs.error());
            }
            auto rhs = expand_expr_impl(binary->right, ctx);
            if (rhs.is_error()) {
                return fail<ExprPtr>(rhs.error());
            }
            return divide_exprs(lhs.value(), rhs.value(), ctx);
        }
        case BinaryOp::Pow: {
            auto exponent = parse_integer_exponent(binary->right);
            if (exponent.is_ok()) {
                return expand_power(binary->left, exponent.value(), ctx);
            }
            if (is_explicit_integer_exponent_expr(binary->right)) {
                return fail<ExprPtr>(exponent.error());
            }

            auto expanded_base = expand_expr_impl(binary->left, ctx);
            if (expanded_base.is_error()) {
                return fail<ExprPtr>(expanded_base.error());
            }
            auto expanded_exponent = expand_expr_impl(binary->right, ctx);
            if (expanded_exponent.is_error()) {
                return fail<ExprPtr>(expanded_exponent.error());
            }
            return simplify_expr(
                ctx.arena().make<Binary>(BinaryOp::Pow, expanded_base.value(), expanded_exponent.value()),
                ctx);
        }
        case BinaryOp::Mod: {
            auto lhs = expand_expr_impl(binary->left, ctx);
            if (lhs.is_error()) {
                return fail<ExprPtr>(lhs.error());
            }
            auto rhs = expand_expr_impl(binary->right, ctx);
            if (rhs.is_error()) {
                return fail<ExprPtr>(rhs.error());
            }
            return simplify_expr(ctx.arena().make<Binary>(BinaryOp::Mod, lhs.value(), rhs.value()), ctx);
        }
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> expanded_terms;
        for (ExprPtr term : sum->terms) {
            auto expanded = expand_expr_impl(term, ctx);
            if (expanded.is_error()) {
                return fail<ExprPtr>(expanded.error());
            }
            append_sum_terms(expanded_terms, expanded.value());
        }
        return make_sum_expr(std::move(expanded_terms), ctx);
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        std::vector<ExprPtr> expanded_factors;
        for (ExprPtr factor : product->factors) {
            auto expanded = expand_expr_impl(factor, ctx);
            if (expanded.is_error()) {
                return fail<ExprPtr>(expanded.error());
            }
            expanded_factors.push_back(expanded.value());
        }
        return distribute_product_terms(expanded_factors, ctx);
    }

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        std::vector<ExprPtr> args;
        args.reserve(call->args.size());
        for (ExprPtr arg : call->args) {
            auto expanded = expand_expr_impl(arg, ctx);
            if (expanded.is_error()) {
                return fail<ExprPtr>(expanded.error());
            }
            args.push_back(expanded.value());
        }
        return simplify_expr(ctx.arena().make<FuncCall>(call->name, std::move(args)), ctx);
    }

    if (const auto* integral = expr_cast<Integral>(expr)) {
        auto integrand = expand_expr_impl(integral->integrand, ctx);
        if (integrand.is_error()) {
            return fail<ExprPtr>(integrand.error());
        }
        std::optional<ExprPtr> lower;
        std::optional<ExprPtr> upper;
        if (integral->lower.has_value()) {
            auto expanded = expand_expr_impl(*integral->lower, ctx);
            if (expanded.is_error()) {
                return fail<ExprPtr>(expanded.error());
            }
            lower = expanded.value();
        }
        if (integral->upper.has_value()) {
            auto expanded = expand_expr_impl(*integral->upper, ctx);
            if (expanded.is_error()) {
                return fail<ExprPtr>(expanded.error());
            }
            upper = expanded.value();
        }
        return simplify_expr(
            ctx.arena().make<Integral>(integrand.value(), integral->variable, std::move(lower), std::move(upper)),
            ctx);
    }

    if (const auto* derivative = expr_cast<Derivative>(expr)) {
        auto expanded = expand_expr_impl(derivative->expression, ctx);
        if (expanded.is_error()) {
            return fail<ExprPtr>(expanded.error());
        }
        return simplify_expr(
            ctx.arena().make<Derivative>(expanded.value(), derivative->variable, derivative->order),
            ctx);
    }

    if (const auto* limit = expr_cast<Limit>(expr)) {
        auto expression = expand_expr_impl(limit->expression, ctx);
        if (expression.is_error()) {
            return fail<ExprPtr>(expression.error());
        }
        auto point = expand_expr_impl(limit->point, ctx);
        if (point.is_error()) {
            return fail<ExprPtr>(point.error());
        }
        return simplify_expr(
            ctx.arena().make<Limit>(expression.value(), limit->variable, point.value(), limit->direction),
            ctx);
    }

    if (const auto* root = expr_cast<RootOf>(expr)) {
        auto polynomial = expand_expr_impl(root->polynomial, ctx);
        if (polynomial.is_error()) {
            return fail<ExprPtr>(polynomial.error());
        }
        return simplify_expr(
            ctx.arena().make<RootOf>(polynomial.value(), root->variable, root->root_index),
            ctx);
    }

    if (const auto* matrix = expr_cast<Matrix>(expr)) {
        std::vector<ExprPtr> elements;
        elements.reserve(matrix->elements.size());
        for (ExprPtr element : matrix->elements) {
            auto expanded = expand_expr_impl(element, ctx);
            if (expanded.is_error()) {
                return fail<ExprPtr>(expanded.error());
            }
            elements.push_back(expanded.value());
        }
        return simplify_expr(ctx.arena().make<Matrix>(matrix->rows, matrix->cols, std::move(elements)), ctx);
    }

    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Tipo di espressione non supportato in expand"));
}

Result<ExprPtr> expand(ExprPtr expr, symbolic::CASContext& ctx) {
    return expand_expr_impl(expr, ctx);
}

Result<ExprPtr> collect(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "collect richiede un'espressione non nulla"));
    }

    auto poly = parse_polynomial(expr, var, ctx);
    if (poly.is_error()) {
        return fail<ExprPtr>(poly.error());
    }

    return polynomial_to_expr(poly.value(), var, ctx);
}

} // namespace cas::algebra
