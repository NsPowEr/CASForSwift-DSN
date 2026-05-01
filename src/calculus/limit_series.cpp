#include "cas/calculus.hpp"
#include "calculus_internal.hpp"
#include "cas/rational.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus {
namespace {

[[nodiscard]] ExprPtr make_integer(AstArena& arena, long long value) {
    return arena.make<IntegerLit>(BigInt(value));
}

[[nodiscard]] ExprPtr make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs) {
    return arena.make<Binary>(op, lhs, rhs);
}

[[nodiscard]] ExprPtr make_sum(AstArena& arena, std::vector<ExprPtr> terms) {
    if (terms.empty()) {
        return make_integer(arena, 0);
    }
    if (terms.size() == 1U) {
        return terms.front();
    }
    return arena.make<Sum>(std::move(terms));
}

[[nodiscard]] ExprPtr make_product(AstArena& arena, std::vector<ExprPtr> factors) {
    if (factors.empty()) {
        return make_integer(arena, 1);
    }
    if (factors.size() == 1U) {
        return factors.front();
    }
    return arena.make<Product>(std::move(factors));
}

[[nodiscard]] ExprPtr make_power(AstArena& arena, ExprPtr base, ExprPtr exponent) {
    return arena.make<Binary>(BinaryOp::Pow, base, exponent);
}

[[nodiscard]] ExprPtr make_function(AstArena& arena, std::string name, std::vector<ExprPtr> args) {
    return arena.make<FuncCall>(std::move(name), std::move(args));
}

[[nodiscard]] bool is_zero(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value == BigInt(0);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator == BigInt(0);
    }
    return false;
}

[[nodiscard]] bool is_one(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value == BigInt(1);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator == BigInt(1) && rational->denominator == BigInt(1);
    }
    return false;
}

[[nodiscard]] bool is_same_symbol(ExprPtr expr, const Symbol& var) {
    const auto* symbol = expr_cast<Symbol>(expr);
    return symbol != nullptr && symbol->name == var.name;
}

[[nodiscard]] BigInt factorial(unsigned int value) {
    BigInt result(1);
    for (unsigned int i = 2; i <= value; ++i) {
        result = result * BigInt(static_cast<long long>(i));
    }
    return result;
}

[[nodiscard]] ExprPtr make_rational_expr(AstArena& arena, const Rational& value) {
    if (value.is_integer()) {
        return arena.make<IntegerLit>(value.numerator());
    }
    return arena.make<RationalLit>(value.numerator(), value.denominator());
}

[[nodiscard]] bool matches_one_plus_variable(ExprPtr expr, const Symbol& var) {
    auto matches_pair = [&](ExprPtr lhs, ExprPtr rhs) {
        return is_one(lhs) && is_same_symbol(rhs, var);
    };

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op != BinaryOp::Add) {
            return false;
        }
        return matches_pair(binary->left, binary->right) || matches_pair(binary->right, binary->left);
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        if (sum->terms.size() != 2U) {
            return false;
        }
        return matches_pair(sum->terms[0], sum->terms[1]) || matches_pair(sum->terms[1], sum->terms[0]);
    }

    return false;
}

[[nodiscard]] bool matches_named_function_of_variable(
    ExprPtr expr,
    const Symbol& var,
    const std::string& primary_name,
    const std::optional<std::string>& alias_name = std::nullopt) {
    const auto* call = expr_cast<FuncCall>(expr);
    if (call == nullptr || call->args.size() != 1U || !is_same_symbol(call->args.front(), var)) {
        return false;
    }
    return call->name == primary_name || (alias_name.has_value() && call->name == alias_name.value());
}

void append_series_term(
    std::vector<ExprPtr>& terms,
    AstArena& arena,
    ExprPtr variable_expr,
    unsigned int degree,
    ExprPtr coefficient_expr) {
    if (degree == 0U) {
        terms.push_back(coefficient_expr);
        return;
    }

    ExprPtr power = degree == 1U
        ? variable_expr
        : make_power(arena, variable_expr, make_integer(arena, static_cast<long long>(degree)));
    terms.push_back(make_product(arena, {coefficient_expr, power}));
}

void append_series_term(
    std::vector<ExprPtr>& terms,
    AstArena& arena,
    ExprPtr variable_expr,
    unsigned int degree,
    const Rational& coefficient) {
    if (coefficient.numerator().is_zero()) {
        return;
    }
    append_series_term(terms, arena, variable_expr, degree, make_rational_expr(arena, coefficient));
}

// Maclaurin series for sin(x)/x (sinc): 1 - x^2/6 + x^4/120 - ...
// Matched on original expr (before simplification) to avoid product-form noise.
[[nodiscard]] static std::optional<TaylorExpansion> try_sinc_series(
    ExprPtr expr, const Symbol& var, unsigned int order, AstArena& arena) {
    const auto* div = expr_cast<Binary>(expr);
    if (!div || div->op != BinaryOp::Div) return std::nullopt;
    if (!matches_named_function_of_variable(div->left, var, "sin")) return std::nullopt;
    if (!is_same_symbol(div->right, var)) return std::nullopt;

    ExprPtr variable_expr = arena.make<Symbol>(var);
    std::vector<ExprPtr> terms;
    terms.reserve(order / 2U + 1U);
    for (unsigned int degree = 0U; degree <= order; degree += 2U) {
        const unsigned int k = degree / 2U;
        const BigInt sign = (k % 2U == 0U) ? BigInt(1) : BigInt(-1);
        append_series_term(terms, arena, variable_expr, degree,
                           Rational(sign, factorial(degree + 1U)));
    }

    return TaylorExpansion{
        .polynomial = make_sum(arena, std::move(terms)),
        .remainder = make_function(arena, "O", {
            make_power(arena, variable_expr,
                make_integer(arena, static_cast<long long>(order + 1U))),
        }),
        .computed_order = order,
    };
}

[[nodiscard]] std::optional<TaylorExpansion> try_standard_maclaurin_series(
    ExprPtr expr,
    const Symbol& var,
    unsigned int order,
    AstArena& arena) {
    ExprPtr variable_expr = arena.make<Symbol>(var);
    std::vector<ExprPtr> terms;
    terms.reserve(order + 1U);

    if (matches_named_function_of_variable(expr, var, "exp")) {
        for (unsigned int degree = 0U; degree <= order; ++degree) {
            append_series_term(terms, arena, variable_expr, degree, Rational(BigInt(1), factorial(degree)));
        }
    } else if (matches_named_function_of_variable(expr, var, "sin")) {
        for (unsigned int degree = 1U; degree <= order; degree += 2U) {
            const unsigned int index = (degree - 1U) / 2U;
            const BigInt sign = (index % 2U) == 0U ? BigInt(1) : BigInt(-1);
            append_series_term(terms, arena, variable_expr, degree, Rational(sign, factorial(degree)));
        }
    } else if (matches_named_function_of_variable(expr, var, "cos")) {
        for (unsigned int degree = 0U; degree <= order; degree += 2U) {
            const unsigned int index = degree / 2U;
            const BigInt sign = (index % 2U) == 0U ? BigInt(1) : BigInt(-1);
            append_series_term(terms, arena, variable_expr, degree, Rational(sign, factorial(degree)));
        }
    } else if (const auto* call = expr_cast<FuncCall>(expr);
               call != nullptr && call->func_id == BuiltinOp::Ln && call->args.size() == 1U &&
               matches_one_plus_variable(call->args.front(), var)) {
        for (unsigned int degree = 1U; degree <= order; ++degree) {
            const BigInt sign = (degree % 2U) == 1U ? BigInt(1) : BigInt(-1);
            append_series_term(terms, arena, variable_expr, degree, Rational(sign, BigInt(static_cast<long long>(degree))));
        }
    } else if (matches_named_function_of_variable(expr, var, "arctan", std::string("atan"))) {
        for (unsigned int degree = 1U; degree <= order; degree += 2U) {
            const unsigned int index = (degree - 1U) / 2U;
            const BigInt sign = (index % 2U) == 0U ? BigInt(1) : BigInt(-1);
            append_series_term(terms, arena, variable_expr, degree, Rational(sign, BigInt(static_cast<long long>(degree))));
        }
    } else if (const auto* power = expr_cast<Binary>(expr);
               power != nullptr &&
               power->op == BinaryOp::Pow &&
               matches_one_plus_variable(power->left, var) &&
               !depends_on(power->right, var)) {
        append_series_term(terms, arena, variable_expr, 0U, make_integer(arena, 1));
        for (unsigned int degree = 1U; degree <= order; ++degree) {
            std::vector<ExprPtr> coefficient_factors;
            coefficient_factors.reserve(degree + 1U);
            coefficient_factors.push_back(make_rational_expr(arena, Rational(BigInt(1), factorial(degree))));

            for (unsigned int index = 0U; index < degree; ++index) {
                coefficient_factors.push_back(index == 0U
                    ? power->right
                    : make_binary(arena, BinaryOp::Sub, power->right, make_integer(arena, static_cast<long long>(index))));
            }
            append_series_term(terms, arena, variable_expr, degree, make_product(arena, std::move(coefficient_factors)));
        }
    } else {
        return std::nullopt;
    }

    return TaylorExpansion{
        .polynomial = make_sum(arena, std::move(terms)),
        .remainder = make_function(arena, "O", {
            make_power(arena, variable_expr, make_integer(arena, static_cast<long long>(order + 1U))),
        }),
        .computed_order = order,
    };
}

}  // namespace

Result<TaylorExpansion> taylor_series(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    unsigned int order,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    auto simplified_expr = ctx.simplify(expr);
    if (simplified_expr.is_error()) {
        return fail<TaylorExpansion>(simplified_expr.error());
    }

    auto simplified_point = ctx.simplify(point);
    if (simplified_point.is_error()) {
        return fail<TaylorExpansion>(simplified_point.error());
    }

    if (is_zero(simplified_point.value())) {
        // Try sinc (sin(x)/x) on the original expression before simplification
        auto sinc_series = try_sinc_series(expr, var, order, arena);
        if (!sinc_series.has_value()) sinc_series = try_sinc_series(simplified_expr.value(), var, order, arena);
        if (sinc_series.has_value()) {
            auto simplified_poly = ctx.simplify(sinc_series->polynomial);
            if (simplified_poly.is_error()) return fail<TaylorExpansion>(simplified_poly.error());
            sinc_series->polynomial = simplified_poly.value();
            auto mat_poly = symbolic::materialize_expr(sinc_series->polynomial, arena);
            if (mat_poly.is_error()) return fail<TaylorExpansion>(mat_poly.error());
            auto mat_rem = symbolic::materialize_expr(sinc_series->remainder, arena);
            if (mat_rem.is_error()) return fail<TaylorExpansion>(mat_rem.error());
            sinc_series->polynomial = mat_poly.value();
            sinc_series->remainder = mat_rem.value();
            return sinc_series.value();
        }

        auto standard_series = try_standard_maclaurin_series(simplified_expr.value(), var, order, arena);
        if (standard_series.has_value()) {
            auto simplified_polynomial = ctx.simplify(standard_series->polynomial);
            if (simplified_polynomial.is_error()) {
                return fail<TaylorExpansion>(simplified_polynomial.error());
            }
            standard_series->polynomial = simplified_polynomial.value();
            auto materialized_polynomial = symbolic::materialize_expr(standard_series->polynomial, arena);
            if (materialized_polynomial.is_error()) {
                return fail<TaylorExpansion>(materialized_polynomial.error());
            }
            auto materialized_remainder = symbolic::materialize_expr(standard_series->remainder, arena);
            if (materialized_remainder.is_error()) {
                return fail<TaylorExpansion>(materialized_remainder.error());
            }
            standard_series->polynomial = materialized_polynomial.value();
            standard_series->remainder = materialized_remainder.value();
            return standard_series.value();
        }
    }

    ExprPtr x = arena.make<Symbol>(var);
    ExprPtr delta = make_binary(arena, BinaryOp::Sub, x, simplified_point.value());
    std::vector<ExprPtr> terms;
    terms.reserve(order + 1U);

    for (unsigned int degree = 0; degree <= order; ++degree) {
        Result<ExprPtr> derivative = degree == 0U ? ok(simplified_expr.value()) : diff(simplified_expr.value(), var, degree, ctx);
        if (derivative.is_error()) {
            return fail<TaylorExpansion>(derivative.error());
        }

        auto coefficient = ctx.substitute(derivative.value(), var, simplified_point.value());
        if (coefficient.is_error()) {
            return fail<TaylorExpansion>(coefficient.error());
        }

        auto simplified_coefficient = ctx.simplify(coefficient.value());
        if (simplified_coefficient.is_error()) {
            return fail<TaylorExpansion>(simplified_coefficient.error());
        }

        ExprPtr term = simplified_coefficient.value();
        if (degree > 0U) {
            ExprPtr power = degree == 1U ? delta : make_power(arena, delta, make_integer(arena, static_cast<long long>(degree)));
            const BigInt denominator = factorial(degree);
            ExprPtr scaled_coefficient = denominator == BigInt(1)
                ? term
                : make_binary(arena, BinaryOp::Div, term, arena.make<RationalLit>(denominator, BigInt(1)));
            term = make_product(arena, {scaled_coefficient, power});
        }

        auto simplified_term = ctx.simplify(term);
        if (simplified_term.is_error()) {
            return fail<TaylorExpansion>(simplified_term.error());
        }
        terms.push_back(simplified_term.value());
    }

    auto polynomial = ctx.simplify(make_sum(arena, std::move(terms)));
    if (polynomial.is_error()) {
        return fail<TaylorExpansion>(polynomial.error());
    }

    ExprPtr remainder = make_function(arena, "O", {
        make_power(arena, delta, make_integer(arena, static_cast<long long>(order + 1U))),
    });

    auto materialized_polynomial = symbolic::materialize_expr(polynomial.value(), arena);
    if (materialized_polynomial.is_error()) {
        return fail<TaylorExpansion>(materialized_polynomial.error());
    }
    auto materialized_remainder = symbolic::materialize_expr(remainder, arena);
    if (materialized_remainder.is_error()) {
        return fail<TaylorExpansion>(materialized_remainder.error());
    }

    return ok(TaylorExpansion{
        .polynomial = materialized_polynomial.value(),
        .remainder = materialized_remainder.value(),
        .computed_order = order,
    });
}

}  // namespace cas::calculus
