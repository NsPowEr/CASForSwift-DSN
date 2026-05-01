#include "calculus_internal.hpp"

#include "cas/error.hpp"
#include "cas/rational.hpp"

#include <algorithm>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus {
namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

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

[[nodiscard]] ExprPtr make_rational_expr(AstArena& arena, const Rational& value) {
    if (value.is_integer()) {
        return arena.make<IntegerLit>(value.numerator());
    }
    return arena.make<RationalLit>(value.numerator(), value.denominator());
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

[[nodiscard]] int infinity_sign(ExprPtr expr) {
    const auto* constant = expr_cast<Constant>(expr);
    if (constant != nullptr && constant->value == MathConstant::Infinity) {
        return 1;
    }

    const auto* unary = expr_cast<Unary>(expr);
    if (unary != nullptr &&
        unary->op == UnaryOp::Neg &&
        expr_is<Constant>(unary->operand) &&
        expr_ref<Constant>(unary->operand).value == MathConstant::Infinity) {
        return -1;
    }

    return 0;
}

[[nodiscard]] ExprPtr make_signed_infinity(AstArena& arena, int sign) {
    ExprPtr infinity = arena.make<Constant>(MathConstant::Infinity);
    if (sign < 0) {
        return arena.make<Unary>(UnaryOp::Neg, infinity);
    }
    return infinity;
}

[[nodiscard]] bool is_same_symbol(ExprPtr expr, const Symbol& var) {
    const auto* symbol = expr_cast<Symbol>(expr);
    return symbol != nullptr && symbol->name == var.name;
}

[[nodiscard]] std::optional<Rational> rational_from_expr(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return Rational(integer->value);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return Rational(rational->numerator, rational->denominator);
    }
    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op != UnaryOp::Neg) {
            return std::nullopt;
        }
        auto operand = rational_from_expr(unary->operand);
        if (!operand.has_value()) {
            return std::nullopt;
        }
        return -operand.value();
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<int> sign_of_variable_power_at_infinity(ExprPtr expr, const Symbol& var, int point_sign) {
    if (point_sign == 0) {
        return std::nullopt;
    }

    if (is_same_symbol(expr, var)) {
        return point_sign;
    }

    const auto* power = expr_cast<Binary>(expr);
    if (power == nullptr || power->op != BinaryOp::Pow || !is_same_symbol(power->left, var)) {
        return std::nullopt;
    }

    auto exponent = rational_from_expr(power->right);
    if (!exponent.has_value() || !exponent->is_integer() || exponent->numerator().is_negative()) {
        return std::nullopt;
    }

    if (point_sign > 0) {
        return 1;
    }

    return (exponent->numerator() % BigInt(2)).is_zero() ? 1 : -1;
}

[[nodiscard]] bool matches_exp_of_variable(ExprPtr expr, const Symbol& var) {
    const auto* call = expr_cast<FuncCall>(expr);
    return call != nullptr && call->name == "exp" && call->args.size() == 1U && is_same_symbol(call->args.front(), var);
}

[[nodiscard]] bool matches_log_of_variable(ExprPtr expr, const Symbol& var) {
    const auto* call = expr_cast<FuncCall>(expr);
    return call != nullptr && call->name == "ln" && call->args.size() == 1U && is_same_symbol(call->args.front(), var);
}

[[nodiscard]] bool matches_inverse_variable_with_constant(ExprPtr expr, const Symbol& var, Rational& coefficient) {
    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Div && is_same_symbol(binary->right, var)) {
            auto constant = rational_from_expr(binary->left);
            if (!constant.has_value()) {
                return false;
            }
            coefficient = constant.value();
            return true;
        }
        if (binary->op == BinaryOp::Pow && is_same_symbol(binary->left, var)) {
            auto exponent = rational_from_expr(binary->right);
            if (exponent.has_value() && exponent.value() == Rational(BigInt(-1))) {
                coefficient = Rational(BigInt(1));
                return true;
            }
        }
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op != UnaryOp::Neg) {
            return false;
        }
        Rational inner(BigInt(0));
        if (!matches_inverse_variable_with_constant(unary->operand, var, inner)) {
            return false;
        }
        coefficient = -inner;
        return true;
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        if (product->factors.size() != 2U) {
            return false;
        }

        Rational constant(BigInt(0));
        if (auto value = rational_from_expr(product->factors[0]); value.has_value() &&
            matches_inverse_variable_with_constant(product->factors[1], var, constant)) {
            coefficient = value.value() * constant;
            return true;
        }
        if (auto value = rational_from_expr(product->factors[1]); value.has_value() &&
            matches_inverse_variable_with_constant(product->factors[0], var, constant)) {
            coefficient = value.value() * constant;
            return true;
        }
    }

    return false;
}

[[nodiscard]] bool matches_one_plus_inverse_variable(ExprPtr expr, const Symbol& var, Rational& coefficient) {
    auto matches_pair = [&](ExprPtr lhs, ExprPtr rhs) {
        return is_one(lhs) && matches_inverse_variable_with_constant(rhs, var, coefficient);
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

}  // namespace

ExprPtr limit_make_integer(AstArena& arena, long long value) {
    return make_integer(arena, value);
}

ExprPtr limit_make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs) {
    return make_binary(arena, op, lhs, rhs);
}

bool limit_is_zero(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value == BigInt(0);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator == BigInt(0);
    }
    return false;
}

bool limit_is_infinity(ExprPtr expr) {
    const auto* constant = expr_cast<Constant>(expr);
    if (constant != nullptr && constant->value == MathConstant::Infinity) {
        return true;
    }

    const auto* unary = expr_cast<Unary>(expr);
    return unary != nullptr &&
        unary->op == UnaryOp::Neg &&
        expr_is<Constant>(unary->operand) &&
        expr_ref<Constant>(unary->operand).value == MathConstant::Infinity;
}

std::optional<QuotientView> extract_quotient_view(ExprPtr expr, AstArena& arena) {
    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Div) {
            return QuotientView{.numerator = binary->left, .denominator = binary->right};
        }

        if (binary->op == BinaryOp::Pow) {
            auto exponent = rational_from_expr(binary->right);
            if (!exponent.has_value() || !exponent->is_integer() || !exponent->numerator().is_negative()) {
                return std::nullopt;
            }

            const BigInt positive_power = -exponent->numerator();
            ExprPtr denominator = positive_power == BigInt(1)
                ? binary->left
                : make_power(arena, binary->left, arena.make<IntegerLit>(positive_power));
            return QuotientView{.numerator = make_integer(arena, 1), .denominator = denominator};
        }
    }

    const auto* product = expr_cast<Product>(expr);
    if (product == nullptr) {
        if (const auto* unary = expr_cast<Unary>(expr)) {
            if (unary->op == UnaryOp::Neg) {
                auto inner = extract_quotient_view(unary->operand, arena);
                if (inner.has_value()) {
                    return QuotientView{
                        .numerator = arena.make<Unary>(UnaryOp::Neg, inner->numerator),
                        .denominator = inner->denominator,
                    };
                }
            }
        }

        if (const auto* sum = expr_cast<Sum>(expr)) {
            std::vector<ExprPtr> numerator_terms;
            ExprPtr common_denominator{};
            numerator_terms.reserve(sum->terms.size());

            for (ExprPtr term : sum->terms) {
                auto term_quotient = extract_quotient_view(term, arena);
                if (!term_quotient.has_value()) {
                    return std::nullopt;
                }

                if (!common_denominator) {
                    common_denominator = term_quotient->denominator;
                } else if (!structural_equal(common_denominator, term_quotient->denominator)) {
                    return std::nullopt;
                }
                numerator_terms.push_back(term_quotient->numerator);
            }

            if (common_denominator) {
                return QuotientView{
                    .numerator = make_sum(arena, std::move(numerator_terms)),
                    .denominator = common_denominator,
                };
            }
        }

        return std::nullopt;
    }

    std::vector<ExprPtr> numerator_factors;
    std::vector<ExprPtr> denominator_factors;
    numerator_factors.reserve(product->factors.size());
    denominator_factors.reserve(product->factors.size());

    for (ExprPtr factor : product->factors) {
        const auto* power = expr_cast<Binary>(factor);
        if (power == nullptr || power->op != BinaryOp::Pow) {
            numerator_factors.push_back(factor);
            continue;
        }

        auto exponent = rational_from_expr(power->right);
        if (!exponent.has_value() || !exponent->is_integer() || !exponent->numerator().is_negative()) {
            numerator_factors.push_back(factor);
            continue;
        }

        const BigInt positive_power = -exponent->numerator();
        denominator_factors.push_back(positive_power == BigInt(1)
            ? power->left
            : make_power(arena, power->left, arena.make<IntegerLit>(positive_power)));
    }

    if (denominator_factors.empty()) {
        return std::nullopt;
    }

    return QuotientView{
        .numerator = make_product(arena, std::move(numerator_factors)),
        .denominator = make_product(arena, std::move(denominator_factors)),
    };
}

Result<ExprPtr> try_infinite_limit(ExprPtr expr, const Symbol& var, ExprPtr point, AstArena& arena) {
    const int point_sign = infinity_sign(point);
    if (point_sign == 0) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Il punto del limite non e' infinito"));
    }

    if (const auto* power = expr_cast<Binary>(expr)) {
        if (power->op == BinaryOp::Pow && is_same_symbol(power->right, var)) {
            Rational coefficient(BigInt(0));
            if (matches_one_plus_inverse_variable(power->left, var, coefficient)) {
                if (coefficient == Rational(BigInt(1))) {
                    return ok(arena.make<Constant>(MathConstant::E));
                }
                return ok(make_function(arena, "exp", {make_rational_expr(arena, coefficient)}));
            }
        }
    }

    auto quotient = extract_quotient_view(expr, arena);
    if (!quotient.has_value()) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Limite all'infinito non riconosciuto"));
    }

    auto left_power_sign = sign_of_variable_power_at_infinity(quotient->numerator, var, point_sign);
    auto right_power_sign = sign_of_variable_power_at_infinity(quotient->denominator, var, point_sign);

    if (left_power_sign.has_value() && matches_exp_of_variable(quotient->denominator, var)) {
        if (point_sign > 0) {
            return ok(make_integer(arena, 0));
        }
        return ok(make_signed_infinity(arena, left_power_sign.value()));
    }

    if (matches_exp_of_variable(quotient->numerator, var) && right_power_sign.has_value()) {
        if (point_sign > 0) {
            return ok(arena.make<Constant>(MathConstant::Infinity));
        }
        return ok(make_integer(arena, 0));
    }

    if (matches_log_of_variable(quotient->numerator, var) && right_power_sign.has_value()) {
        if (point_sign < 0) {
            return fail<ExprPtr>(make_error(CASErrorKind::Undefined, "ln(x) non e' definito per x -> -inf"));
        }
        return ok(make_integer(arena, 0));
    }

    if (left_power_sign.has_value() && matches_log_of_variable(quotient->denominator, var)) {
        if (point_sign < 0) {
            return fail<ExprPtr>(make_error(CASErrorKind::Undefined, "ln(x) non e' definito per x -> -inf"));
        }
        return ok(arena.make<Constant>(MathConstant::Infinity));
    }

    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Limite all'infinito non implementato per questa forma"));
}

}  // namespace cas::calculus
