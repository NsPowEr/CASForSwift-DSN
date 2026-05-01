#include "integrate_engine.hpp"

#include "cas/error.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus::integrate_detail {

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

[[nodiscard]] ExprPtr make_rational(AstArena& arena, long long numerator, long long denominator) {
    return arena.make<RationalLit>(BigInt(numerator), BigInt(denominator));
}

[[nodiscard]] ExprPtr make_rational(AstArena& arena, const Rational& value) {
    if (value.denominator() == BigInt(1)) {
        return arena.make<IntegerLit>(value.numerator());
    }
    return arena.make<RationalLit>(value.numerator(), value.denominator());
}

[[nodiscard]] ExprPtr make_unary(AstArena& arena, UnaryOp op, ExprPtr operand) {
    return arena.make<Unary>(op, operand);
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

[[nodiscard]] ExprPtr make_function(AstArena& arena, std::string name, std::vector<ExprPtr> args) {
    return arena.make<FuncCall>(std::move(name), std::move(args));
}

[[nodiscard]] std::string canonical_function_name(const std::string& name) {
    if (name == "asin") {
        return "arcsin";
    }
    if (name == "acos") {
        return "arccos";
    }
    if (name == "atan") {
        return "arctan";
    }
    return name;
}

[[nodiscard]] bool depends_on(ExprPtr expr, const Symbol& var) {
    if (!expr) {
        return false;
    }

    if (const auto* symbol = expr_cast<Symbol>(expr)) {
        return symbol->name == var.name;
    }
    if (expr_is<IntegerLit>(expr) || expr_is<RationalLit>(expr) || expr_is<Constant>(expr) || expr_is<DecimalLit>(expr)) {
        return false;
    }
    if (const auto* unary = expr_cast<Unary>(expr)) {
        return depends_on(unary->operand, var);
    }
    if (const auto* binary = expr_cast<Binary>(expr)) {
        return depends_on(binary->left, var) || depends_on(binary->right, var);
    }
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        for (ExprPtr arg : call->args) {
            if (depends_on(arg, var)) {
                return true;
            }
        }
        return false;
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr term : sum->terms) {
            if (depends_on(term, var)) {
                return true;
            }
        }
        return false;
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr factor : product->factors) {
            if (depends_on(factor, var)) {
                return true;
            }
        }
        return false;
    }
    if (const auto* derivative = expr_cast<Derivative>(expr)) {
        return depends_on(derivative->expression, var);
    }
    if (const auto* integral = expr_cast<Integral>(expr)) {
        return integral->variable.name != var.name && depends_on(integral->integrand, var);
    }
    if (const auto* limit = expr_cast<Limit>(expr)) {
        return depends_on(limit->expression, var) || depends_on(limit->point, var);
    }
    if (const auto* root = expr_cast<RootOf>(expr)) {
        return depends_on(root->polynomial, var);
    }
    if (const auto* matrix = expr_cast<Matrix>(expr)) {
        for (ExprPtr element : matrix->elements) {
            if (depends_on(element, var)) {
                return true;
            }
        }
        return false;
    }

    return false;
}

[[nodiscard]] bool is_same_symbol(ExprPtr expr, const Symbol& var) {
    const auto* symbol = expr_cast<Symbol>(expr);
    return symbol != nullptr && symbol->name == var.name;
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

[[nodiscard]] bool is_negative_one(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value == BigInt(-1);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator == BigInt(-1) && rational->denominator == BigInt(1);
    }
    return false;
}

[[nodiscard]] std::optional<Rational> exact_scalar_from_expr(ExprPtr expr) {
    if (!expr) {
        return std::nullopt;
    }

    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return Rational(integer->value);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return Rational(rational->numerator, rational->denominator);
    }
    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op == UnaryOp::Neg) {
            if (auto inner = exact_scalar_from_expr(unary->operand); inner.has_value()) {
                return -(*inner);
            }
        }
    }

    return std::nullopt;
}

[[nodiscard]] std::optional<AffineArgument> extract_affine_argument(ExprPtr expr, const Symbol& var) {
    if (!expr) {
        return std::nullopt;
    }

    if (is_same_symbol(expr, var)) {
        return AffineArgument{Rational(BigInt(1)), Rational(BigInt(0))};
    }
    if (const auto exact = exact_scalar_from_expr(expr); exact.has_value()) {
        return AffineArgument{Rational(BigInt(0)), *exact};
    }
    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op == UnaryOp::Neg) {
            if (auto inner = extract_affine_argument(unary->operand, var); inner.has_value()) {
                return AffineArgument{-inner->coefficient, -inner->constant};
            }
        }
        return std::nullopt;
    }
    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Add || binary->op == BinaryOp::Sub) {
            auto left = extract_affine_argument(binary->left, var);
            auto right = extract_affine_argument(binary->right, var);
            if (!left.has_value() || !right.has_value()) {
                return std::nullopt;
            }
            if (binary->op == BinaryOp::Add) {
                return AffineArgument{
                    left->coefficient + right->coefficient,
                    left->constant + right->constant,
                };
            }
            return AffineArgument{
                left->coefficient - right->coefficient,
                left->constant - right->constant,
            };
        }
        if (binary->op == BinaryOp::Mul) {
            const auto left_scalar = exact_scalar_from_expr(binary->left);
            const auto right_scalar = exact_scalar_from_expr(binary->right);
            const auto left_affine = extract_affine_argument(binary->left, var);
            const auto right_affine = extract_affine_argument(binary->right, var);

            if (left_scalar.has_value() && right_affine.has_value()) {
                if (right_affine->constant != Rational(BigInt(0))) {
                    return std::nullopt;
                }
                return AffineArgument{(*left_scalar) * right_affine->coefficient, Rational(BigInt(0))};
            }
            if (right_scalar.has_value() && left_affine.has_value()) {
                if (left_affine->constant != Rational(BigInt(0))) {
                    return std::nullopt;
                }
                return AffineArgument{left_affine->coefficient * (*right_scalar), Rational(BigInt(0))};
            }
            return std::nullopt;
        }
        if (binary->op == BinaryOp::Div) {
            const auto right_scalar = exact_scalar_from_expr(binary->right);
            const auto left_affine = extract_affine_argument(binary->left, var);
            if (!right_scalar.has_value() || !left_affine.has_value() || right_scalar->numerator() == BigInt(0)) {
                return std::nullopt;
            }
            return AffineArgument{
                left_affine->coefficient / (*right_scalar),
                left_affine->constant / (*right_scalar),
            };
        }
        return std::nullopt;
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        Rational coefficient(BigInt(0));
        Rational constant(BigInt(0));
        for (ExprPtr term : sum->terms) {
            auto term_affine = extract_affine_argument(term, var);
            if (!term_affine.has_value()) {
                return std::nullopt;
            }
            coefficient += term_affine->coefficient;
            constant += term_affine->constant;
        }
        return AffineArgument{coefficient, constant};
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        Rational scalar(BigInt(1));
        std::optional<AffineArgument> variable_factor;
        for (ExprPtr factor : product->factors) {
            if (const auto exact = exact_scalar_from_expr(factor); exact.has_value()) {
                scalar *= *exact;
                continue;
            }

            auto factor_affine = extract_affine_argument(factor, var);
            if (!factor_affine.has_value() ||
                factor_affine->coefficient.numerator() == BigInt(0) ||
                variable_factor.has_value()) {
                return std::nullopt;
            }
            variable_factor = *factor_affine;
        }

        if (!variable_factor.has_value()) {
            return AffineArgument{Rational(BigInt(0)), scalar};
        }
        return AffineArgument{
            scalar * variable_factor->coefficient,
            scalar * variable_factor->constant,
        };
    }

    return std::nullopt;
}

[[nodiscard]] bool matches_square_of_variable(ExprPtr expr, const Symbol& var) {
    const auto* power = expr_cast<Binary>(expr);
    if (power == nullptr || power->op != BinaryOp::Pow || !is_same_symbol(power->left, var)) {
        return false;
    }
    return expr_is<IntegerLit>(power->right) && expr_ref<IntegerLit>(power->right).value == BigInt(2);
}

[[nodiscard]] bool matches_negative_square_of_variable(ExprPtr expr, const Symbol& var) {
    if (const auto* unary = expr_cast<Unary>(expr)) {
        return unary->op == UnaryOp::Neg && matches_square_of_variable(unary->operand, var);
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        return product->factors.size() == 2U &&
            ((is_negative_one(product->factors[0]) && matches_square_of_variable(product->factors[1], var)) ||
                (is_negative_one(product->factors[1]) && matches_square_of_variable(product->factors[0], var)));
    }
    return false;
}

[[nodiscard]] bool extract_squared_base(ExprPtr expr, ExprPtr& base) {
    const auto* power = expr_cast<Binary>(expr);
    if (power == nullptr || power->op != BinaryOp::Pow) {
        return false;
    }
    if (!expr_is<IntegerLit>(power->right) || expr_ref<IntegerLit>(power->right).value != BigInt(2)) {
        return false;
    }
    base = power->left;
    return true;
}

[[nodiscard]] bool extract_negative_squared_base(ExprPtr expr, ExprPtr& base) {
    if (const auto* unary = expr_cast<Unary>(expr)) {
        return unary->op == UnaryOp::Neg && extract_squared_base(unary->operand, base);
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        if (product->factors.size() != 2U) {
            return false;
        }
        if (is_negative_one(product->factors[0]) && extract_squared_base(product->factors[1], base)) {
            return true;
        }
        if (is_negative_one(product->factors[1]) && extract_squared_base(product->factors[0], base)) {
            return true;
        }
    }

    return false;
}

[[nodiscard]] bool matches_square_plus_constant_square(ExprPtr expr, const Symbol& var, ExprPtr& constant_base) {
    auto matches_pair = [&](ExprPtr lhs, ExprPtr rhs) {
        ExprPtr candidate{};
        return matches_square_of_variable(lhs, var) &&
            extract_squared_base(rhs, candidate) &&
            !depends_on(candidate, var) &&
            (constant_base = candidate, true);
    };

    if (const auto* binary = expr_cast<Binary>(expr)) {
        return binary->op == BinaryOp::Add &&
            (matches_pair(binary->left, binary->right) || matches_pair(binary->right, binary->left));
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        if (sum->terms.size() != 2U) {
            return false;
        }
        return matches_pair(sum->terms[0], sum->terms[1]) || matches_pair(sum->terms[1], sum->terms[0]);
    }

    return false;
}

[[nodiscard]] bool matches_square_minus_constant_square(ExprPtr expr, const Symbol& var, ExprPtr& constant_base) {
    auto matches_pair = [&](ExprPtr lhs, ExprPtr rhs) {
        ExprPtr candidate{};
        return matches_square_of_variable(lhs, var) &&
            extract_negative_squared_base(rhs, candidate) &&
            !depends_on(candidate, var) &&
            (constant_base = candidate, true);
    };

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Sub) {
            ExprPtr candidate{};
            if (matches_square_of_variable(binary->left, var) &&
                extract_squared_base(binary->right, candidate) &&
                !depends_on(candidate, var)) {
                constant_base = candidate;
                return true;
            }
        }
        if (binary->op == BinaryOp::Add) {
            return matches_pair(binary->left, binary->right) || matches_pair(binary->right, binary->left);
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        if (sum->terms.size() != 2U) {
            return false;
        }
        return matches_pair(sum->terms[0], sum->terms[1]) || matches_pair(sum->terms[1], sum->terms[0]);
    }

    return false;
}

[[nodiscard]] bool matches_one_plus_square(ExprPtr expr, const Symbol& var) {
    if (const auto* sum = expr_cast<Sum>(expr)) {
        bool found_one = false;
        bool found_square = false;
        for (ExprPtr term : sum->terms) {
            found_one = found_one || is_one(term);
            found_square = found_square || matches_square_of_variable(term, var);
        }
        return found_one && found_square && sum->terms.size() == 2U;
    }
    if (const auto* binary = expr_cast<Binary>(expr)) {
        return binary->op == BinaryOp::Add &&
            ((is_one(binary->left) && matches_square_of_variable(binary->right, var)) ||
                (matches_square_of_variable(binary->left, var) && is_one(binary->right)));
    }
    return false;
}

[[nodiscard]] bool matches_one_minus_square(ExprPtr expr, const Symbol& var) {
    if (const auto* sum = expr_cast<Sum>(expr)) {
        bool found_one = false;
        bool found_neg_square = false;
        for (ExprPtr term : sum->terms) {
            if (is_one(term)) {
                found_one = true;
                continue;
            }
            if (matches_negative_square_of_variable(term, var)) {
                found_neg_square = true;
            }
        }
        return found_one && found_neg_square && sum->terms.size() == 2U;
    }
    if (const auto* binary = expr_cast<Binary>(expr)) {
        return binary->op == BinaryOp::Sub && is_one(binary->left) && matches_square_of_variable(binary->right, var);
    }
    return false;
}

[[nodiscard]] bool is_rational_value(ExprPtr expr, long long numerator, long long denominator) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return denominator == 1 && integer->value == BigInt(numerator);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator == BigInt(numerator) && rational->denominator == BigInt(denominator);
    }
    return false;
}

[[nodiscard]] bool matches_constant_square_minus_variable_square(ExprPtr expr, const Symbol& var, ExprPtr& constant_base) {
    auto matches_pair = [&](ExprPtr const_term, ExprPtr var_term) {
        ExprPtr candidate{};
        return extract_squared_base(const_term, candidate) &&
            !depends_on(candidate, var) &&
            matches_negative_square_of_variable(var_term, var) &&
            (constant_base = candidate, true);
    };

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Sub) {
            ExprPtr candidate{};
            if (extract_squared_base(binary->left, candidate) &&
                !depends_on(candidate, var) &&
                matches_square_of_variable(binary->right, var)) {
                constant_base = candidate;
                return true;
            }
        }
        if (binary->op == BinaryOp::Add) {
            return matches_pair(binary->left, binary->right) || matches_pair(binary->right, binary->left);
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        if (sum->terms.size() != 2U) return false;
        return matches_pair(sum->terms[0], sum->terms[1]) || matches_pair(sum->terms[1], sum->terms[0]);
    }

    return false;
}

[[nodiscard]] bool matches_reciprocal_sqrt_one_minus_square(ExprPtr expr, const Symbol& var) {
    if (const auto* quotient = expr_cast<Binary>(expr);
        quotient != nullptr && quotient->op == BinaryOp::Div && is_one(quotient->left)) {
        const auto* sqrt_call = expr_cast<FuncCall>(quotient->right);
        return sqrt_call != nullptr &&
            sqrt_call->name == "sqrt" &&
            sqrt_call->args.size() == 1U &&
            matches_one_minus_square(sqrt_call->args.front(), var);
    }

    const auto* power = expr_cast<Binary>(expr);
    if (power == nullptr || power->op != BinaryOp::Pow) {
        return false;
    }

    if (is_rational_value(power->right, -1, 2)) {
        return matches_one_minus_square(power->left, var);
    }

    if (!is_rational_value(power->right, -1, 1)) {
        return false;
    }

    const auto* sqrt_call = expr_cast<FuncCall>(power->left);
    return sqrt_call != nullptr &&
        sqrt_call->name == "sqrt" &&
        sqrt_call->args.size() == 1U &&
        matches_one_minus_square(sqrt_call->args.front(), var);
}

}  // namespace cas::calculus::integrate_detail
