#pragma once

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "cas/result.hpp"
#include "cas/error.hpp"
#include <string>
#include <vector>
#include <utility>

namespace cas::calculus {

[[nodiscard]] inline CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

[[nodiscard]] inline ExprPtr make_integer(AstArena& arena, long long value) {
    return arena.make<IntegerLit>(BigInt(value));
}

[[nodiscard]] inline ExprPtr make_constant(AstArena& arena, MathConstant value) {
    return arena.make<Constant>(value);
}

[[nodiscard]] inline ExprPtr make_unary(AstArena& arena, UnaryOp op, ExprPtr operand) {
    return arena.make<Unary>(op, operand);
}

[[nodiscard]] inline ExprPtr make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs) {
    return arena.make<Binary>(op, lhs, rhs);
}

[[nodiscard]] inline ExprPtr make_sum(AstArena& arena, std::vector<ExprPtr> terms) {
    if (terms.empty()) {
        return make_integer(arena, 0);
    }
    if (terms.size() == 1U) {
        return terms.front();
    }
    return arena.make<Sum>(std::move(terms));
}

[[nodiscard]] inline ExprPtr make_product(AstArena& arena, std::vector<ExprPtr> factors) {
    if (factors.empty()) {
        return make_integer(arena, 1);
    }
    if (factors.size() == 1U) {
        return factors.front();
    }
    return arena.make<Product>(std::move(factors));
}

[[nodiscard]] inline ExprPtr make_power(AstArena& arena, ExprPtr base, ExprPtr exponent) {
    return make_binary(arena, BinaryOp::Pow, base, exponent);
}

[[nodiscard]] inline ExprPtr make_function(AstArena& arena, std::string name, std::vector<ExprPtr> args) {
    return arena.make<FuncCall>(std::move(name), std::move(args));
}

[[nodiscard]] inline ExprPtr make_matrix(AstArena& arena, std::size_t rows, std::size_t cols, std::vector<ExprPtr> elements) {
    return arena.make<Matrix>(rows, cols, std::move(elements));
}

[[nodiscard]] inline bool is_exact_zero(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value.is_zero();
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator.is_zero();
    }
    return false;
}

class Differentiator {
public:
    explicit Differentiator(symbolic::CASContext& context) noexcept : context_(context), arena_(context.arena()) {}

    [[nodiscard]] Result<ExprPtr> differentiate(ExprPtr expr, const Symbol& var, unsigned int order);

private:
    struct Visitor;

    [[nodiscard]] Result<ExprPtr> differentiate_once(ExprPtr expr, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> differentiate_unary(const Unary& unary, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> differentiate_binary(const Binary& binary, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> differentiate_sum(const Sum& sum, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> differentiate_product(const Product& product, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> differentiate_power(const Binary& power, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> differentiate_function(const FuncCall& call, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> differentiate_integral(const Integral& integral, const Symbol& var);

    symbolic::CASContext& context_;
    AstArena& arena_;
};

}  // namespace cas::calculus
