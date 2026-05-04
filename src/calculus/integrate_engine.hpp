#pragma once

#include "calculus_internal.hpp"
#include "cas/rational.hpp"

#include <optional>
#include <string>
#include <vector>

namespace cas::calculus::integrate_detail {

struct AffineArgument {
    Rational coefficient;
    Rational constant;
};

struct QuadraticArgument {
    Rational quadratic;
    Rational linear;
    Rational constant;
};

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message);
[[nodiscard]] ExprPtr make_integer(AstArena& arena, long long value);
[[nodiscard]] ExprPtr make_rational(AstArena& arena, long long numerator, long long denominator);
[[nodiscard]] ExprPtr make_rational(AstArena& arena, const Rational& value);
[[nodiscard]] ExprPtr make_unary(AstArena& arena, UnaryOp op, ExprPtr operand);
[[nodiscard]] ExprPtr make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs);
[[nodiscard]] ExprPtr make_sum(AstArena& arena, std::vector<ExprPtr> terms);
[[nodiscard]] ExprPtr make_product(AstArena& arena, std::vector<ExprPtr> factors);
[[nodiscard]] ExprPtr make_function(AstArena& arena, std::string name, std::vector<ExprPtr> args);
[[nodiscard]] std::string canonical_function_name(const std::string& name);
[[nodiscard]] bool depends_on(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool is_same_symbol(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool is_one(ExprPtr expr);
[[nodiscard]] bool is_negative_one(ExprPtr expr);
[[nodiscard]] std::optional<Rational> exact_scalar_from_expr(ExprPtr expr);
[[nodiscard]] std::optional<AffineArgument> extract_affine_argument(ExprPtr expr, const Symbol& var);
[[nodiscard]] std::optional<QuadraticArgument> extract_quadratic_argument(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool matches_square_of_variable(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool matches_square_plus_constant_square(ExprPtr expr, const Symbol& var, ExprPtr& constant_base);
[[nodiscard]] bool matches_square_minus_constant_square(ExprPtr expr, const Symbol& var, ExprPtr& constant_base);
[[nodiscard]] bool matches_one_plus_square(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool matches_one_minus_square(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool is_rational_value(ExprPtr expr, long long numerator, long long denominator);
[[nodiscard]] bool matches_constant_square_minus_variable_square(ExprPtr expr, const Symbol& var, ExprPtr& constant_base);
[[nodiscard]] bool matches_reciprocal_sqrt_one_minus_square(ExprPtr expr, const Symbol& var);

class Integrator {
public:
    explicit Integrator(symbolic::CASContext& context) noexcept;
    [[nodiscard]] Result<ExprPtr> integrate(ExprPtr expr, const Symbol& var);

    [[nodiscard]] Result<bool> expressions_match_after_simplify(ExprPtr lhs, ExprPtr rhs);
    [[nodiscard]] Result<ExprPtr> try_u_substitution_for_product(const Product& product, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_rational(ExprPtr expr, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_via_partial_fractions(ExprPtr expr, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_once(ExprPtr expr, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_binary(const Binary& binary, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_linear_over_quadratic(const Binary& quotient, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_sqrt_quadratic(ExprPtr radicand, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_product(const Product& product, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_power(const Binary& power, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_function_direct(const std::string& name, ExprPtr argument);
    [[nodiscard]] Result<ExprPtr> integrate_power_direct(ExprPtr base, ExprPtr exponent, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_function(const FuncCall& call, const Symbol& var);

private:
    symbolic::CASContext& context_;
    AstArena& arena_;
    
    struct DepthGuard {
        std::size_t& depth;
        explicit DepthGuard(std::size_t& d) : depth(d) { ++depth; }
        ~DepthGuard() { --depth; }
    };
    static thread_local std::size_t depth_;
};

[[nodiscard]] Result<ExprPtr> integrate_indefinite_impl(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx);

}  // namespace cas::calculus::integrate_detail
