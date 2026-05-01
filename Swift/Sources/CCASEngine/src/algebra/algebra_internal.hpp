#pragma once

#include "cas/algebra.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include "polynomial_internal.hpp"
#include <optional>
#include <string>
#include <vector>

namespace cas::algebra {

// Utility functions
[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message, std::optional<std::string> hint = std::nullopt);
[[nodiscard]] ExprPtr make_rational_expr(AstArena& arena, const Rational& value);

template <typename T>
[[nodiscard]] Result<T> fail_unimplemented(std::string operation, std::string detail) {
    return fail<T>(make_error(
        CASErrorKind::Unimplemented,
        std::move(operation) + " non e' ancora implementata nel modulo algebra",
        std::move(detail)));
}

[[nodiscard]] ExprPtr make_integer(AstArena& arena, long long value);
[[nodiscard]] bool is_zero_expr(ExprPtr expr);
[[nodiscard]] bool is_one_expr(ExprPtr expr);
[[nodiscard]] bool contains_decimal_literal(ExprPtr expr);

[[nodiscard]] Result<ExprPtr> simplify_expr(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> clone_into_context(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> add_exprs(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> negate_expr(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> multiply_exprs(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> divide_exprs(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> pow_expr(ExprPtr base, std::size_t exponent, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> subtract_exprs(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx);

[[nodiscard]] BigInt pow_bigint_nonnegative(BigInt base, unsigned int exponent);

// Polynomial helpers
[[nodiscard]] Result<BigInt> expr_to_integer_coefficient(ExprPtr expr);
[[nodiscard]] BigInt integer_content(const IntPoly& coefficients);
void divide_integer_coefficients_by_scalar(IntPoly& coefficients, const BigInt& scalar);
[[nodiscard]] PolyExpr integer_coefficients_to_poly(const IntPoly& coefficients, AstArena& arena);
[[nodiscard]] Result<ExprPtr> integer_coefficients_to_expr(const IntPoly& coefficients, const Symbol& var, symbolic::CASContext& ctx);

[[nodiscard]] Result<IntegerExponent> parse_integer_exponent(ExprPtr expr);

// Internal polynomial operations
[[nodiscard]] Result<RationalParts> split_num_den(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> expand_expr_impl(ExprPtr expr, symbolic::CASContext& ctx);

// Factorization helpers
void append_factor_with_multiplicity(std::vector<PolynomialFactor>& factors, ExprPtr factor, unsigned int multiplicity = 1U);
[[nodiscard]] IntPoly primitive_integer_poly_local(IntPoly coefficients);

// Internal multivariate helpers
[[nodiscard]] Result<BigInt> expr_to_integer_value_for_multivariate(ExprPtr value);
[[nodiscard]] Result<ExprPtr> build_multivariate_monomial_expr(const MultivariateTerm& term, symbolic::CASContext& ctx);

// Factorization structures
struct IntegerSquareFreeFactor {
    IntPoly factor;
    unsigned int multiplicity{1U};
};

[[nodiscard]] Result<std::vector<IntegerSquareFreeFactor>> square_free_factorize_integer_poly(const IntPoly& primitive, symbolic::CASContext& ctx);
[[nodiscard]] Result<void> append_integer_factor_component(Factorization& factorization, const IntPoly& component, unsigned int multiplicity, const Symbol& var, symbolic::CASContext& ctx);

} // namespace cas::algebra
