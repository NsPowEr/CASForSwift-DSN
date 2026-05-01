#pragma once

#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas {
namespace symbolic {
class CASContext;
}

namespace algebra {

template <typename Coeff>
class UnivariatePolynomial {
public:
    UnivariatePolynomial() = default;
    explicit UnivariatePolynomial(std::vector<Coeff> coefficients)
        : coefficients_(std::move(coefficients)) {}

    [[nodiscard]] bool empty() const noexcept { return coefficients_.empty(); }
    [[nodiscard]] bool is_zero() const noexcept { return coefficients_.empty(); }
    [[nodiscard]] std::size_t size() const noexcept { return coefficients_.size(); }
    [[nodiscard]] std::size_t degree() const noexcept { return coefficients_.empty() ? 0U : coefficients_.size() - 1U; }
    [[nodiscard]] Coeff constant_term() const { return coefficients_.empty() ? Coeff{} : coefficients_.front(); }
    [[nodiscard]] const Coeff& leading_coeff() const { return coefficients_.back(); }
    [[nodiscard]] const std::vector<Coeff>& coefficients() const noexcept { return coefficients_; }
    [[nodiscard]] std::vector<Coeff>& coefficients() noexcept { return coefficients_; }

    void reserve(std::size_t size) { coefficients_.reserve(size); }
    void resize(std::size_t size, const Coeff& value = Coeff{}) { coefficients_.resize(size, value); }
    void push_back(Coeff value) { coefficients_.push_back(std::move(value)); }

    Coeff& operator[](std::size_t index) { return coefficients_[index]; }
    const Coeff& operator[](std::size_t index) const { return coefficients_[index]; }

    template <typename ZeroPredicate>
    void normalize(ZeroPredicate&& is_zero) {
        while (!coefficients_.empty() && is_zero(coefficients_.back())) {
            coefficients_.pop_back();
        }
    }

private:
    std::vector<Coeff> coefficients_;
};

using IntPoly = UnivariatePolynomial<BigInt>;
using RatPoly = UnivariatePolynomial<Rational>;
using PolyExpr = UnivariatePolynomial<ExprPtr>;

struct IntegerExponent {
    std::size_t magnitude{0U};
    bool negative{false};
};

struct PolyDivisionResult {
    PolyExpr quotient;
    PolyExpr remainder;
};

struct RationalRootCandidate {
    BigInt numerator;
    BigInt denominator;
};

enum class IntegerGcdPath {
    Subresultant,
    PrimitiveFallback,
    PrimitiveFallbackPsi,
    PrimitiveFallbackBeta,
};

struct IntegerGcdResult {
    IntPoly gcd;
    IntegerGcdPath path{IntegerGcdPath::PrimitiveFallback};
};

[[nodiscard]] ExprPtr poly_make_integer(AstArena& arena, long long value);
[[nodiscard]] bool poly_is_zero_expr(ExprPtr expr);
[[nodiscard]] bool poly_is_one_expr(ExprPtr expr);
[[nodiscard]] bool poly_is_minus_one_expr(ExprPtr expr);
[[nodiscard]] bool poly_depends_on(ExprPtr expr, const std::string& variable_name);
[[nodiscard]] bool poly_contains_decimal_literal(ExprPtr expr);
[[nodiscard]] Result<ExprPtr> poly_simplify_expr(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> poly_clone_into_context(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] PolyExpr poly_make_monomial(ExprPtr coefficient, std::size_t degree);
[[nodiscard]] Result<PolyExpr> poly_make_constant_poly(ExprPtr coefficient, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> poly_add(const PolyExpr& lhs, const PolyExpr& rhs, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> poly_negate(const PolyExpr& poly, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> poly_subtract(const PolyExpr& lhs, const PolyExpr& rhs, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> poly_multiply(const PolyExpr& lhs, const PolyExpr& rhs, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> poly_divide_by_scalar(const PolyExpr& poly, ExprPtr scalar, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> poly_pow(PolyExpr base, std::size_t exponent, symbolic::CASContext& ctx);
[[nodiscard]] Result<std::size_t> poly_parse_nonnegative_integer_exponent(ExprPtr expr);

template <typename UInt>
[[nodiscard]] std::optional<UInt> parse_bounded_unsigned_decimal(const std::string& decimal) {
    static_assert(std::numeric_limits<UInt>::is_integer && !std::numeric_limits<UInt>::is_signed);

    UInt value = 0;
    for (char ch : decimal) {
        const unsigned int digit = static_cast<unsigned int>(ch - '0');
        if (value > (std::numeric_limits<UInt>::max() - static_cast<UInt>(digit)) / static_cast<UInt>(10)) {
            return std::nullopt;
        }
        value = static_cast<UInt>(value * static_cast<UInt>(10) + static_cast<UInt>(digit));
    }
    return value;
}

void normalize_poly(PolyExpr& poly);
[[nodiscard]] bool is_zero_poly(const PolyExpr& poly);
[[nodiscard]] std::size_t poly_degree(const PolyExpr& poly);
[[nodiscard]] ExprPtr leading_coefficient(const PolyExpr& poly);
[[nodiscard]] Result<PolyDivisionResult> divide_poly_with_remainder(
    const PolyExpr& dividend,
    const PolyExpr& divisor,
    symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> normalize_poly_monic(const PolyExpr& poly, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> parse_polynomial(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> polynomial_to_expr(const PolyExpr& poly, const Symbol& var, symbolic::CASContext& ctx);

[[nodiscard]] Result<IntPoly> poly_to_integer_poly(const PolyExpr& poly);
[[nodiscard]] Result<std::vector<BigInt>> poly_to_integer_coefficients(const PolyExpr& poly);
[[nodiscard]] bool is_zero_integer_poly(const IntPoly& coefficients);
[[nodiscard]] IntPoly gcd_integer_poly_primitive(IntPoly lhs, IntPoly rhs);
[[nodiscard]] std::optional<IntPoly> gcd_integer_poly_subresultant(IntPoly lhs, IntPoly rhs);
[[nodiscard]] IntegerGcdResult gcd_integer_poly_with_subresultant(IntPoly lhs, IntPoly rhs);
[[nodiscard]] std::optional<RationalRootCandidate> find_rational_root_candidate(const IntPoly& coefficients);
[[nodiscard]] Result<RatPoly> poly_to_rational_poly(const PolyExpr& poly);
[[nodiscard]] Result<std::vector<Rational>> poly_to_rational_coefficients(const PolyExpr& poly);
[[nodiscard]] Rational evaluate_integer_polynomial_at(const IntPoly& coefficients, const Rational& value);
[[nodiscard]] Rational evaluate_rational_polynomial_at(const RatPoly& coefficients, const Rational& value);

}  // namespace algebra
}  // namespace cas
