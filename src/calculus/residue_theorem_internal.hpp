#pragma once

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "cas/result.hpp"
#include "cas/rational.hpp"
#include "cas/algebra.hpp"
#include "../algebra/polynomial_internal.hpp"
#include <vector>
#include <optional>

namespace cas::calculus {

[[nodiscard]] inline ExprPtr make_int(AstArena& arena, long long n) {
    return arena.make<IntegerLit>(BigInt(n));
}
[[nodiscard]] ExprPtr make_rational_expr(AstArena& arena, const Rational& r);
[[nodiscard]] Result<ExprPtr> simplify_or_fail(ExprPtr e, symbolic::CASContext& ctx);
[[nodiscard]] std::optional<Rational> as_rational(ExprPtr e);

[[nodiscard]] Result<ExprPtr> poly_coeff_at_zero(
    ExprPtr f, const Symbol& var, unsigned int order, symbolic::CASContext& ctx);
[[nodiscard]] Result<std::optional<std::size_t>> poly_degree_rational(
    ExprPtr f, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] Result<std::vector<Rational>> extract_rational_coeffs(
    ExprPtr f, const Symbol& var, std::size_t degree, symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> contribution_from_quadratic(
    const Rational& b,
    const Rational& c,
    const Rational& discriminant,
    ExprPtr N_over_D,
    const Symbol& var,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> contribution_from_irreducible_biquadratic(
    const Rational& b,
    const Rational& c,
    ExprPtr N_over_D,
    const Symbol& var,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> numeric_residue_contribution(
    const algebra::PolynomialFactor& pf,
    ExprPtr N,
    ExprPtr D,
    const Symbol& var,
    std::size_t deg_N,
    std::size_t deg_D,
    symbolic::CASContext& ctx);

} // namespace cas::calculus
