// L3-06 internal helpers: factorization on a two-level algebraic tower.
//
// Shared between factorization_tower.cpp (entry point) and
// factorization_tower_helpers.cpp (helper bodies).  Not a public header.

#pragma once

#include "cas/algebra.hpp"
#include "cas/algebraic_number.hpp"
#include "cas/algebraic_tower.hpp"
#include "cas/algebraic_tower_bridge.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include "polynomial_internal.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace cas::algebra::factorization_tower_internal {

[[nodiscard]] CASError tower_error(CASErrorKind kind, std::string msg);
[[nodiscard]] bool is_fatal_inner_error(const CASError& err);
[[nodiscard]] ExprPtr rational_lit(AstArena& arena, const Rational& r);
[[nodiscard]] BigInt bigint_lcm(BigInt lhs, BigInt rhs);

[[nodiscard]] Result<ExprPtr> rational_coeffs_to_expr(
    const AlgebraicNumber::CoeffVec& coeffs,
    const Symbol& y,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> lift_outer_min_poly(
    const std::vector<AlgebraicNumber>& min_poly_2,
    const Symbol& y1,
    const Symbol& y2,
    symbolic::CASContext& ctx);

[[nodiscard]] std::size_t compute_default_shift_bound(
    std::size_t deg_f,
    std::size_t deg_m1,
    std::size_t deg_m2);

[[nodiscard]] ExprPtr build_shift_minus(
    const Symbol& x_sym,
    const Symbol& y1_sym,
    const Symbol& y2_sym,
    std::size_t s1,
    std::size_t s2,
    AstArena& arena);

[[nodiscard]] Result<ExprPtr> clear_denominators_to_integer_poly(
    const RatPoly& rat,
    const Symbol& var,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<bool> is_square_free_over_q(
    ExprPtr poly,
    const Symbol& var,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<std::vector<AlgebraicTowerTwoLevel>> shift_rational_factor_in_tower(
    const RatPoly& Ni_rat,
    std::size_t s1,
    std::size_t s2,
    const TowerGenerators& gens);

[[nodiscard]] std::vector<AlgebraicTowerTwoLevel> rational_poly_to_tower_coeffs(
    const RatPoly& f,
    const TowerGenerators& gens);

[[nodiscard]] Result<std::vector<AlgebraicTowerTwoLevel>> monic_tower(
    std::vector<AlgebraicTowerTwoLevel> coeffs);

[[nodiscard]] Result<ExprPtr> tower_coeffs_to_expr(
    const std::vector<AlgebraicTowerTwoLevel>& coeffs,
    const Symbol& var,
    const TowerGenerators& gens,
    symbolic::CASContext& ctx);

[[nodiscard]] bool tower_min_polys_well_formed(const TowerGenerators& gens);

}  // namespace cas::algebra::factorization_tower_internal
