// L3-06 helper bodies for factor_polynomial_tower.  See
// factorization_tower_internal.hpp for declarations and the entry-point
// algorithm in factorization_tower.cpp.

#include "factorization_tower_internal.hpp"

#include "cas/algebraic_number_bridge.hpp"
#include "cas/algebra.hpp"
#include "cas/calculus.hpp"

#include <algorithm>
#include <cstdint>
#include <utility>

namespace cas::algebra::factorization_tower_internal {

CASError tower_error(CASErrorKind kind, std::string msg) {
    return CASError{.kind = kind, .message = std::move(msg), .hint = std::nullopt};
}

bool is_fatal_inner_error(const CASError& err) {
    return err.kind == CASErrorKind::Timeout
        || err.kind == CASErrorKind::InternalError;
}

ExprPtr rational_lit(AstArena& arena, const Rational& r) {
    return arena.make<RationalLit>(r.numerator(), r.denominator());
}

BigInt bigint_lcm(BigInt lhs, BigInt rhs) {
    lhs = lhs.abs();
    rhs = rhs.abs();
    if (lhs.is_zero() || rhs.is_zero()) return BigInt(0);
    return (lhs / gcd(lhs, rhs)) * rhs;
}

Result<ExprPtr> rational_coeffs_to_expr(
    const AlgebraicNumber::CoeffVec& coeffs,
    const Symbol& y,
    symbolic::CASContext& ctx) {
    PolyExpr poly;
    poly.reserve(coeffs.size());
    for (const Rational& c : coeffs) {
        if (c.numerator().is_zero()) {
            poly.push_back(ExprPtr{});
        } else {
            poly.push_back(rational_lit(ctx.arena(), c));
        }
    }
    normalize_poly(poly);
    return polynomial_to_expr(poly, y, ctx);
}

Result<ExprPtr> lift_outer_min_poly(
    const std::vector<AlgebraicNumber>& min_poly_2,
    const Symbol& y1,
    const Symbol& y2,
    symbolic::CASContext& ctx) {
    PolyExpr poly_in_y2;
    poly_in_y2.reserve(min_poly_2.size());
    ExprPtr y1_expr = ctx.arena().make<Symbol>(y1.name);
    for (const AlgebraicNumber& coeff : min_poly_2) {
        if (coeff.is_zero()) {
            poly_in_y2.push_back(ExprPtr{});
            continue;
        }
        ExprPtr raw = algebraic_number_to_expr_raw(coeff, y1_expr, ctx.arena());
        auto simplified = poly_simplify_expr(raw, ctx);
        if (simplified.is_error()) return fail<ExprPtr>(simplified.error());
        poly_in_y2.push_back(simplified.value());
    }
    normalize_poly(poly_in_y2);
    return polynomial_to_expr(poly_in_y2, y2, ctx);
}

std::size_t compute_default_shift_bound(
    std::size_t deg_f,
    std::size_t deg_m1,
    std::size_t deg_m2) {
    const std::size_t df = std::max<std::size_t>(1U, deg_f);
    const std::size_t d1 = std::max<std::size_t>(1U, deg_m1);
    const std::size_t d2 = std::max<std::size_t>(1U, deg_m2);
    const std::size_t bound = 2U * df * d1 * d2 * (df + d1 + d2);
    return bound + 1U;
}

ExprPtr build_shift_minus(
    const Symbol& x_sym,
    const Symbol& y1_sym,
    const Symbol& y2_sym,
    std::size_t s1,
    std::size_t s2,
    AstArena& arena) {
    ExprPtr x = arena.make<Symbol>(x_sym.name);
    if (s1 == 0U && s2 == 0U) return x;
    std::vector<ExprPtr> terms;
    terms.push_back(x);
    if (s1 != 0U) {
        ExprPtr t = arena.make<Binary>(BinaryOp::Mul,
            arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(s1))),
            arena.make<Symbol>(y1_sym.name));
        terms.push_back(arena.make<Unary>(UnaryOp::Neg, t));
    }
    if (s2 != 0U) {
        ExprPtr t = arena.make<Binary>(BinaryOp::Mul,
            arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(s2))),
            arena.make<Symbol>(y2_sym.name));
        terms.push_back(arena.make<Unary>(UnaryOp::Neg, t));
    }
    if (terms.size() == 1U) return terms.front();
    return arena.make<Sum>(std::move(terms));
}

Result<ExprPtr> clear_denominators_to_integer_poly(
    const RatPoly& rat,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    BigInt scale(1);
    for (const Rational& c : rat.coefficients()) {
        scale = bigint_lcm(scale, c.denominator());
    }
    PolyExpr integer_coeffs;
    integer_coeffs.reserve(rat.size());
    for (const Rational& c : rat.coefficients()) {
        const Rational scaled = c * Rational(scale);
        if (scaled.denominator() != BigInt(1)) {
            return fail<ExprPtr>(tower_error(
                CASErrorKind::InternalError,
                "factor_polynomial_tower: lcm-scaling did not produce integer coefficients"));
        }
        integer_coeffs.push_back(
            scaled.numerator().is_zero()
                ? ExprPtr{}
                : ctx.arena().make<IntegerLit>(scaled.numerator()));
    }
    normalize_poly(integer_coeffs);
    return polynomial_to_expr(integer_coeffs, var, ctx);
}

Result<bool> is_square_free_over_q(
    ExprPtr poly,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    auto deriv = cas::calculus::diff(poly, var, 1U, ctx);
    if (deriv.is_error()) return fail<bool>(deriv.error());
    auto gcd_res = polynomial_gcd(poly, deriv.value(), var, ctx);
    if (gcd_res.is_error()) return fail<bool>(gcd_res.error());
    auto parsed = parse_polynomial(gcd_res.value(), var, ctx);
    if (parsed.is_error()) return fail<bool>(parsed.error());
    return ok(poly_degree(parsed.value()) == 0U);
}

Result<std::vector<AlgebraicTowerTwoLevel>> shift_rational_factor_in_tower(
    const RatPoly& Ni_rat,
    std::size_t s1,
    std::size_t s2,
    const TowerGenerators& gens) {
    if (Ni_rat.is_zero()) return ok(std::vector<AlgebraicTowerTwoLevel>{});

    const AlgebraicNumber::CoeffVec& mp1 = gens.min_poly_1;
    const std::vector<AlgebraicNumber>& mp2 = gens.min_poly_2;

    AlgebraicNumber inner_one({Rational(BigInt(1))}, mp1);
    AlgebraicNumber inner_alpha({Rational(BigInt(0)), Rational(BigInt(1))}, mp1);
    AlgebraicNumber inner_s1({Rational(BigInt(static_cast<std::int64_t>(s1)))}, mp1);
    AlgebraicNumber inner_s2({Rational(BigInt(static_cast<std::int64_t>(s2)))}, mp1);

    AlgebraicTowerTwoLevel beta(
        {inner_s1 * inner_alpha, inner_s2 * inner_one},
        mp2);
    AlgebraicTowerTwoLevel tower_one({inner_one}, mp2);

    const std::size_t n = Ni_rat.degree();
    std::vector<AlgebraicTowerTwoLevel> beta_powers;
    beta_powers.reserve(n + 1U);
    beta_powers.push_back(tower_one);
    for (std::size_t k = 1U; k <= n; ++k) {
        beta_powers.push_back(beta_powers.back() * beta);
    }

    std::vector<std::vector<BigInt>> binom(n + 1U, std::vector<BigInt>(n + 1U, BigInt(0)));
    for (std::size_t k = 0U; k <= n; ++k) {
        binom[k][0] = BigInt(1);
        for (std::size_t j = 1U; j <= k; ++j) {
            binom[k][j] = binom[k - 1U][j - 1U] + binom[k - 1U][j];
        }
    }

    std::vector<AlgebraicTowerTwoLevel> result;
    result.reserve(n + 1U);
    AlgebraicTowerTwoLevel zero_tower({}, mp2);
    for (std::size_t j = 0U; j <= n; ++j) {
        AlgebraicTowerTwoLevel coeff = zero_tower;
        for (std::size_t k = j; k <= n; ++k) {
            const Rational& c_k = Ni_rat[k];
            if (c_k.numerator().is_zero()) continue;
            const Rational c_k_scaled = c_k * Rational(binom[k][j]);
            if (c_k_scaled.numerator().is_zero()) continue;
            AlgebraicNumber inner({c_k_scaled}, mp1);
            AlgebraicTowerTwoLevel term =
                AlgebraicTowerTwoLevel({inner}, mp2) * beta_powers[k - j];
            coeff = coeff + term;
        }
        result.push_back(std::move(coeff));
    }
    return ok(std::move(result));
}

std::vector<AlgebraicTowerTwoLevel> rational_poly_to_tower_coeffs(
    const RatPoly& f,
    const TowerGenerators& gens) {
    std::vector<AlgebraicTowerTwoLevel> coeffs;
    coeffs.reserve(f.size());
    for (const Rational& r : f.coefficients()) {
        if (r.numerator().is_zero()) {
            coeffs.push_back(AlgebraicTowerTwoLevel({}, gens.min_poly_2));
            continue;
        }
        AlgebraicNumber inner({r}, gens.min_poly_1);
        coeffs.push_back(AlgebraicTowerTwoLevel({inner}, gens.min_poly_2));
    }
    return coeffs;
}

Result<std::vector<AlgebraicTowerTwoLevel>> monic_tower(
    std::vector<AlgebraicTowerTwoLevel> coeffs) {
    tower_detail::strip_trailing(coeffs);
    if (coeffs.empty()) return ok(std::move(coeffs));
    auto inv_res = coeffs.back().inverse();
    if (inv_res.is_error()) return fail<std::vector<AlgebraicTowerTwoLevel>>(inv_res.error());
    for (auto& c : coeffs) c = c * inv_res.value();
    tower_detail::strip_trailing(coeffs);
    return ok(std::move(coeffs));
}

Result<ExprPtr> tower_coeffs_to_expr(
    const std::vector<AlgebraicTowerTwoLevel>& coeffs,
    const Symbol& var,
    const TowerGenerators& gens,
    symbolic::CASContext& ctx) {
    PolyExpr poly;
    poly.reserve(coeffs.size());
    for (const auto& c : coeffs) {
        if (c.is_zero()) {
            poly.push_back(ExprPtr{});
            continue;
        }
        ExprPtr e = tower_to_expr(c, gens, ctx.arena());
        auto simplified = poly_simplify_expr(e, ctx);
        if (simplified.is_error()) return fail<ExprPtr>(simplified.error());
        poly.push_back(simplified.value());
    }
    normalize_poly(poly);
    auto poly_expr = polynomial_to_expr(poly, var, ctx);
    if (poly_expr.is_error()) return fail<ExprPtr>(poly_expr.error());
    return ctx.simplify(poly_expr.value());
}

bool tower_min_polys_well_formed(const TowerGenerators& gens) {
    if (gens.min_poly_1.empty() || gens.min_poly_1.back().numerator().is_zero()) return false;
    if (gens.min_poly_2.empty()) return false;
    const AlgebraicNumber& lead = gens.min_poly_2.back();
    if (lead.is_zero()) return false;
    return true;
}

}  // namespace cas::algebra::factorization_tower_internal
