// L3-06: factorization of f(x) in Q[x] over a two-level algebraic tower
// Q(alpha_1, alpha_2) via composite Trager shift + iterated absolute
// resultant.
//
// Algorithm (Trager's norm method extended to a 2-level tower):
//
//   1. Build minimal polynomials in two fresh symbols y1, y2:
//        m1(y1)         in Q[y1]      (rational coefficients = gens.min_poly_1)
//        m2(y1, y2)     in Q[y1, y2]  (lifted from vector<AlgebraicNumber>)
//
//   2. Form the symbolic shift  s = s1*y1 + s2*y2  and the shifted
//      polynomial f_shifted(x, y1, y2) = f(x - s) in Q[x, y1, y2].
//
//   3. Compute the absolute norm in two steps:
//        N1(x, y1) = Res_{y2} ( m2(y1, y2), f_shifted(x, y1, y2) )
//        N (x)     = Res_{y1} ( m1(y1),     N1(x, y1) )
//      N(x) lies in Q[x].
//
//   4. Iterate (s1, s2) until N(x) is square-free in Q[x].  The number of
//      attempts is bounded by the discriminant collision bound (see
//      compute_default_shift_bound below) or by the user-supplied
//      ctx.max_trager_tower_shift_attempts.
//
//   5. Factor N(x) = prod N_i(x) over Q.  Each rational factor N_i lifts
//      to a factor f_i(x) of f(x) over Q(alpha_1, alpha_2) via
//        f_i(x) = gcd_{Q(alpha_1, alpha_2)[x]} ( f(x), N_i(x + s1*alpha_1 + s2*alpha_2) ).
//
//   6. The tower-gcd is computed with AlgebraicElement<AlgebraicNumber>
//      arithmetic (the templated tower_detail::poly_extended_gcd).

#include "cas/algebra.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "cas/algebraic_tower_bridge.hpp"
#include "cas/symbolic.hpp"

#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace cas {
namespace algebra {

namespace {

[[nodiscard]] CASError tower_error(CASErrorKind kind, std::string msg) {
    return CASError{.kind = kind, .message = std::move(msg), .hint = std::nullopt};
}

[[nodiscard]] ExprPtr rational_lit(AstArena& arena, const Rational& r) {
    return arena.make<RationalLit>(r.numerator(), r.denominator());
}

// Build symbolic polynomial in symbol y from a CoeffVec of rationals
// (used for m1 in Q[y1]).
[[nodiscard]] Result<ExprPtr> rational_coeffs_to_expr(
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

// Build symbolic polynomial in (y1, y2) from a vector<AlgebraicNumber>
// representing m2 in Q(alpha_1)[y2].  Each AN encodes a poly in y1.
[[nodiscard]] Result<ExprPtr> lift_outer_min_poly(
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

[[nodiscard]] std::size_t compute_default_shift_bound(
    std::size_t deg_f,
    std::size_t deg_m1,
    std::size_t deg_m2) {
    const std::size_t df = std::max<std::size_t>(1U, deg_f);
    const std::size_t d1 = std::max<std::size_t>(1U, deg_m1);
    const std::size_t d2 = std::max<std::size_t>(1U, deg_m2);
    const std::size_t bound = 2U * df * d1 * d2 * (df + d1 + d2);
    return bound + 1U;
}

// Build (x - s1*y1 - s2*y2) as ExprPtr.
[[nodiscard]] ExprPtr build_shift_minus(
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

// Build (x + s1*alpha_1 + s2*alpha_2) as ExprPtr.
[[nodiscard]] ExprPtr build_shift_plus_alpha(
    const Symbol& x_sym,
    ExprPtr alpha_1,
    ExprPtr alpha_2,
    std::size_t s1,
    std::size_t s2,
    AstArena& arena) {
    ExprPtr x = arena.make<Symbol>(x_sym.name);
    if (s1 == 0U && s2 == 0U) return x;
    std::vector<ExprPtr> terms;
    terms.push_back(x);
    if (s1 != 0U) {
        terms.push_back(arena.make<Binary>(BinaryOp::Mul,
            arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(s1))),
            alpha_1));
    }
    if (s2 != 0U) {
        terms.push_back(arena.make<Binary>(BinaryOp::Mul,
            arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(s2))),
            alpha_2));
    }
    if (terms.size() == 1U) return terms.front();
    return arena.make<Sum>(std::move(terms));
}

// Check that a rational univariate factorization has no factor of
// multiplicity > 1 (i.e. the polynomial is square-free).
[[nodiscard]] bool factorization_is_square_free(const Factorization& fact) {
    for (const auto& pf : fact.factors) {
        if (pf.multiplicity > 1U) return false;
    }
    return true;
}

// Re-shift a rational factor Ni(x) by x -> x + s1*alpha_1 + s2*alpha_2 and
// convert each coefficient into the tower representation
// AlgebraicTowerTwoLevel.  Returns the coefficient vector (lowest degree
// first), suitable for tower_detail::poly_extended_gcd.
[[nodiscard]] Result<std::vector<AlgebraicTowerTwoLevel>> rational_factor_to_tower_coeffs(
    ExprPtr Ni_expr,
    const Symbol& var,
    const TowerGenerators& gens,
    std::size_t s1,
    std::size_t s2,
    symbolic::CASContext& ctx) {
    ExprPtr shift_expr = build_shift_plus_alpha(var, gens.alpha_1, gens.alpha_2, s1, s2, ctx.arena());
    auto substituted = ctx.substitute(Ni_expr, var, shift_expr);
    if (substituted.is_error()) return fail<std::vector<AlgebraicTowerTwoLevel>>(substituted.error());
    auto expanded = expand(substituted.value(), ctx);
    if (expanded.is_error()) return fail<std::vector<AlgebraicTowerTwoLevel>>(expanded.error());

    auto parsed = parse_polynomial(expanded.value(), var, ctx);
    if (parsed.is_error()) return fail<std::vector<AlgebraicTowerTwoLevel>>(parsed.error());

    std::vector<AlgebraicTowerTwoLevel> coeffs;
    coeffs.reserve(parsed.value().coefficients().size());
    for (ExprPtr c : parsed.value().coefficients()) {
        if (!c || poly_is_zero_expr(c)) {
            coeffs.push_back(AlgebraicTowerTwoLevel({}, gens.min_poly_2));
            continue;
        }
        auto tower_coeff = try_express_in_tower_two_level(c, gens, ctx);
        if (tower_coeff.is_error()) return fail<std::vector<AlgebraicTowerTwoLevel>>(tower_coeff.error());
        if (!tower_coeff.value().has_value()) {
            // Coefficient not expressible in tower: caller treats this as
            // "this rational factor does not lift" and skips it.
            return ok(std::vector<AlgebraicTowerTwoLevel>{});
        }
        coeffs.push_back(std::move(tower_coeff.value().value()));
    }
    return ok(std::move(coeffs));
}

// Convert f in Q[x] into a coefficient vector with tower constants.
[[nodiscard]] Result<std::vector<AlgebraicTowerTwoLevel>> rational_poly_to_tower_coeffs(
    const RatPoly& f,
    const TowerGenerators& gens) {
    AlgebraicNumber inner_zero_an({}, gens.min_poly_1);
    std::vector<AlgebraicTowerTwoLevel> coeffs;
    coeffs.reserve(f.size());
    for (const Rational& r : f.coefficients()) {
        if (r.numerator().is_zero()) {
            coeffs.push_back(AlgebraicTowerTwoLevel({}, gens.min_poly_2));
            continue;
        }
        AlgebraicNumber inner({r}, gens.min_poly_1);
        coeffs.push_back(AlgebraicTowerTwoLevel({inner}, gens.min_poly_2));
        (void)inner_zero_an;
    }
    return ok(std::move(coeffs));
}

// Normalise a tower-coefficient polynomial to monic form.
[[nodiscard]] Result<std::vector<AlgebraicTowerTwoLevel>> monic_tower(
    std::vector<AlgebraicTowerTwoLevel> coeffs) {
    tower_detail::strip_trailing(coeffs);
    if (coeffs.empty()) return ok(std::move(coeffs));
    auto inv_res = coeffs.back().inverse();
    if (inv_res.is_error()) return fail<std::vector<AlgebraicTowerTwoLevel>>(inv_res.error());
    for (auto& c : coeffs) c = c * inv_res.value();
    tower_detail::strip_trailing(coeffs);
    return ok(std::move(coeffs));
}

// Build an ExprPtr polynomial in var from a tower-coefficient vector.
[[nodiscard]] Result<ExprPtr> tower_coeffs_to_expr(
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

}  // namespace

Result<Factorization> factor_polynomial_tower(
    ExprPtr poly,
    const Symbol& var,
    const TowerGenerators& gens,
    symbolic::CASContext& ctx) {
    if (!poly) {
        return fail<Factorization>(tower_error(
            CASErrorKind::InvalidArgument,
            "factor_polynomial_tower: null polynomial"));
    }
    if (gens.min_poly_1.size() < 2U || gens.min_poly_2.size() < 2U) {
        return fail<Factorization>(tower_error(
            CASErrorKind::InvalidArgument,
            "factor_polynomial_tower: tower generators must have minimal polynomials of positive degree"));
    }

    // Parse f over Q[var].
    auto parsed = parse_polynomial(poly, var, ctx);
    if (parsed.is_error()) return fail<Factorization>(parsed.error());
    auto f_rat = poly_to_rational_poly(parsed.value());
    if (f_rat.is_error()) {
        return fail<Factorization>(tower_error(
            CASErrorKind::Unimplemented,
            "factor_polynomial_tower: input polynomial must lie in Q[var]"));
    }
    RatPoly f_poly = std::move(f_rat.value());
    normalize_rational_coefficients(f_poly);
    if (f_poly.is_zero()) {
        return fail<Factorization>(tower_error(
            CASErrorKind::InvalidArgument,
            "factor_polynomial_tower: zero polynomial has no canonical factorisation"));
    }
    const std::size_t deg_f = f_poly.degree();
    const std::size_t deg_m1 = gens.min_poly_1.size() - 1U;
    const std::size_t deg_m2 = gens.min_poly_2.size() - 1U;
    if (deg_f == 0U) {
        Factorization trivial;
        trivial.content = rational_lit(ctx.arena(), f_poly.constant_term());
        return ok(std::move(trivial));
    }

    // Fresh symbols y1, y2 for the symbolic resultant computation.
    Symbol y1_sym = ctx.make_fresh_symbol("__l306_y1");
    Symbol y2_sym = ctx.make_fresh_symbol("__l306_y2");

    // Build m1(y1), m2(y1, y2), and f(x) as a single-variable poly in x.
    auto m1_expr_res = rational_coeffs_to_expr(gens.min_poly_1, y1_sym, ctx);
    if (m1_expr_res.is_error()) return fail<Factorization>(m1_expr_res.error());
    auto m2_expr_res = lift_outer_min_poly(gens.min_poly_2, y1_sym, y2_sym, ctx);
    if (m2_expr_res.is_error()) return fail<Factorization>(m2_expr_res.error());

    // f as expression: use the rational coefficients verbatim.
    PolyExpr f_poly_expr;
    f_poly_expr.reserve(f_poly.size());
    for (const Rational& r : f_poly.coefficients()) {
        f_poly_expr.push_back(r.numerator().is_zero() ? ExprPtr{} : rational_lit(ctx.arena(), r));
    }
    normalize_poly(f_poly_expr);
    auto f_expr_res = polynomial_to_expr(f_poly_expr, var, ctx);
    if (f_expr_res.is_error()) return fail<Factorization>(f_expr_res.error());

    // Bound the shift search.
    const std::size_t user_bound = ctx.max_trager_tower_shift_attempts();
    const std::size_t default_bound = compute_default_shift_bound(deg_f, deg_m1, deg_m2);
    const std::size_t shift_bound = (user_bound > 0U) ? user_bound : default_bound;

    // Enumerate (s1, s2) by total weight s1 + s2 = 0, 1, 2, ... up to
    // shift_bound.  For each pair, attempt the composite norm.
    for (std::size_t weight = 0U; weight < shift_bound; ++weight) {
        for (std::size_t s1 = 0U; s1 <= weight; ++s1) {
            const std::size_t s2 = weight - s1;
            // Shift f.
            ExprPtr shift_minus = build_shift_minus(var, y1_sym, y2_sym, s1, s2, ctx.arena());
            auto shifted_f_raw = ctx.substitute(f_expr_res.value(), var, shift_minus);
            if (shifted_f_raw.is_error()) return fail<Factorization>(shifted_f_raw.error());
            auto shifted_f_exp = expand(shifted_f_raw.value(), ctx);
            if (shifted_f_exp.is_error()) return fail<Factorization>(shifted_f_exp.error());

            // First resultant: eliminate y2.  Res_{y2}(m2, f_shifted).
            auto n1_res = polynomial_resultant(m2_expr_res.value(), shifted_f_exp.value(), y2_sym, ctx);
            if (n1_res.is_error()) continue;
            auto n1_exp = expand(n1_res.value(), ctx);
            if (n1_exp.is_error()) continue;

            // Second resultant: eliminate y1.  Res_{y1}(m1, N1).
            auto norm_res = polynomial_resultant(m1_expr_res.value(), n1_exp.value(), y1_sym, ctx);
            if (norm_res.is_error()) continue;
            auto norm_exp = expand(norm_res.value(), ctx);
            if (norm_exp.is_error()) continue;

            // Factor N(x) over Q and require square-free.
            auto norm_factorization = factor_over_integers(norm_exp.value(), var, ctx);
            if (norm_factorization.is_error()) {
                // Try Q[x] route: norm may have non-trivial rational content.
                auto rational_factor = parse_polynomial(norm_exp.value(), var, ctx);
                if (rational_factor.is_error()) continue;
                auto rational_poly_res = poly_to_rational_poly(rational_factor.value());
                if (rational_poly_res.is_error()) continue;
                // Bail: cannot factor non-integer rational norm via factor_over_integers without scaling.
                continue;
            }
            if (!factorization_is_square_free(norm_factorization.value())) continue;

            // Lift each rational factor Ni(x) to a tower factor of f via gcd in Q(alpha_1, alpha_2)[x].
            auto f_tower_res = rational_poly_to_tower_coeffs(f_poly, gens);
            if (f_tower_res.is_error()) return fail<Factorization>(f_tower_res.error());
            auto f_tower = monic_tower(std::move(f_tower_res.value()));
            if (f_tower.is_error()) return fail<Factorization>(f_tower.error());

            Factorization out;
            out.content = make_integer(ctx.arena(), 1);
            for (const auto& pf : norm_factorization.value().factors) {
                auto gi_tower_res = rational_factor_to_tower_coeffs(pf.factor, var, gens, s1, s2, ctx);
                if (gi_tower_res.is_error()) return fail<Factorization>(gi_tower_res.error());
                if (gi_tower_res.value().empty()) {
                    // Rational factor did not lift to tower: skip.
                    continue;
                }
                auto gi_tower = monic_tower(std::move(gi_tower_res.value()));
                if (gi_tower.is_error()) return fail<Factorization>(gi_tower.error());

                auto gcd_res = tower_detail::poly_extended_gcd<AlgebraicTowerTwoLevel>(
                    f_tower.value(), gi_tower.value());
                if (gcd_res.is_error()) {
                    return fail<Factorization>(gcd_res.error());
                }
                auto [g_coeffs, s_coeffs, t_coeffs] = std::move(gcd_res.value());
                (void)s_coeffs;
                (void)t_coeffs;
                tower_detail::strip_trailing(g_coeffs);
                if (g_coeffs.empty() || g_coeffs.size() == 1U) {
                    // gcd is a unit (degree 0) -> factor does not survive tower lift.
                    continue;
                }
                auto monic_g = monic_tower(std::move(g_coeffs));
                if (monic_g.is_error()) return fail<Factorization>(monic_g.error());

                auto factor_expr = tower_coeffs_to_expr(monic_g.value(), var, gens, ctx);
                if (factor_expr.is_error()) return fail<Factorization>(factor_expr.error());
                out.factors.push_back(PolynomialFactor{factor_expr.value(), 1U});
            }
            if (out.factors.empty()) {
                // No factor of f decomposes over the tower beyond the unit.
                // Then f is irreducible over Q(alpha_1, alpha_2) at this shift;
                // return f itself.
                auto f_back = polynomial_to_expr(f_poly_expr, var, ctx);
                if (f_back.is_error()) return fail<Factorization>(f_back.error());
                out.factors.push_back(PolynomialFactor{f_back.value(), 1U});
            }
            return ok(std::move(out));
        }
    }

    return fail<Factorization>(tower_error(
        CASErrorKind::Unimplemented,
        "factor_polynomial_tower: no square-free composite Trager shift found within bound"));
}

}  // namespace algebra
}  // namespace cas
