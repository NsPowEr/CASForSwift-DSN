// L3-06: factorization of f(x) in Q[x] over a two-level algebraic tower
// Q(alpha_1, alpha_2) via composite Trager shift + iterated absolute
// resultant.
//
// Algorithm:
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
//      N(x) lies in Q[x].  Denominators are cleared via lcm-scaling
//      before factor_over_integers (no Categoria-4 bail-out on type).
//
//   4. Iterate (s1, s2) until N(x) is square-free in Q[x].  Square-freeness
//      is verified rigorously via gcd(N, N') == 1 (not by inspecting the
//      multiplicity of factor_over_integers, which would be fragile).
//      The total number of attempts is bounded by the discriminant
//      collision bound (compute_default_shift_bound) or the user-supplied
//      ctx.max_trager_tower_shift_attempts.  Timeout errors from inner
//      symbolic operations are propagated, never swallowed.
//
//   5. Factor N(x) = prod N_i(x) over Q.  Each rational factor N_i lifts
//      to a factor f_i(x) of f(x) over Q(alpha_1, alpha_2) via
//        f_i(x) = gcd_{Q(alpha_1, alpha_2)[x]} ( f(x), N_i(x + s1*alpha_1 + s2*alpha_2) ).
//      The tower-gcd is computed with AlgebraicElement<AlgebraicNumber>
//      arithmetic (the templated tower_detail::poly_extended_gcd).
//
//   6. Post-condition: prod_i deg(f_i) == deg(f).  If lift fails for some
//      N_i (coefficient not expressible in tower), the result would be
//      incomplete; in that case we return Unimplemented with a diagnostic
//      message rather than a presumed-irreducible fallback.
//
//   7. If every gcd is trivial AND every lift succeeded, f is genuinely
//      irreducible over the tower (no Galois conjugate of f's roots lies
//      in Q(alpha_1, alpha_2)), so we return {f} as the single factor.

#include "cas/algebra.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "cas/algebraic_tower_bridge.hpp"
#include "cas/calculus.hpp"
#include "cas/symbolic.hpp"

#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <string>

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

[[nodiscard]] bool is_fatal_inner_error(const CASError& err) {
    return err.kind == CASErrorKind::Timeout
        || err.kind == CASErrorKind::InternalError;
}

[[nodiscard]] ExprPtr rational_lit(AstArena& arena, const Rational& r) {
    return arena.make<RationalLit>(r.numerator(), r.denominator());
}

[[nodiscard]] BigInt bigint_lcm(BigInt lhs, BigInt rhs) {
    lhs = lhs.abs();
    rhs = rhs.abs();
    if (lhs.is_zero() || rhs.is_zero()) return BigInt(0);
    return (lhs / gcd(lhs, rhs)) * rhs;
}

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

// Scale a rational polynomial N(x) by lcm of denominators and return the
// integer-coefficient ExprPtr.  Replaces the previous bail-out on
// non-integer rational coefficients (Categoria 4 compliance).
[[nodiscard]] Result<ExprPtr> clear_denominators_to_integer_poly(
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

// Rigorous square-free check via gcd(p, p') == 1 over Q[var].
[[nodiscard]] Result<bool> is_square_free_over_q(
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

// Shift a rational factor N_i(x) by x -> x + s1*alpha_1 + s2*alpha_2 and
// produce its coefficients in the tower Q(alpha_1, alpha_2) directly via
// algebraic arithmetic on AlgebraicTowerTwoLevel.  This bypasses the
// symbolic substitute/expand pipeline (which is fragile under simplifier
// rewrites such as RootOf <-> sqrt) and computes
//
//   N_i(x + beta) = sum_k c_k (x + beta)^k = sum_j ( sum_{k>=j} c_k * C(k,j) * beta^(k-j) ) x^j
//
// using exact tower arithmetic.
[[nodiscard]] Result<std::vector<AlgebraicTowerTwoLevel>> shift_rational_factor_in_tower(
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

    // beta = s1*alpha_1 + s2*alpha_2; in tower coords its value vector is
    // [s1*alpha_inner_unit, s2*one_inner].
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

[[nodiscard]] std::vector<AlgebraicTowerTwoLevel> rational_poly_to_tower_coeffs(
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

// Validate gens.min_poly_2 has a monic leading AlgebraicNumber (degree
// preserved at construction time of AlgebraicElement, but we check the
// caller's input to fail fast on misuse).
[[nodiscard]] bool tower_min_polys_well_formed(const TowerGenerators& gens) {
    if (gens.min_poly_1.empty() || gens.min_poly_1.back().numerator().is_zero()) return false;
    if (gens.min_poly_2.empty()) return false;
    const AlgebraicNumber& lead = gens.min_poly_2.back();
    if (lead.is_zero()) return false;
    return true;
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
    if (!tower_min_polys_well_formed(gens)) {
        return fail<Factorization>(tower_error(
            CASErrorKind::InvalidArgument,
            "factor_polynomial_tower: tower minimal polynomials must have non-zero leading coefficient"));
    }

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
            "factor_polynomial_tower: zero polynomial has no canonical factorization"));
    }
    const std::size_t deg_f = f_poly.degree();
    const std::size_t deg_m1 = gens.min_poly_1.size() - 1U;
    const std::size_t deg_m2 = gens.min_poly_2.size() - 1U;
    if (deg_f == 0U) {
        Factorization trivial;
        trivial.content = rational_lit(ctx.arena(), f_poly.constant_term());
        return ok(std::move(trivial));
    }

    Symbol y1_sym = ctx.make_fresh_symbol("__l306_y1");
    Symbol y2_sym = ctx.make_fresh_symbol("__l306_y2");

    auto m1_expr_res = rational_coeffs_to_expr(gens.min_poly_1, y1_sym, ctx);
    if (m1_expr_res.is_error()) return fail<Factorization>(m1_expr_res.error());
    auto m2_expr_res = lift_outer_min_poly(gens.min_poly_2, y1_sym, y2_sym, ctx);
    if (m2_expr_res.is_error()) return fail<Factorization>(m2_expr_res.error());

    PolyExpr f_poly_expr;
    f_poly_expr.reserve(f_poly.size());
    for (const Rational& r : f_poly.coefficients()) {
        f_poly_expr.push_back(r.numerator().is_zero() ? ExprPtr{} : rational_lit(ctx.arena(), r));
    }
    normalize_poly(f_poly_expr);
    auto f_expr_res = polynomial_to_expr(f_poly_expr, var, ctx);
    if (f_expr_res.is_error()) return fail<Factorization>(f_expr_res.error());

    const std::size_t user_bound = ctx.max_trager_tower_shift_attempts();
    const std::size_t default_bound = compute_default_shift_bound(deg_f, deg_m1, deg_m2);
    const std::size_t max_attempts = (user_bound > 0U) ? user_bound : default_bound;

    // Enumerate (s1, s2) by total weight in increasing order, counting
    // attempts rather than weight.  This honours the documented
    // "max attempts" semantics of ctx.max_trager_tower_shift_attempts.
    std::size_t attempts = 0U;
    for (std::size_t weight = 0U; attempts < max_attempts; ++weight) {
        for (std::size_t s1 = 0U; s1 <= weight && attempts < max_attempts; ++s1) {
            const std::size_t s2 = weight - s1;
            ++attempts;

            ExprPtr shift_minus = build_shift_minus(var, y1_sym, y2_sym, s1, s2, ctx.arena());
            auto shifted_f_raw = ctx.substitute(f_expr_res.value(), var, shift_minus);
            if (shifted_f_raw.is_error()) {
                if (is_fatal_inner_error(shifted_f_raw.error())) return fail<Factorization>(shifted_f_raw.error());
                continue;
            }
            auto shifted_f_exp = expand(shifted_f_raw.value(), ctx);
            if (shifted_f_exp.is_error()) {
                if (is_fatal_inner_error(shifted_f_exp.error())) return fail<Factorization>(shifted_f_exp.error());
                continue;
            }

            auto n1_res = polynomial_resultant(m2_expr_res.value(), shifted_f_exp.value(), y2_sym, ctx);
            if (n1_res.is_error()) {
                if (is_fatal_inner_error(n1_res.error())) return fail<Factorization>(n1_res.error());
                continue;
            }
            auto n1_exp = expand(n1_res.value(), ctx);
            if (n1_exp.is_error()) {
                if (is_fatal_inner_error(n1_exp.error())) return fail<Factorization>(n1_exp.error());
                continue;
            }

            auto norm_res = polynomial_resultant(m1_expr_res.value(), n1_exp.value(), y1_sym, ctx);
            if (norm_res.is_error()) {
                if (is_fatal_inner_error(norm_res.error())) return fail<Factorization>(norm_res.error());
                continue;
            }
            auto norm_exp = expand(norm_res.value(), ctx);
            if (norm_exp.is_error()) {
                if (is_fatal_inner_error(norm_exp.error())) return fail<Factorization>(norm_exp.error());
                continue;
            }

            // Rigorous square-free check: gcd(N, N') must be a unit.
            auto sf = is_square_free_over_q(norm_exp.value(), var, ctx);
            if (sf.is_error()) {
                if (is_fatal_inner_error(sf.error())) return fail<Factorization>(sf.error());
                continue;
            }
            if (!sf.value()) continue;

            // Parse and clear denominators before integer factorization.
            auto norm_parsed = parse_polynomial(norm_exp.value(), var, ctx);
            if (norm_parsed.is_error()) {
                if (is_fatal_inner_error(norm_parsed.error())) return fail<Factorization>(norm_parsed.error());
                continue;
            }
            auto norm_rat = poly_to_rational_poly(norm_parsed.value());
            if (norm_rat.is_error()) {
                if (is_fatal_inner_error(norm_rat.error())) return fail<Factorization>(norm_rat.error());
                continue;
            }
            auto integer_norm = clear_denominators_to_integer_poly(norm_rat.value(), var, ctx);
            if (integer_norm.is_error()) return fail<Factorization>(integer_norm.error());

            auto norm_factorization = factor_over_integers(integer_norm.value(), var, ctx);
            if (norm_factorization.is_error()) {
                if (is_fatal_inner_error(norm_factorization.error())) return fail<Factorization>(norm_factorization.error());
                continue;
            }

            // Lift each rational factor Ni(x) to a tower factor of f.
            auto f_tower_raw = rational_poly_to_tower_coeffs(f_poly, gens);
            auto f_tower = monic_tower(std::move(f_tower_raw));
            if (f_tower.is_error()) return fail<Factorization>(f_tower.error());

            Factorization out;
            out.content = make_integer(ctx.arena(), 1);
            std::size_t cumulative_degree = 0U;

            for (const auto& pf : norm_factorization.value().factors) {
                // Parse N_i as Q[var] (factor_over_integers returns integer
                // polynomial factors, but we route them through poly_to_rational_poly
                // for a uniform RatPoly representation).
                auto ni_parsed = parse_polynomial(pf.factor, var, ctx);
                if (ni_parsed.is_error()) {
                    if (is_fatal_inner_error(ni_parsed.error())) return fail<Factorization>(ni_parsed.error());
                    continue;
                }
                auto ni_rat = poly_to_rational_poly(ni_parsed.value());
                if (ni_rat.is_error()) return fail<Factorization>(ni_rat.error());

                auto gi_tower_res = shift_rational_factor_in_tower(ni_rat.value(), s1, s2, gens);
                if (gi_tower_res.is_error()) return fail<Factorization>(gi_tower_res.error());
                auto gi_tower = monic_tower(std::move(gi_tower_res.value()));
                if (gi_tower.is_error()) return fail<Factorization>(gi_tower.error());

                auto gcd_res = tower_detail::poly_extended_gcd<AlgebraicTowerTwoLevel>(
                    f_tower.value(), gi_tower.value());
                if (gcd_res.is_error()) return fail<Factorization>(gcd_res.error());
                auto [g_coeffs, s_coeffs, t_coeffs] = std::move(gcd_res.value());
                (void)s_coeffs;
                (void)t_coeffs;
                tower_detail::strip_trailing(g_coeffs);
                if (g_coeffs.size() <= 1U) continue;  // unit -> no factor

                auto monic_g = monic_tower(std::move(g_coeffs));
                if (monic_g.is_error()) return fail<Factorization>(monic_g.error());
                cumulative_degree += monic_g.value().size() - 1U;

                auto factor_expr = tower_coeffs_to_expr(monic_g.value(), var, gens, ctx);
                if (factor_expr.is_error()) return fail<Factorization>(factor_expr.error());
                out.factors.push_back(PolynomialFactor{factor_expr.value(), 1U});
            }

            if (out.factors.empty()) {
                // Every gcd was a unit: no Galois conjugate of f's roots
                // lies in Q(alpha_1, alpha_2) under this shift, so f is
                // irreducible over the tower.
                auto f_back = polynomial_to_expr(f_poly_expr, var, ctx);
                if (f_back.is_error()) return fail<Factorization>(f_back.error());
                out.factors.push_back(PolynomialFactor{f_back.value(), 1U});
                return ok(std::move(out));
            }

            // Post-condition: cumulative lifted degree must equal deg(f).
            // If it doesn't, the lift missed some Galois orbit; this would
            // indicate an internal bug since shift_rational_factor_in_tower
            // is exact and cannot fail to express coefficients.
            if (cumulative_degree != deg_f) {
                return fail<Factorization>(tower_error(
                    CASErrorKind::InternalError,
                    "factor_polynomial_tower: cumulative lifted degree " +
                    std::to_string(cumulative_degree) +
                    " != deg(f) " + std::to_string(deg_f) +
                    "; Trager invariant violated"));
            }
            return ok(std::move(out));
        }
    }

    return fail<Factorization>(tower_error(
        CASErrorKind::Unimplemented,
        "factor_polynomial_tower: no square-free composite Trager shift found within ctx.max_trager_tower_shift_attempts"));
}

}  // namespace algebra
}  // namespace cas
