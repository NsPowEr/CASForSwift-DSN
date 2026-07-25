// A46 — fattorizzazione di un quartico razionale nei suoi due fattori quadratici
// REALI, quando i residui vivono in un'estensione quadratica di Q.
//
// Spec: Symbolic_Integration_I.md (Bronstein) §2.8, LogToReal righe 2381-2395 —
// la resa in forma chiusa richiede i residui, e i residui di un quartico
// generico non sono esprimibili per radicali reali. Lo sono esattamente quando
// il cubico risolvente del quartico ha una radice razionale.
//
// Diviso da partial_fractions_logtoreal.cpp (limite anti-monolito 500 righe):
// qui vive la parte "aritmetica esatta sui razionali" (Sturm + cubico
// risolvente), la' la conversione di Rioboo in forma reale.

#include "partial_fractions_rioboo.hpp"

#include "cas/algebra.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"

#include <utility>
#include <vector>

namespace cas::algebra {

namespace {

[[nodiscard]] Rational rat(long long value) { return Rational(BigInt(value)); }

[[nodiscard]] RatPoly derivative_rational(const RatPoly& poly) {
    RatPoly result;
    if (poly.size() <= 1U) return result;
    result.reserve(poly.size() - 1U);
    for (std::size_t i = 1; i < poly.size(); ++i) {
        result.push_back(poly[i] * rat(static_cast<long long>(i)));
    }
    normalize_rational_coefficients(result);
    return result;
}

[[nodiscard]] RatPoly negate_rational(const RatPoly& poly) {
    RatPoly result;
    result.reserve(poly.size());
    for (const Rational& c : poly.coefficients()) result.push_back(-c);
    normalize_rational_coefficients(result);
    return result;
}

// Un coefficiente radicale non si annulla per sola `simplify` se il prodotto di
// somme non viene distribuito: (p+√d)(p−√d) resta com'e'. `expand` prima del
// confronto e' quindi parte della verifica, non un ornamento.
[[nodiscard]] bool coefficient_is_zero(ExprPtr expr, symbolic::CASContext& ctx) {
    if (!expr) return true;
    auto expanded = expand(expr, ctx);
    ExprPtr candidate = expanded.is_ok() ? expanded.value() : expr;
    auto simplified = ctx.simplify(candidate);
    if (simplified.is_error()) return false;
    return poly_is_zero_expr(simplified.value());
}

}  // namespace

std::size_t count_real_roots_rational(RatPoly poly) {
    normalize_rational_coefficients(poly);
    if (poly.is_zero() || poly.degree() == 0U) return 0U;

    RatPoly derivative = derivative_rational(poly);
    if (derivative.is_zero()) return 0U;

    std::vector<RatPoly> chain;
    chain.push_back(poly);
    chain.push_back(std::move(derivative));

    while (true) {
        const RatPoly& previous = chain[chain.size() - 2U];
        const RatPoly current = chain.back();
        if (current.is_zero() || current.degree() == 0U) break;
        auto [quotient, remainder] = div_rem_rational_poly(previous, current);
        (void)quotient;
        if (remainder.is_zero()) {
            // gcd(p, p') non banale: p non e' squarefree. Sturm conta le radici
            // distinte solo sulla parte squarefree, quindi si riparte da p/gcd.
            auto [square_free, rest] = div_rem_rational_poly(poly, current);
            if (rest.is_zero() && square_free.degree() < poly.degree()) {
                return count_real_roots_rational(square_free);
            }
            break;
        }
        chain.push_back(negate_rational(remainder));
    }

    // Variazioni di segno a ±∞: dipendono solo dal segno del coefficiente di
    // testa e dalla parita' del grado, quindi il conteggio e' esatto.
    auto sign_variations = [&chain](bool at_plus_infinity) {
        int previous_sign = 0;
        std::size_t variations = 0U;
        for (const RatPoly& q : chain) {
            if (q.is_zero()) continue;
            int sign = (q.leading_coeff() < rat(0)) ? -1 : 1;
            if (!at_plus_infinity && (q.degree() % 2U) == 1U) sign = -sign;
            if (previous_sign != 0 && sign != previous_sign) ++variations;
            previous_sign = sign;
        }
        return variations;
    };

    const std::size_t at_minus = sign_variations(false);
    const std::size_t at_plus = sign_variations(true);
    return at_minus > at_plus ? at_minus - at_plus : 0U;
}

Result<std::vector<std::pair<ExprPtr, ExprPtr>>> real_quadratic_factors_of_quartic(
    const PolyExpr& R_z,
    const Symbol& z_var,
    symbolic::CASContext& ctx) {
    using Pairs = std::vector<std::pair<ExprPtr, ExprPtr>>;
    Pairs none;

    auto rational = poly_to_rational_poly(R_z);
    if (rational.is_error()) return ok(none);
    RatPoly R = rational.value();
    normalize_rational_coefficients(R);
    if (R.degree() != 4U) return ok(none);

    // Radici reali di R ⇒ i due fattori quadratici non hanno entrambi residui
    // coniugati e la forma chiusa richiederebbe i radicali del quartico pieno.
    if (count_real_roots_rational(R) != 0U) return ok(none);

    const Rational lc = R.leading_coeff();
    const Rational p = R[3] / lc;
    const Rational q = R[2] / lc;
    const Rational r = R[1] / lc;
    const Rational s = R[0] / lc;

    // Cubico risolvente: y³ − q·y² + (p·r − 4s)·y − (p²s − 4qs + r²).
    // Una sua radice razionale y₀ = f + h esiste esattamente quando il quartico
    // si spezza in due quadratici su un'estensione quadratica di Q.
    RatPoly cubic;
    cubic.push_back(-(p * p * s - rat(4) * q * s + r * r));
    cubic.push_back(p * r - rat(4) * s);
    cubic.push_back(-q);
    cubic.push_back(rat(1));

    BigInt scale(1);
    for (const Rational& c : cubic.coefficients()) scale = scale * c.denominator();
    IntPoly integer_cubic;
    integer_cubic.reserve(cubic.size());
    for (const Rational& c : cubic.coefficients()) {
        Rational scaled = c * Rational(scale);
        if (!(scaled.denominator() == BigInt(1))) return ok(none);
        integer_cubic.push_back(scaled.numerator());
    }
    normalize_integer_poly(integer_cubic);

    auto cubic_expr = integer_coefficients_to_expr(integer_cubic, z_var, ctx);
    if (cubic_expr.is_error()) return ok(none);
    auto factored = factor_over_integers(cubic_expr.value(), z_var, ctx);
    if (factored.is_error()) return ok(none);

    std::vector<Rational> rational_roots;
    for (const auto& factor : factored.value().factors) {
        auto factor_poly = parse_polynomial(factor.factor, z_var, ctx);
        if (factor_poly.is_error()) continue;
        if (poly_degree(factor_poly.value()) != 1U) continue;
        auto factor_rat = poly_to_rational_poly(factor_poly.value());
        if (factor_rat.is_error()) continue;
        const Rational lead = factor_rat.value().leading_coeff();
        if (lead == rat(0)) continue;
        rational_roots.push_back(-factor_rat.value().constant_term() / lead);
    }
    if (rational_roots.empty()) return ok(none);

    AstArena& arena = ctx.arena();
    ExprPtr one = poly_make_integer(arena, 1);
    ExprPtr two = poly_make_integer(arena, 2);
    PolyExpr monic_R({
        make_rational_expr(arena, s),
        make_rational_expr(arena, r),
        make_rational_expr(arena, q),
        make_rational_expr(arena, p),
        one});

    auto half_of = [&](ExprPtr sum_expr) -> Result<ExprPtr> {
        return ctx.simplify(arena.make<Binary>(BinaryOp::Div, sum_expr, two));
    };

    for (const Rational& y0 : rational_roots) {
        // e + g = p,  e·g = q − y₀  ⇒  e, g = (p ± √(p²−4q+4y₀))/2
        // f + h = y₀, f·h = s       ⇒  f, h = (y₀ ± √(y₀²−4s))/2
        // Entrambi i radicandi devono essere ≥ 0: sono RAZIONALI, quindi la
        // decisione e' esatta e non richiede alcun test di segno simbolico.
        const Rational d1 = p * p - rat(4) * q + rat(4) * y0;
        const Rational d2 = y0 * y0 - rat(4) * s;
        if (d1 < rat(0) || d2 < rat(0)) continue;

        auto sqrt_d1 = ctx.simplify(arena.make<FuncCall>(BuiltinOp::Sqrt,
            std::vector<ExprPtr>{make_rational_expr(arena, d1)}));
        auto sqrt_d2 = ctx.simplify(arena.make<FuncCall>(BuiltinOp::Sqrt,
            std::vector<ExprPtr>{make_rational_expr(arena, d2)}));
        if (sqrt_d1.is_error() || sqrt_d2.is_error()) continue;

        ExprPtr p_expr = make_rational_expr(arena, p);
        ExprPtr y0_expr = make_rational_expr(arena, y0);
        auto e = half_of(arena.make<Sum>(std::vector<ExprPtr>{p_expr, sqrt_d1.value()}));
        auto g = half_of(arena.make<Binary>(BinaryOp::Sub, p_expr, sqrt_d1.value()));
        auto f = half_of(arena.make<Sum>(std::vector<ExprPtr>{y0_expr, sqrt_d2.value()}));
        auto h = half_of(arena.make<Binary>(BinaryOp::Sub, y0_expr, sqrt_d2.value()));
        if (e.is_error() || g.is_error() || f.is_error() || h.is_error()) continue;

        // L'accoppiamento (e,f)·(g,h) oppure (e,h)·(g,f) non e' deducibile dai
        // soli e+g, f+h: si verifica per espansione, mai si assume.
        const std::pair<ExprPtr, ExprPtr> pairings[2][2] = {
            {{e.value(), f.value()}, {g.value(), h.value()}},
            {{e.value(), h.value()}, {g.value(), f.value()}},
        };
        for (const auto& pairing : pairings) {
            PolyExpr first({pairing[0].second, pairing[0].first, one});
            PolyExpr second({pairing[1].second, pairing[1].first, one});
            auto product = poly_multiply(first, second, ctx);
            if (product.is_error()) continue;
            auto difference = poly_subtract(product.value(), monic_R, ctx);
            if (difference.is_error()) continue;
            bool identical = true;
            for (const ExprPtr& coefficient : difference.value().coefficients()) {
                if (!coefficient_is_zero(coefficient, ctx)) { identical = false; break; }
            }
            if (!identical) continue;
            Pairs result;
            result.emplace_back(pairing[0].first, pairing[0].second);
            result.emplace_back(pairing[1].first, pairing[1].second);
            return ok(result);
        }
    }

    return ok(none);
}

}  // namespace cas::algebra
