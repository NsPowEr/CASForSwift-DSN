// CAS-L3-18 — Galois group identification for low-degree polynomials.
//
// MVP scope: degree 2 and 3 via discriminant analysis.
//   Deg 2: D = b² - 4ac. If rational square → roots in Q → trivial group;
//          else → C₂ (quadratic extension).
//   Deg 3: D = -4p³ - 27q² (depressed cubic). If polynomial splits in Q
//          (has rational root + remaining factors split): trivial.
//          Else, if D is rational square: A₃ (cyclic of 3).
//          Else: S₃ (full symmetric).
// Degree 4: deferred follow-up (resolvent cubic + discriminant joint).

#include "cas/galois.hpp"

#include "cas/algebra.hpp"
#include "cas/error.hpp"
#include "polynomial_internal.hpp"

namespace cas::algebra {

namespace {

// Try to interpret expr as Rational. Returns nullopt otherwise.
[[nodiscard]] std::optional<Rational> as_rational(ExprPtr e) {
    if (!e) return std::nullopt;
    if (const auto* il = expr_cast<IntegerLit>(e))
        return Rational(il->value, BigInt(1));
    if (const auto* rl = expr_cast<RationalLit>(e))
        return Rational(rl->numerator, rl->denominator);
    if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
        auto inner = as_rational(un->operand);
        if (inner) return -(*inner);
    }
    return std::nullopt;
}

// Returns true if rational is a perfect square in Q (i.e. p/q with
// p, q both perfect squares as BigInt and same sign).
[[nodiscard]] bool is_rational_square(const Rational& r) {
    if (r.numerator().is_negative()) return false;
    if (r.numerator().is_zero()) return true;
    auto isqrt = [](const BigInt& n) -> std::optional<BigInt> {
        if (n.is_zero()) return BigInt(0);
        BigInt x = n;
        BigInt y = (x + BigInt(1)) / BigInt(2);
        while (y < x) {
            x = y;
            y = (x + n / x) / BigInt(2);
        }
        if (x * x == n) return x;
        return std::nullopt;
    };
    auto num_sqrt = isqrt(r.numerator());
    auto den_sqrt = isqrt(r.denominator());
    return num_sqrt.has_value() && den_sqrt.has_value();
}

}  // namespace

Result<std::string> galois_group(ExprPtr poly, const Symbol& var,
                                  symbolic::CASContext& ctx) {
    // Compute discriminant.
    auto disc_res = polynomial_discriminant(poly, var, ctx);
    if (disc_res.is_error()) return fail<std::string>(disc_res.error());
    auto disc_simp = ctx.simplify(disc_res.value());
    if (disc_simp.is_error()) return fail<std::string>(disc_simp.error());
    auto disc_rat = as_rational(disc_simp.value());
    if (!disc_rat) {
        return ok(std::string("unknown"));
    }

    // Factor over Q to detect rational roots / splitting.
    auto fact_res = factor_over_integers(poly, var, ctx);
    if (fact_res.is_error()) return fail<std::string>(fact_res.error());
    const auto& factors = fact_res.value().factors;

    // Count distinct linear factors (degree 1) and total factor count.
    std::size_t linear_count = 0;
    std::size_t total_deg = 0;
    for (const auto& f : factors) {
        auto pp = parse_polynomial(f.factor, var, ctx);
        if (pp.is_error()) continue;
        std::size_t d = poly_degree(pp.value());
        total_deg += d * f.multiplicity;
        if (d == 1U) linear_count += f.multiplicity;
    }

    // Degree analysis.
    if (total_deg == 2U) {
        if (linear_count >= 2U) return ok(std::string("trivial"));
        // Quadratic irreducible → Galois = C₂.
        return ok(std::string("C2"));
    }
    if (total_deg == 3U) {
        if (linear_count >= 3U) return ok(std::string("trivial"));
        if (linear_count == 1U) {
            // Linear · quadratic. Quadratic might be reducible or not.
            // If discriminant of quadratic part is rational square (already
            // captured by factor_over_integers splitting further), this
            // path means quadratic irreducible. Galois = C₂.
            return ok(std::string("C2"));
        }
        // Irreducible cubic. Group depends on disc:
        //   D rational square → A₃ (cyclic order 3)
        //   else → S₃
        if (is_rational_square(*disc_rat)) return ok(std::string("A3"));
        return ok(std::string("S3"));
    }
    if (total_deg == 4U) {
        if (linear_count >= 4U) return ok(std::string("trivial"));
        if (linear_count >= 1U) {
            // Mix lineari + non-linear → fattorizza già parziale, gruppo C₂.
            return ok(std::string("C2"));
        }
        // Irreducible quartic. Build resolvent cubic + factor over Q.
        //   1. Parse poly as RatPoly degree 4.
        //   2. Monicize + depress: y⁴ + p·y² + q·y + r.
        //   3. Resolvent: R(z) = z³ - 2p·z² + (p²-4r)·z + q².
        //   4. Factor R over Q.
        //   5. Dispatch via Buhler-Reverter:
        //        R irreducibile  + disc square → A₄
        //        R irreducibile  + disc NON sq → S₄
        //        R reducibile    + disc square → V₄
        //        R reducibile    + disc NON sq → D₄
        bool d_is_square = is_rational_square(*disc_rat);

        // Parse poly as RatPoly.
        auto poly_parsed = parse_polynomial(poly, var, ctx);
        if (poly_parsed.is_error()) {
            return ok(d_is_square ? std::string("V4") : std::string("S4"));
        }
        auto rat_coeffs = poly_to_rational_coefficients(poly_parsed.value());
        if (rat_coeffs.is_error() || rat_coeffs.value().size() != 5U) {
            return ok(d_is_square ? std::string("V4") : std::string("S4"));
        }
        // coeffs[0]=e, [1]=d, [2]=c, [3]=b, [4]=a (low→high).
        Rational a = rat_coeffs.value()[4];
        Rational b = rat_coeffs.value()[3] / a;
        Rational c = rat_coeffs.value()[2] / a;
        Rational d = rat_coeffs.value()[1] / a;
        Rational e = rat_coeffs.value()[0] / a;
        // Depressed: y = x - b/4
        // p = c - 3b²/8, q = d - bc/2 + b³/8, r = e - bd/4 + b²c/16 - 3b⁴/256
        Rational eight(BigInt(8), BigInt(1));
        Rational half(BigInt(1), BigInt(2));
        Rational quarter(BigInt(1), BigInt(4));
        Rational sixteen(BigInt(16), BigInt(1));
        Rational two_five_six(BigInt(256), BigInt(1));
        Rational three(BigInt(3), BigInt(1));
        Rational b2 = b * b;
        Rational b3 = b2 * b;
        Rational b4 = b2 * b2;
        Rational p = c - (three * b2) / eight;
        Rational q = d - (b * c) * half + b3 / eight;
        Rational r = e - (b * d) * quarter + (b2 * c) / sixteen - (three * b4) / two_five_six;
        // Resolvent: R(z) = z³ - 2p z² + (p²-4r) z + q²
        // Construct as ExprPtr.
        AstArena& arena = ctx.arena();
        auto rat_to_expr = [&](const Rational& r) -> ExprPtr {
            if (r.denominator() == BigInt(1))
                return arena.make<IntegerLit>(r.numerator());
            return arena.make<RationalLit>(r.numerator(), r.denominator());
        };
        Symbol z_sym = ctx.make_fresh_symbol("zgal");
        ExprPtr z_expr = arena.make<Symbol>(z_sym.name);
        ExprPtr z2 = arena.make<Binary>(BinaryOp::Pow, z_expr,
            arena.make<IntegerLit>(BigInt(2)));
        ExprPtr z3 = arena.make<Binary>(BinaryOp::Pow, z_expr,
            arena.make<IntegerLit>(BigInt(3)));
        Rational neg_two_p = -(Rational(BigInt(2), BigInt(1)) * p);
        Rational p_sq_minus_4r = p * p - Rational(BigInt(4), BigInt(1)) * r;
        Rational q_sq = q * q;
        ExprPtr resolvent = arena.make<Sum>(std::vector<ExprPtr>{
            z3,
            arena.make<Product>(std::vector<ExprPtr>{rat_to_expr(neg_two_p), z2}),
            arena.make<Product>(std::vector<ExprPtr>{rat_to_expr(p_sq_minus_4r), z_expr}),
            rat_to_expr(q_sq)
        });
        // Factor resolvent over Z.
        auto resolvent_fact = factor_over_integers(resolvent, z_sym, ctx);
        bool resolvent_reducibile = false;
        if (resolvent_fact.is_ok()) {
            for (const auto& f : resolvent_fact.value().factors) {
                auto fpoly = parse_polynomial(f.factor, z_sym, ctx);
                if (fpoly.is_ok() && poly_degree(fpoly.value()) == 1U) {
                    resolvent_reducibile = true;
                    break;
                }
            }
        }
        // Dispatch
        if (resolvent_reducibile && d_is_square) return ok(std::string("V4"));
        if (resolvent_reducibile && !d_is_square) return ok(std::string("D4"));
        if (!resolvent_reducibile && d_is_square) return ok(std::string("A4"));
        return ok(std::string("S4"));
    }
    return ok(std::string("unknown"));
}

}  // namespace cas::algebra
