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
        // Deg 4 full analysis deferred: requires resolvent cubic Galois
        // recursion. Return unknown for now (anti-furbizia: no guess).
        return ok(std::string("unknown"));
    }
    return ok(std::string("unknown"));
}

}  // namespace cas::algebra
