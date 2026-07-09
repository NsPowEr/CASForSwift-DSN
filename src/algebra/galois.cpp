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
#include "galois_internal.hpp"
#include "polynomial_internal.hpp"

namespace cas::algebra {

std::optional<Rational> as_rational_q(ExprPtr e) {
    if (!e) return std::nullopt;
    if (const auto* il = expr_cast<IntegerLit>(e))
        return Rational(il->value, BigInt(1));
    if (const auto* rl = expr_cast<RationalLit>(e))
        return Rational(rl->numerator, rl->denominator);
    if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
        auto inner = as_rational_q(un->operand);
        if (inner) return -(*inner);
    }
    return std::nullopt;
}

bool is_rational_square_q(const Rational& r) {
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

namespace {

// Legacy in-file aliases preserve the deg-2..4 code below verbatim.
[[nodiscard]] inline std::optional<Rational> as_rational(ExprPtr e) {
    return as_rational_q(e);
}
[[nodiscard]] inline bool is_rational_square(const Rational& r) {
    return is_rational_square_q(r);
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
    if (total_deg == 5U) {
        // Detect irreducible quintic: exactly one non-constant factor and
        // that factor has degree 5 (multiplicity 1, else discriminant=0).
        bool irreducible_quintic = false;
        std::size_t non_constant_factors = 0U;
        std::size_t max_factor_deg = 0U;
        for (const auto& f : factors) {
            auto pp = parse_polynomial(f.factor, var, ctx);
            if (pp.is_error()) continue;
            std::size_t d = poly_degree(pp.value());
            if (d >= 1U) {
                non_constant_factors += 1U;
                if (d > max_factor_deg) max_factor_deg = d;
            }
        }
        if (non_constant_factors == 1U && max_factor_deg == 5U) {
            irreducible_quintic = true;
        }
        if (irreducible_quintic) {
            return galois_group_quintic_irreducible(poly, var, ctx);
        }
        // Reducible quintic: Galois group = direct product of Galois groups
        // of irreducible factors (over the algebraic closure of Q each factor
        // splits independently).  Recurse into each non-linear factor.
        // HC-F36-REDUCIBLE-COARSE closure: fine-grained labels via recursion.
        if (linear_count >= 5U) return ok(std::string("trivial"));
        std::vector<std::string> sub_labels;
        for (const auto& f : factors) {
            auto pp = parse_polynomial(f.factor, var, ctx);
            if (pp.is_error()) continue;
            std::size_t d = poly_degree(pp.value());
            if (d == 0U) continue;
            for (unsigned int mi = 0U; mi < f.multiplicity; ++mi) {
                if (d == 1U) continue;  // linear factor → trivial Galois.
                auto sub = galois_group(f.factor, var, ctx);
                if (sub.is_error()) {
                    return ok(std::string("reducible"));
                }
                if (sub.value() != std::string("trivial")) {
                    sub_labels.push_back(sub.value());
                }
            }
        }
        if (sub_labels.empty()) return ok(std::string("trivial"));
        if (sub_labels.size() == 1U) return ok(sub_labels[0]);
        std::string joined;
        for (std::size_t i = 0U; i < sub_labels.size(); ++i) {
            if (i > 0U) joined += " x ";
            joined += sub_labels[i];
        }
        return ok(joined);
    }
    if (total_deg == 6U || total_deg == 7U) {
        // A6 — an irreducible degree-6/7 input goes through the exact resolvent
        // pipeline (galois_deg6.cpp); reducible ones recurse into factors
        // exactly like the deg-5 path above.
        std::size_t non_constant_factors = 0U;
        std::size_t max_factor_deg = 0U;
        for (const auto& f : factors) {
            auto pp = parse_polynomial(f.factor, var, ctx);
            if (pp.is_error()) continue;
            std::size_t d = poly_degree(pp.value());
            if (d >= 1U) {
                non_constant_factors += 1U;
                if (d > max_factor_deg) max_factor_deg = d;
            }
        }
        if (non_constant_factors == 1U && max_factor_deg == total_deg) {
            return galois_group_irreducible_resolvent(poly, var, ctx);
        }
        if (linear_count >= total_deg) return ok(std::string("trivial"));
        std::vector<std::string> sub_labels;
        for (const auto& f : factors) {
            auto pp = parse_polynomial(f.factor, var, ctx);
            if (pp.is_error()) continue;
            std::size_t d = poly_degree(pp.value());
            if (d <= 1U) continue;
            for (unsigned int mi = 0U; mi < f.multiplicity; ++mi) {
                auto sub = galois_group(f.factor, var, ctx);
                if (sub.is_error()) return ok(std::string("reducible"));
                if (sub.value() != std::string("trivial")) {
                    sub_labels.push_back(sub.value());
                }
            }
        }
        if (sub_labels.empty()) return ok(std::string("trivial"));
        if (sub_labels.size() == 1U) return ok(sub_labels[0]);
        std::string joined;
        for (std::size_t i = 0U; i < sub_labels.size(); ++i) {
            if (i > 0U) joined += " x ";
            joined += sub_labels[i];
        }
        return ok(joined);
    }
    return ok(std::string("unknown"));
}

}  // namespace cas::algebra
