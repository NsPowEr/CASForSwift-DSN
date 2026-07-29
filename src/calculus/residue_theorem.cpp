// Real-axis improper integration of rational functions via the residue
// theorem.  Handles rational integrands with no real poles and degree gap ≥ 2:
//   - irreducible quadratic factors (any multiplicity) via Q(α) reduction;
//   - irreducible biquadratic factors x⁴+bx²+c (any multiplicity) via the
//     Q(√c, √(2√c+b)) tower;
//   - all remaining irreducible factors (degree ≥ 3, general quartics) via the
//     numeric Aberth residue fallback.
// Higher pole orders are resolved by residue()'s Laurent recurrence; the Q(α)→ℝ
// assembly is a fixed linear functional of the residue coordinates, so it is
// independent of multiplicity.  Real (linear) poles are rejected — Cauchy PV is
// out of scope.  See include/cas/residue_theorem.hpp for the public contract.

#include "cas/residue_theorem.hpp"
#include "residue_theorem_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus {


ExprPtr make_rational_expr(AstArena& arena, const Rational& r) {
    if (r.denominator() == BigInt(1)) {
        return arena.make<IntegerLit>(r.numerator());
    }
    return arena.make<RationalLit>(r.numerator(), r.denominator());
}

Result<ExprPtr> simplify_or_fail(ExprPtr e, symbolic::CASContext& ctx) {
    auto s = ctx.simplify(e);
    if (s.is_error()) return fail<ExprPtr>(s.error());
    return s;
}

std::optional<Rational> as_rational(ExprPtr e) {
    if (const auto* lit = expr_cast<IntegerLit>(e)) {
        return Rational(lit->value);
    }
    if (const auto* lit = expr_cast<RationalLit>(e)) {
        return Rational(lit->numerator, lit->denominator);
    }
    if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
        if (auto inner = as_rational(un->operand)) {
            return -*inner;
        }
    }
    return std::nullopt;
}

Result<ExprPtr> poly_coeff_at_zero(
    ExprPtr f, const Symbol& var, unsigned int order, symbolic::CASContext& ctx) {
    ExprPtr derivative = f;
    if (order > 0U) {
        auto d = diff(f, var, order, ctx);
        if (d.is_error()) return fail<ExprPtr>(d.error());
        derivative = d.value();
    }
    auto subbed = ctx.substitute(derivative, var, make_int(ctx.arena(), 0));
    if (subbed.is_error()) return fail<ExprPtr>(subbed.error());
    auto value = simplify_or_fail(subbed.value(), ctx);
    if (value.is_error()) return value;
    if (order <= 1U) return value;

    BigInt fact(1);
    for (unsigned int k = 2U; k <= order; ++k) fact *= BigInt(static_cast<long long>(k));
    ExprPtr div = ctx.arena().make<Binary>(
        BinaryOp::Div,
        value.value(),
        ctx.arena().make<IntegerLit>(fact));
    return simplify_or_fail(div, ctx);
}

Result<std::optional<std::size_t>> poly_degree_rational(
    ExprPtr f, const Symbol& var, symbolic::CASContext& ctx) {
    auto parsed = algebra::parse_polynomial(f, var, ctx);
    if (parsed.is_error()) {
        return ok(std::optional<std::size_t>{});
    }
    const auto& poly = parsed.value();
    if (poly.is_zero()) {
        return ok(std::optional<std::size_t>{static_cast<std::size_t>(0U)});
    }
    for (const ExprPtr& c : poly.coefficients()) {
        auto r = as_rational(c);
        if (!r.has_value()) {
            return ok(std::optional<std::size_t>{});
        }
    }
    return ok(std::optional<std::size_t>{poly.degree()});
}

Result<std::vector<Rational>> extract_rational_coeffs(
    ExprPtr f,
    const Symbol& var,
    std::size_t degree,
    symbolic::CASContext& ctx) {
    std::vector<Rational> coeffs;
    coeffs.reserve(degree + 1U);
    for (std::size_t k = 0U; k <= degree; ++k) {
        auto c = poly_coeff_at_zero(f, var, static_cast<unsigned int>(k), ctx);
        if (c.is_error()) return fail<std::vector<Rational>>(c.error());
        auto r = as_rational(c.value());
        if (!r.has_value()) {
            return fail<std::vector<Rational>>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Residue theorem: polynomial factor has non-rational coefficient"});
        }
        coeffs.push_back(*r);
    }
    return ok(std::move(coeffs));
}

Result<ExprPtr> integrate_rational_full_real_line(
    ExprPtr rational_expr,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    auto parts = algebra::apart_num_den(rational_expr, ctx);
    if (parts.is_error()) return fail<ExprPtr>(parts.error());
    ExprPtr N = parts.value().numerator;
    ExprPtr D = parts.value().denominator;

    auto deg_N_res = poly_degree_rational(N, var, ctx);
    if (deg_N_res.is_error()) return fail<ExprPtr>(deg_N_res.error());
    auto deg_D_res = poly_degree_rational(D, var, ctx);
    if (deg_D_res.is_error()) return fail<ExprPtr>(deg_D_res.error());
    if (!deg_N_res.value().has_value() || !deg_D_res.value().has_value()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Residue theorem: numerator or denominator is not a rational polynomial in the integration variable"});
    }
    const std::size_t deg_N = *deg_N_res.value();
    const std::size_t deg_D = *deg_D_res.value();
    if (deg_D < deg_N + 2U) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Residue theorem: integral does not converge (deg(Q) < deg(P)+2)"});
    }

    auto factorization = algebra::factor_polynomial(D, var, ctx);
    if (factorization.is_error()) return fail<ExprPtr>(factorization.error());

    AstArena& arena = ctx.arena();
    ExprPtr total = make_int(arena, 0);

    for (const auto& pf : factorization.value().factors) {
        auto fdeg_res = poly_degree_rational(pf.factor, var, ctx);
        if (fdeg_res.is_error()) return fail<ExprPtr>(fdeg_res.error());
        if (!fdeg_res.value().has_value()) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Residue theorem: irreducible factor has non-rational coefficients"});
        }
        const std::size_t fdeg = *fdeg_res.value();
        if (fdeg == 0U) {
            continue;
        }
        if (fdeg == 1U) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Residue theorem: denominator has a real pole (linear factor over Q)"});
        }
        if (fdeg == 4U) {
            auto bq_coeffs = extract_rational_coeffs(pf.factor, var, 4U, ctx);
            if (bq_coeffs.is_error()) return fail<ExprPtr>(bq_coeffs.error());
            const auto& q = bq_coeffs.value();
            if (!q[1].numerator().is_zero() || !q[3].numerator().is_zero()) {
                auto numeric_contrib = numeric_residue_contribution(
                    pf, N, D, var, deg_N, deg_D, ctx);
                if (numeric_contrib.is_error()) return fail<ExprPtr>(numeric_contrib.error());
                ExprPtr added_n = arena.make<Binary>(BinaryOp::Add, total, numeric_contrib.value());
                auto simp_n = simplify_or_fail(added_n, ctx);
                if (simp_n.is_error()) return simp_n;
                total = simp_n.value();
                continue;
            }
            if (q[4].numerator().is_zero()) {
                return fail<ExprPtr>(CASError{
                    .kind = CASErrorKind::Unimplemented,
                    .message = "Residue theorem: degenerate quartic factor"});
            }
            const Rational b_norm = q[2] / q[4];
            const Rational c_norm = q[0] / q[4];
            if (c_norm.numerator().is_negative() || c_norm.numerator().is_zero()) {
                return fail<ExprPtr>(CASError{
                    .kind = CASErrorKind::Unimplemented,
                    .message = "Residue theorem: biquadratic factor with non‑positive constant term"});
            }
            const Rational aux_disc = b_norm * b_norm - Rational(BigInt(4)) * c_norm;
            if (!aux_disc.numerator().is_negative()) {
                return fail<ExprPtr>(CASError{
                    .kind = CASErrorKind::Unimplemented,
                    .message = "Residue theorem: biquadratic factor with non‑negative auxiliary discriminant (real roots in u)"});
            }
            // Multiplicity > 1 is handled uniformly: contribution_from_irreducible
            // _biquadratic delegates to residue(), whose Laurent recurrence detects
            // the pole order, and the Q(α)→ℝ assembly is a fixed linear functional
            // of the residue coordinates (independent of pole order). Verified
            // against Maxima on ∫1/(x⁴+1)² = 3π/2^(5/2) and ∫1/(x⁴+x²+1)² = 2π/3^(3/2).
            auto contrib = contribution_from_irreducible_biquadratic(b_norm, c_norm, rational_expr, var, ctx);
            if (contrib.is_error()) return fail<ExprPtr>(contrib.error());
            ExprPtr added = arena.make<Binary>(BinaryOp::Add, total, contrib.value());
            auto simp = simplify_or_fail(added, ctx);
            if (simp.is_error()) return simp;
            total = simp.value();
            continue;
        }
        if (fdeg > 2U) {
            auto numeric_contrib = numeric_residue_contribution(
                pf, N, D, var, deg_N, deg_D, ctx);
            if (numeric_contrib.is_error()) return fail<ExprPtr>(numeric_contrib.error());
            ExprPtr added_n = arena.make<Binary>(BinaryOp::Add, total, numeric_contrib.value());
            auto simp_n = simplify_or_fail(added_n, ctx);
            if (simp_n.is_error()) return simp_n;
            total = simp_n.value();
            continue;
        }

        auto coeffs_res = extract_rational_coeffs(pf.factor, var, 2U, ctx);
        if (coeffs_res.is_error()) return fail<ExprPtr>(coeffs_res.error());
        const auto& coeffs = coeffs_res.value();
        const Rational a = coeffs[2];
        if (a.numerator().is_zero()) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Residue theorem: degenerate quadratic factor"});
        }
        const Rational b = coeffs[1] / a;
        const Rational c = coeffs[0] / a;
        const Rational disc = b * b - Rational(BigInt(4)) * c;
        const bool disc_negative = disc.numerator().is_negative();
        if (!disc_negative) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Residue theorem: quadratic factor has nonnegative discriminant (real roots)"});
        }

        auto contrib = contribution_from_quadratic(b, c, disc, rational_expr, var, ctx);
        if (contrib.is_error()) return fail<ExprPtr>(contrib.error());
        ExprPtr added = arena.make<Binary>(BinaryOp::Add, total, contrib.value());
        auto simp = simplify_or_fail(added, ctx);
        if (simp.is_error()) return simp;
        total = simp.value();
    }

    return simplify_or_fail(total, ctx);
}

}  // namespace cas::calculus
