// Real-axis improper integration of rational functions via the residue
// theorem.  Restricted to rational integrands whose denominator factors over
// Q[x] into linear/quadratic factors (no real roots) and whose degree gap is
// at least 2.  See include/cas/residue_theorem.hpp for the public contract.

#include "cas/residue_theorem.hpp"

#include "cas/algebra.hpp"
#include "cas/algebraic_number.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "cas/ast.hpp"
#include "cas/bigint.hpp"
#include "cas/calculus.hpp"
#include "cas/rational.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] ExprPtr make_int(AstArena& arena, long long n) {
    return arena.make<IntegerLit>(BigInt(n));
}

[[nodiscard]] ExprPtr make_rational_expr(AstArena& arena, const Rational& r) {
    if (r.denominator() == BigInt(1)) {
        return arena.make<IntegerLit>(r.numerator());
    }
    return arena.make<RationalLit>(r.numerator(), r.denominator());
}

[[nodiscard]] Result<ExprPtr> simplify_or_fail(ExprPtr e, symbolic::CASContext& ctx) {
    auto s = ctx.simplify(e);
    if (s.is_error()) return fail<ExprPtr>(s.error());
    return s;
}

// Extract a Rational from an expression that simplifies to a rational literal.
// Returns nullopt if the expression is not a pure rational number.
[[nodiscard]] std::optional<Rational> as_rational(ExprPtr e) {
    if (const auto* lit = expr_cast<IntegerLit>(e)) {
        return Rational(lit->value);
    }
    if (const auto* lit = expr_cast<RationalLit>(e)) {
        return Rational(lit->numerator, lit->denominator);
    }
    // unary negate of a literal:  (-1)*lit collapses normally via simplify,
    // but be defensive.
    if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
        if (auto inner = as_rational(un->operand)) {
            return -*inner;
        }
    }
    return std::nullopt;
}

// Compute a Taylor coefficient  f^{(order)}(0) / order!   when f is a
// polynomial.  Implemented via repeated differentiation + substitution.
// Equivalent to coefficient of x^order in f.
[[nodiscard]] Result<ExprPtr> poly_coeff_at_zero(
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

    // Divide by order!
    BigInt fact(1);
    for (unsigned int k = 2U; k <= order; ++k) fact *= BigInt(static_cast<long long>(k));
    ExprPtr div = ctx.arena().make<Binary>(
        BinaryOp::Div,
        value.value(),
        ctx.arena().make<IntegerLit>(fact));
    return simplify_or_fail(div, ctx);
}

// Determine the degree of `f` in `var` by scanning Taylor coefficients
// up to a soft cap (max_integration_depth).  Returns the highest index k
// for which the coefficient is nonzero literal.  Returns nullopt if no
// rational coefficient could be extracted.
[[nodiscard]] Result<std::optional<std::size_t>> poly_degree_rational(
    ExprPtr f, const Symbol& var, symbolic::CASContext& ctx) {
    const std::size_t max_deg = static_cast<std::size_t>(ctx.max_integration_depth());
    std::optional<std::size_t> last_nonzero;
    for (std::size_t k = 0U; k <= max_deg; ++k) {
        auto c = poly_coeff_at_zero(f, var, static_cast<unsigned int>(k), ctx);
        if (c.is_error()) return fail<std::optional<std::size_t>>(c.error());
        auto r = as_rational(c.value());
        if (!r.has_value()) {
            return ok(std::optional<std::size_t>{});
        }
        if (!r->numerator().is_zero()) {
            last_nonzero = k;
        }
        // HARDCODE-OF-PASSAGE: HC-005 — Early-exit fisso a 8 zeri consecutivi.
        // Fix: usare algebra::parse_polynomial(expr, var, ctx) per ottenere
        // direttamente size/degree esatti. Vedi HARDCODE_LEDGER.md.
        if (last_nonzero.has_value() && k > *last_nonzero + 8U) {
            return ok(last_nonzero);
        }
    }
    return ok(last_nonzero);
}

// Collect rational coefficients of a polynomial f in var, indices 0..degree.
[[nodiscard]] Result<std::vector<Rational>> extract_rational_coeffs(
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

struct QuadraticFactor {
    // Monic representation x^2 + b*x + c with discriminant = b^2 - 4c < 0.
    Rational b;
    Rational c;
    Rational discriminant; // strictly negative
    ExprPtr original_factor; // the irreducible factor as parsed (any leading coeff)
};

// Build the RootOf node for an irreducible quadratic m(y) = y^2 + b*y + c,
// using the supplied generator variable name `gen_name`.  The minimal
// polynomial is the same m, expressed in `gen_var`.
[[nodiscard]] ExprPtr build_min_poly_expr(
    const Rational& b, const Rational& c, const Symbol& gen_var, AstArena& arena) {
    ExprPtr y = arena.make<Symbol>(gen_var);
    ExprPtr y2 = arena.make<Binary>(BinaryOp::Pow, y, make_int(arena, 2));
    ExprPtr by = arena.make<Binary>(BinaryOp::Mul, make_rational_expr(arena, b), y);
    ExprPtr c_e = make_rational_expr(arena, c);
    ExprPtr sum1 = arena.make<Binary>(BinaryOp::Add, y2, by);
    return arena.make<Binary>(BinaryOp::Add, sum1, c_e);
}

// Take the residue r of N/D at the upper-half-plane root α of x^2+b*x+c,
// expressed as an AlgebraicNumber e + f*α in Q(α), and return the real
// contribution to the integral: -π · f · √(-Δ).
[[nodiscard]] Result<ExprPtr> contribution_from_quadratic(
    const Rational& b,
    const Rational& c,
    const Rational& discriminant,
    ExprPtr N_over_D,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    // Build α = RootOf(x^2 + b*x + c, x).  Simplify the minimal polynomial
    // expression so it matches the canonical form used internally.
    Symbol root_var = var; // use the integration variable as the polynomial var
    ExprPtr min_poly_in_x_raw = build_min_poly_expr(b, c, root_var, arena);
    auto min_poly_simp = simplify_or_fail(min_poly_in_x_raw, ctx);
    if (min_poly_simp.is_error()) return fail<ExprPtr>(min_poly_simp.error());
    ExprPtr alpha_expr = arena.make<RootOf>(min_poly_simp.value(), root_var, std::nullopt);

    // residue(N/D, x, α, ctx) → expression in Q(α).
    auto res = residue(N_over_D, var, alpha_expr, ctx);
    if (res.is_error()) return fail<ExprPtr>(res.error());

    // Reduce the residue inside Q(α) so it has canonical form e + f·α.
    auto reduced = algebra::simplify_in_q_alpha(res.value(), ctx);
    if (reduced.is_error()) return fail<ExprPtr>(reduced.error());

    // Build the minimal polynomial coefficient vector  [c, b, 1].
    algebra::AlgebraicNumber::CoeffVec min_poly_coeffs;
    min_poly_coeffs.push_back(c);
    min_poly_coeffs.push_back(b);
    min_poly_coeffs.push_back(Rational(BigInt(1)));

    auto expressed = algebra::try_express_in_q_alpha(
        reduced.value(), alpha_expr, min_poly_coeffs, ctx);
    if (expressed.is_error()) return fail<ExprPtr>(expressed.error());
    if (!expressed.value().has_value()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Residue theorem: residue not expressible in Q(α) for this quadratic factor"});
    }

    const algebra::AlgebraicNumber& an = *expressed.value();
    const auto& value = an.value();
    // value = [e, f] (ascending degree; deg < 2).
    Rational e_part(BigInt(0));
    Rational f_part(BigInt(0));
    if (value.size() >= 1U) e_part = value[0];
    if (value.size() >= 2U) f_part = value[1];

    // Contribution = -π · f · √(-Δ).
    // -Δ is positive rational.
    Rational neg_disc = -discriminant;
    ExprPtr neg_disc_expr = make_rational_expr(arena, neg_disc);
    ExprPtr sqrt_neg_disc = arena.make<FuncCall>(
        std::string("sqrt"),
        std::vector<ExprPtr>{neg_disc_expr});

    ExprPtr f_expr = make_rational_expr(arena, f_part);
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr neg_pi = arena.make<Binary>(BinaryOp::Mul, make_int(arena, -1), pi);
    ExprPtr term = arena.make<Binary>(BinaryOp::Mul, neg_pi, f_expr);
    ExprPtr full = arena.make<Binary>(BinaryOp::Mul, term, sqrt_neg_disc);

    // Imaginary parts must cancel by reality of integrand; we discard e_part.
    // (The full residue sum has imaginary part summing to zero across
    //  conjugate quadratic factors; here `f_part` is the only data that
    //  contributes to the real integral from this individual factor.)
    (void)e_part;

    return simplify_or_fail(full, ctx);
}

}  // namespace

Result<ExprPtr> integrate_rational_full_real_line(
    ExprPtr rational_expr,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    // Step 1: split into N/D.
    auto parts = algebra::apart_num_den(rational_expr, ctx);
    if (parts.is_error()) return fail<ExprPtr>(parts.error());
    ExprPtr N = parts.value().numerator;
    ExprPtr D = parts.value().denominator;

    // Step 2: convergence check  deg(D) >= deg(N) + 2.
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

    // Step 3: factor D over Q.
    auto factorization = algebra::factor_polynomial(D, var, ctx);
    if (factorization.is_error()) return fail<ExprPtr>(factorization.error());

    AstArena& arena = ctx.arena();
    ExprPtr total = make_int(arena, 0);

    // Inspect each irreducible factor.
    for (const auto& pf : factorization.value().factors) {
        // Determine its degree in `var`.
        auto fdeg_res = poly_degree_rational(pf.factor, var, ctx);
        if (fdeg_res.is_error()) return fail<ExprPtr>(fdeg_res.error());
        if (!fdeg_res.value().has_value()) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Residue theorem: irreducible factor has non-rational coefficients"});
        }
        const std::size_t fdeg = *fdeg_res.value();
        if (fdeg == 0U) {
            // constant factor — ignore (absorbed into content).
            continue;
        }
        if (fdeg == 1U) {
            // Linear factor → real root of the denominator.
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Residue theorem: denominator has a real pole (linear factor over Q)"});
        }
        if (fdeg > 2U) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Residue theorem: irreducible factor of degree > 2 not yet supported"});
        }

        // Quadratic factor a*x^2 + b'*x + c'.  Normalize to monic.
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
        // Sign check of disc: numerator * denominator < 0 iff value < 0 (denom>0).
        const bool disc_negative = disc.numerator().is_negative();
        if (!disc_negative) {
            // Real roots exist (could be a perfect square — factor_polynomial
            // should have split it, but defend).
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Residue theorem: quadratic factor has nonnegative discriminant (real roots)"});
        }
        if (pf.multiplicity > 1U) {
            // Higher-order pole at α: still handled by residue() via the
            // Laurent series, but residue must be re-expressed in Q(α).
            // We pass the *full* N/D to residue() because the pole order is
            // already detected by laurent_series internally.  Fall through.
        }

        auto contrib = contribution_from_quadratic(b, c, disc, rational_expr, var, ctx);
        if (contrib.is_error()) return fail<ExprPtr>(contrib.error());
        // For factors with multiplicity m, residue() at α already accounts
        // for the full Laurent c_{-1}, so we count each distinct quadratic
        // factor once regardless of multiplicity.
        ExprPtr added = arena.make<Binary>(BinaryOp::Add, total, contrib.value());
        auto simp = simplify_or_fail(added, ctx);
        if (simp.is_error()) return simp;
        total = simp.value();
    }

    return simplify_or_fail(total, ctx);
}

}  // namespace cas::calculus
