// Real-axis improper integration of rational functions via the residue
// theorem.  Restricted to rational integrands whose denominator factors over
// Q[x] into linear/quadratic factors (no real roots) and whose degree gap is
// at least 2.  See include/cas/residue_theorem.hpp for the public contract.

#include "cas/residue_theorem.hpp"

#include "cas/algebra.hpp"
#include "cas/algebraic_number.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "cas/ast.hpp"
#include "cas/bigfloat.hpp"
#include "cas/bigint.hpp"
#include "cas/calculus.hpp"
#include "cas/numeric/complex_root_isolator.hpp"
#include "cas/rational.hpp"
#include "../algebra/polynomial_internal.hpp"
#include "../numeric/complex_bigfloat_internal.hpp"

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

// Determine the degree of `f` in `var` using the algebra layer's exact
// polynomial parser.  HC-005 closed: replaces the Taylor-coefficient scan
// + "8 trailing zeros" early-exit with a deterministic structural parse.
// Returns nullopt if `f` cannot be parsed as a polynomial in `var` with
// rational coefficients (e.g. coefficients depend on other symbols).
[[nodiscard]] Result<std::optional<std::size_t>> poly_degree_rational(
    ExprPtr f, const Symbol& var, symbolic::CASContext& ctx) {
    auto parsed = algebra::parse_polynomial(f, var, ctx);
    if (parsed.is_error()) {
        return ok(std::optional<std::size_t>{});
    }
    const auto& poly = parsed.value();
    if (poly.is_zero()) {
        return ok(std::optional<std::size_t>{static_cast<std::size_t>(0U)});
    }
    // Verify every coefficient is a rational literal — required by the
    // downstream extract_rational_coeffs contract.
    for (const ExprPtr& c : poly.coefficients()) {
        auto r = as_rational(c);
        if (!r.has_value()) {
            return ok(std::optional<std::size_t>{});
        }
    }
    return ok(std::optional<std::size_t>{poly.degree()});
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

// Build  y⁴ + b·y² + c  as ExprPtr.
[[nodiscard]] ExprPtr build_biquadratic_min_poly_expr(
    const Rational& b, const Rational& c, const Symbol& gen_var, AstArena& arena) {
    ExprPtr y = arena.make<Symbol>(gen_var);
    ExprPtr y2 = arena.make<Binary>(BinaryOp::Pow, y, make_int(arena, 2));
    ExprPtr y4 = arena.make<Binary>(BinaryOp::Pow, y, make_int(arena, 4));
    ExprPtr by2 = arena.make<Binary>(BinaryOp::Mul, make_rational_expr(arena, b), y2);
    ExprPtr c_e = make_rational_expr(arena, c);
    ExprPtr sum1 = arena.make<Binary>(BinaryOp::Add, y4, by2);
    return arena.make<Binary>(BinaryOp::Add, sum1, c_e);
}

// Closure for an irreducible biquadratic factor  D_q(x) = x⁴ + b·x² + c, with
// c > 0 (Rational) and Δ = b² − 4c < 0.  Such a quartic has four complex roots
// on a circle of radius c^{1/4} arranged as conjugate pairs (α, −α, ᾱ, −ᾱ).
// Two of them lie in the upper half plane; call them α₁ and α₂.
//
// Closed‑form derivation, working modulo y⁴ = −b·y² − c.
//
//   α₁·α₂            = −√c                              (sub‑Vieta)
//   α₁ + α₂          =  i·√(2√c + b)                    (geometric)
//   α₁² + α₂²        = −b                               (Vieta on u = α²)
//   α₁³ + α₂³        =  i·√(2√c + b) · (√c − b)
//
// For any residue r expressed in Q(α) as
//      r(α) = c₀ + c₁·α + c₂·α² + c₃·α³,
// the sum at the upper roots collapses to
//      Σ_upper r = (2c₀ − b·c₂) + i·√(2√c + b) · (c₁ + c₃·(√c − b)).
//
// Reality of the integrand forces the real part 2c₀ − b·c₂ to vanish across
// the full sum of contributions, so the real integral coming from THIS
// quartic factor is
//      contribution = 2πi · Σ_upper
//                   = −2π · √(2√c + b) · (c₁ + c₃·(√c − b))
// after combining 2πi · i = −2π.
//
// All sqrt arguments are positive rationals (c > 0 and 2√c + b > 0 because
// Δ < 0 ⇒ b² < 4c ⇒ |b| < 2√c).  The result therefore lives in the tower
// Q(√c, √(2√c + b)), which is exactly the two‑level extension targeted by
// the STEP A tower bridge.
[[nodiscard]] Result<ExprPtr> contribution_from_irreducible_biquadratic(
    const Rational& b,
    const Rational& c,
    ExprPtr N_over_D,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    // Build α = RootOf(x⁴ + b·x² + c).
    Symbol gen = var;
    ExprPtr min_poly_raw = build_biquadratic_min_poly_expr(b, c, gen, arena);
    auto min_poly_simp = simplify_or_fail(min_poly_raw, ctx);
    if (min_poly_simp.is_error()) return fail<ExprPtr>(min_poly_simp.error());
    ExprPtr alpha_expr = arena.make<RootOf>(min_poly_simp.value(), gen, std::nullopt);

    // residue(N/D, x, α).  Returns an expression in Q(α).
    auto res = residue(N_over_D, var, alpha_expr, ctx);
    if (res.is_error()) return fail<ExprPtr>(res.error());

    auto reduced = algebra::simplify_in_q_alpha(res.value(), ctx);
    if (reduced.is_error()) return fail<ExprPtr>(reduced.error());

    algebra::AlgebraicNumber::CoeffVec min_poly_coeffs;
    min_poly_coeffs.push_back(c);                 // x⁰
    min_poly_coeffs.push_back(Rational(BigInt(0)));// x¹
    min_poly_coeffs.push_back(b);                 // x²
    min_poly_coeffs.push_back(Rational(BigInt(0)));// x³
    min_poly_coeffs.push_back(Rational(BigInt(1)));// x⁴

    auto expressed = algebra::try_express_in_q_alpha(
        reduced.value(), alpha_expr, min_poly_coeffs, ctx);
    if (expressed.is_error()) return fail<ExprPtr>(expressed.error());
    if (!expressed.value().has_value()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Residue theorem: residue not expressible in Q(α) for this biquadratic factor"});
    }

    const algebra::AlgebraicNumber& an = *expressed.value();
    const auto& value = an.value();
    Rational c0(BigInt(0)), c1(BigInt(0)), c2(BigInt(0)), c3(BigInt(0));
    if (value.size() >= 1U) c0 = value[0];
    if (value.size() >= 2U) c1 = value[1];
    if (value.size() >= 3U) c2 = value[2];
    if (value.size() >= 4U) c3 = value[3];

    // Real part across the upper pair: 2c₀ − b·c₂.  For this single factor
    // it must be zero in isolation when the integrand has only one quartic
    // factor; in mixed cases the cancellation happens at the sum level.  We
    // therefore keep the real part in the output as a real-valued correction
    // (multiplied by 2πi it would be purely imaginary, so reality of the
    // overall integral guarantees its cancellation).
    //
    // contribution_real = −2π · √(2√c + b) · (c₁ + c₃·(√c − b)).
    ExprPtr c_expr = make_rational_expr(arena, c);
    ExprPtr sqrt_c = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{c_expr});

    ExprPtr two_sqrt_c = arena.make<Binary>(BinaryOp::Mul, make_int(arena, 2), sqrt_c);
    ExprPtr b_expr = make_rational_expr(arena, b);
    ExprPtr radicand = arena.make<Binary>(BinaryOp::Add, two_sqrt_c, b_expr);
    ExprPtr sqrt_radicand = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{radicand});

    ExprPtr sqrt_c_minus_b = arena.make<Binary>(BinaryOp::Sub, sqrt_c, b_expr);
    ExprPtr c3_expr = make_rational_expr(arena, c3);
    ExprPtr c3_times = arena.make<Binary>(BinaryOp::Mul, c3_expr, sqrt_c_minus_b);
    ExprPtr c1_expr = make_rational_expr(arena, c1);
    ExprPtr inner_sum = arena.make<Binary>(BinaryOp::Add, c1_expr, c3_times);

    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr neg_two_pi = arena.make<Binary>(BinaryOp::Mul, make_int(arena, -2), pi);
    ExprPtr part1 = arena.make<Binary>(BinaryOp::Mul, neg_two_pi, sqrt_radicand);
    ExprPtr full = arena.make<Binary>(BinaryOp::Mul, part1, inner_sum);

    (void)c0;  // imaginary‑side coefficient; cancels in the real integral.
    (void)c2;
    return simplify_or_fail(full, ctx);
}

// ── F5.6 sub-block 2: numeric residue contribution via Aberth ───────────────
//
// When a denominator factor escapes the symbolic closed-form catalogue
// (general quartic with a₁ or a₃ ≠ 0, deg ≥ 5 irreducible, …) we still owe
// the caller a real number, not an Unimplemented.  We invoke the Aberth
// root isolator, evaluate the standard simple-pole residue formula
//
//   Res(N/D, z₀) = N(z₀) / D'(z₀)
//
// at each upper-half-plane root of the factor, and apply the residue
// theorem
//
//   ∫_{-∞}^{∞} (N/D) dx = 2πi · Σ_{Im(z_k) > 0} Res(N/D, z_k).
//
// For a real-coefficient integrand the sum is purely imaginary, so the
// final value reduces to  −2π · Σ Im(Res).  The result is emitted as a
// DecimalLit at the working MPFR precision (40 decimal digits by default,
// matching the Aberth options).
[[nodiscard]] numeric::detail::CBF rational_to_cbf(const Rational& q, mpfr_prec_t prec) {
    BigFloat bf = BigFloat::from_rational_parts(
        q.numerator().decimal(), q.denominator().decimal(), prec);
    return numeric::detail::CBF::from_real(std::move(bf));
}

[[nodiscard]] numeric::detail::CBF horner_eval_cbf(
    const std::vector<numeric::detail::CBF>& coeffs,
    const numeric::detail::CBF& z,
    mpfr_prec_t prec) {
    if (coeffs.empty()) return numeric::detail::CBF::zero(prec);
    numeric::detail::CBF acc = coeffs.back();
    for (std::size_t k = coeffs.size() - 1U; k-- > 0;) {
        acc = acc * z + coeffs[k];
    }
    return acc;
}

[[nodiscard]] Result<ExprPtr> numeric_residue_contribution(
    const algebra::PolynomialFactor& pf,
    ExprPtr N,
    ExprPtr D,
    const Symbol& var,
    std::size_t deg_N,
    std::size_t deg_D,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    if (pf.multiplicity > 1U) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Residue theorem (numeric): pole multiplicity > 1 "
                       "requires higher-order residue computation; raise via "
                       "partial fractions or supply a single-pole denominator",
        });
    }

    // Higher precision than the Aberth defaults: residue evaluation involves
    // catastrophic cancellation in the imaginary direction (real-coefficient
    // integrand → sum is purely imaginary), so we ask for extra working
    // precision to keep the final to_double() lossless.  Both knobs are
    // configurable on the CASContext (defaults: 80 decimal digits,
    // 500 iterations) — see `set_residue_aberth_precision_digits`.
    //
    // The convergence tolerance is derived from precision_digits rather than
    // hardcoded: we ask for 10^{−(precision − 10)} so the iteration stops
    // when |Δz| sits 10 decimal orders inside working precision.  This keeps
    // the tolerance coherent under any user-driven precision bump without
    // requiring a second knob.
    const unsigned int prec_digits = ctx.residue_aberth_precision_digits();
    const numeric::AberthOptions opts{
        /* precision_digits = */ prec_digits,
        /* max_iterations   = */ ctx.residue_aberth_max_iterations(),
        /* convergence_tolerance = */ std::pow(
            10.0, -static_cast<double>(prec_digits > 10U ? prec_digits - 10U : 1U)),
    };
    const mpfr_prec_t prec = decimal_digits_to_bits(opts.precision_digits);

    auto roots_res = numeric::aberth_isolate_complex_roots(pf.factor, var.name, ctx, opts);
    if (roots_res.is_error()) return fail<ExprPtr>(roots_res.error());

    // Numerator coefficients over BigFloat.
    auto N_coeffs_res = extract_rational_coeffs(N, var, deg_N, ctx);
    if (N_coeffs_res.is_error()) return fail<ExprPtr>(N_coeffs_res.error());
    std::vector<numeric::detail::CBF> N_cbf;
    N_cbf.reserve(N_coeffs_res.value().size());
    for (const Rational& q : N_coeffs_res.value()) {
        N_cbf.push_back(rational_to_cbf(q, prec));
    }

    // Denominator and its formal derivative — D'(z) = Σ k·d_k·z^{k−1}.
    auto D_coeffs_res = extract_rational_coeffs(D, var, deg_D, ctx);
    if (D_coeffs_res.is_error()) return fail<ExprPtr>(D_coeffs_res.error());
    std::vector<numeric::detail::CBF> Dp_cbf;
    if (deg_D >= 1U) {
        Dp_cbf.reserve(deg_D);
        for (std::size_t k = 1U; k <= deg_D; ++k) {
            const Rational scaled = D_coeffs_res.value()[k] *
                                    Rational(BigInt(static_cast<std::int64_t>(k)));
            Dp_cbf.push_back(rational_to_cbf(scaled, prec));
        }
    }

    BigFloat imag_sum(prec);
    bool any_uhp = false;
    for (const auto& r : roots_res.value()) {
        // Upper half-plane: imag strictly positive.
        if (r.imag.is_negative() || r.imag.is_zero()) continue;
        any_uhp = true;
        numeric::detail::CBF z{r.real, r.imag};
        numeric::detail::CBF N_val = horner_eval_cbf(N_cbf, z, prec);
        numeric::detail::CBF Dp_val = horner_eval_cbf(Dp_cbf, z, prec);
        if (Dp_val.is_zero()) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Residue theorem (numeric): D'(z₀) vanishes at a "
                           "root — indicates a non-simple pole or a "
                           "degenerate factorisation",
            });
        }
        numeric::detail::CBF residue = N_val / Dp_val;
        imag_sum = imag_sum + residue.im;
    }

    if (!any_uhp) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Residue theorem (numeric): factor has no upper-half-"
                       "plane roots; integral contribution is zero only if "
                       "the symbolic analysis confirms no real poles",
        });
    }

    // ∫ = 2πi · Σ Res = −2π · Σ Im(Res).
    BigFloat two_pi = BigFloat::pi(prec) + BigFloat::pi(prec);
    BigFloat result = -(two_pi * imag_sum);
    return ok(arena.make<DecimalLit>(
        result.to_string(static_cast<int>(opts.precision_digits))));
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
        if (fdeg == 4U) {
            // Closure for irreducible biquadratic factor: a₄·x⁴ + a₂·x² + a₀.
            // a₁ and a₃ must be zero; the discriminant of the auxiliary
            // quadratic in u = x² must be strictly negative so that the four
            // roots are complex.
            auto bq_coeffs = extract_rational_coeffs(pf.factor, var, 4U, ctx);
            if (bq_coeffs.is_error()) return fail<ExprPtr>(bq_coeffs.error());
            const auto& q = bq_coeffs.value();
            if (!q[1].numerator().is_zero() || !q[3].numerator().is_zero()) {
                // General quartic — no biquadratic closed form, but the Aberth
                // numeric driver isolates the four complex roots and applies
                // the residue theorem at upper-half-plane simple poles.
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
            if (pf.multiplicity > 1U) {
                return fail<ExprPtr>(CASError{
                    .kind = CASErrorKind::Unimplemented,
                    .message = "Residue theorem: biquadratic factor with multiplicity > 1 not yet supported"});
            }
            auto contrib = contribution_from_irreducible_biquadratic(b_norm, c_norm, rational_expr, var, ctx);
            if (contrib.is_error()) return fail<ExprPtr>(contrib.error());
            ExprPtr added = arena.make<Binary>(BinaryOp::Add, total, contrib.value());
            auto simp = simplify_or_fail(added, ctx);
            if (simp.is_error()) return simp;
            total = simp.value();
            continue;
        }
        if (fdeg > 2U) {
            // Aberth numeric driver: irreducible factor of degree ≥ 5 has no
            // closed-form roots in radicals (generic case), so we fall back to
            // numerical isolation + simple-pole residue summation.
            auto numeric_contrib = numeric_residue_contribution(
                pf, N, D, var, deg_N, deg_D, ctx);
            if (numeric_contrib.is_error()) return fail<ExprPtr>(numeric_contrib.error());
            ExprPtr added_n = arena.make<Binary>(BinaryOp::Add, total, numeric_contrib.value());
            auto simp_n = simplify_or_fail(added_n, ctx);
            if (simp_n.is_error()) return simp_n;
            total = simp_n.value();
            continue;
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
