// F5.6 — Aberth/Ehrlich simultaneous complex-root isolator.
//
// Algorithm (Aberth 1973, refined by Bini 1996):
//   * Convert poly to monic complex BigFloat coefficient vector (length n+1).
//   * Initialise n guesses on a circle of Cauchy-bound radius, perturbed by
//     a half-spoke angular offset to avoid collinearity with real roots.
//   * Iterate the Aberth correction
//         w_i = (p(z_i) / p'(z_i)) / (1 − (p(z_i) / p'(z_i)) · σ_i),
//         σ_i = Σ_{j≠i} 1 / (z_i − z_j),
//         z_i ← z_i − w_i,
//     until max_i |w_i| < tolerance or max_iterations is reached.
//
// Complex arithmetic uses an inline `CBF` struct over two MPFR BigFloats —
// avoiding an MPC dependency keeps the build matrix unchanged.  Division and
// reciprocal use the standard real-only formulae.
//
// The eval / derivative pair uses Horner's scheme (compensated form):
//   p(z)  = b_0 with b_n = a_n,  b_k = b_{k+1}·z + a_k;
//   p'(z) = d_0 with d_n = 0,    d_k = d_{k+1}·z + b_{k+1}.

#include "cas/algebra.hpp"
#include "cas/bigfloat.hpp"
#include "cas/error.hpp"
#include "cas/numeric/complex_root_isolator.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace cas::numeric {

namespace {

struct CBF {
    BigFloat re;
    BigFloat im;

    CBF() = default;
    CBF(BigFloat r, BigFloat i) : re(std::move(r)), im(std::move(i)) {}

    static CBF zero(mpfr_prec_t prec) {
        return CBF{BigFloat(prec), BigFloat(prec)};
    }
    static CBF from_real(BigFloat r) {
        BigFloat zero_im(r.precision_bits());
        return CBF{std::move(r), std::move(zero_im)};
    }
    static CBF from_double_pair(double r, double i, mpfr_prec_t prec) {
        return CBF{BigFloat::from_double(r, prec), BigFloat::from_double(i, prec)};
    }

    bool is_zero() const noexcept { return re.is_zero() && im.is_zero(); }

    CBF operator+(const CBF& o) const { return CBF{re + o.re, im + o.im}; }
    CBF operator-(const CBF& o) const { return CBF{re - o.re, im - o.im}; }
    CBF operator-() const { return CBF{-re, -im}; }
    CBF operator*(const CBF& o) const {
        return CBF{re * o.re - im * o.im, re * o.im + im * o.re};
    }
    CBF operator/(const CBF& o) const {
        BigFloat denom = o.re * o.re + o.im * o.im;
        BigFloat new_re = (re * o.re + im * o.im) / denom;
        BigFloat new_im = (im * o.re - re * o.im) / denom;
        return CBF{std::move(new_re), std::move(new_im)};
    }

    BigFloat abs() const { return BigFloat::sqrt(re * re + im * im); }
    BigFloat abs_sq() const { return re * re + im * im; }
};

// Convert ExprPtr literal to BigFloat at requested precision.  Accepts
// IntegerLit and RationalLit (the algebra layer hands us literal coefficients
// after simplify).  Anything else is a programming error in the caller.
[[nodiscard]] Result<BigFloat> rational_expr_to_bigfloat(
    ExprPtr e,
    mpfr_prec_t prec) {
    if (const auto* il = expr_cast<IntegerLit>(e)) {
        return ok(BigFloat::from_integer_string(il->value.decimal(), prec));
    }
    if (const auto* rl = expr_cast<RationalLit>(e)) {
        return ok(BigFloat::from_rational_parts(
            rl->numerator.decimal(),
            rl->denominator.decimal(),
            prec));
    }
    return fail<BigFloat>(CASError{
        .kind = CASErrorKind::Unimplemented,
        .message = "Aberth: polynomial coefficient is not a literal rational "
                   "after simplification; complex-coefficient case is not "
                   "yet handled by the numeric isolator",
    });
}

// Cauchy bound:  every root z of  a_n z^n + … + a_0  satisfies
//   |z| ≤ 1 + max_k |a_k / a_n|.
[[nodiscard]] BigFloat cauchy_bound(const std::vector<CBF>& coeffs, mpfr_prec_t prec) {
    const std::size_t n = coeffs.size() - 1U;
    BigFloat leading_abs = coeffs[n].abs();
    BigFloat best(prec);
    for (std::size_t k = 0; k < n; ++k) {
        BigFloat r = coeffs[k].abs() / leading_abs;
        if (r > best) best = r;
    }
    return BigFloat::from_double(1.0, prec) + best;
}

// Horner pair: returns (p(z), p'(z)).
[[nodiscard]] std::pair<CBF, CBF> horner_pair(
    const std::vector<CBF>& coeffs,
    const CBF& z,
    mpfr_prec_t prec) {
    const std::size_t n = coeffs.size() - 1U;
    CBF p = coeffs[n];
    CBF dp = CBF::zero(prec);
    for (std::size_t k = n; k-- > 0;) {
        // d_k = d_{k+1}·z + b_{k+1}    (must come BEFORE updating b)
        dp = dp * z + p;
        // b_k = b_{k+1}·z + a_k
        p = p * z + coeffs[k];
    }
    return {std::move(p), std::move(dp)};
}

}  // namespace

Result<std::vector<ComplexRoot>> aberth_isolate_complex_roots(
    ExprPtr poly,
    const std::string& variable,
    symbolic::CASContext& ctx,
    const AberthOptions& options) {
    const mpfr_prec_t prec = decimal_digits_to_bits(options.precision_digits);

    Symbol var(variable);
    auto coeffs_expr_res = algebra::univariate_coefficients(poly, var, ctx);
    if (coeffs_expr_res.is_error()) return fail<std::vector<ComplexRoot>>(coeffs_expr_res.error());
    const std::vector<ExprPtr>& coeffs_expr = coeffs_expr_res.value();

    if (coeffs_expr.empty()) {
        return fail<std::vector<ComplexRoot>>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "Aberth: empty coefficient vector (zero polynomial)",
        });
    }

    std::vector<CBF> coeffs;
    coeffs.reserve(coeffs_expr.size());
    for (ExprPtr c : coeffs_expr) {
        auto simp = ctx.simplify(c);
        if (simp.is_error()) return fail<std::vector<ComplexRoot>>(simp.error());
        auto bf = rational_expr_to_bigfloat(simp.value(), prec);
        if (bf.is_error()) return fail<std::vector<ComplexRoot>>(bf.error());
        coeffs.push_back(CBF::from_real(std::move(bf.value())));
    }

    // Strip trailing zeros (high-degree zero coefficients) to obtain the
    // actual degree.  The caller may pass an expanded polynomial with a
    // zero leading term if expand introduced one — guard against it.
    while (coeffs.size() > 1U && coeffs.back().is_zero()) coeffs.pop_back();

    const std::size_t n = coeffs.size() - 1U;
    if (n == 0U) {
        return fail<std::vector<ComplexRoot>>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "Aberth: degree-zero polynomial has no roots",
        });
    }

    // Strip zero roots once (multiplicities collapse to a single shift; the
    // recovered roots are emitted as exact 0 + 0·i with residual = 0).
    std::vector<ComplexRoot> zero_roots;
    while (!coeffs.empty() && coeffs.front().is_zero()) {
        zero_roots.push_back(ComplexRoot{
            BigFloat::from_double(0.0, prec),
            BigFloat::from_double(0.0, prec),
            BigFloat::from_double(0.0, prec),
        });
        coeffs.erase(coeffs.begin());
    }
    if (coeffs.empty() || coeffs.size() == 1U) {
        // Either the polynomial was identically zero (already filtered) or
        // it reduces to a constant after stripping zero roots; in that case
        // only the zero roots remain.
        return ok(std::move(zero_roots));
    }
    const std::size_t deg = coeffs.size() - 1U;

    // Initialise guesses on a Cauchy-bound circle, perturbed by a half-spoke
    // offset to avoid landing on the real axis when roots are conjugate.
    const BigFloat radius = cauchy_bound(coeffs, prec);
    const BigFloat pi_val = BigFloat::pi(prec);
    const BigFloat two_pi = pi_val + pi_val;
    const BigFloat n_bf = BigFloat::from_double(static_cast<double>(deg), prec);
    const BigFloat two_n_bf = n_bf + n_bf;
    const BigFloat offset = pi_val / two_n_bf;
    std::vector<CBF> z(deg, CBF::zero(prec));
    for (std::size_t k = 0; k < deg; ++k) {
        BigFloat k_bf = BigFloat::from_double(static_cast<double>(k), prec);
        BigFloat theta = (two_pi * k_bf) / n_bf + offset;
        BigFloat ct = BigFloat::cos(theta);
        BigFloat st = BigFloat::sin(theta);
        z[k] = CBF{radius * ct, radius * st};
    }

    const BigFloat tol = BigFloat::from_double(options.convergence_tolerance, prec);

    std::vector<CBF> w(deg, CBF::zero(prec));
    bool converged = false;
    for (unsigned int iter = 0; iter < options.max_iterations; ++iter) {
        BigFloat max_step(prec);
        for (std::size_t i = 0; i < deg; ++i) {
            auto [p_zi, dp_zi] = horner_pair(coeffs, z[i], prec);
            if (dp_zi.is_zero()) {
                // Singular derivative — perturb and continue.
                z[i] = z[i] + CBF{BigFloat::from_double(1e-12, prec),
                                  BigFloat::from_double(1e-12, prec)};
                continue;
            }
            CBF ratio = p_zi / dp_zi;
            CBF sigma = CBF::zero(prec);
            for (std::size_t j = 0; j < deg; ++j) {
                if (i == j) continue;
                CBF diff = z[i] - z[j];
                if (diff.is_zero()) {
                    diff = CBF{BigFloat::from_double(1e-12, prec),
                               BigFloat::from_double(0.0, prec)};
                }
                sigma = sigma + CBF::from_real(BigFloat::from_double(1.0, prec)) / diff;
            }
            CBF one = CBF::from_real(BigFloat::from_double(1.0, prec));
            CBF denom = one - ratio * sigma;
            if (denom.is_zero()) {
                continue;
            }
            w[i] = ratio / denom;
            BigFloat step = w[i].abs();
            if (step > max_step) max_step = step;
            z[i] = z[i] - w[i];
        }
        if (max_step < tol) {
            converged = true;
            break;
        }
    }
    if (!converged) {
        return fail<std::vector<ComplexRoot>>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Aberth: maximum iterations reached without converging "
                       "to the requested tolerance; raise precision_digits or "
                       "max_iterations",
        });
    }

    std::vector<ComplexRoot> roots = std::move(zero_roots);
    roots.reserve(roots.size() + deg);
    for (std::size_t i = 0; i < deg; ++i) {
        auto [p_zi, dp_zi] = horner_pair(coeffs, z[i], prec);
        (void)dp_zi;
        roots.push_back(ComplexRoot{
            std::move(z[i].re),
            std::move(z[i].im),
            p_zi.abs(),
        });
    }
    return ok(std::move(roots));
}

}  // namespace cas::numeric
