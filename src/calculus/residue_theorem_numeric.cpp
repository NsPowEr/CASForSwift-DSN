#include "residue_theorem_internal.hpp"
#include "cas/bigfloat.hpp"
#include "cas/bigint.hpp"
#include "cas/calculus.hpp"
#include "cas/numeric/complex_root_isolator.hpp"
#include "../numeric/complex_bigfloat_internal.hpp"
#include <cmath>
#include <vector>

namespace cas::calculus {

namespace {
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
} // namespace

Result<ExprPtr> numeric_residue_contribution(
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

    const unsigned int prec_digits = ctx.residue_aberth_precision_digits();
    const numeric::AberthOptions opts{
        prec_digits,
        ctx.residue_aberth_max_iterations(),
        std::pow(10.0, -static_cast<double>(prec_digits > 10U ? prec_digits - 10U : 1U)),
    };
    const mpfr_prec_t prec = decimal_digits_to_bits(opts.precision_digits);

    auto roots_res = numeric::aberth_isolate_complex_roots(pf.factor, var.name, ctx, opts);
    if (roots_res.is_error()) return fail<ExprPtr>(roots_res.error());

    auto N_coeffs_res = extract_rational_coeffs(N, var, deg_N, ctx);
    if (N_coeffs_res.is_error()) return fail<ExprPtr>(N_coeffs_res.error());
    std::vector<numeric::detail::CBF> N_cbf;
    N_cbf.reserve(N_coeffs_res.value().size());
    for (const Rational& q : N_coeffs_res.value()) {
        N_cbf.push_back(rational_to_cbf(q, prec));
    }

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

    BigFloat two_pi = BigFloat::pi(prec) + BigFloat::pi(prec);
    BigFloat result = -(two_pi * imag_sum);
    return ok(arena.make<DecimalLit>(
        result.to_string(static_cast<int>(opts.precision_digits))));
}

}  // namespace cas::calculus
