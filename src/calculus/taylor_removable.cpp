// taylor_removable.cpp — Taylor expansion across a removable singularity (A37).
//
// Problem. The generic Taylor path (limit_series.cpp) builds coefficients as
//   c_k = f^(k)(center) / k!
// by differentiating and SUBSTITUTING the center. For a quotient whose
// denominator vanishes there — (e^x-1)/x, sin(x)/x, log(1+x)/x, x/(e^x-1) at
// x = 0 — every such substitution is a 0/0 form, so the arithmetic layer
// reports "ComplexRational: division by zero" and the failure surfaces as if
// the expansion did not exist. It does exist: the singularity is removable and
// the series is perfectly ordinary.
//
// Fix. Do not evaluate the quotient at the center at all. Expand numerator and
// denominator separately (both analytic there) and divide the two truncated
// series — which is precisely what laurent_series_general already implements
// (Taylor coefficients of N and D, leading-exponent detection, geometric-series
// inversion of D). The quotient's Laurent expansion has
//   leading_order = val(N) - val(D),
// so the singularity is removable exactly when leading_order >= 0, and in that
// case the Laurent coefficients ARE the Taylor coefficients.
//
// Depth of the expansion. To get `order` correct coefficients of N/D when D
// vanishes to order v, both series must be carried to order + v. v is not known
// before expanding, so the extra depth doubles (1, 2, 4, ...) until D's leading
// term is visible, bounded by ctx.max_removable_singularity_order(). Past that
// bound the result is Unimplemented — never a silently truncated series.
//
// A genuine pole (leading_order < 0, e.g. sin(x)/x^2) is NOT an error either:
// it is reported as Unimplemented with the pole order, pointing the caller at
// laurent_series, which is the function that can represent it.

#include "cas/calculus.hpp"
#include "cas/error_helpers.hpp"
#include "calculus_internal.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus {
namespace {

[[nodiscard]] bool has_vanishing_denominator_error(const CASError& err) {
    // laurent_series_general reports this when D's expansion is zero to the
    // depth tried — i.e. the depth was too shallow to see val(D).
    return err.message.find("denominator is identically zero") != std::string::npos;
}

}  // namespace

[[nodiscard]] Result<TaylorExpansion> taylor_series_removable(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr center,
    unsigned int order,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    // laurent_series_general splits on a top-level Div node only, but simplify
    // has long since rewritten f/g into the canonical Product·Pow(-1) form, so
    // hand it an explicit quotient recovered by the general splitter (which
    // understands Product, negative Pow and Unary shapes). Without this the
    // whole expression is taken as the numerator over 1 and its coefficients
    // are again evaluated at the singular point.
    auto quotient = extract_quotient_view(expr, arena);
    ExprPtr target = quotient.has_value()
        ? arena.make<Binary>(BinaryOp::Div, quotient->numerator, quotient->denominator)
        : expr;

    const unsigned int max_extra = ctx.max_removable_singularity_order();
    std::optional<LaurentExpansion> expansion;

    // Two independent depths are searched here, and neither may be guessed:
    //   `extra` — how far past the requested order N and D are carried, so
    //             that the denominator's leading term becomes visible;
    //   `pad`   — how many positive powers are requested, because the series
    //             product truncates at the shallower of its two operands and
    //             can therefore return one coefficient short of the request.
    // Rather than pad by a fixed amount, ask deeper and CHECK that the top
    // requested exponent is actually present before accepting the expansion.
    const unsigned int max_pad = ctx.max_removable_singularity_order();
    for (unsigned int pad = 1U; pad <= max_pad && !expansion.has_value(); ++pad) {
        for (unsigned int extra = 1U; extra <= max_extra; extra *= 2U) {
            auto attempt = laurent_series_general(target, var, center,
                                                  order + pad, extra, ctx);
            if (attempt.is_ok()) {
                const auto& lau = attempt.value();
                const long long top = static_cast<long long>(lau.leading_order) +
                                      static_cast<long long>(lau.coefficients.size()) - 1;
                // A pole is a final answer; otherwise require full coverage of
                // the requested order before accepting.
                if (lau.leading_order < 0 || top >= static_cast<long long>(order)) {
                    expansion = lau;
                }
                break;
            }
            if (!has_vanishing_denominator_error(attempt.error())) {
                return fail<TaylorExpansion>(attempt.error());
            }
            // Denominator still invisible at this depth: go deeper.
        }
    }

    if (!expansion.has_value()) {
        return make_unimplemented<TaylorExpansion>(
            "calculus", "taylor_series_removable",
            "denominator vanishes beyond depth " + std::to_string(max_extra),
            "TAYLOR_REMOVABLE_DEPTH",
            "Raise CASContext::set_max_removable_singularity_order if the "
            "denominator legitimately vanishes to a higher order",
            "A37");
    }

    const LaurentExpansion& lau = expansion.value();
    if (lau.leading_order < 0) {
        return make_unimplemented<TaylorExpansion>(
            "calculus", "taylor_series_removable",
            "pole of order " + std::to_string(-lau.leading_order),
            "TAYLOR_POLE_NOT_REMOVABLE",
            "The singularity is not removable — no Taylor expansion exists; "
            "use laurent_series for the negative-power part",
            "A37");
    }

    // leading_order >= 0: rebuild Σ_{k} c_k · (var - center)^k as a plain
    // polynomial in (var - center), the same shape the Taylor path returns.
    ExprPtr x = arena.make<Symbol>(var);
    ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, x, center);

    std::vector<ExprPtr> terms;
    terms.reserve(lau.coefficients.size());
    for (std::size_t i = 0; i < lau.coefficients.size(); ++i) {
        const long long exponent = static_cast<long long>(lau.leading_order) +
                                   static_cast<long long>(i);
        // The expansion was requested deeper than `order` to guarantee the top
        // coefficient exists; the caller asked for `order`.
        if (exponent > static_cast<long long>(order)) break;
        ExprPtr coeff = lau.coefficients[i];
        if (exponent == 0) {
            terms.push_back(coeff);
            continue;
        }
        ExprPtr power = exponent == 1
            ? delta
            : arena.make<Binary>(BinaryOp::Pow, delta,
                                 arena.make<IntegerLit>(BigInt(exponent)));
        terms.push_back(arena.make<Product>(std::vector<ExprPtr>{coeff, power}));
    }

    ExprPtr polynomial = terms.empty()
        ? arena.make<IntegerLit>(BigInt(0))
        : (terms.size() == 1U ? terms.front()
                              : arena.make<Sum>(std::move(terms)));

    auto simplified = ctx.simplify(polynomial);
    if (simplified.is_error()) return fail<TaylorExpansion>(simplified.error());

    ExprPtr remainder = arena.make<FuncCall>("O", std::vector<ExprPtr>{
        arena.make<Binary>(BinaryOp::Pow, delta,
                           arena.make<IntegerLit>(BigInt(static_cast<long long>(order) + 1)))});

    auto materialized_poly = symbolic::materialize_expr(simplified.value(), arena);
    if (materialized_poly.is_error()) return fail<TaylorExpansion>(materialized_poly.error());
    auto materialized_rem = symbolic::materialize_expr(remainder, arena);
    if (materialized_rem.is_error()) return fail<TaylorExpansion>(materialized_rem.error());

    return ok(TaylorExpansion{
        .polynomial = materialized_poly.value(),
        .remainder = materialized_rem.value(),
        .computed_order = order,
    });
}

}  // namespace cas::calculus
