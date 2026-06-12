// F8.0-5.4: AlgebraicNumber / RootOf with isolating bounds.
// Constructs RootOf nodes whose identity is decidable from the triple
//   (polynomial, variable, isolating_bound)
// without recomputing Sturm at evaluation time.
//
// Bound provenance: rigorous Sturm-sequence isolation on the squarefree
// part of `polynomial`, yielding a rational interval [low, high] s.t. the
// polynomial has exactly one real root in that interval.

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/numeric.hpp"

namespace cas::algebra {

[[nodiscard]] Result<ExprPtr> make_rootof_isolated(
    ExprPtr polynomial,
    const Symbol& variable,
    std::size_t root_index,
    symbolic::CASContext& ctx,
    double search_low,
    double search_high,
    double tol) {

    if (!polynomial) {
        return fail<ExprPtr>(CASError{
            .kind    = CASErrorKind::InvalidArgument,
            .message = "make_rootof_isolated: null polynomial",
            .hint    = std::nullopt});
    }
    if (search_low >= search_high) {
        return fail<ExprPtr>(CASError{
            .kind    = CASErrorKind::InvalidArgument,
            .message = "make_rootof_isolated: empty search interval",
            .hint    = std::nullopt});
    }

    // Normalize the polynomial expression so parse_polynomial inside Sturm
    // can consume it. Unsimplified ASTs (e.g. Binary(Sub, Pow(x,2), 2))
    // may bypass the polynomial parser depending on its conventions.
    auto simp_res = ctx.simplify(polynomial);
    ExprPtr normalised = simp_res.is_ok() ? simp_res.value() : polynomial;

    auto bounds_res = numeric::find_polynomial_isolating_intervals(
        normalised, variable.name, ctx, search_low, search_high, tol);
    if (bounds_res.is_error()) return fail<ExprPtr>(bounds_res.error());

    const auto& bounds = bounds_res.value();
    if (root_index >= bounds.size()) {
        return fail<ExprPtr>(CASError{
            .kind    = CASErrorKind::InvalidArgument,
            .message = "make_rootof_isolated: root_index out of range "
                       "(polynomial has " + std::to_string(bounds.size()) +
                       " real roots in the search interval, requested index " +
                       std::to_string(root_index) + ")",
            .hint    = std::nullopt});
    }

    return ok(ctx.arena().make<RootOf>(
        polynomial,
        variable,
        bounds[root_index],
        std::optional<std::size_t>{root_index}));
}

} // namespace cas::algebra
