#include "integrate_engine.hpp"
#include <cas/algebra.hpp>

namespace cas::calculus {

Result<ExprPtr> hermite_reduction(
    ExprPtr /*numerator*/,
    ExprPtr /*denominator*/,
    ExprPtr var,
    symbolic::CASContext& /*ctx*/) {
    if (!expr_is<Symbol>(var)) {
        return fail<ExprPtr>(integrate_detail::make_error(
            CASErrorKind::InvalidArgument,
            "Hermite reduction richiede una variabile simbolica"));
    }
    return fail<ExprPtr>(integrate_detail::make_error(
        CASErrorKind::Unimplemented,
        "Hermite reduction: Bezout step non ancora implementato"));
}

} // namespace cas::calculus
