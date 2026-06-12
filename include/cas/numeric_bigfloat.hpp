#pragma once
// F8.0-5.3: MPFR-precision numeric routines that wrap the Rational-exact
// Sturm pipeline with BigFloat Newton refinement. Separated from the main
// numeric.hpp so that translation units that don't need MPFR don't pay
// the mpfr.h include cost.

#include "cas/ast.hpp"
#include "cas/bigfloat.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <string>
#include <vector>

namespace cas {
namespace numeric {

// MPFR-precision Newton polish of Sturm-isolated real roots.
// Same input contract as find_polynomial_roots_sturm but returns BigFloat
// values polished at the requested precision (in bits). Internally uses
// the exact Sturm sequence on Rational coefficients (sign-variation count
// is rigorous), then refines each isolated interval midpoint via Newton
// iteration in BigFloat arithmetic. Returned values lie inside the
// corresponding rational isolating interval at the requested precision.
//
// Precision argument is the MPFR mantissa precision in bits (≥ 64
// recommended; 128 default). Higher precision narrows the rounding
// envelope of the polished value, NOT of the Sturm isolation itself
// (Sturm is exact by construction).
[[nodiscard]] Result<std::vector<BigFloat>> find_polynomial_roots_sturm_bigfloat(
    ExprPtr expr,
    const std::string& variable,
    symbolic::CASContext& ctx,
    double low,
    double high,
    double tol,
    mpfr_prec_t precision_bits = BigFloat::DEFAULT_PREC);

} // namespace numeric
} // namespace cas
