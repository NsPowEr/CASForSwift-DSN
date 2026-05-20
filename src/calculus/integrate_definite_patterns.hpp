#pragma once

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <vector>

namespace cas::calculus {

// Context passed to every definite-integral pattern matcher.
// `integrand_normalized` is `integrand` after together() + simplify().
struct DefiniteContext {
    ExprPtr integrand;
    ExprPtr integrand_normalized;
    const Symbol& var;
    ExprPtr lower;
    ExprPtr upper;
    symbolic::CASContext& ctx;
};

// A pattern returns:
//   - ok(nullopt)            → pattern not applicable, try next matcher
//   - ok(value)              → pattern matched, value is the integral
//   - fail(error)            → fatal error; aborts the whole dispatch
using DefinitePatternFn = Result<std::optional<ExprPtr>> (*)(const DefiniteContext&);

// Registry. The order matters: first matching pattern wins.  Patterns must
// be conservative — return nullopt unless the structural shape is unambiguous.
[[nodiscard]] const std::vector<DefinitePatternFn>& definite_patterns();

// Individual patterns exported for direct unit testing / registration.
[[nodiscard]] Result<std::optional<ExprPtr>> pattern_gaussian_full_line(const DefiniteContext& dc);
[[nodiscard]] Result<std::optional<ExprPtr>> pattern_rational_full_real_line(const DefiniteContext& dc);
[[nodiscard]] Result<std::optional<ExprPtr>> pattern_bessel_orthogonality(const DefiniteContext& dc);
[[nodiscard]] Result<std::optional<ExprPtr>> pattern_legendre_orthogonality(const DefiniteContext& dc);
[[nodiscard]] Result<std::optional<ExprPtr>> pattern_hermite_h_orthogonality(const DefiniteContext& dc);
[[nodiscard]] Result<std::optional<ExprPtr>> pattern_hermite_he_orthogonality(const DefiniteContext& dc);
[[nodiscard]] Result<std::optional<ExprPtr>> pattern_chebyshev_t_orthogonality(const DefiniteContext& dc);
[[nodiscard]] Result<std::optional<ExprPtr>> pattern_chebyshev_u_orthogonality(const DefiniteContext& dc);

}  // namespace cas::calculus
