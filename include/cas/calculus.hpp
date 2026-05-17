#pragma once

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <vector>

namespace cas::calculus {

struct TaylorExpansion {
    ExprPtr polynomial;
    ExprPtr remainder;
    unsigned int computed_order;
};

// A Laurent expansion around `center` of a function f(x):
//   f(x) = Σ_{k = leading_order .. positive_order} coefficients[k - leading_order] * (x - center)^k
//        + remainder        (= O((x - center)^(positive_order + 1)) )
//
// `leading_order` is negative when the expansion has a true pole at `center`,
// zero or positive when the function is analytic there.
// `coefficients` is stored in ascending power order so that
//   coefficients[0]                          = c_{leading_order}
//   coefficients[-leading_order]             = c_0
//   coefficients[positive_order-leading_order] = c_{positive_order}.
struct LaurentExpansion {
    ExprPtr center;
    int leading_order;
    std::vector<ExprPtr> coefficients;
    unsigned int positive_order;
    ExprPtr remainder;
};

[[nodiscard]] Result<ExprPtr> diff(
    ExprPtr expr,
    const Symbol& var,
    unsigned int order,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> partial_diff(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> implicit_diff(
    ExprPtr relation,
    const Symbol& dependent,
    const Symbol& independent,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<std::vector<ExprPtr>> gradient(
    ExprPtr expr,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> jacobian(
    const std::vector<ExprPtr>& funcs,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> hessian(
    ExprPtr expr,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> integrate(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> definite_integral(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> limit(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    LimitDirection dir,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> sum(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> residue(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr pole,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<TaylorExpansion> taylor_series(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    unsigned int order,
    symbolic::CASContext& ctx);

// Compute the Laurent series expansion of `expr` around `center` with
// `positive_order` positive-power terms beyond the constant.  Detects the
// pole order automatically (returns leading_order = 0 when `expr` is
// analytic at the center).  Currently restricted to expressions that can
// be written as a rational function of `var`.
[[nodiscard]] Result<LaurentExpansion> laurent_series(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr center,
    unsigned int positive_order,
    symbolic::CASContext& ctx);

struct PadeApproximant {
    ExprPtr center;
    ExprPtr numerator;       // polynomial in (var − center) of degree ≤ numerator_order
    ExprPtr denominator;     // polynomial in (var − center) of degree ≤ denominator_order
    unsigned int numerator_order;
    unsigned int denominator_order;
};

// Pade [m/n] approximant of `expr` at `center`.  Builds P, Q with
// deg(P) ≤ m, deg(Q) ≤ n, Q(0) = 1, so that  P(x)/Q(x) ≡ f(x) mod (x − c)^{m+n+1}.
// Internally Taylor-expands `expr` to order m+n at the centre and solves
// the linear Toeplitz system for Q's coefficients (Gauss elimination over
// Rational), then back-substitutes to recover P.  Returns Unimplemented
// if the Toeplitz system is singular, which signals a degenerate Pade
// table entry (defect rank) rather than a silent wrong answer.
[[nodiscard]] Result<PadeApproximant> pade_approximant(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr center,
    unsigned int numerator_order,
    unsigned int denominator_order,
    symbolic::CASContext& ctx);

struct Asymptote {
    enum class Type { Vertical, Horizontal, Slant };
    Type type;
    ExprPtr expression; // x = c for Vertical, y = c for Horizontal, y = mx + q for Slant
};

[[nodiscard]] Result<std::vector<Asymptote>> find_asymptotes(
    ExprPtr f,
    const Symbol& x,
    symbolic::CASContext& ctx);

}  // namespace cas::calculus
