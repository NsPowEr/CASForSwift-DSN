// F5.6 — Aberth/Ehrlich simultaneous complex-root isolator for univariate
// polynomials.  Operates over arbitrary-precision floating-point complex
// arithmetic (MPFR BigFloat).  Returns one ComplexRoot per root of the
// squarefree part of f ∈ Q[x], converged to within the requested tolerance.
//
// The isolator is the numeric backbone for the residue-theorem driver at
// denominator degree ≥ 5 (where roots are not expressible in radicals) and
// for any downstream module that needs a complete set of complex roots
// without solving a closed-form formula.
//
// References
//   * O. Aberth, "Iteration methods for finding all zeros of a polynomial
//     simultaneously" (Math. Comp. 27, 1973).
//   * V. Y. Pan, "Solving a polynomial equation: some history and recent
//     progress" (SIAM Rev. 39, 1997).
//   * D. A. Bini, "Numerical computation of polynomial zeros by means of
//     Aberth's method" (Numer. Algorithms 13, 1996) — initialisation and
//     stopping criteria adopted here.

#pragma once

#include "cas/ast.hpp"
#include "cas/bigfloat.hpp"
#include "cas/result.hpp"

#include <string>
#include <vector>

namespace cas::symbolic { class CASContext; }

namespace cas::numeric {

struct ComplexRoot {
    BigFloat real;
    BigFloat imag;
    BigFloat residual;        // |p(z)| at convergence; a posteriori error proxy.
};

struct AberthOptions {
    unsigned int precision_digits{40};   // working precision, in decimal digits.
    unsigned int max_iterations{200};    // Aberth main-loop cap.
    double convergence_tolerance{1e-30}; // ‖Δz‖_∞ termination threshold.
};

// Isolate all complex roots of `poly` interpreted as a univariate polynomial
// in `variable`.  Coefficients must reduce to rationals via the algebra
// layer.  On success returns exactly deg(poly) ComplexRoot entries (counted
// with multiplicity); on failure (non-polynomial input, zero polynomial,
// budget exhausted) returns a diagnostic Unimplemented / InvalidArgument
// error — never a silent partial result.
[[nodiscard]] Result<std::vector<ComplexRoot>> aberth_isolate_complex_roots(
    ExprPtr poly,
    const std::string& variable,
    symbolic::CASContext& ctx,
    const AberthOptions& options = {});

}  // namespace cas::numeric
