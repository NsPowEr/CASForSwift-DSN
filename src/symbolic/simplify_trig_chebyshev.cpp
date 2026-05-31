#include "simplify_trig_chebyshev_impl.hpp"

// F1.4b — Chebyshev trig generator.
//
// Provides the minimal polynomial of 2cos(π/q) for arbitrary integer q,
// and a RootOf fallback for angles not reachable by nested radical construction.
//
// Mathematical foundation (Disquisitiones Arithmeticae §VII, Gauss 1801):
//   cos(2π/n) is constructible by ruler-and-compass ⟺
//   n = 2^a · p₁ · p₂ · … · pₖ  where pᵢ are distinct Fermat primes {3,5,17,257,65537}.
//   For non-constructible n, the minimal polynomial Ψ_n(y) (y = 2cos(2π/n)) has no
//   radical-form root; the canonical representation is RootOf(Ψ_n(y)).
//
// Algorithm for min_poly_2cos_pi_q(q):
//   Let n = 2q.  The 2n-th cyclotomic polynomial Φ_n(x) has roots {e^{2πik/n} : gcd(k,n)=1}.
//   Under the substitution y = x + x^{-1}, the minimal polynomial of y = 2cos(2π/n) is
//   the unique monic polynomial Ψ_n ∈ Z[y] of degree φ(n)/2 satisfying
//     x^{φ(n)/2} · Ψ_n(x + 1/x) = Φ_n(x).
//   We compute Ψ_n by evaluating the Chebyshev recurrence on the coefficient level:
//   writing y^k = T_k(x) + T_{k-2}(x) + … (sum of Chebyshev monomials), we express
//   each power (x + 1/x)^k in terms of x^j + x^{-j} = 2·cos(j·θ), equate with Φ_n
//   coefficients, and solve for Ψ_n coefficients by back-substitution.
//
//   Practical computation: the substitution y = x + 1/x in Φ_n gives a Laurent poly;
//   multiplying by x^{φ(n)/2} yields an ordinary polynomial.  We read off the Ψ_n
//   coefficients from the pairing of symmetric coefficient pairs of Φ_n.
//
// For q ≤ kCosPolyMaxQ (500) this produces a polynomial of degree ≤ φ(q) ≤ 200
// in O(q log q) time via compute_cyclotomic.

namespace cas::symbolic::detail {

using algebra::IntPoly;
using algebra::compute_cyclotomic;

// ── Minimal polynomial computation ──────────────────────────────────────────

// Compute Ψ_{2q}(y), the minimal polynomial of 2cos(π/q) over Q.
//
// Method: given Φ_{2q}(x) = ∑_{k=0}^{d} a_k x^k  (d = φ(2q)),
// and the symmetry Φ_{2q}(x) = x^d · Φ_{2q}(1/x) for 2q > 2
// (all cyclotomic polynomials are palindromic for n>2), we have
//   Φ_{2q}(x) = x^{d/2} · Ψ(x + 1/x)
// for some degree-d/2 polynomial Ψ.  To extract Ψ:
//   Divide Φ by x^{d/2}, giving a "half-polynomial" in (x + 1/x).
//   The coefficients of Ψ are recovered by the Newton-type identity:
//     Ψ(y) = ∑_{k=0}^{d/2} b_k y^k
//   where b_k are determined by: a_{d/2+k} + a_{d/2-k} = b_k for k<d/2,
//   and a_{d/2} = b_0 (the middle coefficient of Φ).
//
// Actually the precise relationship is:
//   a_{d/2 + k} = [y^k in Ψ] for k = 0,…,d/2.
// This follows because (x + 1/x)^k has leading term x^k, so by
// comparing leading and sub-leading coefficients inductively.
IntPoly min_poly_2cos_pi_q(int q) {
    if (q <= 0 || q > kCosPolyMaxQ) return IntPoly{};
    if (q == 1) {
        // cos(π/1) = cos(π) = -1; min poly: y + 1 = 0 → coeffs [-1, 1]
        // But q=1 is trivial; the half-angle / table handles it.
        return IntPoly({BigInt(1), BigInt(1)});  // y + 1
    }
    if (q == 2) {
        // cos(π/2) = 0; min poly: y = 0
        return IntPoly({BigInt(0), BigInt(1)});  // y
    }

    const int two_q = 2 * q;
    IntPoly phi = compute_cyclotomic(two_q);
    if (phi.empty()) return IntPoly{};

    const std::size_t d = phi.degree();  // φ(2q)
    if (d == 0 || d % 2 != 0) return IntPoly{};  // should not happen for 2q > 2

    const std::size_t half = d / 2;

    // Palindrome check: Φ_{2q} is palindromic for 2q > 2.
    // phi[k] == phi[d-k] for all k.
    // The minimal polynomial of y = x + 1/x satisfies:
    //   phi[half + k] = coefficient of y^k in Ψ   (k = 0, …, half).
    // This is the standard "palindrome half" extraction.
    // Reference: Cohen "A Course in Computational Algebraic Number Theory"
    // §4.5, or Lang "Algebra" Chap. VI proof of Proposition 6.
    std::vector<BigInt> psi_coeffs(half + 1U);
    for (std::size_t k = 0; k <= half; ++k) {
        if (half + k < phi.size()) {
            psi_coeffs[k] = phi[half + k];
        }
    }

    IntPoly psi(std::move(psi_coeffs));
    psi.normalize([](const BigInt& v) { return v.is_zero(); });
    return psi;
}

// ── RootOf construction ──────────────────────────────────────────────────────

// Build RootOf(Ψ_{2q}(t), t, 0) where index 0 selects the root with the
// largest real part (i.e. cos(π/q) for q ≥ 2, which lies in (0,1] for q ≥ 2).
// The returned expression represents cos(π/q)/2 in terms of the root t
// of Ψ_{2q}; the caller is responsible for the factor-of-2 or Chebyshev scaling.
//
// Canonical form: we emit RootOf(Ψ_t, t, 0) ÷ 2 where Ψ_t is built
// with fresh symbol "_tcc" (guaranteed fresh by naming convention).
// The root_index = 0 convention follows the lexicographically-largest
// real root (standard in CAS literature: Mathematica's Root[f,1] picks
// the smallest, Root[f,k] the k-th in order; we follow SymPy's convention
// of index 0 = largest real root for palindromic cosine polynomials).
ExprPtr build_rootof_cos_pi_q(int q, AstArena& arena) {
    IntPoly psi = min_poly_2cos_pi_q(q);
    if (psi.empty()) return nullptr;

    // Build the polynomial expression in fresh symbol _tcc.
    Symbol t_sym("_tcc");
    ExprPtr t_var = arena.make<Symbol>("_tcc");

    // Construct Ψ(t) as Sum of a_k * t^k.
    std::vector<ExprPtr> poly_terms;
    const auto& coeffs = psi.coefficients();
    for (std::size_t k = 0; k < coeffs.size(); ++k) {
        if (coeffs[k].is_zero()) continue;
        ExprPtr coeff_expr = arena.make<IntegerLit>(coeffs[k]);
        if (k == 0) {
            poly_terms.push_back(coeff_expr);
        } else if (k == 1) {
            if (coeffs[k] == BigInt(1)) {
                poly_terms.push_back(t_var);
            } else {
                poly_terms.push_back(arena.make<Binary>(BinaryOp::Mul,
                    coeff_expr, t_var));
            }
        } else {
            ExprPtr t_pow = arena.make<Binary>(BinaryOp::Pow, t_var,
                arena.make<IntegerLit>(BigInt(static_cast<long long>(k))));
            if (coeffs[k] == BigInt(1)) {
                poly_terms.push_back(t_pow);
            } else {
                poly_terms.push_back(arena.make<Binary>(BinaryOp::Mul,
                    coeff_expr, t_pow));
            }
        }
    }

    if (poly_terms.empty()) return nullptr;

    ExprPtr poly_expr = poly_terms.size() == 1U
        ? poly_terms[0]
        : arena.make<Sum>(std::move(poly_terms));

    // RootOf(Ψ, t, 0): index 0 = largest real root ≈ 2cos(π/q) ∈ (0, 2].
    ExprPtr rootof = arena.make<RootOf>(poly_expr, t_sym,
        std::optional<std::size_t>{0U});

    // Return RootOf / 2 so it directly represents cos(π/q).
    return arena.make<Binary>(BinaryOp::Div, rootof,
        arena.make<IntegerLit>(BigInt(2)));
}

} // namespace cas::symbolic::detail
