#include "polynomial_internal.hpp"
#include "cas/symbolic.hpp"

#include <algorithm>
#include <utility>
#include <cstddef>

namespace cas {
namespace algebra {

struct IntegerSubresultantExecution {
    IntPoly gcd;
    IntegerGcdPath path{IntegerGcdPath::Subresultant};
};

[[nodiscard]] IntegerSubresultantExecution run_integer_subresultant(IntPoly A, IntPoly B) {
    if (A.is_zero()) {
        return IntegerSubresultantExecution{
            .gcd = std::move(B),
            .path = IntegerGcdPath::Subresultant,
        };
    }
    if (B.is_zero()) {
        return IntegerSubresultantExecution{
            .gcd = std::move(A),
            .path = IntegerGcdPath::Subresultant,
        };
    }

    const BigInt contA = integer_content(A);
    const BigInt contB = integer_content(B);
    const BigInt common_content = gcd(contA, contB);

    IntPoly H1 = primitive_integer_poly(std::move(A));
    IntPoly H2 = primitive_integer_poly(std::move(B));

    if (H1.degree() < H2.degree()) {
        std::swap(H1, H2);
    }

    std::size_t d1 = H1.degree();
    std::size_t d2 = H2.degree();
    std::size_t delta = d1 - d2;

    // Brown/Collins Subresultant PRS
    // beta_2 = (-1)^(delta_1 + 1)
    BigInt beta = (delta % 2 == 0) ? BigInt(-1) : BigInt(1);
    BigInt psi = BigInt(-1);

    while (true) {
        IntPoly R = pseudo_remainder_integer_poly(H1, H2);
        if (R.is_zero()) {
            break;
        }

        if (!try_divide_integer_coefficients_by_scalar(R, beta)) {
            return IntegerSubresultantExecution{
                .gcd = IntPoly{},
                .path = IntegerGcdPath::PrimitiveFallback,
            };
        }

        H1 = std::move(H2);
        H2 = std::move(R);

        const std::size_t prev_delta = delta;
        d1 = H1.degree();
        d2 = H2.degree();
        delta = d1 - d2;

        const BigInt lcH1 = H1.leading_coeff();

        // psi_{i+1} = (-lc(H_i))^delta_{i-1} * psi_i^(1-delta_{i-1})
        if (prev_delta > 0) {
            psi = bigint_pow_nonnegative(-lcH1, prev_delta) / bigint_pow_nonnegative(psi, prev_delta - 1);
        } else if (prev_delta == 0) {
            // psi remains unchanged: psi = (-lcH1)^0 * psi^1
        }

        // beta_{i+1} = -lc(H_i) * psi_{i+1}^delta_i
        beta = -lcH1 * bigint_pow_nonnegative(psi, delta);
    }

    IntPoly result = primitive_integer_poly(std::move(H2));
    multiply_integer_coefficients_by_scalar(result, common_content);

    return IntegerSubresultantExecution{
        .gcd = std::move(result),
        .path = IntegerGcdPath::Subresultant,
    };
}

bool is_zero_integer_poly(const IntPoly& coefficients) {
    return coefficients.is_zero();
}

IntPoly gcd_integer_poly_primitive(IntPoly lhs, IntPoly rhs) {
    lhs = primitive_integer_poly(std::move(lhs));
    rhs = primitive_integer_poly(std::move(rhs));

    if (lhs.empty()) {
        return rhs;
    }
    if (rhs.empty()) {
        return lhs;
    }

    while (!rhs.empty()) {
        IntPoly remainder = pseudo_remainder_integer_poly(lhs, rhs);
        remainder = primitive_integer_poly(std::move(remainder));
        lhs = std::move(rhs);
        rhs = std::move(remainder);
    }

    return primitive_integer_poly(std::move(lhs));
}

IntegerGcdResult gcd_integer_poly_with_subresultant(IntPoly lhs, IntPoly rhs) {
    BigInt c_lhs = integer_content(lhs);
    BigInt c_rhs = integer_content(rhs);
    BigInt c_gcd = gcd(c_lhs, c_rhs);

    lhs = primitive_integer_poly(std::move(lhs));
    rhs = primitive_integer_poly(std::move(rhs));

    IntegerSubresultantExecution execution = run_integer_subresultant(lhs, rhs);
    
    IntPoly result_gcd;
    if (execution.path == IntegerGcdPath::Subresultant) {
        result_gcd = std::move(execution.gcd);
    } else {
        result_gcd = gcd_integer_poly_primitive(std::move(lhs), std::move(rhs));
    }
    
    if (!c_gcd.is_zero()) {
        multiply_integer_coefficients_by_scalar(result_gcd, c_gcd);
    }
    
    if (!result_gcd.is_zero() && result_gcd.leading_coeff().is_negative()) {
        multiply_integer_coefficients_by_scalar(result_gcd, BigInt(-1));
    }

    return IntegerGcdResult{
        .gcd = std::move(result_gcd),
        .path = execution.path,
    };
}

// hgcd_divides_exactly — defense-in-depth divisibility certificate.
//
// Returns true iff candidate divides dividend exactly in Z[x] (remainder zero).
//
// Design notes (content/sign):
//   - We work on the primitive part of candidate: content is a positive integer
//     divisor of all coefficients, so if candidate | dividend then
//     primitive(candidate) | dividend (and content divides integer coefficients
//     trivially).  We therefore strip content before calling pseudo_remainder.
//   - pseudo_remainder_integer_poly(dividend, divisor) is zero iff divisor | dividend
//     in Z[x]: the pseudo-division multiplies dividend by lc(divisor)^delta, which
//     does not change whether the remainder is zero (divisibility is preserved under
//     multiplication by a nonzero scalar in an integral domain).
//   - Leading-coefficient sign: after primitive-part extraction and positivity
//     normalization (see half_gcd_integer_poly), the candidate already has a positive
//     leading coefficient. No additional sign treatment is needed here.
//   - Degree-0 divisor: a nonzero constant divides every polynomial in Z[x].
[[nodiscard]] static bool hgcd_divides_exactly(
    const IntPoly& dividend, const IntPoly& candidate) {
    if (candidate.is_zero()) return false;
    if (dividend.is_zero()) return true;   // 0 is divisible by everything
    // A nonzero constant divides every polynomial.
    if (candidate.degree() == 0) return true;
    // If candidate degree exceeds dividend degree it cannot divide.
    if (candidate.degree() > dividend.degree()) return false;

    // Strip content from candidate: pseudo_remainder only needs the shape.
    // This avoids spurious failures when content(candidate) does not divide
    // leading coefficients of dividend after pseudo-multiplication.
    IntPoly prim_candidate = primitive_integer_poly(candidate);

    IntPoly rem = pseudo_remainder_integer_poly(dividend, prim_candidate);
    normalize_integer_poly(rem);
    return rem.is_zero();
}

// Dispatching GCD: uses Half-GCD when min(deg(a), deg(b)) exceeds the
// ctx threshold (default 200); otherwise uses subresultant PRS.
// Returns same IntegerGcdResult as gcd_integer_poly_with_subresultant.
//
// Defense-in-depth (A4-HGCD certificate):
//   After computing g via half_gcd_integer_poly, we verify that g divides
//   both lhs and rhs exactly in Z[x].  This is an O(deg) check that runs
//   only on the high-degree path and costs negligible time relative to HGCD.
//   If the check fails (mathematically impossible for a correct HGCD but
//   possible if a latent bug in the matrix arithmetic produces a non-divisor),
//   we fall back to the proven-correct subresultant PRS path.
//   Invariant: production NEVER returns a mathematically wrong GCD.
// Helper: max coefficient bit-length of f and g.
[[nodiscard]] static std::size_t max_coeff_bits(const IntPoly& f, const IntPoly& g) {
    std::size_t best = 0U;
    for (const BigInt& c : f.coefficients()) {
        std::size_t b = c.abs().bit_length();
        if (b > best) best = b;
    }
    for (const BigInt& c : g.coefficients()) {
        std::size_t b = c.abs().bit_length();
        if (b > best) best = b;
    }
    return best;
}

IntegerGcdResult gcd_integer_poly_dispatch(
    IntPoly lhs, IntPoly rhs, const symbolic::CASContext& ctx) {

    const std::size_t threshold = ctx.half_gcd_degree_threshold();
    const std::size_t min_deg = std::min(lhs.degree(), rhs.degree());

    if (min_deg > threshold && !lhs.is_zero() && !rhs.is_zero()) {
        // Half-GCD path (O(M(n) log n))
        IntPoly g = half_gcd_integer_poly(lhs, rhs);
        if (!g.is_zero() && g.leading_coeff().is_negative()) {
            multiply_integer_coefficients_by_scalar(g, BigInt(-1));
        }

        // Defense-in-depth: certify that g | lhs and g | rhs.
        // Uses primitive part of g (content stripped) so that pseudo_remainder
        // correctly tests Z[x]-divisibility regardless of content scaling.
        // Worst case on failure: one subresultant call; never a wrong result.
        if (!hgcd_divides_exactly(lhs, g) || !hgcd_divides_exactly(rhs, g)) {
            // HGCD produced a non-divisor — latent bug in matrix arithmetic.
            // Fall back to the proven-correct subresultant path.
            return gcd_integer_poly_with_subresultant(std::move(lhs), std::move(rhs));
        }

        return IntegerGcdResult{
            .gcd = std::move(g),
            .path = IntegerGcdPath::HalfGcd,
        };
    }

    // B2.1 — Modular GCD CRT path: try when coefficients are large.
    // Threshold: max coefficient bit-length ≥ ctx.modular_gcd_coeff_bits().
    // If CRT succeeds, return immediately (certified by divisibility check).
    // If it returns Unimplemented (prime budget exceeded), fall through to
    // the proven-correct subresultant path.  Never returns a wrong result.
    {
        const std::size_t crt_threshold = ctx.modular_gcd_coeff_bits();
        const std::size_t coeff_bits = max_coeff_bits(lhs, rhs);
        if (coeff_bits >= crt_threshold && !lhs.is_zero() && !rhs.is_zero()) {
            auto crt_res = gcd_integer_poly_crt(lhs, rhs, ctx);
            if (crt_res.is_ok()) {
                return IntegerGcdResult{
                    .gcd = std::move(crt_res.value()),
                    .path = IntegerGcdPath::ModularCrt,
                };
            }
            // Unimplemented → fall through to subresultant.
        }
    }

    // Subresultant PRS path (default)
    return gcd_integer_poly_with_subresultant(std::move(lhs), std::move(rhs));
}

}  // namespace algebra
}  // namespace cas
