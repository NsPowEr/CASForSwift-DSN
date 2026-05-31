// bigint_div_knuth_d.cpp — Knuth Algorithm D for BigInt division.
//
// Algorithm: Knuth TAOCP Vol 2 §4.3.1 Algorithm D (normalized long division).
// Reference: Donald E. Knuth, "The Art of Computer Programming, Volume 2:
//   Seminumerical Algorithms", 3rd edition, §4.3.1 Algorithm D, pp. 272-276.
//
// Also referenced in: GMP library internals (mpn_divrem), and
// Brent-Zimmermann "Modern Computer Arithmetic" §1.4.
//
// Complexity: O(m * n) where m = limb_count(dividend) - limb_count(divisor),
//             n = limb_count(divisor). Each inner loop step is O(n).
// This replaces the previous naive bit-shift loop which was
// O(bit_length(dividend)² / 32) inner-loop iterations.
//
// Implementation follows Knuth §4.3.1 very closely, using uint32 limbs (base B=2^32).
//
// PRECONDITIONS (enforced by asserts):
//   - v (divisor) has >= 2 limbs (single-limb case handled by divide_by_small).
//   - compare_magnitude(u, v) >= 0 (u >= v).
//   - Both u and v are non-negative.

#include "cas/bigint.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace cas {

// ── Knuth Algorithm D ──────────────────────────────────────────────────────
//
// The algorithm works as follows:
//
// D1. Find the normalization shift d: smallest integer such that
//     v[n-1] * 2^d >= B/2 (i.e. the top bit of v[n-1] is set after shifting).
//     Shift both u and v left by d bits. Append a zero limb to u.
//
// D2-D7. For j = m-1 downto 0 (where m = len(u_orig) - len(v)):
//   D3. Compute quotient estimate q̂ = floor((u[j+n]*B + u[j+n-1]) / v[n-1]).
//       Refine: while q̂*v[n-2] > (remainder)*B + u[j+n-2], decrement q̂.
//       Knuth Theorem B guarantees at most 2 decrements.
//   D4. Multiply-subtract: u[j..j+n] -= q̂ * v[0..n-1].
//   D5. If borrow (u went negative), add-back v[0..n-1] once and decrement q̂.
//   D7. Store q[j] = q̂.
//
// D8. Denormalize: remainder = u[0..n-1] >> d.
//
// Notes:
// - "limbs" are stored little-endian (index 0 = least significant).
// - When Knuth writes u[j], we use u_limbs[j].

std::pair<BigInt, BigInt> BigInt::divide_knuth_d(const BigInt& u_in, const BigInt& v_in) {
    const std::size_t m_plus_n = u_in.limb_count();  // n+m in Knuth notation
    const std::size_t n = v_in.limb_count();
    assert(n >= 2U);
    assert(compare_magnitude(u_in, v_in) >= 0);

    const std::size_t m = m_plus_n - n;  // number of quotient digits

    // D1. Normalize.
    // Find shift d such that top bit of v[n-1] is set after left-shift by d.
    const std::uint32_t vn1 = v_in.limb_at(n - 1U);
    assert(vn1 != 0U);
    const std::size_t d = static_cast<std::size_t>(__builtin_clz(vn1));

    // Build normalized v_limbs of length n.
    // v_norm[i] = v shifted left by d bits.
    BigInt v_norm = v_in.shift_left_bits(d);
    // v_norm may have grown a limb; ensure it has exactly n limbs (the top shift
    // bit is absorbed into the normalization — v[n-1] after shift has top bit set).
    // Due to shift, v_norm may produce n+1 limbs if d=0 and top bit was already set;
    // in that case, shrink. For d>0, v_norm has exactly n limbs.
    std::vector<std::uint32_t> v_limbs(n, 0U);
    for (std::size_t i = 0U; i < n && i < v_norm.limb_count(); ++i) {
        v_limbs[i] = v_norm.limb_at(i);
    }

    // Build normalized u_limbs of length m+n+1.
    // u_norm = u shifted left by d bits, with an extra zero limb at position m+n.
    BigInt u_norm = u_in.shift_left_bits(d);
    std::vector<std::uint32_t> u_limbs(m + n + 1U, 0U);
    for (std::size_t i = 0U; i < u_norm.limb_count() && i < m + n + 1U; ++i) {
        u_limbs[i] = u_norm.limb_at(i);
    }
    // u[m+n] may have a carry from the shift; collect it.
    if (u_norm.limb_count() == m + n + 1U) {
        u_limbs[m + n] = u_norm.limb_at(m + n);
    }

    const std::uint64_t B = static_cast<std::uint64_t>(1U) << 32U;
    const std::uint64_t vn1_norm = v_limbs[n - 1U];
    const std::uint64_t vn2_norm = (n >= 2U) ? v_limbs[n - 2U] : 0U;

    // Allocate quotient.
    std::vector<std::uint32_t> q_limbs(m + 1U, 0U);

    // D2/D7. Main loop: j = m downto 0.
    for (std::size_t j = m + 1U; j > 0U; --j) {
        const std::size_t jj = j - 1U;  // 0-indexed quotient position

        // D3. Estimate q̂.
        const std::uint64_t u_jn   = u_limbs[jj + n];
        const std::uint64_t u_jn1  = u_limbs[jj + n - 1U];
        const std::uint64_t u_jn2  = (jj + n >= 2U) ? u_limbs[jj + n - 2U] : 0U;

        std::uint64_t q_hat;
        std::uint64_t r_hat;

        if (u_jn >= vn1_norm) {
            // q̂ saturates to B-1.
            q_hat = B - 1U;
            r_hat = u_jn + u_jn1 - vn1_norm * B;  // approximate; may overflow conceptually
            // But r_hat = u_jn * B + u_jn1 - q_hat * vn1_norm
            //           = u_jn * B + u_jn1 - (B-1)*vn1_norm
            // This can be large; clamp to valid range.
            r_hat = u_jn1 + (u_jn - vn1_norm) * B;
            // Actually, compute properly:
            r_hat = u_jn * B + u_jn1 - q_hat * vn1_norm;
        } else {
            const std::uint64_t num = u_jn * B + u_jn1;
            q_hat = num / vn1_norm;
            r_hat = num % vn1_norm;
        }

        // Refine q̂: at most 2 iterations (Knuth Theorem B ensures 2 corrections suffice).
        while (q_hat >= B || q_hat * vn2_norm > r_hat * B + u_jn2) {
            --q_hat;
            r_hat += vn1_norm;
            if (r_hat >= B) break;  // r_hat overflow: no more refinement
        }

        // D4. Multiply-subtract: u[jj..jj+n] -= q_hat * v[0..n-1].
        // Standard Knuth D4 with separate mul-carry and sub-borrow tracking.
        // q_hat * v[i] < B^2, so we split into low/high 32-bit halves.
        {
            std::uint64_t mul_carry = 0U;   // carry from multiplication
            std::uint64_t sub_borrow = 0U;  // borrow from subtraction

            for (std::size_t i = 0U; i <= n; ++i) {
                const std::uint64_t vi = (i < n) ? v_limbs[i] : 0U;
                // product = q_hat * v[i] + mul_carry (fits in ~64 bits since q_hat<B, vi<B)
                const std::uint64_t prod = q_hat * vi + mul_carry;
                const std::uint64_t prod_lo = prod & 0xFFFFFFFFULL;
                mul_carry = prod >> 32U;
                // Subtract (prod_lo + sub_borrow) from u[jj+i].
                const std::uint64_t ui = u_limbs[jj + i];
                const std::uint64_t total_sub = prod_lo + sub_borrow;
                if (ui >= total_sub) {
                    u_limbs[jj + i] = static_cast<std::uint32_t>(ui - total_sub);
                    sub_borrow = 0U;
                } else {
                    u_limbs[jj + i] = static_cast<std::uint32_t>(B + ui - total_sub);
                    sub_borrow = 1U;
                }
            }
            const bool overall_borrow = (mul_carry > 0U || sub_borrow > 0U);

            // D5. Test remainder (borrow): if negative, add back.
            if (overall_borrow) {
                --q_hat;  // q_hat was 1 too large
                std::uint64_t carry = 0U;
                for (std::size_t i = 0U; i <= n; ++i) {
                    const std::uint64_t vi = (i < n) ? v_limbs[i] : 0U;
                    const std::uint64_t sum = static_cast<std::uint64_t>(u_limbs[jj + i]) + vi + carry;
                    u_limbs[jj + i] = static_cast<std::uint32_t>(sum);
                    carry = sum >> 32U;
                }
                // The carry-out here cancels the borrow (Knuth proof). Ignore it.
            }
        }

        // D7. Store q[jj] = q̂.
        q_limbs[jj] = static_cast<std::uint32_t>(q_hat);
    }

    // D8. Denormalize the remainder.
    // The remainder occupies u_limbs[0..n-1] (shifted left by d).
    // Shift right by d to undo the normalization.
    BigInt rem_shifted = from_limbs_le(
        std::vector<std::uint32_t>(u_limbs.begin(),
                                   u_limbs.begin() + static_cast<std::ptrdiff_t>(n)));
    BigInt remainder = rem_shifted.shift_right_bits(d);

    BigInt quotient = from_limbs_le(std::move(q_limbs));
    return {std::move(quotient), std::move(remainder)};
}

}  // namespace cas
