// bigint_div_burnikel_ziegler.cpp — Recursive division (Burnikel-Ziegler 1998).
//
// References:
//   - Christoph Burnikel and Joachim Ziegler, "Fast Recursive Division",
//     Research Report MPI-I-98-1-022, MPI für Informatik, Saarbrücken, 1998.
//   - Richard P. Brent and Paul Zimmermann, "Modern Computer Arithmetic",
//     Cambridge University Press, 2010, §1.4.3.
//
// Complexity: O(M(n) · log n) where M(n) is the cost of n-limb multiplication.
// With Toom-3 multiplication (M(n) ≈ n^1.46) Burnikel-Ziegler beats Knuth-D
// (O(n²)) once n exceeds the threshold kBzThreshold.
//
// Structure (Brent-Zimmermann §1.4.3):
//   bz_div_2by1(A, B, n):
//     |A| ≤ 2n limbs, |B| = n limbs, A < B·β^n.
//     Returns (Q, R) with A = Q·B + R, 0 ≤ R < B, Q ≤ β^n.
//     Base case: small n or odd n → Knuth Algorithm D.
//     Recursive case: split A into 4 blocks of n/2 limbs and call bz_div_3by2 twice.
//
//   bz_div_3by2(A, B, n):
//     |A| ≤ 3n limbs, |B| = 2n limbs (= B1·β^n + B0).
//     Precondition: A < B·β^n.
//     Returns (Q, R) with A = Q·B + R, Q ≤ β^n.
//     Computes (Q, R1) via bz_div_2by1 on top 2n limbs of A and B1,
//     then subtracts Q·B0 with at most two corrections.
//
//   divide_burnikel_ziegler(u, v):
//     Top-level entry. Normalizes v (top bit set), iterates bz_div_2by1 on
//     blocks of u from high to low. Maintains the invariant that the
//     partial remainder is strictly less than v_norm at every step, so the
//     subsequent block-shift keeps the input to bz_div_2by1 within 2n limbs
//     and A < B·β^n.
//
// HARDCODE NOTE — kBzThreshold = 64 limbs:
//   BigInt arithmetic is context-free (no CASContext access), see HPP-020
//   for the architectural rationale. Threshold 64 derives from Brent-Zimmermann
//   §1.4.3 (BZ pays off once n ≥ 50-100 limbs, coinciding with Toom-3 onset).
//   Falls back to Knuth-D below threshold and for odd n. Correctness is
//   independent of the threshold value.

#include "cas/bigint.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas {

namespace {

constexpr std::size_t kBzThreshold = 64U;

[[nodiscard]] BigInt extract_limbs(const BigInt& x, std::size_t start, std::size_t end) {
    if (end <= start) return BigInt(0);
    std::vector<std::uint32_t> limbs;
    limbs.reserve(end - start);
    for (std::size_t i = start; i < end; ++i) {
        limbs.push_back(i < x.limb_count() ? x.limb_at(i) : 0U);
    }
    return BigInt::from_limbs_le(std::move(limbs));
}

[[nodiscard]] BigInt shift_left_limbs(const BigInt& x, std::size_t k) {
    if (k == 0U || x.is_zero()) return x;
    return x.shift_left_bits(k * 32U);
}

}  // namespace

// bz_div_3by2: |A| ≤ 3n limbs, |B| = 2n limbs (= B1·β^n + B0). A < B·β^n.
std::pair<BigInt, BigInt> BigInt::bz_div_3by2(const BigInt& A, const BigInt& B, std::size_t n) {
    // A = A2·β^(2n) + A1·β^n + A0 (each block n limbs, A2 may be ≤ n limbs).
    const std::size_t a_top_end = std::max(A.limb_count(), 3U * n);
    BigInt A0 = extract_limbs(A, 0U, n);
    BigInt A1 = extract_limbs(A, n, 2U * n);
    BigInt A2 = extract_limbs(A, 2U * n, a_top_end);
    BigInt B1 = extract_limbs(B, n, 2U * n);
    BigInt B0 = extract_limbs(B, 0U, n);

    BigInt Q;
    BigInt R1;
    if (BigInt::compare_magnitude_pub(A2, B1) < 0) {
        // Standard case: divide top 2n limbs of A by B1.
        BigInt top = shift_left_limbs(A2, n) + A1;  // = A2·β^n + A1, ≤ 2n limbs
        auto qr = bz_div_2by1(top, B1, n);
        Q = std::move(qr.first);
        R1 = std::move(qr.second);
    } else {
        // Overflow case: A2 >= B1 implies the true quotient is β^n - 1.
        std::vector<std::uint32_t> q_limbs(n, 0xFFFFFFFFU);
        Q = BigInt::from_limbs_le(std::move(q_limbs));
        // R1 = (A2·β^n + A1) - (β^n - 1)·B1 = (A2 - B1)·β^n + A1 + B1.
        BigInt A2_minus_B1 = A2 - B1;  // ≥ 0 by hypothesis
        R1 = shift_left_limbs(A2_minus_B1, n) + A1 + B1;
    }

    // D = Q · B0; R = R1·β^n + A0 - D.
    BigInt D = Q * B0;
    BigInt R_unsigned = shift_left_limbs(R1, n) + A0;

    BigInt R;
    bool negative;
    if (BigInt::compare_magnitude_pub(R_unsigned, D) >= 0) {
        R = R_unsigned - D;
        negative = false;
    } else {
        R = D - R_unsigned;
        negative = true;
    }

    // Correction loop: while R < 0, R += B and Q -= 1. At most two iterations
    // (Burnikel-Ziegler Lemma 2.2).
    for (int corr = 0; corr < 2 && negative; ++corr) {
        if (BigInt::compare_magnitude_pub(R, B) <= 0) {
            R = B - R;
            negative = false;
        } else {
            R = R - B;
        }
        Q = Q - BigInt(1);
    }
    assert(!negative);

    return {std::move(Q), std::move(R)};
}

// bz_div_2by1: |A| ≤ 2n limbs, |B| = n limbs, A < B·β^n.
std::pair<BigInt, BigInt> BigInt::bz_div_2by1(const BigInt& A, const BigInt& B, std::size_t n) {
    if (BigInt::compare_magnitude_pub(A, B) < 0) {
        return {BigInt(0), A};
    }
    // Base case: small or odd n → Knuth Algorithm D.
    if (n < kBzThreshold || (n & 1U) != 0U) {
        return divide_knuth_d(A, B);
    }
    const std::size_t k = n / 2U;
    // Top 3k limbs of A (= A3·β^(2k) + A2·β^k + A1).
    BigInt A_top = extract_limbs(A, k, 4U * k);
    BigInt A0 = extract_limbs(A, 0U, k);
    auto qr1 = bz_div_3by2(A_top, B, k);
    BigInt Q1 = std::move(qr1.first);
    BigInt R1 = std::move(qr1.second);
    BigInt mid = shift_left_limbs(R1, k) + A0;
    auto qr0 = bz_div_3by2(mid, B, k);
    BigInt Q0 = std::move(qr0.first);
    BigInt R0 = std::move(qr0.second);
    BigInt Q = shift_left_limbs(Q1, k) + Q0;
    return {std::move(Q), std::move(R0)};
}

std::pair<BigInt, BigInt> BigInt::divide_burnikel_ziegler(const BigInt& u, const BigInt& v) {
    assert(!v.is_zero());
    if (u.is_zero() || compare_magnitude(u, v) < 0) {
        return {BigInt(0), u};
    }
    if (v.limb_count() < kBzThreshold) {
        return divide_knuth_d(u, v);
    }
    // Normalize v so its top bit is set. Apply same shift to u.
    const std::uint32_t v_top = v.limb_at(v.limb_count() - 1U);
    assert(v_top != 0U);
    const std::size_t d = static_cast<std::size_t>(__builtin_clz(v_top));
    BigInt v_norm = v.shift_left_bits(d);
    BigInt u_norm = u.shift_left_bits(d);
    const std::size_t n = v_norm.limb_count();  // length of normalized divisor

    // Process u_norm in chunks of n limbs from high to low. Maintain the
    // running remainder R (always < v_norm) and accumulate the quotient.
    const std::size_t u_lc = u_norm.limb_count();
    const std::size_t num_chunks = (u_lc + n - 1U) / n;
    BigInt Q(0);
    BigInt R(0);
    for (std::size_t i = num_chunks; i > 0U; --i) {
        const std::size_t chunk_idx = i - 1U;
        BigInt block = extract_limbs(u_norm, chunk_idx * n, (chunk_idx + 1U) * n);
        // A = R·β^n + block, length ≤ 2n limbs (R < v_norm ⇒ R has ≤ n limbs).
        BigInt A = shift_left_limbs(R, n) + block;
        // A < v_norm·β^n holds because the top bit of v_norm is set:
        //   v_norm ≥ β^n / 2 + 1 ⇒ β^n ≤ 2·v_norm − 2 ⇒ block ≤ β^n − 1 ≤ v_norm,
        //   so A = R·β^n + block ≤ (v_norm − 1)·β^n + (β^n − 1) = v_norm·β^n − 1.
        auto qr = bz_div_2by1(A, v_norm, n);
        BigInt q_i = std::move(qr.first);
        R = std::move(qr.second);
        Q = shift_left_limbs(Q, n) + q_i;
    }
    BigInt remainder = R.shift_right_bits(d);
    return {std::move(Q), std::move(remainder)};
}

}  // namespace cas
