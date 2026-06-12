// bigint_gcd_lehmer.cpp — Binary GCD (Stein) and double-digit Lehmer GCD.
//
// BINARY GCD (Stein's Algorithm):
//   Reference: Knuth TAOCP Vol 2 §4.5.2 Algorithm B.
//   Avoids division entirely; uses halving (shift) and subtraction.
//
// LARGE-INTEGER GCD — Lehmer (Knuth Algorithm L + Jebelean double-digit surrogates):
//   References:
//     - Knuth D.E., The Art of Computer Programming Vol.2 §4.5.2 Algorithm L.
//     - Jebelean T., "A Double-Digit Lehmer-Euclid Algorithm for Finding the GCD
//       of Long Integers", JSC vol.19, 1995, pp.145-157.
//
//   Surrogates: top 64 bits of A and B taken at the same shift s = bit_length(A)-64.
//     â = A >> s   (∈ [2^63, 2^64))
//     b̂ = B >> s   (∈ [0,    2^64))
//   This handles na ≠ nb correctly: b̂ is naturally smaller (or 0) when B has
//   fewer leading bits than A. No bail-out on differing limb counts.
//
//   Simulated Euclidean iteration maintains a signed unimodular cofactor matrix
//   M = [[A_c, B_c], [C_c, D_c]] with alternating-sign entries (Knuth L3).
//   Validity (Knuth L3): a candidate quotient q is accepted iff
//       q == (â + A_c) / (b̂ + C_c)  AND  q == (â + B_c) / (b̂ + D_c)
//   The two divisions bracket the true full-integer quotient; equality
//   guarantees correctness of the simulated step (Knuth, Vol.2 §4.5.2 Thm L).
//
//   When simulation cannot make progress (B_c == 0 after L2), a single full
//   multi-precision Euclidean step A,B ← B, A mod B is performed to unblock.
//   Otherwise the accumulated matrix is applied in one shot:
//       A_new = A_c·A + B_c·B
//       B_new = C_c·A + D_c·B
//   yielding ~64-bit reduction per outer iteration.
//
//   Cofactor overflow protection: every multiplication q·C_c, q·D_c is checked
//   against int64 limits via __builtin_mul_overflow; on overflow the iteration
//   commits the current matrix and exits — correctness preserved.
//
// kLehmerThreshold = 16 limbs (512 bits):
//   Per-step overhead (extract 64-bit surrogate, run simulated Euclid, apply
//   matrix multiply on full integers) pays off once inputs exceed ~16 limbs.
//   HARDCODE HPP-020: BigInt is context-free (no CASContext access), so the
//   threshold cannot be exposed via ctx.* without an architectural change.
//   Documented as a legitimate hardware-safety limit (CLAUDE.md exception #4):
//   using Binary GCD below threshold is correct and only suboptimal.

#include "cas/bigint.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstddef>
#include <limits>
#include <utility>

namespace cas {

// Lehmer threshold (limb count).
static constexpr std::size_t kLehmerThreshold = 16U;

// ── Binary GCD (Stein) ────────────────────────────────────────────────────
BigInt binary_gcd(BigInt u, BigInt v) {
    u = u.abs();
    v = v.abs();

    if (u.is_zero()) return v;
    if (v.is_zero()) return u;

    auto trailing_zero_bits = [](const BigInt& x) -> std::size_t {
        std::size_t tz = 0;
        std::size_t low_limb = 0;
        while (low_limb < x.limb_count() && x.limb_at(low_limb) == 0U) {
            tz += 32U;
            ++low_limb;
        }
        if (low_limb < x.limb_count()) {
            tz += static_cast<std::size_t>(__builtin_ctz(x.limb_at(low_limb)));
        }
        return tz;
    };

    const std::size_t trailing_u = trailing_zero_bits(u);
    const std::size_t trailing_v = trailing_zero_bits(v);

    const std::size_t k = std::min(trailing_u, trailing_v);
    u = u.shift_right_bits(trailing_u);
    v = v.shift_right_bits(trailing_v);

    while (!v.is_zero()) {
        if (BigInt::compare_magnitude_pub(u, v) > 0) {
            std::swap(u, v);
        }
        v = v - u;
        if (v.is_zero()) break;
        v = v.shift_right_bits(trailing_zero_bits(v));
    }

    return u.shift_left_bits(k);
}

// ── Lehmer double-digit GCD (Knuth Algorithm L + Jebelean surrogates) ────

namespace {

// Extract top 64 bits of |x| starting at bit position `shift`.
// Returns floor(|x| / 2^shift) truncated to 64 bits. Returns 0 if x is zero
// or if shift >= bit_length(x).
std::uint64_t top64_at_shift(const BigInt& x, std::size_t shift) {
    if (x.is_zero()) return 0U;
    const std::size_t bl = x.bit_length();
    if (shift >= bl) return 0U;
    const BigInt shifted = x.shift_right_bits(shift);
    // shifted may have up to 65 bits if shift == bl-64 and shift is exact.
    // We want low 64 bits.
    if (shifted.limb_count() == 0U) return 0U;
    const std::uint64_t lo = static_cast<std::uint64_t>(shifted.limb_at(0));
    const std::uint64_t hi = shifted.limb_count() >= 2U
        ? static_cast<std::uint64_t>(shifted.limb_at(1))
        : 0ULL;
    return (hi << 32) | lo;
}

// Signed multiply with overflow detection. Returns false if overflow occurred.
bool safe_mul_i64(std::int64_t a, std::int64_t b, std::int64_t& out) {
    return !__builtin_mul_overflow(a, b, &out);
}

// Signed add with overflow detection. Returns false if overflow occurred.
bool safe_add_i64(std::int64_t a, std::int64_t b, std::int64_t& out) {
    return !__builtin_add_overflow(a, b, &out);
}

// Multiply BigInt by signed int64 coefficient. Returns sign-correct BigInt.
BigInt scale_bigint_i64(const BigInt& x, std::int64_t coeff) {
    if (coeff == 0) return BigInt(0);
    if (coeff == 1) return x;
    if (coeff == -1) return -x;
    const std::uint64_t mag = coeff > 0
        ? static_cast<std::uint64_t>(coeff)
        : static_cast<std::uint64_t>(-(coeff + 1)) + 1ULL;  // safe |INT64_MIN|
    BigInt r = x * BigInt::from_u64(mag);
    if (coeff < 0) r = -r;
    return r;
}

}  // namespace

BigInt lehmer_gcd(BigInt a, BigInt b) {
    a = a.abs();
    b = b.abs();

    if (a.is_zero()) return b;
    if (b.is_zero()) return a;

    // Ensure a >= b.
    if (BigInt::compare_magnitude_pub(a, b) < 0) {
        std::swap(a, b);
    }

    while (b.limb_count() >= 2U && a.bit_length() > 64U) {
        // Step L1: extract 64-bit surrogates aligned to A's top.
        const std::size_t bl_a = a.bit_length();
        const std::size_t shift = bl_a - 64U;

        const std::uint64_t a_hat_u = top64_at_shift(a, shift);
        const std::uint64_t b_hat_u = top64_at_shift(b, shift);

        if (b_hat_u == 0U) {
            // B is too small to give a useful surrogate: one full Euclidean step.
            BigInt rem = a % b;
            a = std::move(b);
            b = std::move(rem);
            if (b.is_zero()) return a;
            continue;
        }

        // Use signed int64 surrogates. Since a_hat_u ∈ [2^63, 2^64), it may exceed
        // INT64_MAX. We work with a/b ratios; cofactors stay bounded so signed
        // arithmetic on (â+A_c) etc. needs care. To stay safe, when a_hat_u or
        // b_hat_u exceeds INT64_MAX we right-shift both by 1 (preserves quotient).
        std::uint64_t a_hat_top = a_hat_u;
        std::uint64_t b_hat_top = b_hat_u;
        while (a_hat_top > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
            a_hat_top >>= 1;
            b_hat_top >>= 1;
        }
        if (b_hat_top == 0U) {
            BigInt rem = a % b;
            a = std::move(b);
            b = std::move(rem);
            if (b.is_zero()) return a;
            continue;
        }

        std::int64_t a_hat = static_cast<std::int64_t>(a_hat_top);
        std::int64_t b_hat = static_cast<std::int64_t>(b_hat_top);

        // Cofactor matrix M = [[A_c, B_c], [C_c, D_c]] with alternating signs.
        std::int64_t A_c = 1, B_c = 0, C_c = 0, D_c = 1;

        // Step L2: simulated Euclidean iterations with Knuth L3 validity check.
        for (;;) {
            // Validity preconditions: denominators (b̂+C_c) and (b̂+D_c) must be
            // non-zero and same sign as b̂. If matrix entries grew so large that
            // these no longer hold, commit current matrix.
            std::int64_t den1, den2;
            if (!safe_add_i64(b_hat, C_c, den1)) break;
            if (!safe_add_i64(b_hat, D_c, den2)) break;
            if (den1 == 0 || den2 == 0) break;
            // Both denominators must be positive (signs match positive b̂).
            if (den1 < 0 || den2 < 0) break;

            std::int64_t num1, num2;
            if (!safe_add_i64(a_hat, A_c, num1)) break;
            if (!safe_add_i64(a_hat, B_c, num2)) break;
            if (num1 < 0 || num2 < 0) break;

            const std::int64_t q1 = num1 / den1;  // upper bracket
            const std::int64_t q2 = num2 / den2;  // lower bracket
            if (q1 != q2 || q1 == 0) break;
            const std::int64_t q = q1;

            // Step L3: apply simulated step to matrix and surrogates.
            //   new A_c = C_c, new B_c = D_c
            //   new C_c = A_c - q·C_c, new D_c = B_c - q·D_c
            //   new â = b̂, new b̂ = â - q·b̂
            std::int64_t qC, qD, qB;
            if (!safe_mul_i64(q, C_c, qC)) break;
            if (!safe_mul_i64(q, D_c, qD)) break;
            if (!safe_mul_i64(q, b_hat, qB)) break;

            std::int64_t new_C, new_D, new_b;
            if (!safe_add_i64(A_c, -qC, new_C)) break;
            if (!safe_add_i64(B_c, -qD, new_D)) break;
            if (!safe_add_i64(a_hat, -qB, new_b)) break;

            A_c = C_c;
            B_c = D_c;
            C_c = new_C;
            D_c = new_D;
            a_hat = b_hat;
            b_hat = new_b;

            if (b_hat == 0) break;
        }

        // Step L4: apply accumulated matrix or fall back to one Euclidean step.
        if (B_c == 0) {
            // No progress from simulation: one full multi-precision Euclidean step.
            BigInt rem = a % b;
            a = std::move(b);
            b = std::move(rem);
            if (b.is_zero()) return a;
            continue;
        }

        // Apply M to (a, b):
        //   a_new = A_c·a + B_c·b
        //   b_new = C_c·a + D_c·b
        // Signs alternate so that |a_new|, |b_new| are positive when matrix is
        // applied to the original (a, b) pair. We compute signed combinations
        // and take absolute values defensively.
        BigInt a_new = scale_bigint_i64(a, A_c) + scale_bigint_i64(b, B_c);
        BigInt b_new = scale_bigint_i64(a, C_c) + scale_bigint_i64(b, D_c);

        a = a_new.abs();
        b = b_new.abs();

        if (BigInt::compare_magnitude_pub(a, b) < 0) {
            std::swap(a, b);
        }
    }

    // Fall through to plain Euclidean for small remainders.
    while (!b.is_zero()) {
        BigInt rem = a % b;
        a = std::move(b);
        b = std::move(rem);
    }
    return a;
}

// Dispatcher: choose algorithm based on input size.
BigInt gcd(BigInt lhs, BigInt rhs) {
    lhs = lhs.abs();
    rhs = rhs.abs();
    if (lhs.is_zero()) return rhs;
    if (rhs.is_zero()) return lhs;

    const std::size_t n = std::max(lhs.limb_count(), rhs.limb_count());

    if (n < 2U) {
        // Single-limb: Euclidean on raw uint64 is fastest.
        std::uint64_t a = lhs.to_u64();
        std::uint64_t b = rhs.to_u64();
        while (b != 0U) { a %= b; std::swap(a, b); }
        return BigInt::from_u64(a);
    }

    if (n >= kLehmerThreshold) {
        return lehmer_gcd(std::move(lhs), std::move(rhs));
    }

    // Binary GCD for [2, kLehmerThreshold) limbs.
    return binary_gcd(std::move(lhs), std::move(rhs));
}

}  // namespace cas
