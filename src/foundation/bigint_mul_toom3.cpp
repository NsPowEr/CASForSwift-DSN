// bigint_mul_toom3.cpp — Toom-Cook 3-way multiplication for BigInt.
//
// Algorithm: Toom-Cook 3 (Toom-3).
// Reference: M. Bodrato, "Towards Optimal Toom-Cook Multiplication for Univariate
//   and Multivariate Polynomials in Characteristic 2 and 0", 2007 (Algorithm 4).
//   Also: Brent-Zimmermann "Modern Computer Arithmetic" §1.3.3 (2010).
//
// Complexity: O(n^{log_3(5)}) ≈ O(n^{1.465}) vs Karatsuba O(n^{1.585}).
//
// Split each operand into 3 parts of equal size k = ceil(n/3):
//   A = A2·B^{2k} + A1·B^k + A0
//   B = B2·B^{2k} + B1·B^k + B0
// 5 evaluation points: {0, 1, -1, 2, ∞}.
// 5 recursive multiplications, then interpolation.
//
// Dispatch:
//   n < kToom3Threshold (64)         → Karatsuba (bigint_arithmetic.cpp)
//   n in [64, kFFTThreshold=4096)   → Toom-3 (this file)
//   n >= kFFTThreshold              → Karatsuba fallback (Schönhage-Strassen deferred)
//
// HARDCODE-OF-PASSAGE HPP-F1.1-MUL:
//   kFFTThreshold = 4096 — falls back to Karatsuba for n >= 4096.
//   Schönhage-Strassen FFT (O(n log n log log n)) deferred as Aperta permanente.
//   See HARDCODE_LEDGER.md.

#include "cas/bigint.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace cas {

// Toom-3 threshold: use Toom-3 when max(n_lhs, n_rhs) >= kToom3Threshold.
static constexpr std::size_t kToom3Threshold = 64U;

// Maximum limbs for Toom-3. Above this, fall back to Karatsuba.
static constexpr std::size_t kFFTThreshold = 8192U;

namespace {

// Signed magnitude BigInt: (magnitude, is_negative).
// Used only within Toom-3 interpolation to track sign.
struct SignedBig {
    BigInt mag;
    bool   neg = false;

    explicit SignedBig(BigInt m, bool n = false) : mag(std::move(m)), neg(n) {}
    explicit SignedBig(int v) : mag(v < 0 ? -v : v), neg(v < 0) {}

    [[nodiscard]] bool is_zero() const {
        return mag == BigInt(0);
    }

    // Negate
    [[nodiscard]] SignedBig operator-() const {
        if (is_zero()) return SignedBig(BigInt(0), false);
        return SignedBig(mag, !neg);
    }
};

[[nodiscard]] static SignedBig sadd(const SignedBig& a, const SignedBig& b) {
    if (a.is_zero()) return b;
    if (b.is_zero()) return a;
    if (a.neg == b.neg) {
        // Same sign: add magnitudes.
        // Use public a.mag + b.mag (both non-negative) = add magnitudes.
        // Create non-negative BigInt temporaries to force unsigned path.
        BigInt sum = a.mag + b.mag;  // both positive, so result is positive
        return SignedBig(std::move(sum), a.neg);
    }
    // Different signs: use public subtraction. We need |a| - |b|.
    // Compare magnitudes via compare_magnitude_pub.
    int cmp = BigInt::compare_magnitude_pub(a.mag, b.mag);
    if (cmp == 0) return SignedBig(BigInt(0), false);
    if (cmp > 0) {
        // |a| > |b|, result has a's sign.
        BigInt diff = a.mag - b.mag;  // positive BigInt minus positive BigInt
        return SignedBig(std::move(diff), a.neg);
    }
    BigInt diff = b.mag - a.mag;
    return SignedBig(std::move(diff), b.neg);
}

[[nodiscard]] static SignedBig ssub(const SignedBig& a, const SignedBig& b) {
    return sadd(a, -b);
}

[[nodiscard]] static SignedBig smul_small(const SignedBig& a, std::uint32_t small_val) {
    // a * small_val  (small_val is an unsigned small constant like 2, 3, 4, 6, 8)
    // Use repeated addition or direct multiply via shift for powers of 2.
    // General case: BigInt multiplication by small unsigned int.
    // We rely on BigInt * BigInt (schoolbook kicks in for small values).
    BigInt result = a.mag * BigInt(small_val);
    return SignedBig(std::move(result), a.neg);
}

[[nodiscard]] static SignedBig sdiv_exact(const SignedBig& a, std::uint32_t divisor) {
    // Exact integer division (no remainder).
    BigInt q = a.mag / BigInt(divisor);
    return SignedBig(std::move(q), a.neg);
}

// Shift left by bits (multiply by 2^bits)
[[nodiscard]] static SignedBig sshift_left(const SignedBig& a, std::size_t bits) {
    return SignedBig(a.mag.shift_left_bits(bits), a.neg);
}

// Shift right by bits (divide by 2^bits) — exact (no fractional part)
[[nodiscard]] static SignedBig sshift_right(const SignedBig& a, std::size_t bits) {
    return SignedBig(a.mag.shift_right_bits(bits), a.neg);
}

// Extract a "slice" of bigint starting at limb offset `start`, length `len`.
[[nodiscard]] static BigInt limb_slice(const BigInt& src, std::size_t start, std::size_t len) {
    if (start >= src.limb_count()) return BigInt(0);
    const std::size_t actual = std::min(len, src.limb_count() - start);
    std::vector<std::uint32_t> limbs(actual);
    for (std::size_t i = 0; i < actual; ++i) {
        limbs[i] = src.limb_at(start + i);
    }
    return BigInt::from_limbs_le(std::move(limbs));
}

} // namespace

// Toom-3 multiplication (non-negative integers).
// Uses Bodrato's "Algorithm 4" interpolation sequence.
BigInt BigInt::multiply_magnitude_toom3(const BigInt& lhs, const BigInt& rhs) {
    const std::size_t n = std::max(lhs.limb_count(), rhs.limb_count());

    // Guard: if below Toom-3 threshold or above max, delegate.
    if (n < kToom3Threshold || n >= kFFTThreshold) {
        return multiply_magnitude(lhs, rhs);
    }

    // k = ceil(n/3) — split size in limbs.
    const std::size_t k = (n + 2U) / 3U;
    const std::size_t k_bits = k * 32U;

    // Split lhs = A2*B^{2k} + A1*B^k + A0
    //      rhs = B2*B^{2k} + B1*B^k + B0
    const BigInt A0 = limb_slice(lhs, 0,     k);
    const BigInt A1 = limb_slice(lhs, k,     k);
    const BigInt A2 = limb_slice(lhs, 2*k,   k);
    const BigInt B0 = limb_slice(rhs, 0,     k);
    const BigInt B1 = limb_slice(rhs, k,     k);
    const BigInt B2 = limb_slice(rhs, 2*k,   k);

    // Evaluate at 5 points using signed arithmetic.
    // P(w) = A0 + A1*w + A2*w²,  Q(w) = B0 + B1*w + B2*w²
    //
    // p0  = A0                  q0  = B0
    // pinf = A2                 qinf = B2
    // p1  = A0 + A1 + A2        q1  = B0 + B1 + B2
    // pm1 = A0 - A1 + A2        qm1 = B0 - B1 + B2
    // p2  = A0 + 2*A1 + 4*A2   q2  = B0 + 2*B1 + 4*B2

    SignedBig p0   {A0,            false};
    SignedBig pinf {A2,            false};
    SignedBig q0   {B0,            false};
    SignedBig qinf {B2,            false};

    // p1 = A0 + A1 + A2
    SignedBig p1   = sadd(sadd(SignedBig{A0,false}, SignedBig{A1,false}), SignedBig{A2,false});
    SignedBig q1   = sadd(sadd(SignedBig{B0,false}, SignedBig{B1,false}), SignedBig{B2,false});

    // pm1 = A0 - A1 + A2
    SignedBig pm1  = sadd(ssub(SignedBig{A0,false}, SignedBig{A1,false}), SignedBig{A2,false});
    SignedBig qm1  = sadd(ssub(SignedBig{B0,false}, SignedBig{B1,false}), SignedBig{B2,false});

    // p2 = A0 + 2*A1 + 4*A2
    SignedBig p2   = sadd(SignedBig{A0,false},
                    sadd(SignedBig{A1.shift_left_bits(1), false},
                         SignedBig{A2.shift_left_bits(2), false}));
    SignedBig q2   = sadd(SignedBig{B0,false},
                    sadd(SignedBig{B1.shift_left_bits(1), false},
                         SignedBig{B2.shift_left_bits(2), false}));

    // 5 recursive multiplications (products of ~k-limb numbers)
    // r(w) = P(w) * Q(w)
    SignedBig r0   {multiply_magnitude(p0.mag,   q0.mag),   p0.neg   != q0.neg};
    SignedBig r1   {multiply_magnitude(p1.mag,   q1.mag),   p1.neg   != q1.neg};
    SignedBig rm1  {multiply_magnitude(pm1.mag,  qm1.mag),  pm1.neg  != qm1.neg};
    SignedBig r2   {multiply_magnitude(p2.mag,   q2.mag),   p2.neg   != q2.neg};
    SignedBig rinf {multiply_magnitude(pinf.mag, qinf.mag), pinf.neg != qinf.neg};

    // Interpolation (Bodrato Algorithm 4 / standard Toom-3 coefficients):
    // We recover C0..C4 such that:
    //   r(w) = C0 + C1*w + C2*w² + C3*w³ + C4*w⁴
    //
    // Step sequence (reference: Bodrato 2007, Algorithm 4):
    //
    //   r0   = C0
    //   rinf = C4
    //
    //   -- Step A: reduce r1 and rm1 --
    //   r1   ← r1   - r0                     [= C1 + C2 + C3 + C4]
    //   rm1  ← rm1  - r0                     [= -C1 + C2 - C3 + C4]
    //   r2   ← r2   - r0                     [= 2C1 + 4C2 + 8C3 + 16C4]
    //
    //   -- Step B: halve --
    //   t1   = (r1 + rm1) / 2                 [= C2 + C4]        (even)
    //   t2   = (r1 - rm1) / 2                 [= C1 + C3]        (even? not necessarily — but C1 + C3 is what we want)
    //
    //   -- Step C: recover C2 --
    //   C2   = t1 - rinf                      [= C2]
    //
    //   -- Step D: recover C3, C1 --
    //   r2   ← r2 - 2*r1                     [= 2C2 + 4C3 + 8C4 - 2C3 - 2C4? No, careful]
    //   Actually use a cleaner known sequence:
    //
    // Cleaner standard sequence from GMP toom33 (gmp-impl.h):
    //   v0   = r(0)    = C0
    //   v1   = r(1)    = C0+C1+C2+C3+C4
    //   vm1  = r(-1)   = C0-C1+C2-C3+C4
    //   v2   = r(2)    = C0+2C1+4C2+8C3+16C4
    //   vinf = r(inf)  = C4
    //
    //   Step 1: vm1 = vm1 + v1       → 2*(C0+C2+C4)   [call it s1]
    //           v1  = v1  - vm1_orig → 2*(C1+C3)       [call it s2]  ← use original vm1!
    //           So we need to do: s2 = v1 - vm1; s1 = v1 + vm1.
    //
    //   s1 = v1 + vm1       = 2*(C0+C2+C4)
    //   s2 = v1 - vm1       = 2*(C1+C3)
    //   s1 = s1/2 - v0      = C2 + C4             [subtract v0 = C0]
    //   C2 = s1 - vinf      = C2
    //   s2 = s2/2           = C1 + C3
    //   v2 = v2 - v0        = 2C1 + 4C2 + 8C3 + 16C4
    //   v2 = (v2 - s2*2 - C2*4 + v0 + vinf*16) / 6
    //      Wait, let's redo:
    //   v2 - v0 = 2C1 + 4C2 + 8C3 + 16C4
    //   v2 - v0 - 2*s2 = 2C1 + 4C2 + 8C3 + 16C4 - 2*(C1+C3)
    //                  = 4C2 + 6C3 + 16C4
    //   (v2-v0-2*s2)/2 = 2C2 + 3C3 + 8C4
    //   (v2-v0-2*s2)/2 - 2*C2 - 8*vinf = 3*C3
    //   C3 = [(v2-v0-2*s2)/2 - 2*C2 - 8*vinf] / 3
    //   C1 = s2 - C3

    // Use r0, r1, rm1, r2, rinf as mutable working values:
    SignedBig s1 = sadd(r1, rm1);       // 2*(C0+C2+C4)
    SignedBig s2 = ssub(r1, rm1);       // 2*(C1+C3)

    s1 = sshift_right(s1, 1);           // C0+C2+C4
    s1 = ssub(s1, r0);                  // C2+C4
    const SignedBig C2 = ssub(s1, rinf);// C2

    s2 = sshift_right(s2, 1);           // C1+C3  (this is exact: r1-rm1 = 2*(C1+C3))

    // v2 - v0 - 2*s2:
    SignedBig tmp = ssub(r2, r0);                  // v2 - v0
    tmp = ssub(tmp, sshift_left(s2, 1));           // v2 - v0 - 2*s2  = 4C2 + 6C3 + 16C4
    tmp = sshift_right(tmp, 1);                    // 2C2 + 3C3 + 8C4
    tmp = ssub(tmp, sshift_left(C2, 1));           // 3C3 + 8C4   (subtract 2*C2)
    tmp = ssub(tmp, smul_small(rinf, 8));           // 3C3          (subtract 8*C4)
    const SignedBig C3 = sdiv_exact(tmp, 3);       // C3

    const SignedBig C1 = ssub(s2, C3);            // C1 = (C1+C3) - C3
    const SignedBig C0 = r0;
    const SignedBig C4 = rinf;

    // Assemble: result = C0 + C1*B^k + C2*B^{2k} + C3*B^{3k} + C4*B^{4k}
    // All coefficients must sum to a positive result for positive inputs.
    // We fold into a running signed accumulator and assert final sign is non-negative.

    auto apply_coeff = [](SignedBig acc, const SignedBig& coeff, std::size_t shift_bits) -> SignedBig {
        return sadd(acc, sshift_left(coeff, shift_bits));
    };

    SignedBig result = apply_coeff(C0,              C1, k_bits);
    result           = apply_coeff(std::move(result), C2, 2*k_bits);
    result           = apply_coeff(std::move(result), C3, 3*k_bits);
    result           = apply_coeff(std::move(result), C4, 4*k_bits);

    // Toom-3 of two positive integers always yields a positive result.
    assert(!result.neg && "Toom-3: result should be non-negative for positive inputs");
    return result.mag;
}

}  // namespace cas
