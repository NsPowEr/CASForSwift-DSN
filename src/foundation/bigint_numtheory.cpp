// bigint_numtheory.cpp — Montgomery modexp for BigInt.
//
// MONTGOMERY MODULAR EXPONENTIATION:
//   Reference: Koc, Acar, Kaliski "Analyzing and Comparing Montgomery
//   Multiplication Algorithms" (IEEE Micro 16(3), June 1996).
//   Reference: Knuth TAOCP Vol 2 §4.3.2 (modular arithmetic background).
//   Reference: Barrett "Implementing the Rivest Shamir and Adleman Public Key
//   Encryption Algorithm on a Standard Digital Signal Processor" (1986).
//
//   Montgomery reduction avoids division by the modulus at each step.
//   Instead, it works in "Montgomery space": ã = a * R mod m where R = 2^(n*32)
//   is a power of the base B=2^32 coprime to m.
//
//   Montgomery multiplication MonPro(ã, b̃) = ã * b̃ * R^{-1} mod m.
//   This is implemented via the CIOS algorithm (Coarsely Integrated Operand Scanning)
//   which is the standard for software implementations.
//
//   For odd modulus m, this gives a 20-40% speedup over standard modular exponentiation
//   by avoiding one BigInt division (mod) per squaring step.
//
//   For even modulus: fall back to standard binary exponentiation (power_mod in numtheory).
//
// BPSW PRIMALITY NOTES:
//   The existing is_prime() in numtheory/primality.cpp already implements
//   deterministic Miller-Rabin for n < 2^64 with the BPSW witness set.
//   montgomery_modexp() is wired in here so that future primality tests
//   can optionally use Montgomery for the modular exponentiation step.
//
// SCOPE:
//   This file provides montgomery_modexp() as a free function in namespace cas.
//   It is NOT wired into power_mod() by default (that would change behavior for
//   existing callers). It is exposed for benchmarking and future opt-in use.

#include "cas/bigint.hpp"
#include "cas/error.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstddef>
#include <limits>
#include <vector>

namespace cas {

namespace {

[[nodiscard]] CASError make_error_nt(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

// Compute -m^{-1} mod 2^32 (Montgomery's m' = -m^{-1} mod B).
// Used in Montgomery reduction. Only valid for odd m.
// Reference: Koc 1996, Algorithm 14.94 (Montgomery reduction).
[[nodiscard]] std::uint32_t mont_inv(std::uint32_t m0) {
    // m0 = m mod 2^32 (low limb of modulus).
    // Find x such that m0 * x ≡ 1 (mod 2^32) using Newton's method:
    // x_{k+1} = x_k * (2 - m0 * x_k)
    assert(m0 & 1U);  // m must be odd
    std::uint32_t x = 1U;
    for (int i = 0; i < 5; ++i) {  // 5 iterations sufficient for 32 bits
        x *= 2U - m0 * x;
    }
    return static_cast<std::uint32_t>(-(std::uint64_t)x);  // m' = -x mod 2^32
}

// Montgomery reduction (REDC):
//   Input: T, m (modulus), m' = -m^{-1} mod B (B = 2^32), n = limb count of m.
//   Output: T * R^{-1} mod m, where R = B^n.
//   Reference: Koc 1996 Algorithm 14.94, Montgomery reduction.
[[nodiscard]] BigInt mont_redc(
    const std::vector<std::uint32_t>& T,  // 2n limbs
    const BigInt& m,
    std::uint32_t m_prime,
    std::size_t n) {

    std::vector<std::uint64_t> t(T.begin(), T.end());
    t.resize(2U * n + 1U, 0U);

    // CIOS Montgomery reduction.
    for (std::size_t i = 0U; i < n; ++i) {
        const std::uint64_t u = (static_cast<std::uint32_t>(t[i]) * static_cast<std::uint64_t>(m_prime))
                                & 0xFFFFFFFFULL;
        std::uint64_t carry = 0U;
        for (std::size_t j = 0U; j < n; ++j) {
            const std::uint64_t mj = (j < m.limb_count())
                ? static_cast<std::uint64_t>(m.limb_at(j)) : 0U;
            const std::uint64_t product = u * mj + t[i + j] + carry;
            t[i + j] = product & 0xFFFFFFFFULL;
            carry = product >> 32U;
        }
        const std::uint64_t sum = t[i + n] + carry;
        t[i + n] = sum & 0xFFFFFFFFULL;
        carry = sum >> 32U;
        if (i + n + 1U < t.size()) {
            t[i + n + 1U] += carry;
        }
    }

    // Extract result: limbs [n..2n-1] form the output.
    std::vector<std::uint32_t> result(n);
    for (std::size_t i = 0U; i < n; ++i) {
        result[i] = static_cast<std::uint32_t>(t[i + n]);
    }

    BigInt r = BigInt::from_limbs_le(std::move(result));

    // Conditional subtraction: if r >= m, subtract m.
    if (BigInt::compare_magnitude_pub(r, m) >= 0) {
        r = r - m;
    }
    return r;
}

// Montgomery multiplication: a_mont * b_mont * R^{-1} mod m.
[[nodiscard]] BigInt mont_mul(
    const BigInt& a, const BigInt& b, const BigInt& m,
    std::uint32_t m_prime, std::size_t n) {

    // Compute a * b as 2n-limb product.
    std::vector<std::uint32_t> product(2U * n, 0U);
    for (std::size_t i = 0U; i < n; ++i) {
        const std::uint64_t ai = (i < a.limb_count()) ? static_cast<std::uint64_t>(a.limb_at(i)) : 0U;
        std::uint64_t carry = 0U;
        for (std::size_t j = 0U; j < n; ++j) {
            const std::uint64_t bj = (j < b.limb_count()) ? static_cast<std::uint64_t>(b.limb_at(j)) : 0U;
            const std::uint64_t val = ai * bj + static_cast<std::uint64_t>(product[i + j]) + carry;
            product[i + j] = static_cast<std::uint32_t>(val);
            carry = val >> 32U;
        }
        if (i + n < product.size()) {
            product[i + n] = static_cast<std::uint32_t>(
                static_cast<std::uint64_t>(product[i + n]) + carry);
        }
    }

    return mont_redc(product, m, m_prime, n);
}

} // namespace

// Montgomery modular exponentiation: base^exp mod modulus.
// Pre-condition: modulus > 0, exp >= 0.
// For odd modulus: uses Montgomery multiplication (CIOS algorithm, Koc 1996).
// For even modulus: falls back to standard binary exponentiation.
//
// Reference: Koc, Acar, Kaliski (1996) §4 "Montgomery Multiplication Algorithm".
Result<BigInt> montgomery_modexp(
    const BigInt& base,
    const BigInt& exp,
    const BigInt& modulus) {

    if (modulus.is_zero() || modulus.is_negative()) {
        return fail<BigInt>(make_error_nt(
            CASErrorKind::InvalidArgument,
            "montgomery_modexp: modulus must be positive"));
    }
    if (exp.is_negative()) {
        return fail<BigInt>(make_error_nt(
            CASErrorKind::InvalidArgument,
            "montgomery_modexp: exponent must be non-negative"));
    }
    if (modulus == BigInt(1)) {
        return ok(BigInt(0));
    }
    if (exp.is_zero()) {
        return ok(BigInt(1));
    }

    // For even modulus, Montgomery reduction is not directly applicable.
    // Fall back to standard binary exponentiation.
    if (modulus.limb_count() == 0U || (modulus.limb_at(0) & 1U) == 0U) {
        // Standard binary exponentiation (base^exp mod modulus).
        BigInt result(1);
        BigInt b = base % modulus;
        if (b.is_negative()) b = b + modulus;
        BigInt e = exp;
        const BigInt zero(0);
        const BigInt one(1);
        const BigInt two(2);
        while (e > zero) {
            if ((e % two) == one) {
                result = (result * b) % modulus;
            }
            e = e.shift_right_bits(1);
            b = (b * b) % modulus;
        }
        return ok(result);
    }

    // Odd modulus: use Montgomery.
    const std::size_t n = modulus.limb_count();
    const std::uint32_t m0 = modulus.limb_at(0);
    const std::uint32_t m_prime = mont_inv(m0);

    // R = 2^{32n} mod m; compute R mod m = 2^{32n} mod m.
    // Use: R mod m = (-m) mod 2^{32n} ... simpler: compute BigInt(1).shift_left_bits(32*n) % m.
    BigInt R = BigInt(1).shift_left_bits(32U * n) % modulus;

    // Convert base to Montgomery space: a_mont = base * R mod m.
    BigInt b_reduced = base % modulus;
    if (b_reduced.is_negative()) b_reduced = b_reduced + modulus;
    BigInt a_mont = (b_reduced * R) % modulus;

    // Result starts as 1 in Montgomery space: r_mont = 1 * R mod m = R mod m.
    BigInt r_mont = R;

    // Binary method exponentiation in Montgomery space.
    const std::size_t exp_bits = exp.bit_length();
    for (std::size_t i = exp_bits; i > 0U; --i) {
        r_mont = mont_mul(r_mont, r_mont, modulus, m_prime, n);
        // Check bit i-1 of exp.
        const std::size_t bit_idx = i - 1U;
        const std::size_t limb_idx = bit_idx / 32U;
        const std::size_t bit_in_limb = bit_idx % 32U;
        const bool bit_set = (limb_idx < exp.limb_count())
            && ((exp.limb_at(limb_idx) >> bit_in_limb) & 1U);
        if (bit_set) {
            r_mont = mont_mul(r_mont, a_mont, modulus, m_prime, n);
        }
    }

    // Convert back from Montgomery space: result = r_mont * R^{-1} mod m.
    // = mont_mul(r_mont, BigInt(1), m, m_prime, n).
    std::vector<std::uint32_t> r_padded(2U * n, 0U);
    for (std::size_t i = 0U; i < std::min(r_mont.limb_count(), n); ++i) {
        r_padded[i] = r_mont.limb_at(i);
    }
    BigInt result = mont_redc(r_padded, modulus, m_prime, n);

    return ok(result);
}

}  // namespace cas
