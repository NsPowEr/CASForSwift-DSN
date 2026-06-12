#pragma once

#include "cas/result.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace cas {

class BigInt {
public:
    static constexpr std::uint32_t API_VERSION = 1;

    BigInt() = default;
    BigInt(std::int64_t value);
    [[nodiscard]] static BigInt from_u64(std::uint64_t value) noexcept;
    [[nodiscard]] static Result<BigInt> parse(std::string decimal);

    [[nodiscard]] const std::string& decimal() const;
    [[nodiscard]] bool is_negative() const noexcept;
    [[nodiscard]] bool is_zero() const noexcept;
    [[nodiscard]] bool is_positive() const noexcept;
    [[nodiscard]] std::size_t hash() const noexcept;

    [[nodiscard]] BigInt operator-() const;
    [[nodiscard]] BigInt abs() const;

    friend bool operator==(const BigInt& lhs, const BigInt& rhs) noexcept;
    friend bool operator!=(const BigInt& lhs, const BigInt& rhs) noexcept;
    friend bool operator<(const BigInt& lhs, const BigInt& rhs) noexcept;
    friend bool operator<=(const BigInt& lhs, const BigInt& rhs) noexcept;
    friend bool operator>(const BigInt& lhs, const BigInt& rhs) noexcept;
    friend bool operator>=(const BigInt& lhs, const BigInt& rhs) noexcept;
    friend Result<std::pair<BigInt, BigInt>> checked_divide_with_remainder(
        const BigInt& dividend,
        const BigInt& divisor);

    friend BigInt operator+(const BigInt& lhs, const BigInt& rhs);
    friend BigInt operator-(const BigInt& lhs, const BigInt& rhs);
    friend BigInt operator*(const BigInt& lhs, const BigInt& rhs);
    friend BigInt operator/(const BigInt& lhs, const BigInt& rhs);
    friend BigInt operator%(const BigInt& lhs, const BigInt& rhs);

    BigInt& operator+=(const BigInt& rhs);
    BigInt& operator-=(const BigInt& rhs);
    BigInt& operator*=(const BigInt& rhs);
    BigInt& operator/=(const BigInt& rhs);
    BigInt& operator%=(const BigInt& rhs);

    [[nodiscard]] std::size_t bit_length() const noexcept;
    [[nodiscard]] BigInt shift_left_bits(std::size_t bits) const;
    [[nodiscard]] BigInt shift_right_bits(std::size_t bits) const;
    [[nodiscard]] std::uint64_t to_u64() const noexcept;
    [[nodiscard]] double to_double() const noexcept;

    // F1.1 Toom-3 helpers: limb access for slice extraction.
    [[nodiscard]] std::size_t limb_count() const noexcept;
    [[nodiscard]] std::uint32_t limb_at(std::size_t index) const noexcept;

    // F1.1 public compare_magnitude for GCD algorithms.
    [[nodiscard]] static int compare_magnitude_pub(const BigInt& a, const BigInt& b) noexcept;

    // F1.1: public factory from little-endian uint32 limbs (for Toom-3, Montgomery).
    // Equivalent to private from_parts(limbs, false) [non-negative].
    [[nodiscard]] static BigInt from_limbs_le(std::vector<std::uint32_t> limbs) noexcept;

    // F7.0-A3.6: hard cap on the number of 32-bit limbs a BigInt result is
    // allowed to occupy. Enforced (pre-allocation) in multiply_magnitude.
    // 0 = unlimited (default). Thread-local — each worker thread carries
    // its own cap. When the bound would be exceeded, the operation
    // produces BigInt(0) AND sets the "exhausted" flag the caller can
    // poll via budget_exhausted().
    //
    // Use case: prevent factorial(10000), pow(2, 10^9), or pathological
    // CRT accumulation from saturating RAM before the OS OOM-killer fires.
    static void set_max_limbs(std::size_t n) noexcept;
    [[nodiscard]] static std::size_t max_limbs() noexcept;
    [[nodiscard]] static bool budget_exhausted() noexcept;
    static void clear_budget_exhausted() noexcept;

    // Public access to the underlying division back-ends — required by the
    // production test suite for cross-checking algorithms against each other.
    // F1.1 Knuth Algorithm D (TAOCP Vol 2 §4.3.1).
    [[nodiscard]] static std::pair<BigInt, BigInt> divide_knuth_d(const BigInt& u, const BigInt& v);
    // F1.2 Burnikel-Ziegler recursive division (HPP-023 closure).
    // Reference: Burnikel-Ziegler 1998, Brent-Zimmermann §1.4.3.
    [[nodiscard]] static std::pair<BigInt, BigInt> divide_burnikel_ziegler(const BigInt& u, const BigInt& v);

private:
    [[nodiscard]] Result<void> assign_decimal_checked(std::string decimal);
    void multiply_by_small(std::uint32_t value) noexcept;
    void add_small(std::uint32_t value) noexcept;
    [[nodiscard]] std::uint32_t divide_by_small(std::uint32_t value) noexcept;
    void normalize() noexcept;
    void invalidate_decimal_cache() noexcept;
    void rebuild_decimal_cache() const;

    [[nodiscard]] static BigInt from_parts(std::vector<std::uint32_t> limbs, bool negative) noexcept;
    [[nodiscard]] static int compare_magnitude(const BigInt& lhs, const BigInt& rhs) noexcept;
    [[nodiscard]] static BigInt add_magnitude(const BigInt& lhs, const BigInt& rhs);
    [[nodiscard]] static BigInt subtract_magnitude(const BigInt& lhs, const BigInt& rhs);
    [[nodiscard]] static BigInt multiply_magnitude(const BigInt& lhs, const BigInt& rhs);
    // F1.1 Toom-3: used by multiply_magnitude for n in [kToom3Threshold, kFFTThreshold).
    // Reference: Brent-Zimmermann "Modern Computer Arithmetic" §1.3.3.
    [[nodiscard]] static BigInt multiply_magnitude_toom3(const BigInt& lhs, const BigInt& rhs);
    // Internal Burnikel-Ziegler primitives (declared here so the implementation
    // can call divide_knuth_d as base case via member-access).
    [[nodiscard]] static std::pair<BigInt, BigInt> bz_div_2by1(const BigInt& A, const BigInt& B, std::size_t n);
    [[nodiscard]] static std::pair<BigInt, BigInt> bz_div_3by2(const BigInt& A, const BigInt& B, std::size_t n);
    [[nodiscard]] static std::pair<BigInt, BigInt> divide_magnitude(const BigInt& dividend, const BigInt& divisor);
    [[nodiscard]] static std::pair<BigInt, BigInt> divide_with_remainder(const BigInt& dividend, const BigInt& divisor);

    std::vector<std::uint32_t> limbs_;
    bool negative_{false};
    mutable std::string decimal_cache_{"0"};
    mutable bool decimal_cache_valid_{true};
};

[[nodiscard]] BigInt gcd(BigInt lhs, BigInt rhs);
[[nodiscard]] Result<std::pair<BigInt, BigInt>> checked_divide_with_remainder(
    const BigInt& dividend,
    const BigInt& divisor);
[[nodiscard]] Result<BigInt> checked_divide(const BigInt& lhs, const BigInt& rhs);
[[nodiscard]] Result<BigInt> checked_mod(const BigInt& lhs, const BigInt& rhs);

// ── F1.1 production-grade algorithms ────────────────────────────────────────
//
// Binary GCD (Stein algorithm).
// Reference: Knuth TAOCP Vol 2 §4.5.2 Algorithm B.
// Faster than Euclidean for large integers due to halving instead of modulo.
// For n < 256 limbs: identical to gcd(). For n >= 256: ~20-30% faster.
[[nodiscard]] BigInt binary_gcd(BigInt lhs, BigInt rhs);

// Lehmer GCD for large integers.
// Reference: Knuth TAOCP Vol 2 §4.5.2 Algorithm L.
// Uses only the leading limbs as surrogates to avoid full modulo.
// Asymptotically O(n²) but with a smaller constant than Euclidean.
// Dispatches: n<256 limbs → Euclidean, n>=256 → Lehmer.
[[nodiscard]] BigInt lehmer_gcd(BigInt a, BigInt b);

// Montgomery modular exponentiation.
// Reference: Koc, Acar, Kaliski "Analyzing and Comparing Montgomery
// Multiplication Algorithms" (IEEE Micro 1996).
// Computes base^exp mod modulus without repeated division.
// For odd modulus only (asserts otherwise). Falls back to binary
// exponentiation for even moduli.
// Pre-condition: modulus > 0, exp >= 0.
[[nodiscard]] Result<BigInt> montgomery_modexp(
    const BigInt& base,
    const BigInt& exp,
    const BigInt& modulus);

// Pollard p-1 factorization.
// Reference: Pollard, J.M. "Theorems on factorization and primality testing"
// (Proc. Cambridge Phil. Soc. 76, 1974).
// Finds a factor of n when n has a prime factor p such that p-1 is B-smooth
// (all prime factors of p-1 are ≤ B). B is derived from ctx_bound parameter
// (default: B = 10^6, covers most "smooth" composites up to ~10^20).
// Returns Unimplemented if no factor found within budget.
[[nodiscard]] Result<BigInt> pollard_p1_factor(const BigInt& n, std::uint64_t bound = 1000000ULL);

}  // namespace cas
