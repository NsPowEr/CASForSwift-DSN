// bigint_factor_pollard.cpp — Pollard p-1 factorization for BigInt.
//
// POLLARD p-1 ALGORITHM:
//   Reference: Pollard, J.M. "Theorems on factorization and primality testing"
//   (Proceedings of the Cambridge Philosophical Society, Vol. 76, 1974, pp. 521-528).
//
//   Idea: if n has a factor p such that p-1 is B-smooth (all prime power factors
//   of p-1 are ≤ B), then p-1 | lcm(1..B) = M(B).
//   By Fermat's little theorem, a^{M(B)} ≡ 1 (mod p) for any a coprime to p.
//   Therefore gcd(a^{M(B)} - 1, n) is divisible by p.
//
//   Stage 1 (basic): compute M_1 = ∏ (p^e for p prime, p^e ≤ B1).
//   Compute g = gcd(a^{M_1} - 1, n). If 1 < g < n, g is a non-trivial factor.
//
//   Stage 2 (optional, not implemented): handles primes p-1 = q * m where
//   q is a single large prime > B1 but ≤ B2. Marked as future work.
//
//   Bound B = 10^6 (default): finds factors p where p-1 is 10^6-smooth.
//   This covers many cryptographically weak composites and numbers arising
//   in CAS integer arithmetic (cyclotomic coefficients, factorial sums, etc.).
//
//   Boundary handling: for very smooth p-1 (p-1 = 2^k * 3^j * ...), even
//   B = 1000 finds the factor. For harder composites, increase B to 10^7.
//
// COMPLEXITY:
//   Time: O(B * log B * log n) using fast sieving + repeated squaring.
//   Space: O(π(B)) for the prime sieve (π(10^6) = 78498 primes).
//
// INTEGRATION:
//   pollard_p1_factor(n, bound) returns a Result<BigInt> factor of n.
//   Returns Unimplemented if no factor is found within the bound.
//   The numtheory::factor_integer() function calls this as a fallback after
//   trial division and before Pollard rho.

#include "cas/bigint.hpp"
#include "cas/error.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace cas {

namespace {

[[nodiscard]] CASError make_error_p1(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

// Simple Sieve of Eratosthenes up to `limit`.
// Returns vector of primes p ≤ limit.
// Used to enumerate the prime base for stage 1.
[[nodiscard]] std::vector<std::uint64_t> sieve_primes(std::uint64_t limit) {
    if (limit < 2U) return {};
    // Cap sieve at practical limit to avoid OOM.
    const std::uint64_t cap = std::min(limit, static_cast<std::uint64_t>(10'000'000ULL));
    std::vector<bool> is_composite(static_cast<std::size_t>(cap + 1U), false);
    std::vector<std::uint64_t> primes;
    primes.reserve(static_cast<std::size_t>(cap / 10U + 100U));
    for (std::uint64_t i = 2U; i <= cap; ++i) {
        if (!is_composite[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
            if (i <= cap / i) {
                for (std::uint64_t j = i * i; j <= cap; j += i) {
                    is_composite[static_cast<std::size_t>(j)] = true;
                }
            }
        }
    }
    return primes;
}

// Modular multiplication of BigInt: (a * b) % m.
// Used in the repeated-squaring phase.
[[nodiscard]] BigInt mulmod(const BigInt& a, const BigInt& b, const BigInt& m) {
    return (a * b) % m;
}

// Modular exponentiation: base^exp mod modulus (simple binary method).
// Used for the p-1 accumulation step.
[[nodiscard]] BigInt powmod(BigInt base, BigInt exp, const BigInt& modulus) {
    BigInt result(1);
    base = base % modulus;
    if (base.is_negative()) base = base + modulus;
    const BigInt zero(0);
    const BigInt two(2);
    const BigInt one(1);
    while (exp > zero) {
        if ((exp % two) == one) {
            result = mulmod(result, base, modulus);
        }
        exp = exp.shift_right_bits(1);
        base = mulmod(base, base, modulus);
    }
    return result;
}

} // namespace

// Pollard p-1 factorization.
// Reference: Pollard 1974, §2 "Algorithm P".
//
// bound: the B-smoothness bound for p-1 (default: 10^6).
// base: starting base a (usually 2; multiple bases tried if needed).
// Returns a non-trivial factor of n, or Unimplemented if none found.
Result<BigInt> pollard_p1_factor(const BigInt& n, std::uint64_t bound) {
    const BigInt zero(0);
    const BigInt one(1);
    const BigInt two(2);

    if (n <= one) {
        return fail<BigInt>(make_error_p1(
            CASErrorKind::InvalidArgument,
            "pollard_p1_factor: n must be > 1"));
    }
    if ((n % two).is_zero()) {
        return ok(two);  // trivial factor
    }

    // Sieve primes up to bound.
    const std::vector<std::uint64_t> primes = sieve_primes(bound);

    // Try several bases: a = 2, 3, 5, 7.
    // Pollard 1974 notes that base 2 works for most composites; multiple
    // bases improve coverage without significant overhead.
    static constexpr std::int64_t kBases[] = {2, 3, 5, 7, 11, 13};
    for (std::int64_t base_val : kBases) {
        const BigInt a_init(base_val);
        if ((n % a_init).is_zero()) {
            // Trivial factor found.
            return ok(a_init);
        }

        BigInt a = a_init;

        // Stage 1: accumulate a ← a^{M(B)} mod n.
        // M(B) = ∏ p^{floor(log_p B)} for all primes p ≤ B.
        for (const std::uint64_t prime : primes) {
            if (a.is_zero()) break;
            // Compute p^e where p^e ≤ bound.
            std::uint64_t q = prime;
            while (q <= bound / prime) {
                q *= prime;
            }
            a = powmod(a, BigInt(static_cast<long long>(q)), n);
        }

        // Check gcd(a - 1, n).
        BigInt g = gcd(a - one, n);
        if (g > one && g < n) {
            return ok(g);
        }
        if (g == n) {
            // Overshot: n divides a^M - 1. Try with smaller bound (halve it).
            // For simplicity, move to next base.
            continue;
        }
        // g == 1: no factor found with this base, try next.
    }

    return fail<BigInt>(make_error_p1(
        CASErrorKind::Unimplemented,
        "pollard_p1_factor: no factor found within bound " +
        std::to_string(bound) + "; increase bound or use ECM "
        "(Aperta permanente CAS-F1.1-ECM in CAS_TASKS.md)"));
}

}  // namespace cas
