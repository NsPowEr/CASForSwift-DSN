// CAS-F1.6-NEW — Gaussian prime factorization of Z>0.
//
// Algorithm reference: Henri Cohen "A Course in Computational Algebraic
// Number Theory" (Springer GTM 138) §4.2 (Gaussian integers).
//
// Classification of rational primes p in Z[i]:
//   p = 2         → 2 = -i·(1+i)²  [ramified;  N(1+i) = 2]
//   p ≡ 1 mod 4  → p = π·π̄  [split;  π = a+bi with a²+b²=p, found via
//                               Hermite-Serret / Tonelli-Shanks + Euclidean]
//   p ≡ 3 mod 4  → p stays prime in Z[i] [inert]
//
// For a general positive integer n with rational factorization ∏ pᵢ^eᵢ:
//   lift each pᵢ^eᵢ to Z[i] according to the classification above.
//
// No hardcoded prime tables — all primality tests via numtheory::is_prime.

#include "cas/gaussian_int.hpp"
#include "cas/numtheory.hpp"
#include "cas/error.hpp"

#include <algorithm>
#include <cassert>

namespace cas {

namespace {

// -------------------------------------------------------------------------
// Hermite-Serret: find a,b ≥ 0 with a² + b² = p for prime p ≡ 1 mod 4.
// Method: find x with x² ≡ -1 (mod p) via Wilson / Euler criterion, then
// run Euclidean algorithm on (p, x+i) in Z[i] until norm(remainder) ≤ p.
// Cohen GTM 138, Proposition 4.2.5.
// -------------------------------------------------------------------------
[[nodiscard]] Result<GaussianInt> split_prime_factor(const BigInt& p) {
    // Step 1: find x such that x² ≡ -1 (mod p).
    // By Fermat's little theorem x = ((p-1)/2)! mod p (Wilson's theorem).
    // Equivalently iterate small bases: for base b, x = b^((p-1)/4) mod p
    // if b is a quadratic non-residue mod p. We try b = 2, 3, 5, ...
    // This is deterministic (not hardcoded) because we iterate bases until
    // we find a QNR — guaranteed to exist with density 1/2.
    BigInt pm1 = p - BigInt(1);
    BigInt exp_qnr = pm1 / BigInt(4);  // (p-1)/4

    BigInt x(0);
    {
        BigInt base(2);
        while (true) {
            auto pm = numtheory::power_mod(base, exp_qnr, p);
            if (!pm.is_ok()) { base = base + BigInt(1); continue; }
            BigInt candidate = pm.value();
            // Check candidate² ≡ -1 mod p
            auto sq = numtheory::power_mod(candidate, BigInt(2), p);
            if (!sq.is_ok()) { base = base + BigInt(1); continue; }
            BigInt sq_val = sq.value();
            // sq_val ≡ -1 (mod p) means sq_val + 1 ≡ 0 (mod p)
            BigInt sum = sq_val + BigInt(1);
            if ((sum % p).is_zero()) {
                x = candidate;
                break;
            }
            base = base + BigInt(1);
        }
    }

    // Step 2: Euclidean GCD of p and (x + i) in Z[i], stopping when norm ≤ p.
    // gcd(p, x+i) in Z[i] will be a Gaussian prime π with N(π) = p.
    // We run the Euclidean algorithm with α=p (a pure real Gaussian int) and
    // β = GaussianInt(x, 1).
    GaussianInt alpha(p, BigInt(0));
    GaussianInt beta(x, BigInt(1));

    // Iterate until norm(alpha) ≤ p: at that point alpha is the factor.
    while (alpha.norm() > p) {
        auto dm = gaussian_divmod(alpha, beta);
        alpha = beta;
        beta  = dm.remainder;
    }

    // Canonicalize: make real > 0 and imag ≥ 0.
    BigInt re = alpha.real().abs();
    BigInt im = alpha.imag().abs();
    // Invariant (Cohen GTM 138 Proposition 4.2.5): the Euclidean GCD
    // algorithm on Z[i] terminates with N(alpha) = p, so re²+im² must
    // equal p. If this fails the Hermite-Serret step produced an
    // incorrect intermediate — swap (re, im) as the last-resort recovery.
    //
    // Note: all four unit-rotations (±re, ±im) have the same norm, so
    // swapping cannot fix a genuine norm mismatch; a swap only helps if
    // the abs() above discarded sign information that affected which of
    // re/im ended up as which component. A true norm != p is an algorithm
    // bug in split_prime_factor, not a caller-input error.
    if (re * re + im * im != p) {
        return fail<GaussianInt>(CASError{
            .kind    = CASErrorKind::InternalError,
            .message = "HPP-018 Hermite-Serret invariant violated: "
                       "N(alpha) != p after Z[i] Euclidean GCD",
            .hint    = std::nullopt,
        });
    }
    if (re < im) {
        std::swap(re, im);
    }
    return ok(GaussianInt(re, im));
}

// -------------------------------------------------------------------------
// Lift a single rational prime power p^e to a list of Gaussian factors.
// -------------------------------------------------------------------------
[[nodiscard]] Result<std::vector<GaussianFactor>> lift_prime_power(
    const BigInt& p, unsigned int e)
{
    std::vector<GaussianFactor> result;

    // p = 2: ramified. 2 = -i · (1+i)². The Gaussian prime is (1+i).
    // Each factor of 2 in the rational integer contributes (1+i)² in Z[i].
    // So p^e = (-i)^e · (1+i)^{2e}.  Unit factor handled by caller.
    if (p == BigInt(2)) {
        // Gaussian prime: 1+i (norm 2, canonical: real>0, imag>0 → already fine)
        result.push_back(GaussianFactor{GaussianInt(BigInt(1), BigInt(1)), 2U * e});
        return ok(std::move(result));
    }

    // Determine p mod 4 via p % 4.
    BigInt mod4 = p % BigInt(4);

    if (mod4 == BigInt(3)) {
        // p ≡ 3 mod 4: inert, p stays prime in Z[i].
        result.push_back(GaussianFactor{GaussianInt(p, BigInt(0)), e});
        return ok(std::move(result));
    }

    // p ≡ 1 mod 4: split. p = π · π̄ with N(π) = p.
    assert(mod4 == BigInt(1));
    auto pi_result = split_prime_factor(p);
    if (!pi_result.is_ok()) {
        return fail<std::vector<GaussianFactor>>(pi_result.error());
    }
    GaussianInt pi = pi_result.value();
    GaussianInt pi_conj = pi.conjugate();

    // Canonicalize pi_conj so its real > 0.
    if (pi_conj.real().is_negative()) {
        pi_conj = -pi_conj;
    }

    // Are π and π̄ distinct (i.e. not unit multiples of each other)?
    // They are distinct iff im(π) ≠ 0.
    bool distinct = !pi.imag().is_zero();

    if (distinct) {
        result.push_back(GaussianFactor{pi,      e});
        result.push_back(GaussianFactor{pi_conj, e});
    } else {
        // Degenerate: pi is real (norm = pi² → shouldn't happen for p prime ≡ 1 mod 4).
        result.push_back(GaussianFactor{pi, e});
    }
    return ok(std::move(result));
}

// -------------------------------------------------------------------------
// Canonical form: (real > 0) or (real = 0 and imag > 0).
// -------------------------------------------------------------------------
[[nodiscard]] GaussianInt canonicalize(GaussianInt g) {
    if (g.real().is_negative() ||
        (g.real().is_zero() && g.imag().is_negative())) {
        return -g;
    }
    return g;
}

}  // namespace

// -------------------------------------------------------------------------
// Public API: factor_gaussian(n), n ∈ Z>0.
// -------------------------------------------------------------------------
Result<GaussianFactorization> factor_gaussian(const BigInt& n) {
    if (n.is_zero() || n.is_negative()) {
        return fail<GaussianFactorization>(CASError{
            .kind    = CASErrorKind::InvalidArgument,
            .message = "factor_gaussian requires a positive integer",
            .hint    = std::nullopt,
        });
    }

    // Rational factorization first.
    auto rat_fact = numtheory::factor_integer(n);
    if (!rat_fact.is_ok()) {
        return fail<GaussianFactorization>(rat_fact.error());
    }

    GaussianFactorization result;
    result.unit = GaussianInt(BigInt(1), BigInt(0));  // start with unit 1

    // For p = 2: each 2^e contributes (-i)^e as a unit adjustment and
    // (1+i)^{2e} as the prime power.
    // For p ≡ 3 mod 4: each p^e contributes p^e (inert prime) directly.
    // For p ≡ 1 mod 4: each p^e contributes π^e · π̄^e (split).
    for (auto& [prime, exp] : rat_fact.value().prime_factors) {
        if (prime == BigInt(2)) {
            // The unit adjustment: 2 = -i·(1+i)², so 2^e = (-i)^e · (1+i)^{2e}.
            // (-i)^1 = -i, (-i)^2 = -1, (-i)^3 = i, (-i)^4 = 1 (period 4).
            static const GaussianInt unit_powers[4] = {
                GaussianInt(BigInt(1), BigInt(0)),   // (-i)^0 = 1
                GaussianInt(BigInt(0), BigInt(-1)),  // (-i)^1 = -i
                GaussianInt(BigInt(-1), BigInt(0)),  // (-i)^2 = -1
                GaussianInt(BigInt(0), BigInt(1)),   // (-i)^3 = i
            };
            result.unit = result.unit * unit_powers[exp % 4U];
        }
        auto lifted = lift_prime_power(prime, exp);
        if (!lifted.is_ok()) {
            return fail<GaussianFactorization>(lifted.error());
        }
        for (auto& f : lifted.value()) {
            f.prime = canonicalize(f.prime);
            result.factors.push_back(std::move(f));
        }
    }

    // Sort factors by norm ascending for deterministic output.
    std::sort(result.factors.begin(), result.factors.end(),
        [](const GaussianFactor& a, const GaussianFactor& b) {
            return a.prime.norm() < b.prime.norm();
        });

    return ok(std::move(result));
}

}  // namespace cas
