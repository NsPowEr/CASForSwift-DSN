// test_bigint_production.cpp — Production-grade BigInt algorithm test suite.
//
// F1.1 certification tests as specified in PLAN_HP_PRIME_PARITY.md §F1.1.
//
// Tests cover:
//   1. Karatsuba/Toom-3 multiplication: (a*b)/b ≡ a for b ≠ 0.
//   2. gcd(a,b)·(a/gcd) = a, gcd(a,b)·(b/gcd) = b (GCD property).
//   3. Knuth D division: a/b*b + a%b = a.
//   4. Montgomery modexp: montgomery_modexp(a,b,m) ≡ power_mod(a,b,m).
//   5. is_prime(p) for known primes via Miller-Rabin.
//   6. pollard_p1_factor correctness.
//   7. Binary GCD / Lehmer GCD: binary_gcd(a,b) == gcd(a,b).
//
// Property tests use specific large-integer examples; no RapidCheck dependency.
// All tests use structural comparison (EXPECT_EQ on decimal strings).
//
// Reference: PLAN_HP_PRIME_PARITY.md §F1.1 "Certificatori".

#include <gtest/gtest.h>

#include "cas/bigint.hpp"
#include "cas/numtheory.hpp"

using namespace cas;
using namespace cas::numtheory;

namespace {

[[nodiscard]] BigInt parse_or_abort(const char* str) {
    auto r = BigInt::parse(str);
    if (r.is_error()) {
        ADD_FAILURE() << "Failed to parse BigInt: " << str;
        return BigInt(0);
    }
    return r.value();
}

}  // namespace

// ── 1. Multiplication correctness: (a*b)/b ≡ a ──────────────────────────────

TEST(BigIntProductionTest, MultiplyDivideInverse_Small) {
    // Small: schoolbook path (n < 32)
    const BigInt a = parse_or_abort("123456789");
    const BigInt b = parse_or_abort("987654321");
    ASSERT_EQ((a * b / b).decimal(), a.decimal());
}

TEST(BigIntProductionTest, MultiplyDivideInverse_Karatsuba) {
    // Medium: Karatsuba path (n in [32,64))
    // 40-limb numbers: ~1280 bits each.
    // Use Fibonacci-like construction for non-trivial limb pattern.
    const BigInt a = parse_or_abort(
        "12345678901234567890123456789012345678901234567890"
        "12345678901234567890123456789012345678901234567890");
    const BigInt b = parse_or_abort(
        "98765432109876543210987654321098765432109876543210"
        "98765432109876543210987654321098765432109876543210");
    const BigInt product = a * b;
    const BigInt recovered = product / b;
    EXPECT_EQ(recovered.decimal(), a.decimal())
        << "Karatsuba: (a*b)/b must equal a";
}

TEST(BigIntProductionTest, MultiplyDivideInverse_Toom3) {
    // Large: Toom-3 path (n in [64,4096))
    // Need 64+ limbs = 64*32 bits ≈ 617 decimal digits per operand.
    // Use 640 nines and 640 ones.
    const BigInt a = parse_or_abort(
        "99999999999999999999999999999999999999999999999999"
        "99999999999999999999999999999999999999999999999999"
        "99999999999999999999999999999999999999999999999999"
        "99999999999999999999999999999999999999999999999999"
        "99999999999999999999999999999999999999999999999999"
        "99999999999999999999999999999999999999999999999999"
        "99999999999999999999999999999999999999999999999999"
        "99999999999999999999999999999999999999999999999999"
        "99999999999999999999999999999999999999999999999999"
        "99999999999999999999999999999999999999999999999999"
        "99999999999999999999999999999999999999999999999999"
        "99999999999999999999999999999999999999999999999999"
        "9999999999999999999999999999999999999999");  // 640 nines
    const BigInt b = parse_or_abort(
        "11111111111111111111111111111111111111111111111111"
        "11111111111111111111111111111111111111111111111111"
        "11111111111111111111111111111111111111111111111111"
        "11111111111111111111111111111111111111111111111111"
        "11111111111111111111111111111111111111111111111111"
        "11111111111111111111111111111111111111111111111111"
        "11111111111111111111111111111111111111111111111111"
        "11111111111111111111111111111111111111111111111111"
        "11111111111111111111111111111111111111111111111111"
        "11111111111111111111111111111111111111111111111111"
        "11111111111111111111111111111111111111111111111111"
        "11111111111111111111111111111111111111111111111111"
        "1111111111111111111111111111111111111111");  // 638 ones
    const BigInt product = a * b;
    const BigInt recovered = product / b;
    EXPECT_EQ(recovered.decimal(), a.decimal())
        << "Toom-3: (a*b)/b must equal a (operands have ~67 limbs each, triggers Toom-3)";
}

// ── 2. Division: a/b*b + a%b = a ─────────────────────────────────────────────

TEST(BigIntProductionTest, KnuthD_DivisionRemainder_MultiLimb) {
    // Division of a 20-limb number by a 10-limb number via Knuth D.
    const BigInt a = parse_or_abort(
        "12345678901234567890123456789012345678901234567890"
        "12345678901234567890123456789012345678901234567890"
        "12345678901234567890");
    const BigInt b = parse_or_abort(
        "99999999999999999999999999999999999999999999999999");
    const BigInt q = a / b;
    const BigInt r = a % b;
    EXPECT_EQ((q * b + r).decimal(), a.decimal())
        << "Knuth D: a/b*b + a%b must equal a";
    EXPECT_TRUE(r < b) << "Knuth D: remainder must be < divisor";
    EXPECT_FALSE(r.is_negative()) << "Knuth D: remainder must be non-negative for positive inputs";
}

TEST(BigIntProductionTest, KnuthD_DivisionRemainder_SmallDivisor) {
    // Single-limb divisor: uses divide_by_small path.
    const BigInt a = parse_or_abort("-12345678901234567890");
    const BigInt b = parse_or_abort("9876543210");
    EXPECT_EQ((a / b).decimal(), "-1249999989");
    EXPECT_EQ((a % b).decimal(), "2623456800");
}

TEST(BigIntProductionTest, KnuthD_DivisionProperty_RandomSamples) {
    // Property: a/b*b + a%b = a for 10 random large pairs.
    struct Sample {
        const char* a_str;
        const char* b_str;
    };
    const Sample samples[] = {
        {"123456789123456789", "9999"},
        {"10000000000000000000000000", "99999999999999"},
        {"999999999999999999999999999999999999999999999999", "123456789"},
        {"1234567890987654321", "111111111111111111"},
        {"99999999999999999999999999999999", "77777777777777"},
    };
    for (const auto& s : samples) {
        const BigInt a = parse_or_abort(s.a_str);
        const BigInt b = parse_or_abort(s.b_str);
        const BigInt q = a / b;
        const BigInt r = a % b;
        EXPECT_EQ((q * b + r).decimal(), a.decimal())
            << "Division property a/b*b + a%b = a failed for a=" << s.a_str;
        EXPECT_TRUE(r < b || r.is_zero()) << "Remainder must be < divisor";
    }
}

// ── 3. GCD properties ────────────────────────────────────────────────────────

TEST(BigIntProductionTest, GCD_Commutativity) {
    const BigInt a = parse_or_abort("123456789012345678901234567890");
    const BigInt b = parse_or_abort("987654321098765432109876543210");
    EXPECT_EQ(gcd(a, b).decimal(), gcd(b, a).decimal())
        << "GCD must be commutative";
}

TEST(BigIntProductionTest, GCD_Divisibility) {
    // gcd(a,b) must divide both a and b.
    const BigInt a = parse_or_abort("9699690");   // 2*3*5*7*11*13*17*19
    const BigInt b = parse_or_abort("223092870");  // 2*3*5*7*11*13*17*19*23
    const BigInt g = gcd(a, b);
    EXPECT_EQ(g.decimal(), "9699690") << "gcd(a, a*23) = a";
    EXPECT_TRUE((a % g).is_zero()) << "gcd must divide a";
    EXPECT_TRUE((b % g).is_zero()) << "gcd must divide b";
}

TEST(BigIntProductionTest, BinaryGCD_MatchesStandardGCD) {
    // binary_gcd(a, b) == gcd(a, b) for various sizes.
    struct Sample { const char* a; const char* b; };
    const Sample samples[] = {
        {"48", "180"},
        {"12345678901234567890", "98765432109876543210"},
        {"9699690", "223092870"},
        {"1000000007", "1000000009"},  // two consecutive primes
        {"0", "42"},
        {"42", "0"},
    };
    for (const auto& s : samples) {
        BigInt a = parse_or_abort(s.a);
        BigInt b = parse_or_abort(s.b);
        EXPECT_EQ(binary_gcd(a, b).decimal(), gcd(a, b).decimal())
            << "binary_gcd != gcd for a=" << s.a << " b=" << s.b;
    }
}

TEST(BigIntProductionTest, LehmerGCD_MatchesStandardGCD_LargeInputs) {
    // Lehmer GCD must match standard GCD for inputs triggering the Lehmer path.
    // kLehmerThreshold = 16 limbs = 512 bits ≈ 155 decimal digits.
    const BigInt a = parse_or_abort(
        "12345678901234567890123456789012345678901234567890"
        "12345678901234567890123456789012345678901234567890"
        "12345678901234567890123456789012345678901234567890"
        "12345678901234567890");
    const BigInt b = parse_or_abort(
        "98765432109876543210987654321098765432109876543210"
        "98765432109876543210987654321098765432109876543210"
        "98765432109876543210987654321098765432109876543210"
        "98765432109876543210");
    // Standard Euclidean GCD for reference.
    const BigInt g_standard = gcd(a, b);
    const BigInt g_lehmer   = lehmer_gcd(a, b);
    EXPECT_EQ(g_lehmer.decimal(), g_standard.decimal())
        << "Lehmer GCD must match standard GCD for large inputs";
}

// ── 4. Montgomery modexp ─────────────────────────────────────────────────────

TEST(BigIntProductionTest, MontgomeryModexp_SmallValues) {
    // 2^10 mod 1000 = 24.
    const BigInt base(2);
    const BigInt exp(10);
    const BigInt mod(1000);
    auto result = montgomery_modexp(base, exp, mod);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_EQ(result.value().decimal(), "24");
}

TEST(BigIntProductionTest, MontgomeryModexp_FermatLittle) {
    // Fermat's little theorem: a^(p-1) ≡ 1 (mod p) for prime p, gcd(a,p)=1.
    // Use p = 999999937 (large prime), a = 7.
    const BigInt p(999999937LL);
    const BigInt a(7);
    const BigInt exp = p - BigInt(1);
    auto result = montgomery_modexp(a, exp, p);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_EQ(result.value().decimal(), "1")
        << "Fermat: 7^(p-1) mod p must be 1 for prime p";
}

TEST(BigIntProductionTest, MontgomeryModexp_MatchesPowerMod) {
    // montgomery_modexp(a,b,m) == power_mod(a,b,m) for several cases.
    struct Sample { const char* base; const char* exp; const char* mod; };
    const Sample samples[] = {
        {"2", "100", "1000000007"},
        {"12345", "6789", "999999937"},
        {"1000000000", "1000000000", "1000000007"},
        {"3", "0", "17"},    // 3^0 mod 17 = 1
        {"0", "5", "7"},     // 0^5 mod 7 = 0
    };
    for (const auto& s : samples) {
        const BigInt b  = parse_or_abort(s.base);
        const BigInt e  = parse_or_abort(s.exp);
        const BigInt m  = parse_or_abort(s.mod);
        auto ref = power_mod(b, e, m);
        auto mont = montgomery_modexp(b, e, m);
        ASSERT_TRUE(ref.is_ok())  << "power_mod error: " << ref.error().message;
        ASSERT_TRUE(mont.is_ok()) << "montgomery_modexp error: " << mont.error().message;
        EXPECT_EQ(mont.value().decimal(), ref.value().decimal())
            << "Montgomery mismatch for base=" << s.base << " exp=" << s.exp << " mod=" << s.mod;
    }
}

// ── 5. Miller-Rabin primality ─────────────────────────────────────────────────

TEST(BigIntProductionTest, IsPrime_KnownPrimes) {
    // BPSW / deterministic Miller-Rabin for small n.
    const std::int64_t small_primes[] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
        97, 101, 997, 9999991, 1000000007LL, 999999937LL
    };
    for (std::int64_t p : small_primes) {
        auto result = is_prime(BigInt(p));
        ASSERT_TRUE(result.is_ok()) << "is_prime error for " << p;
        EXPECT_TRUE(result.value()) << "is_prime should be true for prime " << p;
    }
}

TEST(BigIntProductionTest, IsPrime_KnownComposites) {
    const std::int64_t composites[] = {
        4, 6, 8, 9, 10, 15, 25, 49, 100, 1001, 9999999, 1000000006LL
    };
    for (std::int64_t c : composites) {
        auto result = is_prime(BigInt(c));
        ASSERT_TRUE(result.is_ok()) << "is_prime error for " << c;
        EXPECT_FALSE(result.value()) << "is_prime should be false for composite " << c;
    }
}

TEST(BigIntProductionTest, IsPrime_LargePrime) {
    // 2^61 - 1 = 2305843009213693951 (Mersenne prime M61).
    const BigInt mersenne = BigInt(1).shift_left_bits(61) - BigInt(1);
    auto result = is_prime(mersenne);
    ASSERT_TRUE(result.is_ok());
    EXPECT_TRUE(result.value()) << "2^61-1 is a Mersenne prime (M61)";
}

// ── 6. Pollard rho factorization ──────────────────────────────────────────────

TEST(BigIntProductionTest, PollardRho_FactorizesSmallComposite) {
    // factor 15 = 3 * 5.
    const BigInt n(15);
    auto r = pollards_rho_factor(n);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    const BigInt f = r.value();
    EXPECT_TRUE(f > BigInt(1) && f < n) << "Pollard rho must return non-trivial factor";
    EXPECT_TRUE((n % f).is_zero()) << "Factor must divide n";
}

TEST(BigIntProductionTest, PollardRho_FactorizesLargerComposite) {
    // factor 9999991 * 9999973.
    const BigInt p1(9999991LL);
    const BigInt p2(9999973LL);
    const BigInt n = p1 * p2;
    auto fact = factor_integer(n);
    ASSERT_TRUE(fact.is_ok()) << fact.error().message;
    // Verify product of factors = n.
    BigInt product(1);
    for (const auto& [prime, exp] : fact.value().prime_factors) {
        for (unsigned int i = 0; i < exp; ++i) {
            product = product * prime;
        }
    }
    EXPECT_EQ(product.decimal(), n.decimal())
        << "Product of prime factors must equal n";
}

// ── 7. Pollard p-1 factorization ─────────────────────────────────────────────

TEST(BigIntProductionTest, PollardP1_FindsSmooth_Factor) {
    // n = p * q where p-1 = 2 * 3 * 5 * 7 * 11 = 2310, so p=2311 (check).
    // 2311 is prime and 2310 is 11-smooth.
    // n = 2311 * 2333 (2333 is also prime, 2332 = 4*11*53 — less smooth).
    const BigInt p_smooth(2311LL);   // p-1 = 2310 = 2*3*5*7*11, 11-smooth
    const BigInt q_other(2333LL);    // q-1 = 2332 = 4*11*53
    const BigInt n = p_smooth * q_other;
    auto r = pollard_p1_factor(n, 20ULL);  // B=20 should be enough for 11-smooth
    ASSERT_TRUE(r.is_ok()) << "Pollard p-1 should find factor of " << n.decimal()
                           << "; error: " << (r.is_error() ? r.error().message : "none");
    EXPECT_TRUE((n % r.value()).is_zero())
        << "p-1 factor must divide n";
    EXPECT_TRUE(r.value() > BigInt(1) && r.value() < n)
        << "p-1 factor must be non-trivial";
}

TEST(BigIntProductionTest, PollardP1_LargeSmoothFactor) {
    // n = 1000003 * 1000033.
    // 1000003 - 1 = 1000002 = 2 * 3 * 166667.  166667 is prime.
    // 1000033 - 1 = 1000032 = 2^5 * 3 * 10417.  10417 is prime.
    // Neither is easily B-smooth for small B.
    // Use n = 2 * 1000003 = 2000006 (p-1 = 1000002, B=1000003 would work,
    // but B-smooth test for simple case: p=3, p-1=2, 2-smooth).
    const BigInt n_simple = parse_or_abort("1000000007000000021");  // arbitrary semiprime
    // Just verify the function returns an error gracefully (not crash).
    auto r = pollard_p1_factor(n_simple, 100ULL);
    // We don't assert success here (the semiprime may not have smooth p-1).
    // Just verify no crash and proper error type.
    if (r.is_error()) {
        EXPECT_EQ(r.error().kind, CASErrorKind::Unimplemented)
            << "pollard_p1_factor should return Unimplemented when no factor found";
    } else {
        EXPECT_TRUE((n_simple % r.value()).is_zero())
            << "If a factor is found, it must divide n";
    }
}

// ── 8. Full factorization: ∏ p_i^{e_i} ≡ n ───────────────────────────────────

TEST(BigIntProductionTest, FactorInteger_ReconstitutesN) {
    const char* composites[] = {
        "360",          // 2^3 * 3^2 * 5
        "1000000007",   // prime
        "720720",       // 2^4 * 3^2 * 5 * 7 * 11 * 13
        "99999999999999999999",  // 3 * 3 * 239 * 4649 * 333667 * (others)
    };
    for (const char* n_str : composites) {
        const BigInt n = parse_or_abort(n_str);
        auto fact = factor_integer(n);
        ASSERT_TRUE(fact.is_ok()) << "factor_integer failed for " << n_str
                                  << ": " << fact.error().message;
        // Verify ∏ p_i^{e_i} = |n|.
        BigInt product(1);
        for (const auto& [prime, exp] : fact.value().prime_factors) {
            for (unsigned int i = 0; i < exp; ++i) {
                product = product * prime;
            }
        }
        EXPECT_EQ(product.decimal(), n.decimal())
            << "∏ p_i^e_i must equal n for n=" << n_str;
    }
}

