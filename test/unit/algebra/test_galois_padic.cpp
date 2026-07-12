// A6 Brick 3a — p-adic splitting engine tests (GR(p^k, L), root lifting,
// Frobenius). Everything is validated against mathematical identities
// computed in exact arithmetic — never against transcribed values:
//   • ring axioms and inverses in GF(9);
//   • rational roots lift to their exact integer residues;
//   • ∏(y − r_i) re-expands to f's coefficients inside the ring;
//   • Frobenius cycle type = factor-degree multiset mod p (Dedekind);
//   • precision raise refines roots in place (agreement mod p^k, same
//     Frobenius).

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "../../../src/algebra/galois_padic_internal.hpp"
#include "../../../src/algebra/polynomial_internal.hpp"
#include "cas/bigint.hpp"
#include "cas/error.hpp"

using namespace cas;
using namespace cas::algebra;
using namespace cas::algebra::galois_padic;

namespace {

[[nodiscard]] IntPoly make_poly(std::initializer_list<long long> coeffs) {
    std::vector<BigInt> c;
    c.reserve(coeffs.size());
    for (const long long v : coeffs) c.emplace_back(v);
    IntPoly p(std::move(c));
    p.normalize([](const BigInt& v) { return v.is_zero(); });
    return p;
}

// ∏(y − r_i) expanded over the ring must reproduce f's coefficients: the
// strongest end-to-end identity a splitting can satisfy.
void expect_reconstructs(const PadicSplitting& s) {
    const PadicRing& R = s.ring;
    RingPoly acc{R.one()};
    for (const auto& r : s.roots) {
        const RingPoly lin{R.neg(r), R.one()};
        acc = rp_mul(R, acc, lin);
    }
    ASSERT_EQ(acc.size(), s.f.size());
    for (std::size_t i = 0U; i < acc.size(); ++i) {
        EXPECT_TRUE(R.equal(acc[i], R.from_int(s.f.coefficients()[i])))
            << "coefficient " << i << " mismatch";
    }
}

[[nodiscard]] std::vector<std::size_t> perm_cycle_type(
    const permgrp::Perm& p) {
    std::vector<std::size_t> cyc;
    std::vector<bool> seen(p.size(), false);
    for (std::size_t i = 0U; i < p.size(); ++i) {
        if (seen[i]) continue;
        std::size_t len = 0U;
        std::size_t j = i;
        while (!seen[j]) {
            seen[j] = true;
            j = p[j];
            ++len;
        }
        cyc.push_back(len);
    }
    std::sort(cyc.begin(), cyc.end());
    return cyc;
}

// ── ring arithmetic in GF(9) = F_3[t]/(t² + 1) ──────────────────────────────

TEST(GaloisPadicRing, FieldAxiomsGF9) {
    auto ring = PadicRing::make(BigInt(3), 1U, make_poly({1, 0, 1}));
    ASSERT_TRUE(ring.is_ok());
    const PadicRing& F = ring.value();
    EXPECT_EQ(F.ext_degree(), 2U);
    const RingElem t = F.basis_power(1U);
    // t² = −1 in F_9.
    EXPECT_TRUE(F.equal(F.mul(t, t), F.from_int(BigInt(-1))));
    // Multiplicative order of every unit divides q − 1 = 8.
    EXPECT_TRUE(F.equal(F.pow(t, BigInt(8)), F.one()));
    // Inverse: t · t⁻¹ = 1, and (t+2)(t+2)⁻¹ = 1.
    auto tinv = F.inv(t);
    ASSERT_TRUE(tinv.is_ok());
    EXPECT_TRUE(F.equal(F.mul(t, tinv.value()), F.one()));
    const RingElem u = F.add(t, F.from_int(BigInt(2)));
    auto uinv = F.inv(u);
    ASSERT_TRUE(uinv.is_ok());
    EXPECT_TRUE(F.equal(F.mul(u, uinv.value()), F.one()));
    // Zero is not a unit.
    EXPECT_TRUE(F.inv(F.zero()).is_error());
}

TEST(GaloisPadicRing, InverseNewtonLiftsPrecision) {
    // GR(5^6, 2): inverse must hold at full precision, not just mod 5.
    auto ring = PadicRing::make(BigInt(5), 6U, make_poly({2, 0, 1}));
    ASSERT_TRUE(ring.is_ok());
    const PadicRing& R = ring.value();
    RingElem a = R.add(R.basis_power(1U), R.from_int(BigInt(7)));
    auto ainv = R.inv(a);
    ASSERT_TRUE(ainv.is_ok());
    EXPECT_TRUE(R.equal(R.mul(a, ainv.value()), R.one()));
}

// ── rational roots: exact integer residues ─────────────────────────────────

TEST(GaloisPadicSplit, RationalRootsLiftToIntegers) {
    // f = (x−1)(x−2)(x−3): the engine must recover exactly {1, 2, 3}.
    const IntPoly f = make_poly({-6, 11, -6, 1});
    auto sp = choose_splitting_prime(f, 10U, nullptr);
    ASSERT_TRUE(sp.is_ok());
    EXPECT_EQ(sp.value().ext_degree, 1U);  // splits completely at some p
    auto s = build_padic_splitting(f, sp.value(), 10U, nullptr);
    ASSERT_TRUE(s.is_ok());
    std::vector<BigInt> found;
    for (const auto& r : s.value().roots) {
        auto res = s.value().ring.integer_residue(r);
        ASSERT_TRUE(res.has_value());
        found.push_back(*res);
    }
    std::sort(found.begin(), found.end(),
              [](const BigInt& a, const BigInt& b) { return a < b; });
    ASSERT_EQ(found.size(), 3U);
    EXPECT_TRUE(found[0] == BigInt(1));
    EXPECT_TRUE(found[1] == BigInt(2));
    EXPECT_TRUE(found[2] == BigInt(3));
    // All factors linear ⇒ Frobenius is the identity.
    EXPECT_EQ(perm_cycle_type(s.value().frobenius),
              (std::vector<std::size_t>{1U, 1U, 1U}));
    expect_reconstructs(s.value());
}

// ── inert prime: Frobenius orbit path (Φ = the unique factor) ──────────────

TEST(GaloisPadicSplit, IrreducibleCubicFrobeniusOrbit) {
    // x³ − 2 is irreducible mod 7 (2 is not a cube: 2² = 4 ≠ 1 = 2^((7−1)/3)
    // would be needed). Forcing p = 7 exercises the Frobenius-orbit root
    // path with Φ = f mod p itself.
    const IntPoly f = make_poly({-2, 0, 0, 1});
    auto fact = factor_polynomial_mod_p(f, BigInt(7), nullptr);
    ASSERT_TRUE(fact.is_ok());
    ASSERT_EQ(fact.value().size(), 1U);
    ASSERT_EQ(fact.value()[0].degree(), 3U);
    const SplittingPrime sp{BigInt(7), fact.value(), 3U};
    auto s = build_padic_splitting(f, sp, 8U, nullptr);
    ASSERT_TRUE(s.is_ok());
    const PadicRing& R = s.value().ring;
    ASSERT_EQ(s.value().roots.size(), 3U);
    // Frobenius must be a 3-cycle; no root lies in Z_7 (f has no rational
    // root), and e_2 = coefficient of x must vanish.
    EXPECT_EQ(perm_cycle_type(s.value().frobenius),
              (std::vector<std::size_t>{3U}));
    RingElem e2 = R.zero();
    for (std::size_t i = 0U; i < 3U; ++i) {
        EXPECT_FALSE(R.integer_residue(s.value().roots[i]).has_value());
        for (std::size_t j = i + 1U; j < 3U; ++j) {
            e2 = R.add(e2,
                       R.mul(s.value().roots[i], s.value().roots[j]));
        }
    }
    EXPECT_TRUE(R.is_zero(e2));
    expect_reconstructs(s.value());
}

// ── two quadratic factors: Cantor-Zassenhaus split inside GF(p²) ────────────

TEST(GaloisPadicSplit, TwoQuadraticFactorsCZ) {
    // f = (x²+1)(x²−2) = x⁴ − x² − 2; mod 11 both quadratics are
    // irreducible (−1 and 2 are non-residues mod 11), so one factor is Φ
    // and the other must be split by Cantor-Zassenhaus inside GF(121).
    const IntPoly f = make_poly({-2, 0, -1, 0, 1});
    auto fact = factor_polynomial_mod_p(f, BigInt(11), nullptr);
    ASSERT_TRUE(fact.is_ok());
    ASSERT_EQ(fact.value().size(), 2U);
    ASSERT_EQ(fact.value()[0].degree(), 2U);
    ASSERT_EQ(fact.value()[1].degree(), 2U);
    const SplittingPrime sp{BigInt(11), fact.value(), 2U};
    auto s = build_padic_splitting(f, sp, 6U, nullptr);
    ASSERT_TRUE(s.is_ok());
    EXPECT_EQ(perm_cycle_type(s.value().frobenius),
              (std::vector<std::size_t>{2U, 2U}));
    for (const auto& r : s.value().roots) {
        EXPECT_FALSE(s.value().ring.integer_residue(r).has_value());
    }
    expect_reconstructs(s.value());
}

// ── p = 2: trace splitter + Rabin sweep for Φ ───────────────────────────────

TEST(GaloisPadicSplit, TwoAdicTraceSplitAndRabinPhi) {
    // f = (x²+x+1)(x³−x+1) = x⁵ + x⁴ + 1. Mod 2 the factors have degrees
    // {2, 3}, so L = 6 and no factor supplies Φ: the exhaustive Rabin sweep
    // must find an irreducible sextic over F_2, and both factors' roots are
    // found by the GF(2^L) trace splitter.
    const IntPoly f = make_poly({1, 0, 0, 0, 1, 1});
    auto fact = factor_polynomial_mod_p(f, BigInt(2), nullptr);
    ASSERT_TRUE(fact.is_ok());
    std::size_t total = 0U;
    for (const auto& g : fact.value()) total += g.degree();
    ASSERT_EQ(total, 5U);  // squarefree mod 2
    const SplittingPrime sp{BigInt(2), fact.value(), 6U};
    auto s = build_padic_splitting(f, sp, 6U, nullptr);
    ASSERT_TRUE(s.is_ok());
    EXPECT_EQ(s.value().ring.ext_degree(), 6U);
    EXPECT_EQ(perm_cycle_type(s.value().frobenius),
              (std::vector<std::size_t>{2U, 3U}));
    expect_reconstructs(s.value());
}

// ── precision raise: refinement in place ────────────────────────────────────

TEST(GaloisPadicSplit, PrecisionRaiseRefinesInPlace) {
    const IntPoly f = make_poly({-2, 0, -1, 0, 1});
    auto fact = factor_polynomial_mod_p(f, BigInt(11), nullptr);
    ASSERT_TRUE(fact.is_ok());
    const SplittingPrime sp{BigInt(11), fact.value(), 2U};
    auto s4 = build_padic_splitting(f, sp, 4U, nullptr);
    ASSERT_TRUE(s4.is_ok());
    auto s12 = raise_splitting_precision(s4.value(), 12U, nullptr);
    ASSERT_TRUE(s12.is_ok());
    EXPECT_EQ(s12.value().ring.precision(), 12U);
    EXPECT_EQ(s12.value().frobenius, s4.value().frobenius);
    // Same roots mod 11⁴ (refinement never relocates a root).
    const BigInt p4 = s4.value().ring.modulus();
    ASSERT_EQ(s4.value().roots.size(), s12.value().roots.size());
    for (std::size_t i = 0U; i < s4.value().roots.size(); ++i) {
        for (std::size_t j = 0U; j < 2U; ++j) {
            BigInt d =
                (s12.value().roots[i][j] - s4.value().roots[i][j]) % p4;
            if (d.is_negative()) d += p4;
            EXPECT_TRUE(d.is_zero());
        }
    }
    expect_reconstructs(s12.value());
    // No-op raise returns the same splitting.
    auto same = raise_splitting_precision(s12.value(), 5U, nullptr);
    ASSERT_TRUE(same.is_ok());
    EXPECT_EQ(same.value().ring.precision(), 12U);
}

// ── budgets and structured failures ─────────────────────────────────────────

TEST(GaloisPadicSplit, PrimeBudgetIsStructured) {
    const IntPoly f = make_poly({-6, 11, -6, 1});
    // Zero budget: structured Unimplemented, never a guess.
    auto none = choose_splitting_prime(f, 0U, nullptr);
    ASSERT_TRUE(none.is_error());
    EXPECT_EQ(none.error().kind, CASErrorKind::Unimplemented);
    // Budget 1 tries only p = 2, which is ramified for this f
    // (f ≡ x(x+1)² mod 2): still a structured failure.
    auto one = choose_splitting_prime(f, 1U, nullptr);
    ASSERT_TRUE(one.is_error());
    EXPECT_EQ(one.error().kind, CASErrorKind::Unimplemented);
}

}  // namespace
