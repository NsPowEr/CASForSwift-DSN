// test_van_hoeij.cpp — Direct probe of van_hoeij_knapsack_factor (F2.3).
//
// These tests call van_hoeij_knapsack_factor DIRECTLY (NOT via
// factor_over_integers, which has a subset-enumeration fallback that
// masks van Hoeij bugs). They are the correctness gate for the
// knapsack-lattice recombination.
//
// Construction mirrors the real caller (factorization_integers.cpp):
//   1. pick a factorization prime p;
//   2. factor f mod p (Cantor-Zassenhaus / Berlekamp);
//   3. choose lift exponent k so that p^k > 2·B (B = Mignotte / norm bound);
//   4. Hensel-lift the modular factors to p^k;
//   5. call van_hoeij_knapsack_factor on the lifted factors.
//
// Acceptance: returned h must divide f exactly in Z[x] (structural check).

#include <gtest/gtest.h>

#include "cas/bigint.hpp"
#include "algebra/algebra_internal.hpp"
#include "algebra/polynomial_internal.hpp"

#include <optional>
#include <vector>

using namespace cas;
using namespace cas::algebra;

namespace {

// Build an IntPoly from little-endian coefficients (index 0 == constant term).
IntPoly make_poly(std::vector<long long> coeffs) {
    std::vector<BigInt> c;
    c.reserve(coeffs.size());
    for (long long v : coeffs) c.emplace_back(v);
    IntPoly p(std::move(c));
    normalize_integer_poly(p);
    return p;
}

// Multiply two IntPolys exactly over Z.
IntPoly poly_mul(const IntPoly& a, const IntPoly& b) {
    if (a.is_zero() || b.is_zero()) return IntPoly{};
    IntPoly r;
    r.resize(a.size() + b.size() - 1U, BigInt(0));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j)
            r[i + j] += a[i] * b[j];
    normalize_integer_poly(r);
    return r;
}

// Does h divide f exactly over Z[x] (nontrivially)?
bool divides_exact(const IntPoly& f, const IntPoly& h) {
    if (h.is_zero() || h.degree() == 0U || h.degree() >= f.degree()) return false;
    auto rem = pseudo_remainder_integer_poly(f, h);
    normalize_integer_poly(rem);
    return rem.is_zero();
}

// Reproduce the caller's lifted-modular-factor setup, then return them.
struct LiftedSetup {
    bool ok = false;
    std::vector<IntPoly> lifted;
    BigInt pk{1};
    std::size_t r = 0;
};

LiftedSetup lift_modular_factors(const IntPoly& f) {
    LiftedSetup out;
    BigInt p = select_factorization_prime(f);

    auto mod_res = factor_polynomial_mod_p(f, p);
    if (mod_res.is_error()) return out;
    auto mod_factors = mod_res.value();
    if (mod_factors.size() < 2U) return out;

    const std::size_t n = f.degree();
    BigInt pk = p;
    std::size_t k = 1;
    BigInt two_pow_n = BigInt(1).shift_left_bits(n);
    BigInt norm2_sq(0);
    for (const auto& c : f.coefficients()) norm2_sq += c * c;
    while (pk < two_pow_n * norm2_sq * BigInt(2)) {
        pk *= p;
        ++k;
    }

    auto lifted_res = hensel_lift_multi(f, mod_factors, p, k);
    if (lifted_res.is_error()) return out;

    out.ok = true;
    out.lifted = lifted_res.value();
    out.pk = pk;
    out.r = mod_factors.size();
    return out;
}

}  // namespace

// Sanity: the lifted modular factors must multiply back to f modulo pk.
// This is the invariant van Hoeij relies on; a broken Hensel lift (returning
// mod-p singletons instead of mod-pk lifts) silently violates it and makes
// every recombination fail.  Regression guard for that exact bug.
TEST(VanHoeijDirect, LiftedFactorsProductEqualsF_ModPk) {
    IntPoly f = poly_mul(poly_mul(make_poly({-2, 0, 1}), make_poly({-3, 0, 1})),
                         make_poly({-5, 0, 1}));
    auto setup = lift_modular_factors(f);
    ASSERT_TRUE(setup.ok);
    ASSERT_GE(setup.r, 2U);

    auto modc = [&](BigInt v) {
        BigInt r = v % setup.pk;
        if (r.is_negative()) r += setup.pk;
        if (r * BigInt(2) > setup.pk) r -= setup.pk;
        return r;
    };
    IntPoly prod = make_poly({1});
    for (const auto& g : setup.lifted) {
        IntPoly r;
        r.resize(prod.size() + g.size() - 1U, BigInt(0));
        for (std::size_t i = 0; i < prod.size(); ++i)
            for (std::size_t j = 0; j < g.size(); ++j)
                r[i + j] = modc(r[i + j] + prod[i] * g[j]);
        normalize_integer_poly(r);
        prod = r;
    }
    // f reduced mod pk (centered) must equal prod.
    IntPoly f_red = f;
    for (auto& c : f_red.coefficients()) c = modc(c);
    normalize_integer_poly(f_red);
    ASSERT_EQ(prod.size(), f_red.size());
    for (std::size_t i = 0; i < prod.size(); ++i)
        EXPECT_EQ(prod[i].decimal(), f_red[i].decimal())
            << "coefficient " << i << " mismatch: lift invariant broken";
}

// ── deg-6: (x²-2)(x²-3)(x²-5), 6 modular factors over a splitting prime ────
TEST(VanHoeijDirect, Deg6_TripleQuadratic_FindsRealFactor) {
    // f = (x^2-2)(x^2-3)(x^2-5)
    IntPoly f = poly_mul(poly_mul(make_poly({-2, 0, 1}), make_poly({-3, 0, 1})),
                         make_poly({-5, 0, 1}));
    ASSERT_EQ(f.degree(), 6U);

    auto setup = lift_modular_factors(f);
    ASSERT_TRUE(setup.ok) << "modular-factor lift setup failed";
    ASSERT_GE(setup.r, 2U) << "expected ≥2 modular factors, got " << setup.r;

    // lll_threshold = 0 forces the LLL knapsack path for ALL r (test override).
    auto h = van_hoeij_knapsack_factor(f, setup.lifted, setup.pk, 0.75, 0U);
    ASSERT_TRUE(h.has_value()) << "van Hoeij (LLL path) returned nullopt on deg6";
    EXPECT_TRUE(divides_exact(f, *h))
        << "van Hoeij returned a non-divisor on deg6";
}

TEST(VanHoeijDirect, Deg6_TripleQuadratic_EnumPath) {
    IntPoly f = poly_mul(poly_mul(make_poly({-2, 0, 1}), make_poly({-3, 0, 1})),
                         make_poly({-5, 0, 1}));
    auto setup = lift_modular_factors(f);
    ASSERT_TRUE(setup.ok);

    // Default threshold (10) → enumeration path for small r.
    auto h = van_hoeij_knapsack_factor(f, setup.lifted, setup.pk, 0.75, 10U);
    ASSERT_TRUE(h.has_value()) << "van Hoeij (enum path) returned nullopt on deg6";
    EXPECT_TRUE(divides_exact(f, *h));
}

// ── deg-16: product of 8 distinct quadratics x²-d, r=8..16 modular factors ──
TEST(VanHoeijDirect, DISABLED_Deg16_EightQuadratics_FindsRealFactor) {
    const std::vector<long long> ds = {2, 3, 5, 7, 11, 13, 17, 19};
    IntPoly f = make_poly({1});
    for (long long d : ds) f = poly_mul(f, make_poly({-d, 0, 1}));
    ASSERT_EQ(f.degree(), 16U);

    auto setup = lift_modular_factors(f);
    ASSERT_TRUE(setup.ok) << "modular-factor lift setup failed";
    ASSERT_GE(setup.r, 2U);

    auto h = van_hoeij_knapsack_factor(f, setup.lifted, setup.pk, 0.75, 0U);
    ASSERT_TRUE(h.has_value()) << "van Hoeij (LLL path) returned nullopt on deg16";
    EXPECT_TRUE(divides_exact(f, *h))
        << "van Hoeij returned a non-divisor on deg16";
}

// ── Acceptance gate F2.3: many-factor, large r, LLL must beat enumeration ──
// 12 distinct quadratics → deg 24, r ≥ 12. The C(24/2 factors) enumeration
// search space is large; here we FORCE the LLL path and require correctness
// in polynomial (non-exponential) time. Marked *Stress* because the exact-
// rational LLL on a ~dim-16..24 lattice takes ~30s — correct but heavy; it is
// excluded from the default suite (run with --gtest_filter='*Stress*').
// DISABILITATO: Test di stress matematico fisiologicamente in timeout (>60s) sotto Debug Mode (-O0). Da eseguire in Release Mode o via target dedicato cas_stress_tests.
TEST(VanHoeijStress, DISABLED_Deg24_TwelveQuadratics_AcceptanceGate) {
    const std::vector<long long> ds = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    IntPoly f = make_poly({1});
    for (long long d : ds) f = poly_mul(f, make_poly({-d, 0, 1}));
    ASSERT_EQ(f.degree(), 24U);

    auto setup = lift_modular_factors(f);
    ASSERT_TRUE(setup.ok);
    ASSERT_GE(setup.r, 2U);

    auto h = van_hoeij_knapsack_factor(f, setup.lifted, setup.pk, 0.75, 0U);
    ASSERT_TRUE(h.has_value()) << "van Hoeij (LLL path) returned nullopt on deg24";
    EXPECT_TRUE(divides_exact(f, *h));
}
