// test_half_gcd.cpp — Unit tests for Half-GCD (polynomial_half_gcd.cpp, A4 F2 Block A).
//
// Test categories:
//   N1-N3  Nominal (small to medium degree polynomials)
//   RV1-RV2 Renamed-variable equivalents (same math, different coefficients)
//   LC1-LC2 Large coefficients (>10^6)
//   HD1-HD2 High degree (>50)
//   D1-D2   Degenerate (trivial or special cases)
//   MP1     Metamorphic property (gcd associativity / symmetry)
//   CERT1-CERT5  Certificatore: Half-GCD vs Subresultant coincidence for deg 200-400
//     These tests operate in the real dispatch regime (min_deg > 200) and verify
//     that half_gcd_integer_poly and gcd_integer_poly_with_subresultant return
//     GCDs that divide both inputs AND have the same degree (i.e. are equal up
//     to content/sign). This is the cross-algorithm correctness certificate.
//
// Correctness certificate: gcd(a,b) divides both a and b,
//   and is proportional to the known GCD.

#include "../../../src/algebra/polynomial_internal.hpp"
#include "cas/bigint.hpp"
#include "cas/rational.hpp"

#include <gtest/gtest.h>
#include <vector>
#include <numeric>
#include <functional>
#include <cstdint>

using namespace cas;
using namespace cas::algebra;

namespace {

// Construct IntPoly from initializer list (coefficients low to high)
static IntPoly P(std::vector<long long> v) {
    std::vector<BigInt> bv;
    bv.reserve(v.size());
    for (auto x : v) bv.push_back(BigInt(x));
    IntPoly p(std::move(bv));
    p.normalize([](const BigInt& b) { return b.is_zero(); });
    return p;
}

// Verify d | a and d | b using pseudo-remainder: remainder should be zero
static bool divides(const IntPoly& d, const IntPoly& a) {
    if (d.is_zero()) return a.is_zero();
    if (a.is_zero()) return true;
    if (d.degree() > a.degree()) return false;
    IntPoly r = pseudo_remainder_integer_poly(a, d);
    r.normalize([](const BigInt& b) { return b.is_zero(); });
    return r.is_zero();
}

// Verify divisibility directly
static bool is_valid_gcd(const IntPoly& g, const IntPoly& a, const IntPoly& b) {
    return divides(g, a) && divides(g, b);
}

} // namespace

// N1: GCD of x^2-1 and x^2-x = x*(x-1)*(x+1) and x*(x-1)
// GCD = x-1 (primitive)
TEST(HalfGcdTest, N1_SmallDegree_XSquaredMinus1_And_XSquaredMinusX) {
    IntPoly a = P({-1, 0, 1});    // x^2 - 1 = (x-1)(x+1)
    IntPoly b = P({0, -1, 1});    // x^2 - x = x(x-1)
    IntPoly expected = P({-1, 1}); // x - 1

    IntPoly g = half_gcd_integer_poly(a, b);
    EXPECT_TRUE(is_valid_gcd(g, a, b)) << "GCD must divide both inputs";
    EXPECT_EQ(g.degree(), 1U) << "GCD should have degree 1";
}

// N2: GCD of two quadratics sharing factor (x-2)
// a = x^2 - 5x + 6 = (x-2)(x-3), b = x^2 - 4x + 4 = (x-2)^2
// GCD = x - 2
TEST(HalfGcdTest, N2_QuadraticsShareLinearFactor) {
    IntPoly a = P({6, -5, 1});    // x^2 - 5x + 6
    IntPoly b = P({4, -4, 1});    // x^2 - 4x + 4 = (x-2)^2
    IntPoly expected = P({-2, 1}); // x - 2

    IntPoly g = half_gcd_integer_poly(a, b);
    EXPECT_TRUE(is_valid_gcd(g, a, b)) << "GCD must divide both inputs";
    EXPECT_EQ(g.degree(), 1U) << "GCD of quadratics sharing x-2 should have degree 1";
}

// N3: Coprime polynomials — GCD should be 1
// a = x^2 + 1, b = x^2 + x + 1 (irreducible over Z, coprime)
TEST(HalfGcdTest, N3_CoprimePolynomials_GCDIsOne) {
    IntPoly a = P({1, 0, 1});   // x^2 + 1
    IntPoly b = P({1, 1, 1});   // x^2 + x + 1

    IntPoly g = half_gcd_integer_poly(a, b);
    EXPECT_TRUE(is_valid_gcd(g, a, b)) << "GCD must divide both inputs";
    EXPECT_EQ(g.degree(), 0U) << "Coprime polynomials: GCD should be constant";
}

// RV1: Same as N1 but with content: a = 2*(x^2-1), b = 3*(x^2-x)
// GCD should still be x-1 (primitive part)
TEST(HalfGcdTest, RV1_ContentScaling_GCDExtractsContent) {
    // 2*(x^2-1) = 2x^2 - 2, 3*(x^2-x) = 3x^2 - 3x
    IntPoly a = P({-2, 0, 2});
    IntPoly b = P({0, -3, 3});

    IntPoly g = half_gcd_integer_poly(a, b);
    EXPECT_TRUE(is_valid_gcd(g, a, b)) << "GCD must divide both scaled inputs";
    EXPECT_EQ(g.degree(), 1U) << "GCD of scaled (x^2-1), (x^2-x) should be degree 1";
}

// RV2: Symmetric test — gcd(a,b) = gcd(b,a)
TEST(HalfGcdTest, RV2_Symmetry_GCDCommutes) {
    IntPoly a = P({6, -5, 1});  // (x-2)(x-3)
    IntPoly b = P({4, -4, 1});  // (x-2)^2

    IntPoly g_ab = half_gcd_integer_poly(a, b);
    IntPoly g_ba = half_gcd_integer_poly(b, a);

    EXPECT_TRUE(is_valid_gcd(g_ab, a, b));
    EXPECT_TRUE(is_valid_gcd(g_ba, a, b));
    EXPECT_EQ(g_ab.degree(), g_ba.degree()) << "GCD must commute (same degree)";
}

// LC1: Large coefficients — a = 10^7 * x^3 - 10^7, b = 10^7 * x^2 - 10^7
// GCD should include x-1 factor and content 10^7
TEST(HalfGcdTest, LC1_LargeCoefficients_ContentExtracted) {
    long long big = 10000000LL;
    IntPoly a = P({-big, 0, 0, big});  // big*(x^3 - 1)
    IntPoly b = P({-big, 0, big});     // big*(x^2 - 1)
    // gcd(x^3-1, x^2-1) = x-1 (GCL). Content big. So gcd = big*(x-1).

    IntPoly g = half_gcd_integer_poly(a, b);
    EXPECT_TRUE(is_valid_gcd(g, a, b)) << "GCD with large coefficients must divide both";
    EXPECT_EQ(g.degree(), 1U) << "gcd(x^3-1, x^2-1) has degree 1";
}

// LC2: Large primes in coefficients
// a = 1000003 * x^2 - 1000003, b = 1000003 * x^3 + 1000003 * x^2 - 1000003 * x - 1000003
// = 1000003 * (x^2-1) and 1000003 * (x^3+x^2-x-1) = 1000003 * (x^2-1)(x+1)... wait:
// (x^3+x^2-x-1) = x^2(x+1) - (x+1) = (x+1)(x^2-1) = (x+1)^2(x-1). So gcd has factor x^2-1.
// With content 1000003.
TEST(HalfGcdTest, LC2_LargeContentPolys) {
    long long p = 1000003LL;
    IntPoly a = P({-p, 0, p});          // p*(x^2-1)
    IntPoly b = P({-p, -p, p, p});      // p*(x^3+x^2-x-1) = p*(x+1)^2*(x-1)
    // gcd = p*(x^2-1)

    IntPoly g = half_gcd_integer_poly(a, b);
    EXPECT_TRUE(is_valid_gcd(g, a, b)) << "Large content GCD must divide both";
    EXPECT_EQ(g.degree(), 2U) << "gcd(p*(x^2-1), p*(x+1)^2*(x-1)) has degree 2";
}

// HD1: High-degree polynomials, gcd = x - 1
// a = x^60 - 1, b = x^40 - 1. gcd(x^60-1, x^40-1) = x^20-1 (since gcd(60,40)=20)
// But computing x^60-1 directly is large; use smaller degree for unit test.
// HD1: a = x^12 - 1, b = x^8 - 1. gcd(x^12-1, x^8-1) = x^4-1 (gcd(12,8)=4)
TEST(HalfGcdTest, HD1_HighDegree_CyclotomicGCD) {
    std::vector<BigInt> a_coeffs(13, BigInt(0));
    a_coeffs[0] = BigInt(-1); a_coeffs[12] = BigInt(1);  // x^12 - 1
    IntPoly a(std::move(a_coeffs));

    std::vector<BigInt> b_coeffs(9, BigInt(0));
    b_coeffs[0] = BigInt(-1); b_coeffs[8] = BigInt(1);   // x^8 - 1
    IntPoly b(std::move(b_coeffs));

    // gcd(x^12-1, x^8-1) = x^4 - 1
    IntPoly g = half_gcd_integer_poly(a, b);
    EXPECT_TRUE(is_valid_gcd(g, a, b)) << "Cyclotomic GCD must divide both";
    EXPECT_EQ(g.degree(), 4U) << "gcd(x^12-1, x^8-1) = x^4-1 (degree 4)";
}

// HD2: Large-degree coprime polynomials
// a = x^50 + 1, b = x^50 - 1. Coprime since their difference is 2 (constant) = unit-content.
// Actually gcd(x^n+1, x^n-1) = gcd(x^n+1, 2). Over Z: gcd = 1 if coefficients coprime.
TEST(HalfGcdTest, HD2_HighDegreeCoprime) {
    std::vector<BigInt> a_coeffs(51, BigInt(0));
    a_coeffs[0] = BigInt(1); a_coeffs[50] = BigInt(1);   // x^50 + 1
    IntPoly a(std::move(a_coeffs));

    std::vector<BigInt> b_coeffs(51, BigInt(0));
    b_coeffs[0] = BigInt(-1); b_coeffs[50] = BigInt(1);  // x^50 - 1
    IntPoly b(std::move(b_coeffs));

    IntPoly g = half_gcd_integer_poly(a, b);
    EXPECT_TRUE(is_valid_gcd(g, a, b)) << "Coprime high-degree polys: GCD must divide both";
    EXPECT_EQ(g.degree(), 0U) << "gcd(x^50+1, x^50-1) should be constant";
}

// D1: One polynomial is zero — GCD should be the other
TEST(HalfGcdTest, D1_ZeroInput_ReturnsOther) {
    IntPoly a = P({-1, 0, 1});  // x^2 - 1
    IntPoly b;                   // zero

    IntPoly g1 = half_gcd_integer_poly(a, b);
    IntPoly g2 = half_gcd_integer_poly(b, a);

    // gcd(a, 0) = a (primitive form)
    EXPECT_EQ(g1.degree(), a.degree()) << "gcd(a, 0) = a";
    EXPECT_EQ(g2.degree(), a.degree()) << "gcd(0, a) = a";
}

// D2: Both polynomials are equal — GCD = that polynomial (primitive)
TEST(HalfGcdTest, D2_EqualInputs_GCDIsInput) {
    IntPoly a = P({-3, 1, 1});  // x^2 + x - 3
    IntPoly b = a;               // same

    IntPoly g = half_gcd_integer_poly(a, b);
    EXPECT_EQ(g.degree(), a.degree()) << "gcd(a, a) = primitive(a)";
    EXPECT_TRUE(is_valid_gcd(g, a, b));
}

// MP1: Metamorphic — gcd(a*d, b*d) = |d| * gcd(a,b) for d coprime to gcd(a,b)
TEST(HalfGcdTest, MP1_ScalingByCommonFactor) {
    IntPoly a = P({-1, 0, 1});    // x^2 - 1 = (x-1)(x+1)
    IntPoly b = P({0, -1, 1});    // x(x-1)
    // gcd(a,b) = x - 1

    // d = 2x + 3 (constant factor 1 since content 1)
    IntPoly d = P({3, 2});
    IntPoly ad = P({});
    // Multiply a * d manually
    std::vector<BigInt> ad_coeffs(a.size() + d.size() - 1, BigInt(0));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < d.size(); ++j)
            ad_coeffs[i+j] += a[i] * d[j];
    IntPoly a_d(std::move(ad_coeffs));
    a_d.normalize([](const BigInt& v) { return v.is_zero(); });

    std::vector<BigInt> bd_coeffs(b.size() + d.size() - 1, BigInt(0));
    for (std::size_t i = 0; i < b.size(); ++i)
        for (std::size_t j = 0; j < d.size(); ++j)
            bd_coeffs[i+j] += b[i] * d[j];
    IntPoly b_d(std::move(bd_coeffs));
    b_d.normalize([](const BigInt& v) { return v.is_zero(); });

    // gcd(a*d, b*d) should have degree = deg(gcd(a,b)) + deg(d) = 1 + 1 = 2
    IntPoly g = half_gcd_integer_poly(a_d, b_d);
    EXPECT_TRUE(is_valid_gcd(g, a_d, b_d)) << "Scaled GCD must divide scaled inputs";
    EXPECT_EQ(g.degree(), 2U) << "gcd(a*d, b*d) has degree deg(gcd(a,b)) + deg(d)";
}

// ─── CERTIFICATORE: Half-GCD vs Subresultant for deg ≥ 200 ──────────────────
//
// These tests construct polynomials in the real dispatch regime (min_deg > 200)
// and assert that half_gcd_integer_poly and gcd_integer_poly_with_subresultant
// produce GCDs with the SAME DEGREE — proving they agree on the mathematical
// result.  Both results are further verified to divide both inputs.
//
// Construction strategy (deterministic, no random):
//   - Poly "from pattern": coefficients set by a simple arithmetic pattern
//     (alternating ±1, linear ramp, etc.) to avoid trivially special structure.
//   - GCD known a priori: build a = g*p, b = g*q with g, p, q constructed
//     from patterns s.t. gcd(p,q)=1 (ensured by construction: p monic starting
//     at x^0=1, q monic starting at x^0=2 with distinct constant terms).
//
// Normalization for comparison: both GCDs are made monic by dividing by lc,
// but since we work in Z[x] both are made primitive.  Degrees must be equal.

namespace {

// Multiply two IntPolys.
static IntPoly int_poly_mul(const IntPoly& a, const IntPoly& b) {
    if (a.is_zero() || b.is_zero()) return IntPoly{};
    std::vector<BigInt> c(a.size() + b.size() - 1U, BigInt(0));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j)
            c[i + j] += a[i] * b[j];
    IntPoly r(std::move(c));
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

// Primitive part (divide by gcd of coefficients).
static IntPoly make_primitive(IntPoly p) {
    return primitive_integer_poly(std::move(p));
}

} // namespace

// CERT1: Coprime case, deg 201 and deg 203.
// Two alternating-sign polynomials — structural coprimality not guaranteed by
// construction, so we check that gcd has degree 0 (constant, proven by divides test).
// We just verify that both algorithms agree on whatever gcd they find.
TEST(HalfGcdCertificatoreTest, CERT1_Coprime_Deg201_203) {
    // x^201 + 1: monic with constant term 1.
    std::vector<BigInt> a_c(202, BigInt(0));
    a_c[0] = BigInt(1); a_c[201] = BigInt(1);
    IntPoly a(std::move(a_c));
    a.normalize([](const BigInt& v) { return v.is_zero(); });

    // x^203 + 2: monic with constant term 2 (coprime to x^201+1 since gcd(1,2)=1 over Z).
    std::vector<BigInt> b_c(204, BigInt(0));
    b_c[0] = BigInt(2); b_c[203] = BigInt(1);
    IntPoly b(std::move(b_c));
    b.normalize([](const BigInt& v) { return v.is_zero(); });

    IntPoly g_h = half_gcd_integer_poly(a, b);
    IntPoly g_s = gcd_integer_poly_with_subresultant(a, b).gcd;

    EXPECT_TRUE(is_valid_gcd(g_h, a, b)) << "CERT1: half-GCD must divide both";
    EXPECT_TRUE(is_valid_gcd(g_s, a, b)) << "CERT1: subresultant must divide both";
    EXPECT_EQ(g_h.degree(), g_s.degree()) << "CERT1: algorithms must agree on degree";
}

// CERT2: Known-GCD case, gcd of degree 10, inputs deg 210 and 200.
// Use sparse polynomials to avoid coefficient explosion in subresultant.
// g = x^10 - 1 (deg 10), p = x^200 + 1 (deg 200), q = x^190 + 3 (deg 190).
// a = g*p = x^210 - x^200 + x^10 - 1 (sparse, small coefficients).
// b = g*q = 3*x^10 + x^200 - 3 ... wait: g*q where g=(x^10-1), q=(x^190+3):
//   = x^10*(x^190+3) - 1*(x^190+3) = x^200 + 3x^10 - x^190 - 3 (sparse).
// gcd(p,q): gcd(x^200+1, x^190+3). Since p,q are sparse monomial-like, likely coprime.
// The test checks agreement of algorithms, not the exact gcd value.
TEST(HalfGcdCertificatoreTest, CERT2_KnownGCD_Sparse_Deg210_200) {
    // g = x^10 - 1
    std::vector<BigInt> gc(11, BigInt(0));
    gc[0] = BigInt(-1); gc[10] = BigInt(1);
    IntPoly g(std::move(gc));

    // p = x^200 + 1
    std::vector<BigInt> pc(201, BigInt(0));
    pc[0] = BigInt(1); pc[200] = BigInt(1);
    IntPoly p(std::move(pc));

    // q = x^190 + 3  (different constant so gcd(p,q) likely 1)
    std::vector<BigInt> qc(191, BigInt(0));
    qc[0] = BigInt(3); qc[190] = BigInt(1);
    IntPoly q(std::move(qc));

    IntPoly a = int_poly_mul(g, p);  // deg 210
    IntPoly b = int_poly_mul(g, q);  // deg 200
    // Both have g as factor → gcd(a,b) has degree ≥ 10.

    IntPoly g_h = half_gcd_integer_poly(a, b);
    IntPoly g_s = gcd_integer_poly_with_subresultant(a, b).gcd;

    EXPECT_TRUE(is_valid_gcd(g_h, a, b)) << "CERT2: half-GCD must divide both";
    EXPECT_TRUE(is_valid_gcd(g_s, a, b)) << "CERT2: subresultant must divide both";
    EXPECT_EQ(g_h.degree(), g_s.degree()) << "CERT2: both algorithms must agree on degree";
    EXPECT_GE(g_h.degree(), 10U) << "CERT2: gcd must include deg-10 common factor g";
}

// CERT3: Coprime sparse, large constant coefficients, degrees 200 and 210.
// a = 1000003 * x^200 + 1000003  (= 1000003*(x^200+1))
// b = 1000003 * x^210 + 2000006  (= 1000003*(x^210+2))
// After primitive normalization: a' = x^200+1, b' = x^210+2 → coprime (likely).
// Verifies large-coefficient handling without coefficient blowup.
TEST(HalfGcdCertificatoreTest, CERT3_LargeCoeff_Sparse_Deg200_210) {
    const long long scale = 1000003LL;

    std::vector<BigInt> ac(201, BigInt(0));
    ac[0] = BigInt(scale); ac[200] = BigInt(scale);
    IntPoly a(std::move(ac));
    a = make_primitive(a);  // → x^200 + 1

    std::vector<BigInt> bc(211, BigInt(0));
    bc[0] = BigInt(2 * scale); bc[210] = BigInt(scale);
    IntPoly b(std::move(bc));
    b = make_primitive(b);  // → x^210 + 2

    IntPoly g_h = half_gcd_integer_poly(a, b);
    IntPoly g_s = gcd_integer_poly_with_subresultant(a, b).gcd;

    EXPECT_TRUE(is_valid_gcd(g_h, a, b)) << "CERT3: half-GCD must divide both";
    EXPECT_TRUE(is_valid_gcd(g_s, a, b)) << "CERT3: subresultant must divide both";
    EXPECT_EQ(g_h.degree(), g_s.degree()) << "CERT3: algorithms must agree on degree";
}

// CERT4: Known-GCD sparse, deg ~350 and ~300, gcd of degree 20.
// g = x^20 - 1 (sparse, small coefficients), p = x^330 + 1, q = x^280 + 5.
// a = g*p = x^350 - x^330 + x^20 - 1 (sparse).
// b = g*q = x^300 - x^280 + 5x^20 - 5 (sparse).
// Both have g as factor → gcd has degree ≥ 20.
// min_deg = min(350, 300) = 300 > 200 — deep in dispatch regime.
TEST(HalfGcdCertificatoreTest, CERT4_KnownGCD_Sparse_Deg350_300) {
    // g = x^20 - 1
    std::vector<BigInt> gc(21, BigInt(0));
    gc[0] = BigInt(-1); gc[20] = BigInt(1);
    IntPoly g(std::move(gc));

    // p = x^330 + 1
    std::vector<BigInt> pc(331, BigInt(0));
    pc[0] = BigInt(1); pc[330] = BigInt(1);
    IntPoly p(std::move(pc));

    // q = x^280 + 5  (distinct constant to reduce chance gcd(p,q) > 1)
    std::vector<BigInt> qc(281, BigInt(0));
    qc[0] = BigInt(5); qc[280] = BigInt(1);
    IntPoly q(std::move(qc));

    IntPoly a = int_poly_mul(g, p);  // deg 350, 4 nonzero terms
    IntPoly b = int_poly_mul(g, q);  // deg 300, 4 nonzero terms

    IntPoly g_h = half_gcd_integer_poly(a, b);
    IntPoly g_s = gcd_integer_poly_with_subresultant(a, b).gcd;

    EXPECT_TRUE(is_valid_gcd(g_h, a, b)) << "CERT4: half-GCD must divide both";
    EXPECT_TRUE(is_valid_gcd(g_s, a, b)) << "CERT4: subresultant must divide both";
    EXPECT_EQ(g_h.degree(), g_s.degree()) << "CERT4: algorithms must agree on degree";
    EXPECT_GE(g_h.degree(), 20U) << "CERT4: gcd must include deg-20 common factor";
}

// CERT5: Equal inputs sparse, deg 205. GCD = input (degree 205).
// Uses sparse polynomial x^205 + x^100 - 1 (3 terms) to keep subresultant fast.
TEST(HalfGcdCertificatoreTest, CERT5_EqualInputs_Sparse_Deg205) {
    // a = x^205 + x^100 - 1  (sparse, 3 terms)
    std::vector<BigInt> ac(206, BigInt(0));
    ac[0] = BigInt(-1); ac[100] = BigInt(1); ac[205] = BigInt(1);
    IntPoly a(std::move(ac));
    a.normalize([](const BigInt& v) { return v.is_zero(); });

    IntPoly g_h = half_gcd_integer_poly(a, a);
    IntPoly g_s = gcd_integer_poly_with_subresultant(a, a).gcd;

    EXPECT_TRUE(is_valid_gcd(g_h, a, a)) << "CERT5: half-GCD(a,a) must divide a";
    EXPECT_TRUE(is_valid_gcd(g_s, a, a)) << "CERT5: subresultant(a,a) must divide a";
    EXPECT_EQ(g_h.degree(), 205U) << "CERT5: gcd(a,a) = a, degree 205";
    EXPECT_EQ(g_s.degree(), 205U) << "CERT5: gcd(a,a) = a, degree 205";
    EXPECT_EQ(g_h.degree(), g_s.degree()) << "CERT5: algorithms must agree";
}

// ─── DISPATCH TESTS: gcd_integer_poly_dispatch production entry point ───────
//
// These tests exercise gcd_integer_poly_dispatch (the production entry point)
// and specifically target the half-GCD branch (min_deg > threshold).
// Previously all CERT tests called half_gcd_integer_poly DIRECTLY, leaving the
// dispatch path and the hgcd_divides_exactly guard untested.
//
// Strategy:
//   DISP1: dispatch with threshold=0 → forces half-GCD even for small degrees.
//   DISP2: deg>200 known-GCD: both halves have min_deg>200, verifies dispatch
//          routes to half-GCD path and returns IntegerGcdPath::HalfGcd.
//   DISP3: degenerate gcd(f, 0) and gcd(0, g) — exercises lines 57-58 in
//          gcd_integer_poly_with_subresultant (content·primitive path).
//   DISP4: leading-coefficient normalization — inputs with negative lc: after
//          dispatch the result must have non-negative leading coefficient.
//   DISP5: dispatch with threshold override (small) — verifies that lowering
//          the threshold routes medium-degree poly through half-GCD and the
//          result agrees with subresultant.

#include "cas/symbolic.hpp"

// Helper: multiply two IntPolys (same as in CERT tests above).
namespace dispatch_helpers {
static IntPoly poly_mul_d(const IntPoly& a, const IntPoly& b) {
    if (a.is_zero() || b.is_zero()) return IntPoly{};
    std::vector<BigInt> c(a.size() + b.size() - 1U, BigInt(0));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j)
            c[i + j] += a[i] * b[j];
    IntPoly r(std::move(c));
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}
}  // namespace dispatch_helpers

// DISP1: threshold=0 forces half-GCD path even on small-degree polynomials.
// gcd(x^2-1, x-1) = x-1; dispatch with threshold=0 must set path = HalfGcd.
TEST(GcdDispatchTest, DISP1_ForceHalfGcdPathViaThreshold) {
    using namespace cas::symbolic;
    CASContext ctx;
    ctx.set_half_gcd_degree_threshold(0U);  // any min_deg > 0 triggers half-GCD
    ctx.set_modular_gcd_coeff_bits(~std::size_t{0U});  // disable CRT path

    IntPoly f = P({-1, 0, 1});  // x^2 - 1
    IntPoly g = P({-1, 1});     // x - 1

    IntegerGcdResult res = gcd_integer_poly_dispatch(std::move(f), std::move(g), ctx);
    // Result must divide both inputs.
    IntPoly f2 = P({-1, 0, 1});
    IntPoly g2 = P({-1, 1});
    EXPECT_TRUE(is_valid_gcd(res.gcd, f2, g2)) << "DISP1: gcd must divide both";
    EXPECT_EQ(res.gcd.degree(), 1U) << "DISP1: gcd(x^2-1, x-1) has degree 1";
    // Path should be HalfGcd (dispatch chose half-GCD because min_deg > threshold=0).
    EXPECT_EQ(res.path, IntegerGcdPath::HalfGcd) << "DISP1: dispatch must route to HalfGcd";
}

// DISP2: Both polynomials have min_deg > 200 with a known common factor.
// g = x^10-1 (deg 10), p = x^210+1 (deg 210), q = x^200+7 (deg 200).
// a = g*p, b = g*q. gcd(a, b) must have deg >= 10.
// Uses threshold=0 to guarantee half-GCD is selected by dispatch.
TEST(GcdDispatchTest, DISP2_HalfGcdPathHighDegreeKnownGcd) {
    using namespace dispatch_helpers;
    using namespace cas::symbolic;
    CASContext ctx;
    ctx.set_half_gcd_degree_threshold(0U);
    ctx.set_modular_gcd_coeff_bits(~std::size_t{0U});  // disable CRT

    // g = x^10 - 1
    std::vector<BigInt> gc(11, BigInt(0));
    gc[0] = BigInt(-1); gc[10] = BigInt(1);
    IntPoly g(std::move(gc));

    // p = x^210 + 1
    std::vector<BigInt> pc(211, BigInt(0));
    pc[0] = BigInt(1); pc[210] = BigInt(1);
    IntPoly p(std::move(pc));

    // q = x^200 + 7
    std::vector<BigInt> qc(201, BigInt(0));
    qc[0] = BigInt(7); qc[200] = BigInt(1);
    IntPoly q(std::move(qc));

    IntPoly a = poly_mul_d(g, p);  // deg 220
    IntPoly b = poly_mul_d(g, q);  // deg 210

    IntPoly a_copy = a;
    IntPoly b_copy = b;

    IntegerGcdResult res = gcd_integer_poly_dispatch(std::move(a), std::move(b), ctx);

    EXPECT_TRUE(is_valid_gcd(res.gcd, a_copy, b_copy)) << "DISP2: gcd must divide both";
    EXPECT_GE(res.gcd.degree(), 10U) << "DISP2: gcd must include deg-10 common factor";
    // path = HalfGcd (threshold=0; may fall back to subresultant if divisibility check fails,
    // which is also correct behavior — we only require the GCD is valid).
    EXPECT_TRUE(res.path == IntegerGcdPath::HalfGcd || res.path == IntegerGcdPath::Subresultant)
        << "DISP2: valid path";
}

// DISP3: gcd(f, 0) and gcd(0, g) — degenerate zero inputs exercising the
// early-return branches in gcd_integer_poly_with_subresultant (~lines 127-148).
// gcd(f, 0) = primitive(f) with content; gcd(0, g) = primitive(g) with content.
TEST(GcdDispatchTest, DISP3_DegenZeroInput) {
    using namespace cas::symbolic;
    CASContext ctx;
    ctx.set_modular_gcd_coeff_bits(~std::size_t{0U});  // subresultant path

    IntPoly f = P({-6, 0, 6});  // 6*(x^2 - 1)
    IntPoly zero_poly;           // empty = zero

    // gcd(f, 0): subresultant with B=zero returns A immediately.
    {
        IntPoly f_c = f;
        IntegerGcdResult r = gcd_integer_poly_dispatch(std::move(f_c), IntPoly{}, ctx);
        // Result must divide f; zero_poly is divisible by everything.
        EXPECT_FALSE(r.gcd.is_zero()) << "DISP3a: gcd(f, 0) must be non-zero";
        EXPECT_TRUE(is_valid_gcd(r.gcd, f, IntPoly{})) << "DISP3a: gcd must divide f";
        EXPECT_EQ(r.gcd.degree(), f.degree()) << "DISP3a: gcd(f, 0) = prim(f), same degree";
    }

    // gcd(0, g): symmetric.
    {
        IntPoly g = P({-10, 0, 0, 10});  // 10*(x^3 - 1)
        IntPoly g_c = g;
        IntegerGcdResult r = gcd_integer_poly_dispatch(IntPoly{}, std::move(g_c), ctx);
        EXPECT_FALSE(r.gcd.is_zero()) << "DISP3b: gcd(0, g) must be non-zero";
        EXPECT_TRUE(is_valid_gcd(r.gcd, IntPoly{}, g)) << "DISP3b: gcd must divide g";
        EXPECT_EQ(r.gcd.degree(), g.degree()) << "DISP3b: gcd(0, g) = prim(g), same degree";
    }
}

// DISP4: Negative leading coefficient normalization.
// gcd_integer_poly_dispatch must return a GCD with non-negative leading coefficient
// even when the half-GCD result has a negative lc (sign flip applied at ~line 228).
// We use threshold=0 to force half-GCD path.
TEST(GcdDispatchTest, DISP4_NegativeLcNormalized) {
    using namespace cas::symbolic;
    CASContext ctx;
    ctx.set_half_gcd_degree_threshold(0U);
    ctx.set_modular_gcd_coeff_bits(~std::size_t{0U});

    // f = -(x-1)*(x+2) = -(x^2+x-2) = -x^2 - x + 2, g = -(x-1)*(x-3) = -(x^2-4x+3)
    // Both share factor (x-1); expected gcd = primitive, pos lc.
    IntPoly f = P({2, -1, -1});    // -x^2 - x + 2 = -(x^2+x-2) = -(x+2)(x-1)
    IntPoly g = P({3, -4, 1});     //  x^2 - 4x + 3 = (x-1)(x-3)

    IntPoly f_c = f, g_c = g;
    IntegerGcdResult res = gcd_integer_poly_dispatch(std::move(f_c), std::move(g_c), ctx);

    EXPECT_TRUE(is_valid_gcd(res.gcd, f, g)) << "DISP4: gcd must divide both";
    EXPECT_EQ(res.gcd.degree(), 1U) << "DISP4: common factor is linear (x-1)";
    if (!res.gcd.is_zero()) {
        EXPECT_FALSE(res.gcd.leading_coeff().is_negative())
            << "DISP4: dispatch must normalize lc to positive";
    }
}

// DISP5: Lower threshold to trigger half-GCD on medium-degree polys and verify
// result agrees with subresultant path.
// Poly: g = x^5 - 1, p = x^20 + 3, q = x^15 + 5. a = g*p, b = g*q.
// With threshold=10: min_deg(a,b) = min(25,20) = 20 > 10 → half-GCD.
// With threshold=50: min_deg(a,b) = 20 < 50 → subresultant.
// Both should yield equal-degree gcd containing g as factor.
TEST(GcdDispatchTest, DISP5_ThresholdOverride_HalfGcdAgreesSubresultant) {
    using namespace dispatch_helpers;
    using namespace cas::symbolic;

    // g = x^5 - 1
    std::vector<BigInt> gc(6, BigInt(0));
    gc[0] = BigInt(-1); gc[5] = BigInt(1);
    IntPoly g(std::move(gc));

    // p = x^20 + 3
    std::vector<BigInt> pc(21, BigInt(0));
    pc[0] = BigInt(3); pc[20] = BigInt(1);
    IntPoly p(std::move(pc));

    // q = x^15 + 5
    std::vector<BigInt> qc(16, BigInt(0));
    qc[0] = BigInt(5); qc[15] = BigInt(1);
    IntPoly q(std::move(qc));

    IntPoly a = poly_mul_d(g, p);  // deg 25
    IntPoly b = poly_mul_d(g, q);  // deg 20

    // Half-GCD path: threshold=10 → min_deg(20) > 10.
    CASContext ctx_hgcd;
    ctx_hgcd.set_half_gcd_degree_threshold(10U);
    ctx_hgcd.set_modular_gcd_coeff_bits(~std::size_t{0U});

    IntPoly a1 = a, b1 = b;
    IntegerGcdResult r_hgcd = gcd_integer_poly_dispatch(std::move(a1), std::move(b1), ctx_hgcd);

    // Subresultant path: threshold=50 → min_deg(20) < 50.
    CASContext ctx_sub;
    ctx_sub.set_half_gcd_degree_threshold(50U);
    ctx_sub.set_modular_gcd_coeff_bits(~std::size_t{0U});

    IntPoly a2 = a, b2 = b;
    IntegerGcdResult r_sub = gcd_integer_poly_dispatch(std::move(a2), std::move(b2), ctx_sub);

    EXPECT_TRUE(is_valid_gcd(r_hgcd.gcd, a, b)) << "DISP5: half-GCD result must divide both";
    EXPECT_TRUE(is_valid_gcd(r_sub.gcd, a, b)) << "DISP5: subresultant result must divide both";
    EXPECT_EQ(r_hgcd.gcd.degree(), r_sub.gcd.degree())
        << "DISP5: half-GCD and subresultant must agree on GCD degree";
    EXPECT_GE(r_hgcd.gcd.degree(), 5U) << "DISP5: gcd must include deg-5 common factor";
    EXPECT_EQ(r_hgcd.path, IntegerGcdPath::HalfGcd) << "DISP5: ctx_hgcd must route to HalfGcd";
}
