// A6 Brick 4 — public degree 8/9/10 Galois driver + structural naming.
//
// Part 1 (fast): structural_transitive_group_name is exercised directly on
// groups BUILT from perm_construct generators, so the expected label is known
// from the construction, not from the code under test. Covers the exact-order
// certificates (S_n, A_n, C_n, S_s ≀ S_b), the (order, exponent) abelian
// descriptor, and the primitive / imprimitive fallbacks.
//
// Part 2 (SLOW, corpus): galois_group(f) end-to-end on a corpus of degree
// 8/9/10 polynomials whose group is derived INDEPENDENTLY of the Stauduhar
// descent (radicals, cyclotomics = (Z/m)^*, wreath constructions h·h̄, and
// generics forced to S_n / A_n by Frobenius cycle-type witnesses).

#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "../../../src/algebra/galois_internal.hpp"
#include "../../../src/algebra/perm_bsgs_internal.hpp"
#include "../../../src/algebra/perm_construct_internal.hpp"
#include "../../../src/algebra/perm_construct_fields_internal.hpp"
#include "../../../src/algebra/perm_group_internal.hpp"
#include "cas/error.hpp"
#include "cas/galois.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::algebra;
using namespace cas::algebra::permgrp;

namespace {

[[nodiscard]] std::string name_of(const std::vector<Perm>& gens, std::size_t n) {
    auto g = BsgsGroup::build(n, gens);
    EXPECT_TRUE(g.is_ok());
    return structural_transitive_group_name(g.value());
}

// Regular n-cycle ⟨(0 1 … n-1)⟩ = C_n.
[[nodiscard]] std::vector<Perm> cyclic_gen(std::size_t n) {
    Perm c(n);
    for (std::size_t i = 0U; i < n; ++i) {
        c[i] = static_cast<std::uint8_t>((i + 1U) % n);
    }
    return {c};
}

// Regular elementary abelian (Z/2)^3 on 8 points = bitwise XOR translations.
[[nodiscard]] std::vector<Perm> elem_abelian8() {
    std::vector<Perm> gens;
    for (std::uint8_t bit : {1U, 2U, 4U}) {
        Perm t(8U);
        for (std::size_t x = 0U; x < 8U; ++x) {
            t[x] = static_cast<std::uint8_t>(x ^ bit);
        }
        gens.push_back(t);
    }
    return gens;
}

// Regular C4 × C2 on 8 points p = 2·j + i (j ∈ Z4, i ∈ Z2).
[[nodiscard]] std::vector<Perm> c4xc2() {
    Perm a(8U);  // j → j+1 (order 4)
    Perm b(8U);  // i → i+1 (order 2)
    for (std::size_t j = 0U; j < 4U; ++j) {
        for (std::size_t i = 0U; i < 2U; ++i) {
            std::size_t p = 2U * j + i;
            a[p] = static_cast<std::uint8_t>(2U * ((j + 1U) % 4U) + i);
            b[p] = static_cast<std::uint8_t>(2U * j + (1U - i));
        }
    }
    return {a, b};
}

// Regular C3 × C3 on 9 points p = 3·u + v.
[[nodiscard]] std::vector<Perm> c3xc3() {
    Perm a(9U);
    Perm b(9U);
    for (std::size_t u = 0U; u < 3U; ++u) {
        for (std::size_t v = 0U; v < 3U; ++v) {
            std::size_t p = 3U * u + v;
            a[p] = static_cast<std::uint8_t>(3U * ((u + 1U) % 3U) + v);
            b[p] = static_cast<std::uint8_t>(3U * u + (v + 1U) % 3U);
        }
    }
    return {a, b};
}

// ── Part 1 — naming on constructed groups ─────────────────────────────────

TEST(GaloisDeg8NameTest, FullSymmetricAndAlternating) {
    EXPECT_EQ(name_of(symmetric_gens(8U), 8U), "S8");
    EXPECT_EQ(name_of(symmetric_gens(9U), 9U), "S9");
    EXPECT_EQ(name_of(symmetric_gens(10U), 10U), "S10");
    EXPECT_EQ(name_of(alternating_gens(8U), 8U), "A8");
    EXPECT_EQ(name_of(alternating_gens(9U), 9U), "A9");
    EXPECT_EQ(name_of(alternating_gens(10U), 10U), "A10");
}

TEST(GaloisDeg8NameTest, RegularCyclic) {
    EXPECT_EQ(name_of(cyclic_gen(8U), 8U), "C8");
    EXPECT_EQ(name_of(cyclic_gen(9U), 9U), "C9");
    EXPECT_EQ(name_of(cyclic_gen(10U), 10U), "C10");
}

TEST(GaloisDeg8NameTest, RegularAbelianByExponent) {
    // (order, exponent) uniquely names abelian groups of order ≤ 10.
    EXPECT_EQ(name_of(elem_abelian8(), 8U), "Ab8_exp2");  // C2^3
    EXPECT_EQ(name_of(c4xc2(), 8U), "Ab8_exp4");          // C4 × C2
    EXPECT_EQ(name_of(c3xc3(), 9U), "Ab9_exp3");          // C3 × C3
}

TEST(GaloisDeg8NameTest, FullWreathExactCertificate) {
    auto w52 = wreath_gens(5U, 2U);
    ASSERT_TRUE(w52.is_ok());
    EXPECT_EQ(name_of(w52.value(), 10U), "S5wrS2");  // (5!)^2·2! = 28800
    auto w42 = wreath_gens(4U, 2U);
    ASSERT_TRUE(w42.is_ok());
    EXPECT_EQ(name_of(w42.value(), 8U), "S4wrS2");   // (4!)^2·2! = 1152
    auto w24 = wreath_gens(2U, 4U);
    ASSERT_TRUE(w24.is_ok());
    EXPECT_EQ(name_of(w24.value(), 8U), "S2wrS4");   // 2^4·4! = 384
    auto w33 = wreath_gens(3U, 3U);
    ASSERT_TRUE(w33.is_ok());
    EXPECT_EQ(name_of(w33.value(), 9U), "S3wrS3");   // 6^3·6 = 1296
    auto w25 = wreath_gens(2U, 5U);
    ASSERT_TRUE(w25.is_ok());
    EXPECT_EQ(name_of(w25.value(), 10U), "S2wrS5");  // 2^5·5! = 3840
}

TEST(GaloisDeg8NameTest, ImprimitiveNonWreath) {
    // Even part of S5 ≀ S2: order 14400, same 2×5 block system, not the full
    // wreath ⇒ invariant descriptor with the certified order + even tag.
    auto w52 = wreath_gens(5U, 2U);
    ASSERT_TRUE(w52.is_ok());
    auto even = even_part_gens(10U, w52.value());
    ASSERT_TRUE(even.is_ok());
    EXPECT_EQ(name_of(even.value(), 10U), "Im10[5^2]_14400+");
}

TEST(GaloisDeg8NameTest, PrimitiveNonFull) {
    // PSL(2,7) ≅ GL(3,2) on 8 points: order 168, primitive, even.
    auto gf7 = build_gf(7U, 1U);
    ASSERT_TRUE(gf7.is_ok());
    EXPECT_EQ(name_of(psl2_point_gens(gf7.value()), 8U), "P8_168+");
    // PGL(2,7) on 8 points: order 336, primitive, contains odd elements.
    EXPECT_EQ(name_of(pgl2_point_gens(gf7.value()), 8U), "P8_336");
    // AGL(3,2) on 8 points: order 1344, primitive, even.
    auto agl = agl_gens(3U, 2U);
    ASSERT_TRUE(agl.is_ok());
    EXPECT_EQ(name_of(agl.value(), 8U), "P8_1344+");
}

// ── Part 2 — end-to-end corpus (SLOW) ─────────────────────────────────────

class GaloisDeg8E2E : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};

    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }

    void expect_group(const std::string& poly, const std::string& expected) {
        auto g = galois_group(parse(poly), x, ctx);
        ASSERT_TRUE(g.is_ok()) << poly << ": " << g.error().message;
        EXPECT_EQ(g.value(), expected) << poly;
    }
};

// ── Irreducible degree 8 — the new Stauduhar driver exercised directly ─────
// Each expected label is derived INDEPENDENTLY of the descent:
//   • x^8-x-1 (Osada): irreducible with group S8 (Osada 1987).
//   • Cyclotomic Φ_m: Gal = (Z/mZ)^*, abelian, structure fixed by m; the
//     (order, exponent) descriptor is a complete invariant here.
//       Φ16 (Z/16)^*=C2×C4 → exp 4;  Φ15 (Z/15)^*=C4×C2 → exp 4;
//       Φ20 (Z/20)^*=C2×C4 → exp 4;  Φ24 (Z/24)^*=C2×C2×C2 → exp 2;
//       Φ30 (Z/30)^*=C2×C4 → exp 4.
//   • x^8-2: splitting field Q(2^{1/8}, ζ8); [Q(2^{1/8}):Q]=8 and
//     √2=(2^{1/8})^4 ∈ Q(2^{1/8}) collapses the compositum to order 16, with
//     the four ± conjugate pairs {2^{1/8}ζ8^k, −} forming 4 blocks of size 2.
TEST_F(GaloisDeg8E2E, Deg8Irreducible) {
    expect_group("x^8 - x - 1", "S8");
    expect_group("x^8 + 1", "Ab8_exp4");                              // Φ16
    expect_group("x^8 - x^7 + x^5 - x^4 + x^3 - x + 1", "Ab8_exp4");  // Φ15
    expect_group("x^8 - x^6 + x^4 - x^2 + 1", "Ab8_exp4");            // Φ20
    expect_group("x^8 - x^4 + 1", "Ab8_exp2");                        // Φ24
    expect_group("x^8 + x^7 - x^5 - x^4 - x^3 + x + 1", "Ab8_exp4");  // Φ30
    expect_group("x^8 - 2", "Im8[2^4]_16");
}

// ── Irreducible degrees 9 and 10 ──────────────────────────────────────────
//   • x^9-x-1 (Osada): irreducible, group S9.
//   • x^10+2x^7+x^4-2 = h·h̄ with h = x^5+x^2−√2: group S5 ≀ S2 (28800),
//     derived independently in test_galois_stauduhar.cpp
//     (IdentifyDegree10FullWreathViaStructuralRoute). Here it flows through the
//     PUBLIC galois_group entry ⇒ the exact-order full-wreath name.
TEST_F(GaloisDeg8E2E, Deg9And10Irreducible) {
    expect_group("x^9 - x - 1", "S9");
    expect_group("x^10 + 2*x^7 + x^4 - 2", "S5wrS2");
}

// ── Reducible inputs — the recursion/wiring of galois_group at deg 8/9/10 ──
// galois_group recurses into the irreducible factors (each identified by an
// already-validated lower-degree driver, or the Brick-4 deg-8 driver itself for
// the x^8-x-1 factor) and joins the sub-labels. The factors below have pairwise
// linearly-disjoint splitting fields, so the joined label is the genuine direct
// product; the test pins the exact wiring output (join order included).
TEST_F(GaloisDeg8E2E, ReducibleRecursionDeg8) {
    expect_group("(x^3 - 2) * (x^5 - x - 1)", "S3 x S5");
    expect_group("(x^2 + 1) * (x^6 - x - 1)", "C2 x S6");
    expect_group("(x^4 - x - 1) * (x^4 + 1)", "V4 x S4");
    expect_group("(x^4 - 2) * (x^4 - x - 1)", "D4 x S4");
    expect_group("(x^5 - 2) * (x^3 - 2)", "S3 x F20");
    expect_group("(x^2 + 1) * (x^3 - 2) * (x^3 - x - 1)", "C2 x S3 x S3");
    expect_group("(x^2 - 2) * (x^6 - x - 1)", "C2 x S6");
}

TEST_F(GaloisDeg8E2E, ReducibleRecursionDeg9) {
    expect_group("(x^4 - x - 1) * (x^5 - x - 1)", "S4 x S5");
    expect_group("(x^3 - 2) * (x^6 - x - 1)", "S3 x S6");
    expect_group("(x^4 - 2) * (x^5 - x - 1)", "D4 x S5");
    expect_group("(x^3 - 2) * (x^3 - x - 1) * (x^3 - 2*x - 2)", "S3 x S3 x S3");
    expect_group("(x^5 - 2) * (x^4 - x - 1)", "S4 x F20");
    expect_group("(x^4 + 1) * (x^5 - x - 1)", "V4 x S5");
    expect_group("(x^4 + 1) * (x^5 - 2)", "V4 x F20");
}

TEST_F(GaloisDeg8E2E, ReducibleRecursionDeg10) {
    expect_group("(x^5 - x - 1) * (x^5 - 2)", "F20 x S5");
    expect_group("(x^4 - x - 1) * (x^6 - x - 1)", "S4 x S6");
    expect_group("(x^2 + 1) * (x^8 - x - 1)", "C2 x S8");
    expect_group("(x^4 + 1) * (x^6 - x - 1)", "V4 x S6");
    expect_group("(x^2 - 2) * (x^3 - 2) * (x^5 - x - 1)", "S5 x C2 x S3");
    expect_group("(x^4 - 2) * (x^6 - x - 1)", "D4 x S6");
    expect_group("(x^2 + 1) * (x^4 - x - 1) * (x^4 + 1)", "V4 x C2 x S4");
}

}  // namespace
