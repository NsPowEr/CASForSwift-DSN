// A6 — Constructor-family tests, Brick 2 of the deg ≥ 8 Stauduhar closure.
//
// Every group is BUILT by the engine from mathematical structure; the tests
// validate it against independently-derivable order formulas
// (|S_a≀S_b| = (a!)^b·b!, |AGL(d,p)| = p^d·∏(p^d−p^i), |PGL(2,q)| = q³−q,
// |PΓL(2,q)| = e·(q³−q), …), parity facts and primitivity — the literature
// appears here only as an oracle, never inside the engine (REGOLA 0.1).

#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "../../../src/algebra/perm_blocks_internal.hpp"
#include "../../../src/algebra/perm_bsgs_internal.hpp"
#include "../../../src/algebra/perm_construct_fields_internal.hpp"
#include "../../../src/algebra/perm_construct_internal.hpp"
#include "../../../src/algebra/perm_group_internal.hpp"
#include "cas/error.hpp"

using namespace cas;
using namespace cas::algebra::permgrp;

namespace {

[[nodiscard]] BsgsGroup group_of(std::size_t n, std::vector<Perm> gens) {
    auto g = BsgsGroup::build(n, std::move(gens));
    EXPECT_TRUE(g.is_ok());
    return std::move(g.value());
}

[[nodiscard]] std::uint64_t agl_order(std::size_t d, std::uint64_t p) {
    std::uint64_t pd = 1U;
    for (std::size_t i = 0U; i < d; ++i) pd *= p;
    std::uint64_t o = pd;  // translations
    std::uint64_t pi = 1U;
    for (std::size_t i = 0U; i < d; ++i) {
        o *= (pd - pi);
        pi *= p;
    }
    return o;
}

[[nodiscard]] std::uint64_t wreath_order(std::size_t a, std::size_t b) {
    std::uint64_t o = factorial_u64(b);
    for (std::size_t j = 0U; j < b; ++j) o *= factorial_u64(a);
    return o;
}

}  // namespace

// ── combinatorial families ───────────────────────────────────────────────────

TEST(PermConstructTest, SymmetricAndAlternatingOrders) {
    for (std::size_t n = 3U; n <= 10U; ++n) {
        EXPECT_EQ(group_of(n, symmetric_gens(n)).order(), factorial_u64(n));
        EXPECT_EQ(group_of(n, alternating_gens(n)).order(),
                  factorial_u64(n) / 2U);
        EXPECT_FALSE(group_of(n, alternating_gens(n)).has_odd_element());
    }
}

TEST(PermConstructTest, WreathOrdersAndTransitivity) {
    const std::vector<std::pair<std::size_t, std::size_t>> cases{
        {2U, 2U}, {2U, 3U}, {3U, 2U}, {2U, 4U},
        {4U, 2U}, {3U, 3U}, {2U, 5U}, {5U, 2U}};
    for (const auto& [a, b] : cases) {
        auto gens = wreath_gens(a, b);
        ASSERT_TRUE(gens.is_ok());
        const BsgsGroup w = group_of(a * b, gens.value());
        EXPECT_EQ(w.order(), wreath_order(a, b)) << "S_" << a << " wr S_" << b;
        EXPECT_TRUE(w.is_transitive());
        EXPECT_TRUE(w.has_odd_element());  // an in-block transposition is odd
    }
}

TEST(PermConstructTest, EvenPartHalvesOddGroupsAndFixesEvenOnes) {
    // S_5 ∩ A_5 = A_5.
    auto ev_s5 = even_part_gens(5U, symmetric_gens(5U));
    ASSERT_TRUE(ev_s5.is_ok());
    const BsgsGroup a5 = group_of(5U, ev_s5.value());
    EXPECT_EQ(a5.order(), 60U);
    EXPECT_FALSE(a5.has_odd_element());
    // (S_2 ≀ S_4) ∩ A_8 has order 384/2 = 192 and stays transitive.
    auto w = wreath_gens(2U, 4U);
    ASSERT_TRUE(w.is_ok());
    auto ev_w = even_part_gens(8U, w.value());
    ASSERT_TRUE(ev_w.is_ok());
    const BsgsGroup ew = group_of(8U, ev_w.value());
    EXPECT_EQ(ew.order(), 192U);
    EXPECT_FALSE(ew.has_odd_element());
    EXPECT_TRUE(ew.is_transitive());
    // An even group is returned unchanged (as a generating set).
    auto ev_a5 = even_part_gens(5U, alternating_gens(5U));
    ASSERT_TRUE(ev_a5.is_ok());
    EXPECT_EQ(group_of(5U, ev_a5.value()).order(), 60U);
}

TEST(PermConstructTest, KSubsetActionPetersenAndImprimitiveCase) {
    EXPECT_EQ(binomial_u64(5U, 2U), 10U);
    EXPECT_EQ(binomial_u64(10U, 5U), 252U);
    EXPECT_EQ(binomial_u64(20U, 10U), 184756U);

    // S_5 on the 10 pairs (Petersen action): order 120, primitive, odd.
    auto s5p = ksubset_action_gens(5U, 2U, symmetric_gens(5U));
    ASSERT_TRUE(s5p.is_ok());
    const BsgsGroup g = group_of(10U, s5p.value());
    EXPECT_EQ(g.order(), 120U);
    EXPECT_TRUE(g.is_transitive());
    EXPECT_TRUE(g.has_odd_element());
    auto prim = is_primitive(10U, s5p.value());
    ASSERT_TRUE(prim.is_ok());
    EXPECT_TRUE(prim.value());

    // A_5 on the pairs: order 60, even.
    auto a5p = ksubset_action_gens(5U, 2U, alternating_gens(5U));
    ASSERT_TRUE(a5p.is_ok());
    EXPECT_EQ(group_of(10U, a5p.value()).order(), 60U);
    EXPECT_FALSE(group_of(10U, a5p.value()).has_odd_element());

    // S_4 on the 6 pairs is imprimitive (complementary-pair blocks) and even.
    auto s4p = ksubset_action_gens(4U, 2U, symmetric_gens(4U));
    ASSERT_TRUE(s4p.is_ok());
    EXPECT_EQ(group_of(6U, s4p.value()).order(), 24U);
    EXPECT_FALSE(group_of(6U, s4p.value()).has_odd_element());
    auto prim4 = is_primitive(6U, s4p.value());
    ASSERT_TRUE(prim4.is_ok());
    EXPECT_FALSE(prim4.value());

    // Representation bound: C(12,6) = 924 > 255 fails structured.
    auto too_big = ksubset_action_gens(12U, 6U, symmetric_gens(12U));
    ASSERT_FALSE(too_big.is_ok());
    EXPECT_EQ(too_big.error().kind, CASErrorKind::Unimplemented);
}

// ── affine and projective families over prime fields ─────────────────────────

TEST(PermConstructTest, AffineGroupOrdersParityPrimitivity) {
    struct Case {
        std::size_t d;
        std::uint32_t p;
        bool odd;
    };
    // AGL(1,5) = F_20 (odd), AGL(1,7) = F_42 (odd), AGL(3,2) ⊆ A_8 (even),
    // AGL(2,3) ⊄ A_9 (odd) — parity read off the constructed generators.
    const std::vector<Case> cases{
        {1U, 5U, true}, {1U, 7U, true}, {3U, 2U, false}, {2U, 3U, true}};
    for (const auto& c : cases) {
        auto gens = agl_gens(c.d, c.p);
        ASSERT_TRUE(gens.is_ok()) << gens.error().message;
        std::size_t npts = 1U;
        for (std::size_t i = 0U; i < c.d; ++i) npts *= c.p;
        const BsgsGroup g = group_of(npts, gens.value());
        EXPECT_EQ(g.order(), agl_order(c.d, c.p))
            << "AGL(" << c.d << "," << c.p << ")";
        EXPECT_TRUE(g.is_transitive());
        EXPECT_EQ(g.has_odd_element(), c.odd);
        auto prim = is_primitive(npts, gens.value());
        ASSERT_TRUE(prim.is_ok());
        EXPECT_TRUE(prim.value());  // affine groups are 2-transitive
    }
    EXPECT_FALSE(agl_gens(2U, 4U).is_ok());   // 4 is not prime
    EXPECT_FALSE(agl_gens(9U, 2U).is_ok());   // 2^9 exceeds the image bound
}

TEST(PermConstructTest, ProjectiveGlOnP2F2AndP2F3) {
    // PGL(3,2) = PSL(3,2) on the 7 points of P²(F_2): order 168, even.
    auto p32 = projective_gl_gens(3U, 2U);
    ASSERT_TRUE(p32.is_ok()) << p32.error().message;
    const BsgsGroup g = group_of(7U, p32.value());
    EXPECT_EQ(g.order(), 168U);
    EXPECT_TRUE(g.is_transitive());
    EXPECT_FALSE(g.has_odd_element());
    auto prim = is_primitive(7U, p32.value());
    ASSERT_TRUE(prim.is_ok());
    EXPECT_TRUE(prim.value());

    // PGL(3,3) on the 13 points of P²(F_3): |GL(3,3)|/(3−1) = 5616.
    auto p33 = projective_gl_gens(3U, 3U);
    ASSERT_TRUE(p33.is_ok()) << p33.error().message;
    EXPECT_EQ(group_of(13U, p33.value()).order(), 5616U);
}

// ── GF(q) construction ────────────────────────────────────────────────────────

namespace {

void check_field_axioms(const Gf& gf) {
    const std::uint32_t q = gf.q;
    auto mul = [&](std::uint32_t a, std::uint32_t b) {
        return gf.mul[static_cast<std::size_t>(a) * q + b];
    };
    for (std::uint32_t a = 0U; a < q; ++a) {
        EXPECT_EQ(gf.add(a, 0U), a);
        EXPECT_EQ(mul(a, 1U), a);
        EXPECT_EQ(gf.add(a, gf.neg[a]), 0U);
        if (a != 0U) EXPECT_EQ(mul(a, gf.inv[a]), 1U);
        for (std::uint32_t b = 0U; b < q; ++b) {
            EXPECT_EQ(gf.add(a, b), gf.add(b, a));
            EXPECT_EQ(mul(a, b), mul(b, a));
            // Frobenius is a ring homomorphism.
            EXPECT_EQ(gf.frob[gf.add(a, b)], gf.add(gf.frob[a], gf.frob[b]));
            EXPECT_EQ(gf.frob[mul(a, b)], mul(gf.frob[a], gf.frob[b]));
            for (std::uint32_t c = 0U; c < q; ++c) {
                EXPECT_EQ(mul(a, gf.add(b, c)),
                          gf.add(mul(a, b), mul(a, c)));
                EXPECT_EQ(mul(mul(a, b), c), mul(a, mul(b, c)));
            }
        }
    }
    // γ generates F_q^*.
    std::uint32_t x = gf.gamma;
    std::uint32_t ord = 1U;
    while (x != 1U) {
        x = gf.mul[static_cast<std::size_t>(x) * q + gf.gamma];
        ++ord;
    }
    EXPECT_EQ(ord, q - 1U);
}

}  // namespace

TEST(PermConstructTest, GfFieldAxiomsExhaustive) {
    for (const auto& [p, e] : std::vector<std::pair<std::uint32_t, std::size_t>>{
             {2U, 2U}, {2U, 3U}, {3U, 2U}, {5U, 1U}}) {
        auto gf = build_gf(p, e);
        ASSERT_TRUE(gf.is_ok()) << gf.error().message;
        check_field_axioms(gf.value());
    }
    EXPECT_FALSE(build_gf(4U, 2U).is_ok());   // 4 not prime
    EXPECT_FALSE(build_gf(2U, 8U).is_ok());   // q+1 > 255
}

// ── Möbius / semilinear families on P¹(F_q) ──────────────────────────────────

TEST(PermConstructTest, MoebiusOrdersAcrossFields) {
    struct Case {
        std::uint32_t p;
        std::size_t e;
        std::uint64_t pgl;
        std::uint64_t psl;
    };
    // |PGL(2,q)| = q³−q; |PSL(2,q)| = (q³−q)/gcd(2,q−1).
    const std::vector<Case> cases{{5U, 1U, 120U, 60U},
                                  {7U, 1U, 336U, 168U},
                                  {2U, 3U, 504U, 504U},
                                  {3U, 2U, 720U, 360U}};
    for (const auto& c : cases) {
        auto gf = build_gf(c.p, c.e);
        ASSERT_TRUE(gf.is_ok());
        const std::size_t n = gf.value().q + 1U;
        const BsgsGroup pgl = group_of(n, pgl2_point_gens(gf.value()));
        EXPECT_EQ(pgl.order(), c.pgl) << "PGL(2," << gf.value().q << ")";
        EXPECT_TRUE(pgl.is_transitive());
        const BsgsGroup psl = group_of(n, psl2_point_gens(gf.value()));
        EXPECT_EQ(psl.order(), c.psl) << "PSL(2," << gf.value().q << ")";
        EXPECT_FALSE(psl.has_odd_element());  // PSL(2,q) ⊆ A_{q+1} always
    }
}

TEST(PermConstructTest, SemilinearLayerOverF9) {
    auto gf_r = build_gf(3U, 2U);
    ASSERT_TRUE(gf_r.is_ok());
    const Gf& gf = gf_r.value();
    std::vector<Perm> psl = psl2_point_gens(gf);
    const Perm gmul = gamma_mult_point_perm(gf);
    const Perm frob = frobenius_point_perm(gf);

    // The three index-2 subgroups above PSL(2,9), all of order 720, pairwise
    // distinct: ⟨PSL, γ⟩ = PGL (odd), ⟨PSL, φ⟩ (odd), ⟨PSL, γφ⟩ (even —
    // the point-stabiliser-sharply-3-transitive one), and PΓL(2,9) = 1440.
    std::vector<Perm> pgl = psl;
    pgl.push_back(gmul);
    std::vector<Perm> psigma = psl;
    psigma.push_back(frob);
    std::vector<Perm> twisted = psl;
    twisted.push_back(compose(gmul, frob));
    std::vector<Perm> pgamma = psl;
    pgamma.push_back(gmul);
    pgamma.push_back(frob);

    const BsgsGroup g_pgl = group_of(10U, pgl);
    const BsgsGroup g_psigma = group_of(10U, psigma);
    const BsgsGroup g_twist = group_of(10U, twisted);
    const BsgsGroup g_pgamma = group_of(10U, pgamma);
    EXPECT_EQ(g_pgl.order(), 720U);
    EXPECT_EQ(g_psigma.order(), 720U);
    EXPECT_EQ(g_twist.order(), 720U);
    EXPECT_EQ(g_pgamma.order(), 1440U);
    EXPECT_TRUE(g_pgl.has_odd_element());
    EXPECT_TRUE(g_psigma.has_odd_element());
    EXPECT_FALSE(g_twist.has_odd_element());
    // Pairwise distinct as subgroups: each contains its defining coset rep,
    // not the others'.
    EXPECT_FALSE(g_pgl.contains(frob));
    EXPECT_FALSE(g_psigma.contains(gmul));
    EXPECT_FALSE(g_twist.contains(gmul));
    EXPECT_FALSE(g_twist.contains(frob));
    // PΓL(2,8) = ⟨PSL(2,8), φ⟩, order 3·504 = 1512, even (⊆ A_9).
    auto gf8_r = build_gf(2U, 3U);
    ASSERT_TRUE(gf8_r.is_ok());
    std::vector<Perm> pgamma8 = psl2_point_gens(gf8_r.value());
    pgamma8.push_back(frobenius_point_perm(gf8_r.value()));
    const BsgsGroup g8 = group_of(9U, pgamma8);
    EXPECT_EQ(g8.order(), 1512U);
    EXPECT_FALSE(g8.has_odd_element());
}
