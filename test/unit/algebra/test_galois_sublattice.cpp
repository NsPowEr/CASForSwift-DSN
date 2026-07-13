// A6 Brick 3.5 — maximal transitive subgroup classes inside a small node.
//
// Oracles are elementary group theory, never transcribed tables:
//   • S₄: the proper transitive subgroups are A₄, D₄ (3 conjugates), C₄
//     and the normal V₄; C₄ < D₄ and V₄ < A₄, so the maximal transitive
//     classes are exactly {A₄ (12), D₄ (8)}.
//   • D₄ = ⟨(0123),(02)⟩: transitive proper subgroups are C₄ and the
//     normal Klein group (the point-stabilising Klein group {e,(02),(13),
//     (02)(13)} has orbit {0,2} — intransitive); both have index 2, so
//     both are maximal. C₄ contains the odd 4-cycle, the Klein group is
//     even — a parity distinguisher with no table.
//   • S₅: A₅ and F₂₀ = AGL(1,5) are the maximal transitive classes
//     (C₅ < D₅ < F₂₀); coverage is additionally cross-checked against the
//     INDEPENDENT dense lattice engine (transitive_subgroup_classes):
//     every proper transitive class of S₅ must lie in a conjugate of a
//     returned class — the exact contract the Stauduhar walk relies on.

#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <vector>

#include "../../../src/algebra/galois_sublattice_internal.hpp"
#include "../../../src/algebra/perm_bsgs_internal.hpp"
#include "../../../src/algebra/perm_group_internal.hpp"
#include "cas/error.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::algebra;
using namespace cas::algebra::permgrp;

namespace {

class GaloisSublatticeTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] static Perm make_perm(std::initializer_list<int> img) {
        Perm p;
        p.reserve(img.size());
        for (const int v : img) p.push_back(static_cast<std::uint8_t>(v));
        return p;
    }

    [[nodiscard]] static BsgsGroup symmetric(std::size_t n) {
        Perm ncyc = identity(n);
        for (std::size_t i = 0U; i < n; ++i) {
            ncyc[i] = static_cast<std::uint8_t>((i + 1U) % n);
        }
        Perm swap01 = identity(n);
        swap01[0] = 1U;
        swap01[1] = 0U;
        auto g = BsgsGroup::build(n, {swap01, ncyc});
        EXPECT_TRUE(g.is_ok());
        return g.value();
    }

    [[nodiscard]] std::vector<BsgsGroup> classes_of(const BsgsGroup& h) {
        auto r = transitive_subgroup_classes_in(
            h, ctx.galois_lattice_max_ops(),
            ctx.galois_sublattice_max_bytes(), &ctx);
        EXPECT_TRUE(r.is_ok()) << (r.is_error() ? r.error().message : "");
        return r.is_ok() ? r.value() : std::vector<BsgsGroup>{};
    }
};

TEST_F(GaloisSublatticeTest, S4MaximalTransitiveClasses) {
    const auto classes = classes_of(symmetric(4U));
    ASSERT_EQ(classes.size(), 2U);
    EXPECT_EQ(classes[0].order(), 12U);  // A4
    EXPECT_EQ(classes[1].order(), 8U);   // D4
    for (const auto& c : classes) EXPECT_TRUE(c.is_transitive());
}

TEST_F(GaloisSublatticeTest, D4MaximalTransitiveClasses) {
    auto d4 = BsgsGroup::build(
        4U, {make_perm({1, 2, 3, 0}), make_perm({2, 1, 0, 3})});
    ASSERT_TRUE(d4.is_ok());
    ASSERT_EQ(d4.value().order(), 8U);
    const auto classes = classes_of(d4.value());
    ASSERT_EQ(classes.size(), 2U);
    EXPECT_EQ(classes[0].order(), 4U);
    EXPECT_EQ(classes[1].order(), 4U);
    for (const auto& c : classes) EXPECT_TRUE(c.is_transitive());
    // C4 has the odd 4-cycle; the transitive Klein group is even.
    EXPECT_NE(classes[0].has_odd_element(), classes[1].has_odd_element());
}

TEST_F(GaloisSublatticeTest, S5ClassesAndDenseLatticeCoverage) {
    const auto classes = classes_of(symmetric(5U));
    ASSERT_EQ(classes.size(), 2U);
    EXPECT_EQ(classes[0].order(), 60U);  // A5
    EXPECT_EQ(classes[1].order(), 20U);  // F20 = AGL(1,5)
    // Independent engine: every proper transitive class of S5 must be
    // contained in an S5-conjugate of a returned maximal class.
    auto dense = transitive_subgroup_classes(
        5U, LatticeBudget{.max_degree = 7U,
                          .max_ops = ctx.galois_lattice_max_ops()});
    ASSERT_TRUE(dense.is_ok());
    const std::uint64_t s5_order = factorial_u64(5U);
    for (const auto& k : dense.value()) {
        if (k.order() == s5_order) continue;  // S5 itself
        bool covered = false;
        for (const auto& c : classes) {
            if (k.order() > c.order() || c.order() % k.order() != 0U) {
                continue;
            }
            for (std::uint32_t r = 0U;
                 r < static_cast<std::uint32_t>(s5_order) && !covered; ++r) {
                const Perm s = lehmer_unrank(r, 5U);
                const Perm si = inverse(s);
                bool all_in = true;
                for (const auto& g : k.generators()) {
                    if (!c.contains(compose(si, compose(g, s)))) {
                        all_in = false;
                        break;
                    }
                }
                covered = all_in;
            }
            if (covered) break;
        }
        EXPECT_TRUE(covered) << "dense class of order " << k.order()
                             << " not covered by any maximal class";
    }
}

TEST_F(GaloisSublatticeTest, PrimeCycleHasNoProperTransitiveSubgroup) {
    auto c5 = BsgsGroup::build(5U, {make_perm({1, 2, 3, 4, 0})});
    ASSERT_TRUE(c5.is_ok());
    ASSERT_EQ(c5.value().order(), 5U);
    // Terminal certificate of the walk: an empty class list.
    EXPECT_TRUE(classes_of(c5.value()).empty());
}

TEST_F(GaloisSublatticeTest, StructuredFailures) {
    const auto s4 = symmetric(4U);
    // Zero budgets are a caller bug.
    auto zero_ops = transitive_subgroup_classes_in(s4, 0U, 1024U, &ctx);
    ASSERT_TRUE(zero_ops.is_error());
    EXPECT_EQ(zero_ops.error().kind, CASErrorKind::InvalidArgument);
    auto zero_bytes = transitive_subgroup_classes_in(s4, 1024U, 0U, &ctx);
    ASSERT_TRUE(zero_bytes.is_error());
    EXPECT_EQ(zero_bytes.error().kind, CASErrorKind::InvalidArgument);
    // Exhaustion is structured, never a truncated answer.
    auto tiny_ops = transitive_subgroup_classes_in(
        s4, 1U, ctx.galois_sublattice_max_bytes(), &ctx);
    ASSERT_TRUE(tiny_ops.is_error());
    EXPECT_EQ(tiny_ops.error().kind, CASErrorKind::Unimplemented);
    auto tiny_bytes = transitive_subgroup_classes_in(
        s4, ctx.galois_lattice_max_ops(), 8U, &ctx);
    ASSERT_TRUE(tiny_bytes.is_error());
    EXPECT_EQ(tiny_bytes.error().kind, CASErrorKind::Unimplemented);
    // Beyond the u16-index universe (|S9| = 362880 > 2^16): structured.
    auto big = transitive_subgroup_classes_in(
        symmetric(9U), ctx.galois_lattice_max_ops(),
        ctx.galois_sublattice_max_bytes(), &ctx);
    ASSERT_TRUE(big.is_error());
    EXPECT_EQ(big.error().kind, CASErrorKind::Unimplemented);
}

}  // namespace
