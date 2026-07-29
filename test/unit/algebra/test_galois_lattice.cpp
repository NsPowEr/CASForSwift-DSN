// A6 / CAS-L3-18 — Exhaustive transitive-subgroup lattice tests.
//
// The engine *generates* the lattice; these tests cross-check it against
// classical, independently-published facts (transitive-class counts and
// orders — Butler & McKay 1983). The literature is used only here, as an
// oracle — never inside the engine.

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "../../../src/algebra/perm_group_internal.hpp"
#include "cas/error.hpp"

using namespace cas;
using namespace cas::algebra::permgrp;

namespace {

[[nodiscard]] std::vector<std::uint64_t> orders_of(
    const std::vector<PermGroup>& v) {
    std::vector<std::uint64_t> o;
    o.reserve(v.size());
    for (const auto& g : v) o.push_back(g.order());
    return o;
}

TEST(GaloisLatticeTest, Degree4_FiveTransitiveClasses) {
    // Classical: C4, V4, D4, A4, S4 — orders 4, 4, 8, 12, 24.
    auto r = transitive_subgroup_classes(
        4U, LatticeBudget{.max_degree = 7U, .max_ops = 100'000'000ULL});
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    const auto ords = orders_of(r.value());
    EXPECT_EQ(ords, (std::vector<std::uint64_t>{4U, 4U, 8U, 12U, 24U}));
    for (const auto& g : r.value()) EXPECT_TRUE(g.is_transitive());
    // C4 and V4 share order 4 but are non-conjugate (distinct cycle types).
    EXPECT_FALSE(conjugate_in_sn(r.value()[0], r.value()[1]));
}

TEST(GaloisLatticeTest, Degree5_FiveTransitiveClasses) {
    // Classical: C5, D5, F20, A5, S5 — orders 5, 10, 20, 60, 120.
    auto r = transitive_subgroup_classes(
        5U, LatticeBudget{.max_degree = 7U, .max_ops = 400'000'000ULL});
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    const auto ords = orders_of(r.value());
    EXPECT_EQ(ords, (std::vector<std::uint64_t>{5U, 10U, 20U, 60U, 120U}));
}

TEST(GaloisLatticeTest, Degree6_SixteenTransitiveClasses) {
    // Classical count (Butler-McKay): 16 transitive classes in S_6.
    auto r = transitive_subgroup_classes(
        6U, LatticeBudget{.max_degree = 7U, .max_ops = 4'000'000'000ULL});
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value().size(), 16U);
    // Contains both transitive order-120 classes: PGL(2,5) (exotic S_5)
    // and S_5 does not act transitively on 6 points in its natural action,
    // so exactly one order-120 class here is PGL(2,5); and S_6 (720) tops.
    const auto ords = orders_of(r.value());
    EXPECT_EQ(std::count(ords.begin(), ords.end(), 120U), 1);
    EXPECT_EQ(ords.back(), 720U);
    EXPECT_EQ(std::count(ords.begin(), ords.end(), 360U), 1);  // A_6
}

TEST(GaloisLatticeTest, Degree7_SevenTransitiveClasses) {
    // Classical (Butler-McKay): exactly 7 transitive classes in S_7, with
    // pairwise-distinct orders 7, 14, 21, 42, 168, 2520, 5040 — i.e.
    // C7, D7, F21, F42, PSL(3,2), A7, S7. Used here as an independent oracle;
    // the engine derives the lattice from first principles.
    auto r = transitive_subgroup_classes(
        7U, LatticeBudget{.max_degree = 7U, .max_ops = 40'000'000'000ULL});
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value().size(), 7U);
    EXPECT_EQ(orders_of(r.value()),
              (std::vector<std::uint64_t>{7U, 14U, 21U, 42U, 168U, 2520U,
                                          5040U}));
    for (const auto& g : r.value()) EXPECT_TRUE(g.is_transitive());
}

TEST(GaloisLatticeTest, DegreeBeyondCapIsStructuredUnimplemented) {
    auto r = transitive_subgroup_classes(
        9U, LatticeBudget{.max_degree = 7U, .max_ops = 1'000'000ULL});
    ASSERT_TRUE(r.is_error());
    EXPECT_EQ(r.error().kind, CASErrorKind::Unimplemented);
}

TEST(GaloisLatticeTest, ZeroOpsBudgetRejected) {
    auto r = transitive_subgroup_classes(
        4U, LatticeBudget{.max_degree = 7U, .max_ops = 0ULL});
    ASSERT_TRUE(r.is_error());
    EXPECT_EQ(r.error().kind, CASErrorKind::InvalidArgument);
}

TEST(GaloisLatticeTest, OpsBudgetExhaustionIsStructured) {
    auto r = transitive_subgroup_classes(
        5U, LatticeBudget{.max_degree = 7U, .max_ops = 10ULL});
    ASSERT_TRUE(r.is_error());
    EXPECT_EQ(r.error().kind, CASErrorKind::Unimplemented);
}

}  // namespace
