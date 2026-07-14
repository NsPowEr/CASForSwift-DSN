// A6 Brick 3.75 — structural wreath-preimage maximal generator:
//   • detection of the wreath-preimage shape (standard, even part,
//     relabelled coordinates);
//   • the FA ∪ FB ∪ FD candidate lists of the three over-budget (5,2)
//     nodes, pinned against the hand-derived group theory of
//     galois_wreath_maximal_internal.hpp;
//   • machine verification of the two theorem steps on dense ground
//     truth (Scott's lemma on A₅ × A₅; the N_{S₅}(D) shell coverage);
//   • the new dense API variants (maximal / all classes).
// The end-to-end degree-10 walk through the structural route lives in
// test_galois_stauduhar.cpp (Degree10 tests).

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "../../../src/algebra/galois_sublattice_internal.hpp"
#include "../../../src/algebra/galois_wreath_maximal_internal.hpp"
#include "../../../src/algebra/perm_construct_internal.hpp"
#include "../../../src/algebra/perm_group_internal.hpp"
#include "../../../src/algebra/perm_maximal_internal.hpp"
#include "cas/error.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::algebra;
using namespace cas::algebra::permgrp;

namespace {

constexpr std::uint64_t kOps = 1ULL << 40U;
constexpr std::uint64_t kBytes = 1ULL << 28U;

[[nodiscard]] Perm embed5(const Perm& g, std::size_t block) {
    Perm w = identity(10U);
    for (std::size_t x = 0U; x < 5U; ++x) {
        w[block * 5U + x] = static_cast<std::uint8_t>(block * 5U + g[x]);
    }
    return w;
}

[[nodiscard]] Perm t01_5() {
    Perm t = identity(5U);
    t[0] = 1U;
    t[1] = 0U;
    return t;
}

[[nodiscard]] Perm block_swap10() {
    Perm s(10U);
    for (std::size_t x = 0U; x < 5U; ++x) {
        s[x] = static_cast<std::uint8_t>(x + 5U);
        s[x + 5U] = static_cast<std::uint8_t>(x);
    }
    return s;
}

[[nodiscard]] BsgsGroup must_build(std::size_t n, std::vector<Perm> gens) {
    auto b = BsgsGroup::build(n, std::move(gens));
    EXPECT_TRUE(b.is_ok());
    return std::move(b.value());
}

// The full standard wreath W = S₅ ≀ S₂ on 10 points (order 28800).
[[nodiscard]] BsgsGroup full_wreath52() {
    auto w = wreath_gens(5U, 2U);
    EXPECT_TRUE(w.is_ok());
    return must_build(10U, std::move(w.value()));
}

// ker χ₁ = {(b₁, b₂)σ : sgn b₁ · sgn b₂ = 1} ⋊ swap (order 14400):
// A₅ in both blocks, a both-odd base pair, the plain block swap.
[[nodiscard]] BsgsGroup ker_chi1() {
    std::vector<Perm> gens;
    for (const Perm& g : alternating_gens(5U)) {
        gens.push_back(embed5(g, 0U));
        gens.push_back(embed5(g, 1U));
    }
    gens.push_back(compose(embed5(t01_5(), 0U), embed5(t01_5(), 1U)));
    gens.push_back(block_swap10());
    return must_build(10U, std::move(gens));
}

// W ∩ A₁₀ = ker(χ₁χ₂) for d = 5 (order 14400): A₅ in both blocks, a
// both-odd base pair, a one-block-odd decorated swap (even on Ω).
[[nodiscard]] BsgsGroup even_wreath52() {
    std::vector<Perm> gens;
    for (const Perm& g : alternating_gens(5U)) {
        gens.push_back(embed5(g, 0U));
        gens.push_back(embed5(g, 1U));
    }
    gens.push_back(compose(embed5(t01_5(), 0U), embed5(t01_5(), 1U)));
    gens.push_back(compose(embed5(t01_5(), 1U), block_swap10()));
    return must_build(10U, std::move(gens));
}

[[nodiscard]] std::vector<std::uint64_t> orders_of(
    const std::vector<BsgsGroup>& v) {
    std::vector<std::uint64_t> o;
    o.reserve(v.size());
    for (const BsgsGroup& g : v) o.push_back(g.order());
    return o;
}

void check_candidates_shape(const BsgsGroup& h,
                            const std::vector<BsgsGroup>& cands) {
    for (const BsgsGroup& c : cands) {
        EXPECT_TRUE(c.is_transitive());
        EXPECT_LT(c.order(), h.order());
        EXPECT_EQ(h.order() % c.order(), 0ULL);  // Lagrange
        for (const Perm& g : c.generators()) {
            EXPECT_TRUE(h.contains(g));
        }
    }
    for (std::size_t i = 1U; i < cands.size(); ++i) {
        EXPECT_GE(cands[i - 1U].order(), cands[i].order());  // sorted
    }
}

TEST(GaloisWreathMaximalTest, DenseMaximalClassesOfS4) {
    // Classical pin: the maximal subgroups of S₄ are A₄ (12), D₄ (8),
    // S₃ (6) — one class each.
    symbolic::CASContext ctx;
    auto s4 = must_build(4U, symmetric_gens(4U));
    auto r = maximal_subgroup_classes_in(s4, kOps, kBytes, &ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(orders_of(r.value()),
              (std::vector<std::uint64_t>{12U, 8U, 6U}));
}

TEST(GaloisWreathMaximalTest, DenseAllClassesOfS3) {
    // S₃ has exactly four subgroup classes: S₃, A₃, ⟨(01)⟩, 1.
    symbolic::CASContext ctx;
    auto s3 = must_build(3U, symmetric_gens(3U));
    auto r = all_subgroup_classes_in(s3, kOps, kBytes, &ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(orders_of(r.value()),
              (std::vector<std::uint64_t>{6U, 3U, 2U, 1U}));
}

TEST(GaloisWreathMaximalTest, DetectsStandardEvenAndRelabelled) {
    symbolic::CASContext ctx;
    const BsgsGroup w = full_wreath52();
    ASSERT_EQ(w.order(), 28800ULL);
    auto sw = detect_wreath_preimage(w);
    ASSERT_TRUE(sw.is_ok()) << sw.error().message;
    ASSERT_TRUE(sw.value().has_value());
    EXPECT_EQ(sw.value()->d, 5U);
    EXPECT_EQ(sw.value()->k, 2U);
    EXPECT_EQ(sw.value()->q_image.order(), 8ULL);

    const BsgsGroup even = even_wreath52();
    ASSERT_EQ(even.order(), 14400ULL);
    EXPECT_FALSE(even.has_odd_element());
    auto se = detect_wreath_preimage(even);
    ASSERT_TRUE(se.is_ok()) << se.error().message;
    ASSERT_TRUE(se.value().has_value());
    EXPECT_EQ(se.value()->q_image.order(), 4ULL);

    // Relabelled copy: conjugate by ρ(i) = 3i + 1 (mod 10) — blocks are
    // scattered; detection must recover the shape through `align`.
    Perm rho(10U);
    for (std::size_t i = 0U; i < 10U; ++i) {
        rho[i] = static_cast<std::uint8_t>((3U * i + 1U) % 10U);
    }
    std::vector<Perm> conj_gens;
    for (const Perm& g : w.generators()) {
        conj_gens.push_back(compose(rho, compose(g, inverse(rho))));
    }
    const BsgsGroup wr = must_build(10U, std::move(conj_gens));
    ASSERT_EQ(wr.order(), 28800ULL);
    auto sr = detect_wreath_preimage(wr);
    ASSERT_TRUE(sr.is_ok()) << sr.error().message;
    ASSERT_TRUE(sr.value().has_value());
    auto cands = wreath_preimage_maximal_transitive(wr, *sr.value(), kOps,
                                                    kBytes, &ctx);
    ASSERT_TRUE(cands.is_ok()) << cands.error().message;
    EXPECT_EQ(orders_of(cands.value()),
              (std::vector<std::uint64_t>{14400U, 14400U, 800U}));
    check_candidates_shape(wr, cands.value());

    // A primitive group is no wreath preimage.
    auto sa = detect_wreath_preimage(must_build(5U, alternating_gens(5U)));
    ASSERT_TRUE(sa.is_ok());
    EXPECT_FALSE(sa.value().has_value());
}

TEST(GaloisWreathMaximalTest, FullWreath52Candidates) {
    // Theory (header of galois_wreath_maximal_internal.hpp): exactly the
    // two transitive index-2 kernels (14400 = ker χ₁ and ker χ₁χ₂) and
    // F₂₀ ≀ S₂ (800). A₅ ≀ S₂ (7200) and N_W(Δ) (240) are pruned as
    // direct subgroups of ker χ₁; the base S₅ × S₅ is intransitive.
    symbolic::CASContext ctx;
    const BsgsGroup w = full_wreath52();
    auto sh = detect_wreath_preimage(w);
    ASSERT_TRUE(sh.is_ok());
    ASSERT_TRUE(sh.value().has_value());
    auto cands = wreath_preimage_maximal_transitive(w, *sh.value(), kOps,
                                                    kBytes, &ctx);
    ASSERT_TRUE(cands.is_ok()) << cands.error().message;
    check_candidates_shape(w, cands.value());
    EXPECT_EQ(orders_of(cands.value()),
              (std::vector<std::uint64_t>{14400U, 14400U, 800U}));
    // The two 14400 classes are DISTINCT subgroups (χ₁ vs χ₁χ₂ kernels).
    bool same = true;
    for (const Perm& g : cands.value()[0].generators()) {
        same = same && cands.value()[1].contains(g);
    }
    EXPECT_FALSE(same);
}

TEST(GaloisWreathMaximalTest, KerChi1Candidates) {
    // Theory: two 7200 preimage lifts (A₅² ⋊ decorated swaps), the
    // F₂₀-shell intersections of order 400, and the TWO diagonal classes
    // N_H(Δ) ≅ S₅ × C₂ of order 240 (the H-class of Δ splits — the twin
    // transversal is essential here).
    symbolic::CASContext ctx;
    const BsgsGroup h = ker_chi1();
    ASSERT_EQ(h.order(), 14400ULL);
    EXPECT_TRUE(h.has_odd_element());  // (1,1)·swap is odd on Ω
    auto sh = detect_wreath_preimage(h);
    ASSERT_TRUE(sh.is_ok());
    ASSERT_TRUE(sh.value().has_value());
    auto cands = wreath_preimage_maximal_transitive(h, *sh.value(), kOps,
                                                    kBytes, &ctx);
    ASSERT_TRUE(cands.is_ok()) << cands.error().message;
    check_candidates_shape(h, cands.value());
    const auto ord = orders_of(cands.value());
    EXPECT_EQ(std::count(ord.begin(), ord.end(), 7200U), 2);
    EXPECT_GE(std::count(ord.begin(), ord.end(), 400U), 1);
    EXPECT_EQ(std::count(ord.begin(), ord.end(), 240U), 2);
    EXPECT_EQ(ord.front(), 7200U);
    // The two 240 diagonal classes are distinct subgroups.
    std::vector<const BsgsGroup*> diag;
    for (const BsgsGroup& c : cands.value()) {
        if (c.order() == 240U) diag.push_back(&c);
    }
    ASSERT_EQ(diag.size(), 2U);
    bool same = true;
    for (const Perm& g : diag[0]->generators()) {
        same = same && diag[1]->contains(g);
    }
    EXPECT_FALSE(same);
}

TEST(GaloisWreathMaximalTest, EvenWreath52Candidates) {
    // Theory: W ∩ A₁₀ has NO transitive index-2 subgroup (the unique
    // maximal of its C₄ sign image lifts to a block-preserving,
    // intransitive group), the A₅-shell intersection collapses to the
    // intransitive A₅², and the diagonal normalizer is intransitive —
    // so every candidate is of F₂₀-shell type, order 400.
    symbolic::CASContext ctx;
    const BsgsGroup h = even_wreath52();
    ASSERT_EQ(h.order(), 14400ULL);
    auto sh = detect_wreath_preimage(h);
    ASSERT_TRUE(sh.is_ok());
    ASSERT_TRUE(sh.value().has_value());
    EXPECT_EQ(sh.value()->q_image.order(), 4ULL);
    auto cands = wreath_preimage_maximal_transitive(h, *sh.value(), kOps,
                                                    kBytes, &ctx);
    ASSERT_TRUE(cands.is_ok()) << cands.error().message;
    check_candidates_shape(h, cands.value());
    ASSERT_FALSE(cands.value().empty());
    for (const BsgsGroup& c : cands.value()) {
        EXPECT_EQ(c.order(), 400ULL);
    }
}

TEST(GaloisWreathMaximalTest, HypothesisGuardsAreStructured) {
    // d = 2 (A_d not non-abelian simple): S₂ ≀ S₂ = D₄ on 4 points is a
    // wreath preimage (trivial kernel) but the structural generator must
    // refuse it — in production such nodes always fit the dense route.
    symbolic::CASContext ctx;
    auto w22 = wreath_gens(2U, 2U);
    ASSERT_TRUE(w22.is_ok());
    const BsgsGroup d4 = must_build(4U, std::move(w22.value()));
    auto sh = detect_wreath_preimage(d4);
    ASSERT_TRUE(sh.is_ok());
    ASSERT_TRUE(sh.value().has_value());
    auto r = wreath_preimage_maximal_transitive(d4, *sh.value(), kOps,
                                                kBytes, &ctx);
    ASSERT_TRUE(r.is_error());
    EXPECT_EQ(r.error().kind, CASErrorKind::Unimplemented);
}

TEST(GaloisWreathMaximalTest, NodeDispatchDenseAndStructural) {
    // Small node → identical to the dense route; over-budget wreath
    // node → structural; over-budget non-wreath node → structured
    // Unimplemented.
    symbolic::CASContext ctx;
    const BsgsGroup s4 = must_build(4U, symmetric_gens(4U));
    auto dense = node_maximal_transitive_classes(s4, kOps, kBytes, &ctx);
    ASSERT_TRUE(dense.is_ok()) << dense.error().message;
    auto direct = transitive_subgroup_classes_in(s4, kOps, kBytes, &ctx);
    ASSERT_TRUE(direct.is_ok());
    EXPECT_EQ(orders_of(dense.value()), orders_of(direct.value()));

    const BsgsGroup w = full_wreath52();
    // Byte budget below the dense need of a 28800-node forces the
    // structural route (the production default is exactly in this
    // regime: 2·28800² ≈ 1.66 GiB ≫ 256 MiB).
    ASSERT_GT(dense_sublattice_min_bytes(w.order(), w.degree()), kBytes);
    auto structural = node_maximal_transitive_classes(w, kOps, kBytes, &ctx);
    ASSERT_TRUE(structural.is_ok()) << structural.error().message;
    EXPECT_EQ(orders_of(structural.value()),
              (std::vector<std::uint64_t>{14400U, 14400U, 800U}));

    // A₁₀ (order 1814400) is primitive: no wreath shape, must fail
    // structured (never silently wrong).
    const BsgsGroup a10 = must_build(10U, alternating_gens(10U));
    auto bad = node_maximal_transitive_classes(a10, kOps, kBytes, &ctx);
    ASSERT_TRUE(bad.is_error());
    EXPECT_EQ(bad.error().kind, CASErrorKind::Unimplemented);
}

TEST(GaloisWreathMaximalTest, ScottLemmaGroundTruthOnA5xA5) {
    // Machine check of the Scott step used by family FD: every MAXIMAL
    // class of A₅ × A₅ whose two block projections are BOTH all of A₅
    // must be a full diagonal (order 60) — dense exhaustive ground
    // truth, no transcription.
    symbolic::CASContext ctx;
    std::vector<Perm> gens;
    for (const Perm& g : alternating_gens(5U)) {
        gens.push_back(embed5(g, 0U));
        gens.push_back(embed5(g, 1U));
    }
    const BsgsGroup a5a5 = must_build(10U, std::move(gens));
    ASSERT_EQ(a5a5.order(), 3600ULL);
    auto r = maximal_subgroup_classes_in(a5a5, kOps, 1ULL << 30U, &ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    ASSERT_FALSE(r.value().empty());
    bool saw_diagonal = false;
    for (const BsgsGroup& m : r.value()) {
        std::uint64_t proj_order[2] = {0U, 0U};
        for (std::size_t b = 0U; b < 2U; ++b) {
            std::vector<Perm> proj;
            for (const Perm& g : m.generators()) {
                Perm p(5U);
                for (std::size_t x = 0U; x < 5U; ++x) {
                    p[x] = static_cast<std::uint8_t>(g[b * 5U + x] -
                                                     b * 5U);
                }
                proj.push_back(std::move(p));
            }
            proj_order[b] = must_build(5U, std::move(proj)).order();
        }
        if (proj_order[0] == 60ULL && proj_order[1] == 60ULL) {
            EXPECT_EQ(m.order(), 60ULL);  // Scott: full diagonal
            saw_diagonal = true;
        }
    }
    EXPECT_TRUE(saw_diagonal);
}

TEST(GaloisWreathMaximalTest, NormalizerShellGroundTruthOnA5) {
    // Machine check of the FB step: for every proper nontrivial subgroup
    // class D of A₅, N_{S₅}(D) is proper, and when it is transitive it
    // lies (up to S₅-conjugacy) inside a Brick-2 maximal transitive
    // shell of S₅ — so the K ≀ S_k shells cover case (i) of the theorem.
    symbolic::CASContext ctx;
    const BsgsGroup a5 = must_build(5U, alternating_gens(5U));
    auto classes = all_subgroup_classes_in(a5, kOps, kBytes, &ctx);
    ASSERT_TRUE(classes.is_ok()) << classes.error().message;
    // Enumerate S₅'s 120 elements once.
    std::vector<Perm> s5_elems{identity(5U)};
    for (std::size_t i = 0U; i < s5_elems.size(); ++i) {
        for (const Perm& g : symmetric_gens(5U)) {
            Perm nxt = compose(g, s5_elems[i]);
            if (std::find(s5_elems.begin(), s5_elems.end(), nxt) ==
                s5_elems.end()) {
                s5_elems.push_back(std::move(nxt));
            }
        }
    }
    ASSERT_EQ(s5_elems.size(), 120U);
    auto shells = maximal_transitive_candidates(AmbientGroup::Symmetric, 5U);
    ASSERT_TRUE(shells.is_ok());
    for (const BsgsGroup& d : classes.value()) {
        if (d.order() == 1ULL || d.order() == 60ULL) continue;
        std::vector<Perm> normalizers;
        for (const Perm& s : s5_elems) {
            bool norm = true;
            for (const Perm& g : d.generators()) {
                if (!d.contains(compose(s, compose(g, inverse(s))))) {
                    norm = false;
                    break;
                }
            }
            if (norm) normalizers.push_back(s);
        }
        const BsgsGroup n = must_build(5U, std::move(normalizers));
        EXPECT_LT(n.order(), 120ULL) << "D order " << d.order();
        if (!n.is_transitive()) continue;
        bool covered = false;
        for (const auto& shell : shells.value()) {
            for (const Perm& s : s5_elems) {
                bool inside = true;
                for (const Perm& g : n.generators()) {
                    if (!shell.group.contains(
                            compose(s, compose(g, inverse(s))))) {
                        inside = false;
                        break;
                    }
                }
                if (inside) {
                    covered = true;
                    break;
                }
            }
            if (covered) break;
        }
        EXPECT_TRUE(covered) << "N(D) order " << n.order()
                             << " uncovered, D order " << d.order();
    }
}

}  // namespace
