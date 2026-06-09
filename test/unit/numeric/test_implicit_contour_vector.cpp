// F7.0-B1 — Implicit / Contour / VectorField sampler tests.

#include <gtest/gtest.h>

#include "cas/numeric/sampler.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

#include <cmath>

using namespace cas;
using namespace cas::numeric;
using namespace cas::symbolic;

namespace {

ExprPtr add_(AstArena& ar, ExprPtr lhs, ExprPtr rhs) {
    return ar.make<Sum>(std::vector<ExprPtr>{lhs, rhs});
}
ExprPtr pow_(AstArena& ar, ExprPtr base, long long e) {
    return ar.make<Binary>(BinaryOp::Pow, base, ar.make<IntegerLit>(BigInt(e)));
}
ExprPtr neg_(AstArena& ar, ExprPtr e) {
    return ar.make<Unary>(UnaryOp::Neg, e);
}
ExprPtr int_(AstArena& ar, long long v) {
    return ar.make<IntegerLit>(BigInt(v));
}
ExprPtr sym_(AstArena& ar, const std::string& s) {
    return ar.make<Symbol>(s);
}

TEST(ImplicitSampler, UnitCircleProducesSegments) {
    CASContext ctx;
    auto& ar = ctx.arena();
    auto expr = add_(ar, add_(ar, pow_(ar, sym_(ar, "x"), 2),
                                  pow_(ar, sym_(ar, "y"), 2)),
                         int_(ar, -1));
    ImplicitSampler::Options o;
    o.nx = 64;
    o.ny = 64;
    o.level = 0.0;
    ImplicitSampler s(expr, "x", "y", o);
    auto segs = s.sample(-1.5, 1.5, -1.5, 1.5);
    ASSERT_TRUE(segs.is_ok());
    EXPECT_GT(segs.value().size(), 50U);
    for (const auto& sg : segs.value()) {
        const double mx = 0.5 * (sg.a.x + sg.b.x);
        const double my = 0.5 * (sg.a.y + sg.b.y);
        const double r = std::sqrt(mx * mx + my * my);
        EXPECT_NEAR(r, 1.0, 0.1);
    }
}

TEST(ImplicitSampler, EmptyDomainNoSegments) {
    CASContext ctx;
    auto& ar = ctx.arena();
    auto expr = add_(ar, add_(ar, pow_(ar, sym_(ar, "x"), 2),
                                  pow_(ar, sym_(ar, "y"), 2)),
                         int_(ar, -100));
    ImplicitSampler s(expr, "x", "y");
    auto segs = s.sample(-1.0, 1.0, -1.0, 1.0);
    ASSERT_TRUE(segs.is_ok());
    EXPECT_EQ(segs.value().size(), 0U);
}

TEST(ContourSampler, ThreeLevelsThreeSegmentSets) {
    CASContext ctx;
    auto& ar = ctx.arena();
    auto expr = add_(ar, pow_(ar, sym_(ar, "x"), 2),
                         pow_(ar, sym_(ar, "y"), 2));
    ContourSampler cs(expr, "x", "y", {0.25, 1.0, 2.25});
    auto res = cs.sample(-2.0, 2.0, -2.0, 2.0);
    ASSERT_TRUE(res.is_ok());
    ASSERT_EQ(res.value().size(), 3U);
    for (const auto& lvl_segs : res.value()) {
        EXPECT_GT(lvl_segs.size(), 10U);
    }
}

TEST(VectorFieldSampler, ProducesGridArrows) {
    CASContext ctx;
    auto& ar = ctx.arena();
    auto fx = neg_(ar, sym_(ar, "y"));
    auto fy = sym_(ar, "x");
    VectorFieldSampler::Options o;
    o.nx = 10;
    o.ny = 10;
    VectorFieldSampler vfs(fx, fy, "x", "y", o);
    auto arrows = vfs.sample(-1.0, 1.0, -1.0, 1.0);
    ASSERT_TRUE(arrows.is_ok());
    EXPECT_EQ(arrows.value().size(), 100U);
    const auto& arr = arrows.value().front();
    EXPECT_NEAR(arr.dir.x, -arr.base.y, 1e-12);
    EXPECT_NEAR(arr.dir.y,  arr.base.x, 1e-12);
}

TEST(ImplicitSampler, DegenerateDomainEmpty) {
    CASContext ctx;
    auto& ar = ctx.arena();
    auto expr = add_(ar, sym_(ar, "x"), sym_(ar, "y"));
    ImplicitSampler s(expr, "x", "y");
    auto res = s.sample(1.0, 1.0, 0.0, 1.0);  // empty x range
    ASSERT_TRUE(res.is_ok());
    EXPECT_EQ(res.value().size(), 0U);
}

}  // namespace
