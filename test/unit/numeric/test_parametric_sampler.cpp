// F7.1-T1 — Parametric plot sampler unit tests.

#include "cas/lexer.hpp"
#include "cas/numeric.hpp"
#include "cas/numeric/sampler.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

namespace cas::numeric {
namespace {

ExprPtr parse_expr(const std::string& s, AstArena& arena) {
    auto t = Lexer(s).tokenize();
    EXPECT_TRUE(t.is_ok());
    Parser p(t.value(), arena);
    auto r = p.parse();
    EXPECT_TRUE(r.is_ok());
    return r.value();
}

class ParametricSamplerTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

TEST_F(ParametricSamplerTest, UnitCircleClosesNumerically) {
    // (cos t, sin t) on [0, 2π].
    auto x = parse_expr("cos(t)", ctx.arena());
    auto y = parse_expr("sin(t)", ctx.arena());
    ParametricSampler sampler(x, y, "t");
    auto res = sampler.sample(0.0, 2.0 * std::numbers::pi);
    ASSERT_TRUE(res.is_ok());
    ASSERT_GT(res.value().size(), 8U);
    // Every sample sits on the unit circle x² + y² = 1.
    for (const auto& p : res.value()) {
        EXPECT_NEAR(p.x * p.x + p.y * p.y, 1.0, 1e-10);
    }
}

TEST_F(ParametricSamplerTest, StraightLineNeedsNoRefinement) {
    // (t, 2t + 1).
    auto x = parse_expr("t", ctx.arena());
    auto y = parse_expr("2*t + 1", ctx.arena());
    ParametricSampler::Options opts;
    opts.angle_tolerance = 0.01;
    opts.max_depth = 8;
    opts.initial_samples = 5;
    ParametricSampler sampler(x, y, "t", opts);
    auto res = sampler.sample(0.0, 1.0);
    ASSERT_TRUE(res.is_ok());
    // Refinement should not trigger on a straight line.
    EXPECT_EQ(res.value().size(), opts.initial_samples + 1U);
}

TEST_F(ParametricSamplerTest, EmptyRangeReturnsEmpty) {
    auto x = parse_expr("t", ctx.arena());
    auto y = parse_expr("t", ctx.arena());
    ParametricSampler sampler(x, y, "t");
    auto res = sampler.sample(1.0, 1.0);
    ASSERT_TRUE(res.is_ok());
    EXPECT_TRUE(res.value().empty());
}

}  // namespace
}  // namespace cas::numeric
