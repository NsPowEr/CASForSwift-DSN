#include <gtest/gtest.h>

#include "CasGuiSession.hpp"

#include <algorithm>

namespace {

using cas::gui::CasGuiSession;

TEST(CasGuiSessionTest, SimplifyReturnsStructuredRepresentations) {
    CasGuiSession session;
    const auto result = session.simplify("1/2");

    ASSERT_TRUE(result.ok);
    EXPECT_EQ(result.latex, "\\frac{1}{2}");
    ASSERT_GE(result.representations.size(), 2U);

    const auto numeric = std::find_if(
        result.representations.begin(),
        result.representations.end(),
        [](const auto& item) { return item.id == "numeric"; });
    ASSERT_NE(numeric, result.representations.end());
    EXPECT_EQ(numeric->value, "0.5");
}

TEST(CasGuiSessionTest, SimplifyCapturesRewriteTraceForInspector) {
    CasGuiSession session;
    const auto result = session.simplify("x+0");

    ASSERT_TRUE(result.ok);
    EXPECT_EQ(result.latex, "x");
    ASSERT_FALSE(result.steps.empty());
    EXPECT_EQ(result.steps.back().after_latex, "x");
}

TEST(CasGuiSessionTest, SessionDefinitionsAffectLaterCells) {
    CasGuiSession session;

    const auto define = session.simplify("a = 2");
    ASSERT_TRUE(define.ok);
    EXPECT_EQ(define.latex, "a = 2");

    const auto value = session.simplify("a + 3");
    ASSERT_TRUE(value.ok);
    EXPECT_EQ(value.latex, "5");

    const auto variables = session.list_variables();
    ASSERT_EQ(variables.size(), 1U);
    EXPECT_EQ(variables.front().first, "a");
    EXPECT_EQ(variables.front().second, "2");
}

TEST(CasGuiSessionTest, DefinitionsRoundTripAcrossSnapshotRestore) {
    CasGuiSession source;
    ASSERT_TRUE(source.simplify("a = 2").ok);
    ASSERT_TRUE(source.simplify("b = a + 5").ok);

    const auto snapshot = source.snapshot_definitions();
    ASSERT_EQ(snapshot.size(), 2U);

    CasGuiSession restored;
    const auto restore = restored.restore_definitions(snapshot);
    ASSERT_TRUE(restore.is_ok());

    const auto result = restored.simplify("b + 1");
    ASSERT_TRUE(result.ok);
    EXPECT_EQ(result.latex, "8");
}

TEST(CasGuiSessionTest, PlotVariableShadowsStoredDefinition) {
    CasGuiSession session;
    ASSERT_TRUE(session.simplify("x = 2").ok);

    const auto sampled = session.sample_2d("x", "x", -1.0, 1.0);
    ASSERT_TRUE(sampled.is_ok());
    ASSERT_FALSE(sampled.value().empty());
    EXPECT_LT(sampled.value().front().x, sampled.value().back().x);
    EXPECT_NE(sampled.value().front().y, sampled.value().back().y);
}

} // namespace
