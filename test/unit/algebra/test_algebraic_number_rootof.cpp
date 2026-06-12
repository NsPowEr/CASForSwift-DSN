// F8.0-5.4: tests for RootOf with rigorous rational isolating bounds.
// Verifies that make_rootof_isolated:
//   1. Produces a RootOf node carrying an IsolatingBound.
//   2. The bound brackets the requested root and contains no other real
//      roots of the polynomial.
//   3. structural_equal distinguishes RootOf nodes by their bound.
//   4. Bound is preserved through clone_into_arena.

#include "cas/algebra.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "cas/lexer.hpp"
#include "cas/numeric.hpp"
#include "cas/parser.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>
#include <unordered_map>

namespace cas {
namespace {

// Parse a polynomial via the CAS Parser+Simplifier so it lands in the same
// canonical form parse_polynomial expects.
[[nodiscard]] ExprPtr parse_poly(symbolic::CASContext& ctx, const std::string& s) {
    auto tokens = Lexer(s).tokenize();
    if (tokens.is_error()) return nullptr;
    Parser p(tokens.value(), ctx.arena());
    auto e = p.parse();
    if (e.is_error()) return nullptr;
    auto simp = ctx.simplify(e.value());
    return simp.is_ok() ? simp.value() : e.value();
}

TEST(AlgebraicNumberRootOf, Sqrt2_IsolatedBound_Brackets_Sqrt2) {
    symbolic::CASContext ctx;
    Symbol x("x");
    ExprPtr poly = parse_poly(ctx, "x^2 - 2");
    ASSERT_NE(poly, nullptr);

    // In [-5, 5] there are two real roots; root 0 = -√2, root 1 = +√2.
    auto root_res = algebra::make_rootof_isolated(poly, x, 1U, ctx, -5.0, 5.0);
    ASSERT_TRUE(root_res.is_ok()) << root_res.error().message;

    const auto* rootof = expr_cast<RootOf>(root_res.value());
    ASSERT_NE(rootof, nullptr) << "expected a RootOf node";
    ASSERT_TRUE(rootof->isolating_bound.has_value())
        << "rigorous constructor must populate isolating_bound";

    const auto& b = *rootof->isolating_bound;
    Rational low(b.low_num, b.low_den);
    Rational high(b.high_num, b.high_den);

    // √2 = 1.4142135623730951... — the bound must straddle it.
    constexpr double sqrt2 = 1.41421356237309515;
    EXPECT_LT(low.to_double(), sqrt2)
        << "low=" << low.to_double() << " must be < √2";
    EXPECT_GT(high.to_double(), sqrt2)
        << "high=" << high.to_double() << " must be > √2";

    // The bound must be sharper than the search interval [-5, 5] and
    // must isolate the positive root (away from the negative root −√2).
    EXPECT_GT(low.to_double(), 0.0);
    EXPECT_LT(high.to_double(), 5.0);
}

TEST(AlgebraicNumberRootOf, CubicThreeRoots_AllIsolated) {
    symbolic::CASContext ctx;
    Symbol x("x");
    ExprPtr poly = parse_poly(ctx, "x^3 - x");
    ASSERT_NE(poly, nullptr);

    // Roots: -1, 0, +1 in ascending order.
    auto r0 = algebra::make_rootof_isolated(poly, x, 0U, ctx, -2.0, 2.0);
    auto r1 = algebra::make_rootof_isolated(poly, x, 1U, ctx, -2.0, 2.0);
    auto r2 = algebra::make_rootof_isolated(poly, x, 2U, ctx, -2.0, 2.0);
    ASSERT_TRUE(r0.is_ok()) << r0.error().message;
    ASSERT_TRUE(r1.is_ok()) << r1.error().message;
    ASSERT_TRUE(r2.is_ok()) << r2.error().message;

    auto bound_of = [](ExprPtr e) {
        const auto* rf = expr_cast<RootOf>(e);
        return Rational(rf->isolating_bound->low_num, rf->isolating_bound->low_den).to_double();
    };

    // Sturm sorts by midpoint; verify ascending order.
    EXPECT_LT(bound_of(r0.value()), bound_of(r1.value()));
    EXPECT_LT(bound_of(r1.value()), bound_of(r2.value()));

    // Out-of-range index must fail diagnostically.
    auto r3 = algebra::make_rootof_isolated(poly, x, 3U, ctx, -2.0, 2.0);
    ASSERT_TRUE(r3.is_error());
    EXPECT_NE(r3.error().message.find("out of range"), std::string::npos);
}

TEST(AlgebraicNumberRootOf, BoundDistinguishes_DifferentRoots) {
    symbolic::CASContext ctx;
    Symbol x("x");
    ExprPtr poly = parse_poly(ctx, "x^3 - x");
    ASSERT_NE(poly, nullptr);

    auto r0_res = algebra::make_rootof_isolated(poly, x, 0U, ctx, -2.0, 2.0);
    auto r2_res = algebra::make_rootof_isolated(poly, x, 2U, ctx, -2.0, 2.0);
    ASSERT_TRUE(r0_res.is_ok() && r2_res.is_ok());

    // structural_equal must distinguish RootOf nodes by their bounds even
    // though they share the same polynomial+variable.
    EXPECT_FALSE(structural_equal(r0_res.value(), r2_res.value()));
}

TEST(AlgebraicNumberRootOf, BoundSurvives_CloneIntoArena) {
    symbolic::CASContext ctx;
    Symbol x("x");
    ExprPtr poly = parse_poly(ctx, "x^2 - 2");
    ASSERT_NE(poly, nullptr);
    auto rootof_res = algebra::make_rootof_isolated(poly, x, 1U, ctx, -5.0, 5.0);
    ASSERT_TRUE(rootof_res.is_ok()) << rootof_res.error().message;

    AstArena target;
    std::unordered_map<ExprPtr, ExprPtr> cache;
    auto cloned = clone_into_arena(rootof_res.value(), target, cache);
    const auto* rf = expr_cast<RootOf>(cloned);
    ASSERT_NE(rf, nullptr);
    ASSERT_TRUE(rf->isolating_bound.has_value())
        << "clone_into_arena must preserve isolating_bound";

    const auto& orig_bound = *expr_cast<RootOf>(rootof_res.value())->isolating_bound;
    const auto& clone_bound = *rf->isolating_bound;
    EXPECT_EQ(orig_bound.low_num, clone_bound.low_num);
    EXPECT_EQ(orig_bound.low_den, clone_bound.low_den);
    EXPECT_EQ(orig_bound.high_num, clone_bound.high_num);
    EXPECT_EQ(orig_bound.high_den, clone_bound.high_den);
}

TEST(AlgebraicNumberRootOf, NumericEval_UsesBound_NoSiblingDrift) {
    // Numeric evaluator: a RootOf carrying an isolating bound for +√2 must
    // return +√2 (not -√2), even when both roots lie inside the legacy
    // ±1e3 sweep window. Verifies the bound short-circuits the search.
    symbolic::CASContext ctx;
    Symbol x("x");
    ExprPtr poly = parse_poly(ctx, "x^2 - 2");
    ASSERT_NE(poly, nullptr);

    auto rootof_res = algebra::make_rootof_isolated(poly, x, 1U, ctx, -5.0, 5.0);
    ASSERT_TRUE(rootof_res.is_ok()) << rootof_res.error().message;

    cas::numeric::NumericEvaluator eval;
    auto val = eval.evaluate(rootof_res.value());
    ASSERT_TRUE(val.is_ok()) << val.error().message;
    EXPECT_NEAR(val.value(), 1.41421356237309515, 1e-9);
    EXPECT_GT(val.value(), 0.0) << "bound should pin to +√2, not −√2";
}

TEST(AlgebraicNumberRootOf, RoundTrip_PreservesBound) {
    // Print → parse → bound must survive intact.
    symbolic::CASContext ctx;
    Symbol x("x");
    ExprPtr poly = parse_poly(ctx, "x^2 - 2");
    ASSERT_NE(poly, nullptr);

    auto root_res = algebra::make_rootof_isolated(poly, x, 1U, ctx, -5.0, 5.0);
    ASSERT_TRUE(root_res.is_ok());

    const auto* orig = expr_cast<RootOf>(root_res.value());
    ASSERT_TRUE(orig->isolating_bound.has_value());
    const auto& orig_bound = *orig->isolating_bound;

    // Print using the round-trip printer.
    auto printed = to_round_trip_text(root_res.value());
    ASSERT_TRUE(printed.is_ok()) << printed.error().message;

    // Re-parse.
    auto tokens = Lexer(printed.value()).tokenize();
    ASSERT_TRUE(tokens.is_ok());
    Parser p(tokens.value(), ctx.arena());
    auto reparsed = p.parse();
    ASSERT_TRUE(reparsed.is_ok()) << reparsed.error().message;

    const auto* rf2 = expr_cast<RootOf>(reparsed.value());
    ASSERT_NE(rf2, nullptr);
    ASSERT_TRUE(rf2->isolating_bound.has_value())
        << "round-trip lost the isolating bound; printed: " << printed.value();

    const auto& parsed_bound = *rf2->isolating_bound;
    EXPECT_EQ(orig_bound.low_num, parsed_bound.low_num);
    EXPECT_EQ(orig_bound.low_den, parsed_bound.low_den);
    EXPECT_EQ(orig_bound.high_num, parsed_bound.high_num);
    EXPECT_EQ(orig_bound.high_den, parsed_bound.high_den);
}

TEST(AlgebraicNumberRootOf, QAlphaArithmetic_Interop_With_IsolatedRootOf) {
    // Bound-enriched RootOf integrates with Q(α) arithmetic: building the
    // AlgebraicNumber α from a RootOf carrying an IsolatingBound must
    // recover the same minimal polynomial regardless of the bound.
    symbolic::CASContext ctx;
    Symbol x("x");
    ExprPtr poly = parse_poly(ctx, "x^2 - 2");
    ASSERT_NE(poly, nullptr);

    auto root_res = algebra::make_rootof_isolated(poly, x, 1U, ctx, -5.0, 5.0);
    ASSERT_TRUE(root_res.is_ok());
    const auto* root_node = expr_cast<RootOf>(root_res.value());
    ASSERT_NE(root_node, nullptr);
    ASSERT_TRUE(root_node->isolating_bound.has_value());

    // alpha_from_rootof must succeed regardless of the IsolatingBound payload.
    auto alpha_res = algebra::alpha_from_rootof(*root_node, ctx);
    ASSERT_TRUE(alpha_res.is_ok()) << alpha_res.error().message;

    const auto& alpha = alpha_res.value();
    // α² − 2 = 0  ⟹  the minimal polynomial in monic ascending form is
    // [-2, 0, 1]  (i.e. -2 + 0·x + 1·x²).
    const auto& mp = alpha.min_poly();
    ASSERT_EQ(mp.size(), 3U);
    EXPECT_EQ(mp[0], Rational(BigInt(-2)));
    EXPECT_EQ(mp[1], Rational(BigInt(0)));
    EXPECT_EQ(mp[2], Rational(BigInt(1)));

    // α + α  must reduce to 2α  (value vector = [0, 2]).
    auto sum = alpha + alpha;
    ASSERT_EQ(sum.value().size(), 2U);
    EXPECT_EQ(sum.value()[0], Rational(BigInt(0)));
    EXPECT_EQ(sum.value()[1], Rational(BigInt(2)));

    // α · α  must reduce via the minimal polynomial to 2  (constant).
    auto sq = alpha * alpha;
    // Either size 1 (just [2]) or size 2/3 with non-trivial constant only.
    ASSERT_GE(sq.value().size(), 1U);
    EXPECT_EQ(sq.value()[0], Rational(BigInt(2)));
    for (std::size_t k = 1; k < sq.value().size(); ++k)
        EXPECT_EQ(sq.value()[k], Rational(BigInt(0)))
            << "α² must reduce to a pure rational, but coeff[" << k
            << "] = " << sq.value()[k].numerator().decimal();
}

TEST(AlgebraicNumberRootOf, LegacyRootOf_NoBound_BackwardCompatible) {
    // The legacy constructor (no bound) must remain functional and
    // produce a RootOf with isolating_bound == nullopt.
    symbolic::CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr poly = parse_poly(ctx, "x^2 - 2");
    ASSERT_NE(poly, nullptr);

    ExprPtr legacy = arena.make<RootOf>(poly, Symbol("x"),
        std::optional<std::size_t>{1U});
    const auto* rf = expr_cast<RootOf>(legacy);
    ASSERT_NE(rf, nullptr);
    EXPECT_FALSE(rf->isolating_bound.has_value());
    EXPECT_TRUE(rf->root_index.has_value());
    EXPECT_EQ(*rf->root_index, 1U);
}

} // namespace
} // namespace cas
