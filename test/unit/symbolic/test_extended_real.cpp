// F7.5.F1 — Extended-Real AST tests.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Extended_Real_AST.md

#include "cas/ast.hpp"
#include "cas/extended_real.hpp"

#include <gtest/gtest.h>

namespace cas {
namespace {

TEST(ExtendedRealAST, EnumValuesAreDistinct) {
    EXPECT_NE(static_cast<int>(MathConstant::Infinity),
              static_cast<int>(MathConstant::NegInfinity));
    EXPECT_NE(static_cast<int>(MathConstant::Infinity),
              static_cast<int>(MathConstant::ComplexInfinity));
    EXPECT_NE(static_cast<int>(MathConstant::ComplexInfinity),
              static_cast<int>(MathConstant::Indeterminate));
}

TEST(ExtendedRealAST, IsPosInfinityDetectsCanonical) {
    AstArena arena;
    EXPECT_TRUE(is_pos_infinity(make_pos_infinity(arena)));
    EXPECT_FALSE(is_pos_infinity(make_neg_infinity(arena)));
    EXPECT_FALSE(is_pos_infinity(make_complex_infinity(arena)));
    EXPECT_FALSE(is_pos_infinity(make_indeterminate(arena)));
}

TEST(ExtendedRealAST, IsNegInfinityAcceptsCanonical) {
    AstArena arena;
    EXPECT_TRUE(is_neg_infinity(make_neg_infinity(arena)));
    EXPECT_FALSE(is_neg_infinity(make_pos_infinity(arena)));
    EXPECT_FALSE(is_neg_infinity(make_complex_infinity(arena)));
}

TEST(ExtendedRealAST, IsNegInfinityAcceptsLegacyUnaryNeg) {
    AstArena arena;
    ExprPtr legacy = arena.make<Unary>(
        UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity));
    EXPECT_TRUE(is_neg_infinity(legacy))
        << "Transition: predicate must accept legacy Unary(Neg, Infinity)";
}

TEST(ExtendedRealAST, IsComplexInfinityIsolated) {
    AstArena arena;
    EXPECT_TRUE(is_complex_infinity(make_complex_infinity(arena)));
    EXPECT_FALSE(is_complex_infinity(make_pos_infinity(arena)));
    EXPECT_FALSE(is_complex_infinity(make_neg_infinity(arena)));
    EXPECT_FALSE(is_complex_infinity(make_indeterminate(arena)));
}

TEST(ExtendedRealAST, IsIndeterminateIsolated) {
    AstArena arena;
    EXPECT_TRUE(is_indeterminate(make_indeterminate(arena)));
    EXPECT_FALSE(is_indeterminate(make_pos_infinity(arena)));
    EXPECT_FALSE(is_indeterminate(make_neg_infinity(arena)));
    EXPECT_FALSE(is_indeterminate(make_complex_infinity(arena)));
}

TEST(ExtendedRealAST, IsSignedInfinityCovers) {
    AstArena arena;
    EXPECT_TRUE(is_signed_infinity(make_pos_infinity(arena)));
    EXPECT_TRUE(is_signed_infinity(make_neg_infinity(arena)));
    EXPECT_FALSE(is_signed_infinity(make_complex_infinity(arena)));
    EXPECT_FALSE(is_signed_infinity(make_indeterminate(arena)));
}

TEST(ExtendedRealAST, IsAnyInfinityCoversAllInfinities) {
    AstArena arena;
    EXPECT_TRUE(is_any_infinity(make_pos_infinity(arena)));
    EXPECT_TRUE(is_any_infinity(make_neg_infinity(arena)));
    EXPECT_TRUE(is_any_infinity(make_complex_infinity(arena)));
    EXPECT_FALSE(is_any_infinity(make_indeterminate(arena)));
}

TEST(ExtendedRealAST, NullPtrPredicatesReturnFalse) {
    EXPECT_FALSE(is_pos_infinity(nullptr));
    EXPECT_FALSE(is_neg_infinity(nullptr));
    EXPECT_FALSE(is_complex_infinity(nullptr));
    EXPECT_FALSE(is_indeterminate(nullptr));
    EXPECT_FALSE(is_signed_infinity(nullptr));
    EXPECT_FALSE(is_any_infinity(nullptr));
}

TEST(ExtendedRealAST, NonConstantNodesReturnFalse) {
    AstArena arena;
    ExprPtr sym = arena.make<Symbol>("x");
    EXPECT_FALSE(is_pos_infinity(sym));
    EXPECT_FALSE(is_neg_infinity(sym));
    EXPECT_FALSE(is_complex_infinity(sym));
    EXPECT_FALSE(is_indeterminate(sym));
}

TEST(ExtendedRealAST, FactoryRoundTripPreservesValue) {
    AstArena arena;
    ExprPtr c1 = make_complex_infinity(arena);
    ExprPtr c2 = make_complex_infinity(arena);
    const auto* k1 = expr_cast<Constant>(c1);
    const auto* k2 = expr_cast<Constant>(c2);
    ASSERT_NE(k1, nullptr);
    ASSERT_NE(k2, nullptr);
    EXPECT_EQ(k1->value, MathConstant::ComplexInfinity);
    EXPECT_EQ(k2->value, MathConstant::ComplexInfinity);
}

TEST(ExtendedRealAST, LegacyUnaryNegInfinityNotConfusedWithPosInfinity) {
    AstArena arena;
    ExprPtr legacy = arena.make<Unary>(
        UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity));
    EXPECT_FALSE(is_pos_infinity(legacy));
}

TEST(ExtendedRealAST, UnaryNegOfNonInfinityIsNotNegInfinity) {
    AstArena arena;
    ExprPtr neg_pi = arena.make<Unary>(
        UnaryOp::Neg, arena.make<Constant>(MathConstant::Pi));
    EXPECT_FALSE(is_neg_infinity(neg_pi));
}

}  // namespace
}  // namespace cas
