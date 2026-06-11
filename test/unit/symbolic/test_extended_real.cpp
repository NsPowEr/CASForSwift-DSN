// F7.5.F1 — Extended-Real AST tests.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Extended_Real_AST.md

#include "cas/ast.hpp"
#include "cas/extended_real.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

namespace cas {
namespace {

// ── Phase 2 arithmetic test fixture ─────────────────────────────────────
// Drives the simplifier end-to-end and asserts that operations involving
// the new extended-real constants (NegInfinity, ComplexInfinity,
// Indeterminate) propagate per Extended_Real_AST.md §"Aritmetica
// extended-real".
class ExtendedRealArithmeticTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] ExprPtr pos_inf() { return ctx.arena().make<Constant>(MathConstant::Infinity); }
    [[nodiscard]] ExprPtr neg_inf() { return ctx.arena().make<Constant>(MathConstant::NegInfinity); }
    [[nodiscard]] ExprPtr cplx_inf() { return ctx.arena().make<Constant>(MathConstant::ComplexInfinity); }
    [[nodiscard]] ExprPtr indet() { return ctx.arena().make<Constant>(MathConstant::Indeterminate); }
    [[nodiscard]] ExprPtr lit(long long v) {
        return ctx.arena().make<IntegerLit>(BigInt(v));
    }

    [[nodiscard]] ExprPtr simp(ExprPtr e) {
        auto r = ctx.simplify(e);
        if (!r.is_ok()) return nullptr;
        return r.value();
    }

    [[nodiscard]] ExprPtr sum(std::vector<ExprPtr> terms) {
        return ctx.arena().make<Sum>(std::move(terms));
    }
    [[nodiscard]] ExprPtr prod(std::vector<ExprPtr> factors) {
        return ctx.arena().make<Product>(std::move(factors));
    }
    [[nodiscard]] ExprPtr pow(ExprPtr b, ExprPtr e) {
        return ctx.arena().make<Binary>(BinaryOp::Pow, b, e);
    }
};

TEST_F(ExtendedRealArithmeticTest, NegInfPlusNegInfYieldsNegInf) {
    EXPECT_TRUE(is_neg_infinity(simp(sum({neg_inf(), neg_inf()}))));
}

TEST_F(ExtendedRealArithmeticTest, PosInfPlusNegInfIsIndeterminate) {
    EXPECT_TRUE(is_indeterminate(simp(sum({pos_inf(), neg_inf()}))));
}

TEST_F(ExtendedRealArithmeticTest, ComplexInfPlusPosInfIsIndeterminate) {
    EXPECT_TRUE(is_indeterminate(simp(sum({cplx_inf(), pos_inf()}))));
}

TEST_F(ExtendedRealArithmeticTest, ComplexInfPlusFiniteIsComplexInf) {
    EXPECT_TRUE(is_complex_infinity(simp(sum({cplx_inf(), lit(7)}))));
}

TEST_F(ExtendedRealArithmeticTest, NegInfPlusFiniteIsNegInf) {
    EXPECT_TRUE(is_neg_infinity(simp(sum({neg_inf(), lit(1000)}))));
}

TEST_F(ExtendedRealArithmeticTest, IndeterminateAbsorbsUnderAddition) {
    EXPECT_TRUE(is_indeterminate(simp(sum({indet(), lit(5), pos_inf()}))));
}

TEST_F(ExtendedRealArithmeticTest, ZeroTimesNegInfIsIndeterminate) {
    EXPECT_TRUE(is_indeterminate(simp(prod({lit(0), neg_inf()}))));
}

TEST_F(ExtendedRealArithmeticTest, ZeroTimesComplexInfIsIndeterminate) {
    EXPECT_TRUE(is_indeterminate(simp(prod({lit(0), cplx_inf()}))));
}

TEST_F(ExtendedRealArithmeticTest, NegInfTimesNegInfIsPosInf) {
    EXPECT_TRUE(is_pos_infinity(simp(prod({neg_inf(), neg_inf()}))));
}

TEST_F(ExtendedRealArithmeticTest, PosInfTimesNegInfIsNegInf) {
    EXPECT_TRUE(is_neg_infinity(simp(prod({pos_inf(), neg_inf()}))));
}

TEST_F(ExtendedRealArithmeticTest, NegInfTimesNegativeLiteralIsPosInf) {
    EXPECT_TRUE(is_pos_infinity(simp(prod({neg_inf(), lit(-3)}))));
}

TEST_F(ExtendedRealArithmeticTest, ComplexInfTimesNonzeroFiniteIsComplexInf) {
    EXPECT_TRUE(is_complex_infinity(simp(prod({cplx_inf(), lit(5)}))));
}

TEST_F(ExtendedRealArithmeticTest, IndeterminateAbsorbsUnderProduct) {
    EXPECT_TRUE(is_indeterminate(simp(prod({indet(), lit(2)}))));
}

TEST_F(ExtendedRealArithmeticTest, PowNegInfToPositiveOddIsNegInf) {
    EXPECT_TRUE(is_neg_infinity(simp(pow(neg_inf(), lit(3)))));
}

TEST_F(ExtendedRealArithmeticTest, PowNegInfToPositiveEvenIsPosInf) {
    EXPECT_TRUE(is_pos_infinity(simp(pow(neg_inf(), lit(4)))));
}

TEST_F(ExtendedRealArithmeticTest, PowNegInfToZeroIsIndeterminate) {
    EXPECT_TRUE(is_indeterminate(simp(pow(neg_inf(), lit(0)))));
}

TEST_F(ExtendedRealArithmeticTest, PowComplexInfToPositiveIsComplexInf) {
    EXPECT_TRUE(is_complex_infinity(simp(pow(cplx_inf(), lit(2)))));
}

TEST_F(ExtendedRealArithmeticTest, PowComplexInfToNegativeIsZero) {
    ExprPtr r = simp(pow(cplx_inf(), lit(-1)));
    ASSERT_NE(r, nullptr);
    const auto* il = expr_cast<IntegerLit>(r);
    ASSERT_NE(il, nullptr);
    EXPECT_TRUE(il->value.is_zero());
}

TEST_F(ExtendedRealArithmeticTest, PowOneToComplexInfIsIndeterminate) {
    EXPECT_TRUE(is_indeterminate(simp(pow(lit(1), cplx_inf()))));
}

TEST_F(ExtendedRealArithmeticTest, NegInfTimesZeroChainIsIndeterminate) {
    // Verifies that the helper triggers even when one finite-zero operand
    // is buried in a longer product factor list.
    EXPECT_TRUE(is_indeterminate(simp(prod({lit(7), neg_inf(), lit(0)}))));
}

// ────────────────────────────────────────────────────────────────────────

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
