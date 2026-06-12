// CAS-L2-21 — Branch-cut policy global toggle test.
//
// Verifies that strict_branch_cuts() flag exists, defaults false,
// roundtrips via setter. Future integration: simplifier rules should
// consult this flag to refuse positivity-dependent identities under
// strict mode.

#include <gtest/gtest.h>

#include "cas/symbolic.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"

using namespace cas;

namespace {

TEST(BranchCutsGlobalTest, DefaultIsLenient) {
    symbolic::CASContext ctx;
    EXPECT_FALSE(ctx.strict_branch_cuts());
}

TEST(BranchCutsGlobalTest, SetStrictRoundtrips) {
    symbolic::CASContext ctx;
    ctx.set_strict_branch_cuts(true);
    EXPECT_TRUE(ctx.strict_branch_cuts());
    ctx.set_strict_branch_cuts(false);
    EXPECT_FALSE(ctx.strict_branch_cuts());
}

// F4.K — sqrt(x²) branch-cut gating tests (Branch_Cut_Propagation.md §3.6).

TEST(BranchCutsGlobalTest, SqrtOfSquare_KnownPositive_ReducesToX) {
    symbolic::CASContext ctx;
    ExprPtr x = ctx.arena().make<Symbol>("x");
    ctx.assumptions().assume_positive(*expr_cast<Symbol>(x));
    ExprPtr x_sq = ctx.arena().make<Binary>(BinaryOp::Pow, x,
        ctx.arena().make<IntegerLit>(BigInt(2)));
    ExprPtr expr = ctx.arena().make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{x_sq});
    auto r = ctx.simplify(expr);
    ASSERT_TRUE(r.is_ok());
    const auto* sym = expr_cast<Symbol>(r.value());
    EXPECT_TRUE(sym != nullptr && sym->name == "x");
}

TEST(BranchCutsGlobalTest, SqrtOfSquare_KnownReal_ReducesToAbs) {
    symbolic::CASContext ctx;
    ExprPtr x = ctx.arena().make<Symbol>("x");
    ctx.assumptions().assume_real(*expr_cast<Symbol>(x));
    ExprPtr x_sq = ctx.arena().make<Binary>(BinaryOp::Pow, x,
        ctx.arena().make<IntegerLit>(BigInt(2)));
    ExprPtr expr = ctx.arena().make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{x_sq});
    auto r = ctx.simplify(expr);
    ASSERT_TRUE(r.is_ok());
    const auto* fc = expr_cast<FuncCall>(r.value());
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->func_id, BuiltinOp::Abs);
}

TEST(BranchCutsGlobalTest, SqrtOfSquare_StrictMode_ComplexGeneric_PreservesStructure) {
    symbolic::CASContext ctx;
    ctx.set_strict_branch_cuts(true);
    ExprPtr x = ctx.arena().make<Symbol>("x");  // no assumption → may be complex
    ExprPtr x_sq = ctx.arena().make<Binary>(BinaryOp::Pow, x,
        ctx.arena().make<IntegerLit>(BigInt(2)));
    ExprPtr expr = ctx.arena().make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{x_sq});
    auto r = ctx.simplify(expr);
    ASSERT_TRUE(r.is_ok());
    // Must stay as sqrt(x^2) — no abs, no reduction to x.
    const auto* fc = expr_cast<FuncCall>(r.value());
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->func_id, BuiltinOp::Sqrt)
        << "strict mode must preserve sqrt(x^2) for complex generic x";
}

TEST(BranchCutsGlobalTest, SqrtOfSquare_LegacyMode_ComplexGeneric_EmitsAbs) {
    symbolic::CASContext ctx;
    // strict_branch_cuts left false (legacy default).
    ExprPtr x = ctx.arena().make<Symbol>("x");
    ExprPtr x_sq = ctx.arena().make<Binary>(BinaryOp::Pow, x,
        ctx.arena().make<IntegerLit>(BigInt(2)));
    ExprPtr expr = ctx.arena().make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{x_sq});
    auto r = ctx.simplify(expr);
    ASSERT_TRUE(r.is_ok());
    const auto* fc = expr_cast<FuncCall>(r.value());
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->func_id, BuiltinOp::Abs);
}

}  // namespace
