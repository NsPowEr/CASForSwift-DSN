// CAS-L2-21 — Branch-cut policy global toggle test.
//
// Verifies that strict_branch_cuts() flag exists, defaults false,
// roundtrips via setter. Future integration: simplifier rules should
// consult this flag to refuse positivity-dependent identities under
// strict mode.

#include <gtest/gtest.h>
#include <functional>
#include <vector>

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

TEST(BranchCutsGlobalTest, SqrtOfSquare_StrictMode_ComplexGeneric_EmitsUnwindingCorrection) {
    // F8.0-6.2 / Task 20 BC-1 (Branch_Cut_Propagation.md §2 rule 1):
    //   sqrt(z²) = z · (-1)^K(2·ln(z))
    // strict_branch_cuts must produce the explicit K(·) correction so the
    // identity stays algebraically exact in the complex plane.
    symbolic::CASContext ctx;
    ctx.set_strict_branch_cuts(true);
    ExprPtr x = ctx.arena().make<Symbol>("x");  // no assumption → may be complex
    ExprPtr x_sq = ctx.arena().make<Binary>(BinaryOp::Pow, x,
        ctx.arena().make<IntegerLit>(BigInt(2)));
    ExprPtr expr = ctx.arena().make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{x_sq});
    auto r = ctx.simplify(expr);
    ASSERT_TRUE(r.is_ok());
    // Strict mode must NOT reduce to plain Symbol(x) or Abs(x); the result
    // must mention UnwindingNumber to evidence the correction.
    bool found_K = false;
    std::function<void(ExprPtr)> walk = [&](ExprPtr e) {
        if (!e) return;
        if (const auto* fc = expr_cast<FuncCall>(e);
            fc && fc->func_id == BuiltinOp::UnwindingNumber) {
            found_K = true;
            return;
        }
        if (const auto* bin = expr_cast<Binary>(e)) { walk(bin->left); walk(bin->right); }
        if (const auto* un = expr_cast<Unary>(e)) walk(un->operand);
        if (const auto* prod = expr_cast<Product>(e))
            for (auto f : prod->factors) walk(f);
        if (const auto* sum = expr_cast<Sum>(e))
            for (auto t : sum->terms) walk(t);
        if (const auto* fc = expr_cast<FuncCall>(e))
            for (auto a : fc->args) walk(a);
    };
    walk(r.value());
    EXPECT_TRUE(found_K)
        << "strict mode must emit K(·) correction for sqrt(x^2) with complex generic x";
}

// F8.0-6.2 / Task 20 BC-2 — (z^a)^b unwinding correction under strict mode.
TEST(BranchCutsGlobalTest, PowOfPow_StrictMode_ComplexGeneric_EmitsUnwindingCorrection) {
    // (z^2)^(1/3) with strict_branch_cuts: must NOT silently flatten to
    // z^(2/3) — that drops the K(·) correction needed for complex generic z.
    symbolic::CASContext ctx;
    ctx.set_strict_branch_cuts(true);
    ExprPtr z = ctx.arena().make<Symbol>("z");  // generic complex
    ExprPtr two = ctx.arena().make<IntegerLit>(BigInt(2));
    ExprPtr one_third = ctx.arena().make<RationalLit>(BigInt(1), BigInt(3));
    ExprPtr inner = ctx.arena().make<Binary>(BinaryOp::Pow, z, two);
    ExprPtr outer = ctx.arena().make<Binary>(BinaryOp::Pow, inner, one_third);
    auto r = ctx.simplify(outer);
    ASSERT_TRUE(r.is_ok());
    // Result must mention UnwindingNumber.
    bool found_K = false;
    std::function<void(ExprPtr)> walk = [&](ExprPtr e) {
        if (!e) return;
        if (const auto* fc = expr_cast<FuncCall>(e);
            fc && fc->func_id == BuiltinOp::UnwindingNumber) { found_K = true; return; }
        if (const auto* bin = expr_cast<Binary>(e)) { walk(bin->left); walk(bin->right); }
        if (const auto* un = expr_cast<Unary>(e)) walk(un->operand);
        if (const auto* prod = expr_cast<Product>(e))
            for (auto f : prod->factors) walk(f);
        if (const auto* sum = expr_cast<Sum>(e))
            for (auto t : sum->terms) walk(t);
        if (const auto* fc = expr_cast<FuncCall>(e))
            for (auto a : fc->args) walk(a);
    };
    walk(r.value());
    EXPECT_TRUE(found_K)
        << "strict mode must emit K(·) correction for (z^a)^b with complex generic z";
}

// F8.0-6.2 / Task 20 BC-2 legacy — without strict, (z^(1/2))^(1/2) still flattens.
TEST(BranchCutsGlobalTest, PowOfPow_LegacyMode_ComplexGeneric_FlattensNaively) {
    symbolic::CASContext ctx;
    // strict_branch_cuts left false (legacy)
    ExprPtr z = ctx.arena().make<Symbol>("z");
    ExprPtr two = ctx.arena().make<IntegerLit>(BigInt(2));
    ExprPtr one_third = ctx.arena().make<RationalLit>(BigInt(1), BigInt(3));
    ExprPtr inner = ctx.arena().make<Binary>(BinaryOp::Pow, z, two);
    ExprPtr outer = ctx.arena().make<Binary>(BinaryOp::Pow, inner, one_third);
    auto r = ctx.simplify(outer);
    ASSERT_TRUE(r.is_ok());
    // Legacy must NOT emit K(·).
    bool found_K = false;
    std::function<void(ExprPtr)> walk = [&](ExprPtr e) {
        if (!e) return;
        if (const auto* fc = expr_cast<FuncCall>(e);
            fc && fc->func_id == BuiltinOp::UnwindingNumber) { found_K = true; return; }
        if (const auto* bin = expr_cast<Binary>(e)) { walk(bin->left); walk(bin->right); }
        if (const auto* un = expr_cast<Unary>(e)) walk(un->operand);
        if (const auto* prod = expr_cast<Product>(e))
            for (auto f : prod->factors) walk(f);
        if (const auto* sum = expr_cast<Sum>(e))
            for (auto t : sum->terms) walk(t);
        if (const auto* fc = expr_cast<FuncCall>(e))
            for (auto a : fc->args) walk(a);
    };
    walk(r.value());
    EXPECT_FALSE(found_K)
        << "legacy mode must not introduce K(·)";
}

// F8.0-6.2 / Task 20 BC-3 — ln(z1·z2) unwinding correction under strict mode.
TEST(BranchCutsGlobalTest, LnOfProduct_StrictMode_ComplexGeneric_EmitsUnwindingCorrection) {
    // ln(z1·z2) = ln(z1) + ln(z2) - 2πi·K(ln(z1)+ln(z2)) for complex generic z.
    symbolic::CASContext ctx;
    ctx.set_strict_branch_cuts(true);
    ExprPtr z1 = ctx.arena().make<Symbol>("z1");
    ExprPtr z2 = ctx.arena().make<Symbol>("z2");
    ExprPtr prod = ctx.arena().make<Product>(std::vector<ExprPtr>{z1, z2});
    ExprPtr expr = ctx.arena().make<FuncCall>(BuiltinOp::Ln,
        std::vector<ExprPtr>{prod});
    auto r = ctx.simplify(expr);
    ASSERT_TRUE(r.is_ok());
    bool found_K = false;
    std::function<void(ExprPtr)> walk = [&](ExprPtr e) {
        if (!e) return;
        if (const auto* fc = expr_cast<FuncCall>(e);
            fc && fc->func_id == BuiltinOp::UnwindingNumber) { found_K = true; return; }
        if (const auto* bin = expr_cast<Binary>(e)) { walk(bin->left); walk(bin->right); }
        if (const auto* un = expr_cast<Unary>(e)) walk(un->operand);
        if (const auto* prod_n = expr_cast<Product>(e))
            for (auto f : prod_n->factors) walk(f);
        if (const auto* sum = expr_cast<Sum>(e))
            for (auto t : sum->terms) walk(t);
        if (const auto* fc = expr_cast<FuncCall>(e))
            for (auto a : fc->args) walk(a);
    };
    walk(r.value());
    EXPECT_TRUE(found_K)
        << "strict mode must emit K(·) correction for ln(z1·z2) with complex generic z";
}

// F8.0-6.2 / Task 20 BC-3 — ln(z1·z2) with provably positive factors must
// still reduce to ln(z1) + ln(z2) even under strict mode (no spurious K(·)).
TEST(BranchCutsGlobalTest, LnOfProduct_StrictMode_AllPositive_ReducesCleanly) {
    symbolic::CASContext ctx;
    ctx.set_strict_branch_cuts(true);
    ExprPtr a = ctx.arena().make<Symbol>("a");
    ExprPtr b = ctx.arena().make<Symbol>("b");
    ctx.assumptions().assume_positive(*expr_cast<Symbol>(a));
    ctx.assumptions().assume_positive(*expr_cast<Symbol>(b));
    ExprPtr prod = ctx.arena().make<Product>(std::vector<ExprPtr>{a, b});
    ExprPtr expr = ctx.arena().make<FuncCall>(BuiltinOp::Ln,
        std::vector<ExprPtr>{prod});
    auto r = ctx.simplify(expr);
    ASSERT_TRUE(r.is_ok());
    bool found_K = false;
    std::function<void(ExprPtr)> walk = [&](ExprPtr e) {
        if (!e) return;
        if (const auto* fc = expr_cast<FuncCall>(e);
            fc && fc->func_id == BuiltinOp::UnwindingNumber) { found_K = true; return; }
        if (const auto* bin = expr_cast<Binary>(e)) { walk(bin->left); walk(bin->right); }
        if (const auto* un = expr_cast<Unary>(e)) walk(un->operand);
        if (const auto* prod_n = expr_cast<Product>(e))
            for (auto f : prod_n->factors) walk(f);
        if (const auto* sum = expr_cast<Sum>(e))
            for (auto t : sum->terms) walk(t);
        if (const auto* fc = expr_cast<FuncCall>(e))
            for (auto a : fc->args) walk(a);
    };
    walk(r.value());
    EXPECT_FALSE(found_K)
        << "all-positive factors must not require K(·) correction";
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
