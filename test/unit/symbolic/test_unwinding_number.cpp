// F8.0-6.1: tests for the UnwindingNumber K(z) builtin.
// K(z) = ⌈(Im(z) − π) / (2π)⌉  encodes the branch-cut wrap so that
//   ln(exp(z)) = z + 2πi · K(z).
// At this stage the node is registered as a builtin and parsed; arithmetic
// propagation through log/exp composition is reserved for Task 6.2.

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <functional>
#include <string>

namespace cas {
namespace {

[[nodiscard]] ExprPtr parse_expr(symbolic::CASContext& ctx, const std::string& s) {
    auto tokens = Lexer(s).tokenize();
    if (tokens.is_error()) return nullptr;
    Parser p(tokens.value(), ctx.arena());
    auto e = p.parse();
    return e.is_ok() ? e.value() : nullptr;
}

TEST(UnwindingNumber, ParsesShortName_K) {
    symbolic::CASContext ctx;
    ExprPtr e = parse_expr(ctx, "K(z)");
    ASSERT_NE(e, nullptr);
    const auto* call = expr_cast<FuncCall>(e);
    ASSERT_NE(call, nullptr);
    EXPECT_EQ(call->func_id, BuiltinOp::UnwindingNumber);
    ASSERT_EQ(call->args.size(), 1U);
    const auto* z = expr_cast<Symbol>(call->args[0]);
    ASSERT_NE(z, nullptr);
    EXPECT_EQ(z->name, "z");
}

TEST(UnwindingNumber, ParsesLongName_UnwindingNumber) {
    symbolic::CASContext ctx;
    ExprPtr e = parse_expr(ctx, "UnwindingNumber(z)");
    ASSERT_NE(e, nullptr);
    const auto* call = expr_cast<FuncCall>(e);
    ASSERT_NE(call, nullptr);
    EXPECT_EQ(call->func_id, BuiltinOp::UnwindingNumber);
}

TEST(UnwindingNumber, SimplifyPreservesIdentity) {
    // The simplifier must not collapse K(z) to anything else — branch-cut
    // bookkeeping is non-arithmetic and must survive verbatim until Task 6.2
    // implements log/exp branch-aware propagation.
    symbolic::CASContext ctx;
    ExprPtr e = parse_expr(ctx, "K(z)");
    ASSERT_NE(e, nullptr);

    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok()) << s.error().message;
    const auto* call = expr_cast<FuncCall>(s.value());
    ASSERT_NE(call, nullptr);
    EXPECT_EQ(call->func_id, BuiltinOp::UnwindingNumber);
}

TEST(UnwindingNumber, BuiltinNameRoundTrip) {
    // The string-name mapping must round-trip: name → enum → name.
    EXPECT_EQ(get_builtin_op("K"), BuiltinOp::UnwindingNumber);
    EXPECT_EQ(get_builtin_op("UnwindingNumber"), BuiltinOp::UnwindingNumber);
    EXPECT_EQ(builtin_op_name(BuiltinOp::UnwindingNumber), "UnwindingNumber");
}

TEST(UnwindingNumber, BranchCutAware_LnExp_NonRealZ_PropagatesK) {
    // Task 6.2 opt-in: with branch_cut_aware_logexp() = true and z not
    // declared real, ln(exp(z)) → z + 2πi·K(z).
    symbolic::CASContext ctx;
    ctx.set_branch_cut_aware_logexp(true);
    ExprPtr e = parse_expr(ctx, "ln(exp(z))");
    ASSERT_NE(e, nullptr);
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok()) << s.error().message;

    // The simplified form must contain K(z) as a sub-expression.
    bool found_K = false;
    std::function<void(ExprPtr)> walk = [&](ExprPtr x) {
        if (!x) return;
        if (const auto* call = expr_cast<FuncCall>(x)) {
            if (call->func_id == BuiltinOp::UnwindingNumber) found_K = true;
            for (auto a : call->args) walk(a);
        }
        if (const auto* bin = expr_cast<Binary>(x)) {
            walk(bin->left); walk(bin->right);
        }
        if (const auto* sum = expr_cast<Sum>(x))
            for (auto t : sum->terms) walk(t);
        if (const auto* prod = expr_cast<Product>(x))
            for (auto f : prod->factors) walk(f);
        if (const auto* un = expr_cast<Unary>(x)) walk(un->operand);
    };
    walk(s.value());
    EXPECT_TRUE(found_K)
        << "expected ln(exp(z)) → z + 2πi·K(z) in branch-aware mode";
}

TEST(UnwindingNumber, BranchCutAware_LnExp_RealZ_StaysZ) {
    // Realness assumption suppresses K(z): ln(exp(x)) → x for real x.
    symbolic::CASContext ctx;
    ctx.set_branch_cut_aware_logexp(true);
    ctx.assumptions().assume_real(Symbol("x"));

    ExprPtr e = parse_expr(ctx, "ln(exp(x))");
    ASSERT_NE(e, nullptr);
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok()) << s.error().message;

    // K(x) must NOT appear.
    bool found_K = false;
    std::function<void(ExprPtr)> walk = [&](ExprPtr x) {
        if (!x) return;
        if (const auto* call = expr_cast<FuncCall>(x)) {
            if (call->func_id == BuiltinOp::UnwindingNumber) found_K = true;
            for (auto a : call->args) walk(a);
        }
        if (const auto* bin = expr_cast<Binary>(x)) {
            walk(bin->left); walk(bin->right);
        }
        if (const auto* sum = expr_cast<Sum>(x))
            for (auto t : sum->terms) walk(t);
        if (const auto* prod = expr_cast<Product>(x))
            for (auto f : prod->factors) walk(f);
    };
    walk(s.value());
    EXPECT_FALSE(found_K) << "real-x assumption must suppress K(x)";
    // And result must be x.
    const auto* sym = expr_cast<Symbol>(s.value());
    ASSERT_NE(sym, nullptr) << "expected the result to collapse to x";
    EXPECT_EQ(sym->name, "x");
}

TEST(UnwindingNumber, BranchCutAware_Default_Off_LegacyBehaviour) {
    // Default = off: ln(exp(z)) → z (legacy), regardless of z realness.
    symbolic::CASContext ctx;
    ExprPtr e = parse_expr(ctx, "ln(exp(z))");
    ASSERT_NE(e, nullptr);
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok()) << s.error().message;

    const auto* sym = expr_cast<Symbol>(s.value());
    ASSERT_NE(sym, nullptr);
    EXPECT_EQ(sym->name, "z");
}

TEST(UnwindingNumber, StructuralEqualityWithDistinctArguments) {
    symbolic::CASContext ctx;
    ExprPtr a = parse_expr(ctx, "K(z)");
    ExprPtr b = parse_expr(ctx, "K(z)");
    ExprPtr c = parse_expr(ctx, "K(w)");
    ASSERT_NE(a, nullptr);
    ASSERT_NE(b, nullptr);
    ASSERT_NE(c, nullptr);
    EXPECT_TRUE(structural_equal(a, b));
    EXPECT_FALSE(structural_equal(a, c));
}

} // namespace
} // namespace cas
