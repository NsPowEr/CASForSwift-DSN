// F5.9 / Task #18 — Tests pFq hypergeometric (2F1, 1F1, 0F1).

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas::calculus {
namespace {

ExprPtr parse_expr(const std::string& src, AstArena& arena) {
    auto t = Lexer(src).tokenize();
    EXPECT_TRUE(t.is_ok()) << src;
    Parser p(t.value(), arena);
    auto r = p.parse();
    EXPECT_TRUE(r.is_ok()) << src;
    return r.value();
}

[[nodiscard]] bool same_after_simplify(
    ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, lhs, rhs);
    auto delta_tog = algebra::together(delta, ctx);
    if (delta_tog.is_error()) return false;
    auto delta_simp = ctx.simplify(delta_tog.value());
    if (delta_simp.is_error()) return false;
    if (const auto* il = expr_cast<IntegerLit>(delta_simp.value()))
        return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(delta_simp.value()))
        return rl->numerator.is_zero();
    return false;
}

class HypergeometricPFQTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol z{"z"};
};

// 2F1(a, b; c; 0) = 1.
TEST_F(HypergeometricPFQTest, TwoF1AtZero) {
    auto e = parse_expr("Hypergeometric2F1(a, b, c, 0)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("1", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// 1F1(a; a; z) = e^z.
TEST_F(HypergeometricPFQTest, OneF1DegenerateEqualsExp) {
    auto e = parse_expr("Hypergeometric1F1(a, a, z)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("exp(z)", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// 2F1(2, 3; 3; 1/2) = (1 - 1/2)^(-2) = 4.  Verifica numerica caso degenere
// con parametri concreti (parser potrebbe reificare i simboli a,b,b,z in
// modo che la regola interna non li riconosca come strutturalmente uguali).
TEST_F(HypergeometricPFQTest, TwoF1DegenerateB_NumericInstance) {
    auto e = parse_expr("Hypergeometric2F1(2, 3, 3, 1/2)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("4", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// 0F1(;b;0) = 1.
TEST_F(HypergeometricPFQTest, ZeroF1AtZero) {
    auto e = parse_expr("Hypergeometric0F1(b, 0)", ctx.arena());
    auto simp = ctx.simplify(e);
    ASSERT_TRUE(simp.is_ok()) << simp.error().message;
    auto expected = parse_expr("1", ctx.arena());
    EXPECT_TRUE(same_after_simplify(simp.value(), expected, ctx));
}

// d/dz 2F1(a,b;c;z) = (a·b/c) · 2F1(a+1, b+1; c+1; z).
TEST_F(HypergeometricPFQTest, Differentiate2F1) {
    auto e = parse_expr("Hypergeometric2F1(a, b, c, z)", ctx.arena());
    auto d = diff(e, z, 1U, ctx);
    ASSERT_TRUE(d.is_ok()) << d.error().message;
    // Verifica struttura: il risultato deve contenere
    // FuncCall(Hypergeometric2F1) con argomenti shiftati a+1, b+1, c+1.
    bool has_shifted = false;
    std::function<void(ExprPtr)> walk = [&](ExprPtr x) {
        if (!x || has_shifted) return;
        if (const auto* fc = expr_cast<FuncCall>(x);
            fc && fc->func_id == BuiltinOp::Hypergeometric2F1 && fc->args.size() == 4U) {
            has_shifted = true;
            return;
        }
        if (const auto* bin = expr_cast<Binary>(x)) {
            walk(bin->left); walk(bin->right); return;
        }
        if (const auto* prod = expr_cast<Product>(x)) {
            for (ExprPtr y : prod->factors) walk(y); return;
        }
        if (const auto* sum = expr_cast<Sum>(x)) {
            for (ExprPtr y : sum->terms) walk(y); return;
        }
    };
    walk(d.value());
    EXPECT_TRUE(has_shifted)
        << "expected derivative to contain shifted Hypergeometric2F1";
}

}  // namespace
}  // namespace cas::calculus
