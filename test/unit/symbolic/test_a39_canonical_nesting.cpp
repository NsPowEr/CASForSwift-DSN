// A39 — simplify must not leave nested Product/Sum nodes in its result.
//
// The engine has a strict-canonicity checker (is_strictly_canonical,
// simplify_utils.cpp) wired to a debug-only canary in CASContext::simplify.
// It used to fire 17 times across the Mellin/MeijerG suites, from two distinct
// producers; these tests pin both so the shapes cannot come back.
//
// Correctness was never at stake (results compare equal via
// mathematically_equal) — canonical form was, and a non-canonical result
// blocks structural comparison and defeats cancellation: u * u^-1 survives
// when the second u sits inside a nested Product.

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class A39CanonicalNestingTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }

    // True if any Product has a Product factor, or any Sum has a Sum term.
    [[nodiscard]] static bool has_nested_same_kind(ExprPtr e) {
        if (!e) return false;
        if (const auto* pr = expr_cast<Product>(e)) {
            for (ExprPtr f : pr->factors) {
                if (expr_is<Product>(f)) return true;
                if (has_nested_same_kind(f)) return true;
            }
            return false;
        }
        if (const auto* sm = expr_cast<Sum>(e)) {
            for (ExprPtr t : sm->terms) {
                if (expr_is<Sum>(t)) return true;
                if (has_nested_same_kind(t)) return true;
            }
            return false;
        }
        if (const auto* un = expr_cast<Unary>(e)) return has_nested_same_kind(un->operand);
        if (const auto* bi = expr_cast<Binary>(e))
            return has_nested_same_kind(bi->left) || has_nested_same_kind(bi->right);
        if (const auto* fc = expr_cast<FuncCall>(e)) {
            for (ExprPtr a : fc->args) if (has_nested_same_kind(a)) return true;
            return false;
        }
        return false;
    }
};

// Producer 1: the Product flatten did not descend through Unary(Neg, ...), so
// Neg(Product(...)) landed whole in the factor list and its inner Product was
// never flattened. Repro: (-P)*pi simplified to -(pi*P) with P still nested.
TEST_F(A39CanonicalNestingTest, NegatedProductFactorIsFlattened) {
    AstArena& a = ctx.arena();
    ExprPtr x = a.make<Symbol>("x");
    ExprPtr inner = a.make<Product>(std::vector<ExprPtr>{
        a.make<Binary>(BinaryOp::Pow, x, a.make<IntegerLit>(BigInt(2))),
        a.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{x})});
    ExprPtr e = a.make<Product>(std::vector<ExprPtr>{
        a.make<Unary>(UnaryOp::Neg, inner), a.make<Constant>(MathConstant::Pi)});
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok()) << s.error().message;
    EXPECT_FALSE(has_nested_same_kind(s.value()));
}

// The nesting also defeated cancellation: with u^-1 * (u * x) the two u's never
// met, so they did not cancel. Sign must still come out right.
TEST_F(A39CanonicalNestingTest, NestingUnderNegStillCancels) {
    AstArena& a = ctx.arena();
    ExprPtr u = a.make<Symbol>("u");
    ExprPtr x = a.make<Symbol>("x");
    ExprPtr inner = a.make<Product>(std::vector<ExprPtr>{u, x});
    ExprPtr e = a.make<Product>(std::vector<ExprPtr>{
        a.make<Binary>(BinaryOp::Pow, u, a.make<IntegerLit>(BigInt(-1))),
        a.make<Unary>(UnaryOp::Neg, inner)});
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok()) << s.error().message;
    EXPECT_FALSE(has_nested_same_kind(s.value()));
    auto eq = symbolic::mathematically_equal(s.value(), parse("-x"), ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

// Producer 2: collecting like terms over a shared factor built Sum([qi, qj])
// without splicing a quotient that was itself a Sum.
// (1 - 2x^2)*e^x - e^x  ->  the coefficient sum must stay flat.
TEST_F(A39CanonicalNestingTest, LikeTermCollectionSplicesSumQuotients) {
    auto e = parse("(1 - 2*x^2) * exp(x) - exp(x)");
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok()) << s.error().message;
    EXPECT_FALSE(has_nested_same_kind(s.value()));
    auto eq = symbolic::mathematically_equal(s.value(), parse("-2*x^2*exp(x)"), ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

// Idempotency: simplifying an already-simplified result must not reintroduce
// nesting nor change the value.
TEST_F(A39CanonicalNestingTest, ResimplifyIsStable) {
    for (const char* src : {"(1 - 2*x^2) * exp(x) - exp(x)",
                            "(3*x - 1) * sin(x) - 2 * sin(x)",
                            "-(x^2 * exp(x)) * 5"}) {
        auto once = ctx.simplify(parse(src));
        ASSERT_TRUE(once.is_ok()) << src;
        auto twice = ctx.simplify(once.value());
        ASSERT_TRUE(twice.is_ok()) << src;
        EXPECT_FALSE(has_nested_same_kind(twice.value())) << src;
        auto eq = symbolic::mathematically_equal(once.value(), twice.value(), ctx);
        ASSERT_TRUE(eq.is_ok()) << src;
        EXPECT_TRUE(eq.value()) << src;
    }
}

}  // namespace
