// CAS-L2-26 — Piecewise expression evaluation MVP.
//
// piecewise(c1, e1, c2, e2, ..., default):
//   sequentially evaluate conditions; first that's exact True (IntegerLit≠0)
//   selects its branch; exact False (IntegerLit 0) drops the pair;
//   undetermined conditions keep their pair; if all dropped → default.

#include <gtest/gtest.h>

#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class PiecewiseTest : public ::testing::Test {
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
    [[nodiscard]] ExprPtr make_pw(std::vector<ExprPtr> args) {
        return ctx.arena().make<FuncCall>(BuiltinOp::Piecewise, std::move(args));
    }
};

TEST_F(PiecewiseTest, FirstTrueConditionSelectsBranch) {
    // piecewise(1, expr1, _, expr2, default) → expr1
    auto expr1 = parse("x + 1");
    auto expr2 = parse("x - 1");
    auto deflt = parse("0");
    auto cond_true = ctx.arena().make<IntegerLit>(BigInt(1));
    auto cond_false = ctx.arena().make<IntegerLit>(BigInt(0));
    auto pw = make_pw({cond_true, expr1, cond_false, expr2, deflt});
    auto s = ctx.simplify(pw);
    ASSERT_TRUE(s.is_ok());
    EXPECT_TRUE(structural_equal(s.value(), expr1));
}

TEST_F(PiecewiseTest, AllFalseReturnsDefault) {
    // piecewise(0, e1, 0, e2, default) → default
    auto e1 = parse("x + 1");
    auto e2 = parse("x - 1");
    auto deflt = parse("42");
    auto cond_false = ctx.arena().make<IntegerLit>(BigInt(0));
    auto pw = make_pw({cond_false, e1, cond_false, e2, deflt});
    auto s = ctx.simplify(pw);
    ASSERT_TRUE(s.is_ok());
    EXPECT_TRUE(structural_equal(s.value(), deflt));
}

TEST_F(PiecewiseTest, UndecidedConditionKept) {
    // piecewise(x>0, expr1, default) — x>0 is Binary, no decision → keep
    auto cond = parse("x > 0");  // depending on parser; might be Binary(Greater,...)
    auto expr1 = parse("x");
    auto deflt = parse("-x");
    auto pw = make_pw({cond, expr1, deflt});
    auto s = ctx.simplify(pw);
    ASSERT_TRUE(s.is_ok());
    // Should remain a Piecewise FuncCall.
    auto* fc = expr_cast<FuncCall>(s.value());
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->func_id, BuiltinOp::Piecewise);
}

TEST_F(PiecewiseTest, FalseBranchesDroppedKeepUndecided) {
    // piecewise(0, e_skip, x>0, e_keep, default)
    auto cond_false = ctx.arena().make<IntegerLit>(BigInt(0));
    auto e_skip = parse("999");
    auto cond_undec = parse("x > 0");
    auto e_keep = parse("x");
    auto deflt = parse("0");
    auto pw = make_pw({cond_false, e_skip, cond_undec, e_keep, deflt});
    auto s = ctx.simplify(pw);
    ASSERT_TRUE(s.is_ok());
    auto* fc = expr_cast<FuncCall>(s.value());
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->func_id, BuiltinOp::Piecewise);
    // Should have 3 args: undecided cond + branch + default.
    EXPECT_EQ(fc->args.size(), 3U);
}

TEST_F(PiecewiseTest, AntiHardcodeOnlyDefaultReturned) {
    // piecewise() with only default arg (1 arg, odd) returns default.
    auto deflt = parse("7");
    auto pw = make_pw({deflt});
    auto s = ctx.simplify(pw);
    ASSERT_TRUE(s.is_ok());
    EXPECT_TRUE(structural_equal(s.value(), deflt));
}

}  // namespace
