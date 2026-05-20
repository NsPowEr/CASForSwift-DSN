// CAS-L3-04 — LambertW principal branch W₀ tests.

#include <gtest/gtest.h>

#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class LambertWTest : public ::testing::Test {
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
};

TEST_F(LambertWTest, WZeroIsZero) {
    auto e = parse("LambertW(0)");
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    auto* il = expr_cast<IntegerLit>(r.value());
    ASSERT_NE(il, nullptr);
    EXPECT_TRUE(il->value.is_zero());
}

TEST_F(LambertWTest, WeIsOne) {
    auto e = parse("LambertW(e)");
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    auto* il = expr_cast<IntegerLit>(r.value());
    ASSERT_NE(il, nullptr);
    EXPECT_EQ(il->value, BigInt(1));
}

TEST_F(LambertWTest, AntiHardcodeSymbolicReturnsUnchanged) {
    // LambertW(x) symbolic → should NOT silently simplify.
    auto e = parse("LambertW(x)");
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    auto* fc = expr_cast<FuncCall>(r.value());
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->func_id, BuiltinOp::LambertW);
}

}  // namespace
