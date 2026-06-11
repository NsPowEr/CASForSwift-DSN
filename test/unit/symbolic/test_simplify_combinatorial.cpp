// F7.5 follow-up — exact-value identities for factorial, binomial,
// and erfc. Closes the bulk of the special_fn corpus regressions.

#include <gtest/gtest.h>

#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class CombinatorialTest : public ::testing::Test {
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

    [[nodiscard]] ExprPtr simp(const std::string& s) {
        auto r = ctx.simplify(parse(s));
        EXPECT_TRUE(r.is_ok()) << s;
        return r.is_ok() ? r.value() : nullptr;
    }

    [[nodiscard]] bool is_int(ExprPtr e, long long n) {
        if (!e) return false;
        const auto* lit = expr_cast<IntegerLit>(e);
        return lit != nullptr && lit->value == BigInt(n);
    }
};

TEST_F(CombinatorialTest, FactorialZeroAndOne) {
    EXPECT_TRUE(is_int(simp("factorial(0)"), 1));
    EXPECT_TRUE(is_int(simp("factorial(1)"), 1));
}

TEST_F(CombinatorialTest, FactorialSmallPositiveIntegers) {
    EXPECT_TRUE(is_int(simp("factorial(5)"), 120));
    EXPECT_TRUE(is_int(simp("factorial(7)"), 5040));
    EXPECT_TRUE(is_int(simp("factorial(10)"), 3628800));
}

TEST_F(CombinatorialTest, FactorialLargeIntegerStillExact) {
    // 15! = 1307674368000 — well within BigInt exact range.
    EXPECT_TRUE(is_int(simp("factorial(15)"), 1307674368000LL));
}

TEST_F(CombinatorialTest, FactorialNegativeIsComplexInfinity) {
    auto r = simp("factorial(-3)");
    ASSERT_NE(r, nullptr);
    const auto* c = expr_cast<Constant>(r);
    ASSERT_NE(c, nullptr);
    EXPECT_EQ(c->value, MathConstant::ComplexInfinity);
}

TEST_F(CombinatorialTest, BinomialBaseCases) {
    EXPECT_TRUE(is_int(simp("binomial(5, 0)"), 1));
    EXPECT_TRUE(is_int(simp("binomial(5, 5)"), 1));
    EXPECT_TRUE(is_int(simp("binomial(10, 1)"), 10));
}

TEST_F(CombinatorialTest, BinomialIntegerComputation) {
    EXPECT_TRUE(is_int(simp("binomial(5, 2)"), 10));
    EXPECT_TRUE(is_int(simp("binomial(10, 3)"), 120));
    EXPECT_TRUE(is_int(simp("binomial(8, 4)"), 70));
    EXPECT_TRUE(is_int(simp("binomial(10, 5)"), 252));
}

TEST_F(CombinatorialTest, BinomialOutOfRangeIsZero) {
    EXPECT_TRUE(is_int(simp("binomial(3, 5)"), 0));
}

TEST_F(CombinatorialTest, BinomialSymbolicNZeroK) {
    EXPECT_TRUE(is_int(simp("binomial(n, 0)"), 1));
}

TEST_F(CombinatorialTest, BinomialSymbolicNOneK) {
    // binomial(n, 1) = n (any n)
    auto r = simp("binomial(n, 1)");
    auto expected = parse("n");
    auto eq = cas::symbolic::mathematically_equal(r, expected, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

TEST_F(CombinatorialTest, BinomialSymbolicNAndKEqualsTwo) {
    // binomial(n, 2) = n(n-1)/2.
    auto r = simp("binomial(n, 2)");
    auto expected = simp("n*(n-1)/2");
    auto eq = cas::symbolic::mathematically_equal(r, expected, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

TEST_F(CombinatorialTest, ErfcZeroEqualsOne) {
    EXPECT_TRUE(is_int(simp("erfc(0)"), 1));
}

TEST_F(CombinatorialTest, ErfcAtInfinity) {
    auto inf = ctx.arena().make<Constant>(MathConstant::Infinity);
    auto erfc_inf = ctx.arena().make<FuncCall>(
        BuiltinOp::Erfc, std::vector<ExprPtr>{inf});
    auto r = ctx.simplify(erfc_inf);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(is_int(r.value(), 0));
}

TEST_F(CombinatorialTest, ErfcAtNegativeInfinity) {
    auto neg_inf = ctx.arena().make<Unary>(UnaryOp::Neg,
        ctx.arena().make<Constant>(MathConstant::Infinity));
    auto erfc_neg = ctx.arena().make<FuncCall>(
        BuiltinOp::Erfc, std::vector<ExprPtr>{neg_inf});
    auto r = ctx.simplify(erfc_neg);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(is_int(r.value(), 2));
}

}  // namespace
