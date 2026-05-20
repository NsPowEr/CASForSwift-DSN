#include <gtest/gtest.h>
#include "cas/symbolic.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"

using namespace cas;

namespace {

Result<ExprPtr> parse(const std::string& s, AstArena& arena) {
    auto tok = Lexer(s).tokenize();
    if (tok.is_error()) return fail<ExprPtr>(tok.error());
    return Parser(tok.value(), arena).parse();
}

bool math_eq(ExprPtr a, ExprPtr b, symbolic::CASContext& ctx) {
    auto r = symbolic::mathematically_equal(a, b, ctx);
    return r.is_ok() && r.value();
}

} // namespace

TEST(CalculusLimitTest, SimplifyFractionEquivalence) {
    symbolic::CASContext ctx;

    auto e1 = parse("(1/(1+1/n) * (-1*n^-2)) / (-1*n^-2)", ctx.arena());
    ASSERT_TRUE(e1.is_ok()) << e1.error().message;
    auto s1 = ctx.simplify(e1.value());
    ASSERT_TRUE(s1.is_ok()) << s1.error().message;

    auto e2 = parse("n/(n+1)", ctx.arena());
    ASSERT_TRUE(e2.is_ok()) << e2.error().message;
    auto s2 = ctx.simplify(e2.value());
    ASSERT_TRUE(s2.is_ok()) << s2.error().message;

    EXPECT_TRUE(math_eq(s1.value(), s2.value(), ctx))
        << "(1/(1+1/n)*(-n^-2))/(-n^-2) should equal n/(n+1)";
}

TEST(CalculusLimitTest, LimitOneOverNAtInfinity) {
    symbolic::CASContext ctx;

    auto expr = parse("1/n", ctx.arena());
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    ExprPtr inf = ctx.arena().make<Constant>(MathConstant::Infinity);
    auto lim = calculus::limit(expr.value(), Symbol{"n"}, inf, LimitDirection::Both, ctx);
    ASSERT_TRUE(lim.is_ok()) << lim.error().message;

    auto zero = parse("0", ctx.arena());
    ASSERT_TRUE(zero.is_ok());

    EXPECT_TRUE(math_eq(lim.value(), zero.value(), ctx))
        << "lim(n->inf, 1/n) should equal 0";
}

TEST(CalculusLimitTest, LimitSinXOverXAtZero) {
    symbolic::CASContext ctx;

    auto expr = parse("sin(x)/x", ctx.arena());
    ASSERT_TRUE(expr.is_ok()) << expr.error().message;

    auto zero = ctx.arena().make<IntegerLit>(BigInt(0));
    auto lim = calculus::limit(expr.value(), Symbol{"x"}, zero, LimitDirection::Both, ctx);
    ASSERT_TRUE(lim.is_ok()) << lim.error().message;

    auto one = parse("1", ctx.arena());
    ASSERT_TRUE(one.is_ok());

    EXPECT_TRUE(math_eq(lim.value(), one.value(), ctx))
        << "lim(x->0, sin(x)/x) should equal 1";
}
