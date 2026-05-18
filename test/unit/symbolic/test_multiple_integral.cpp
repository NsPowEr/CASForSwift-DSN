#include "cas/calculus.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

using namespace cas;
using namespace cas::calculus;

namespace {

class MultipleIntegralTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto toks = Lexer(s).tokenize();
        EXPECT_TRUE(toks.is_ok()) << s;
        Parser p(toks.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }

    [[nodiscard]] ExprPtr make_int(int n) {
        return ctx.arena().make<IntegerLit>(BigInt(n));
    }

    // Simplify and check structural equality
    [[nodiscard]] bool equal_after_simplify(ExprPtr a, ExprPtr b) {
        auto sa = ctx.simplify(a);
        auto sb = ctx.simplify(b);
        if (!sa.is_ok() || !sb.is_ok()) return false;
        auto eq = symbolic::mathematically_equal(sa.value(), sb.value(), ctx);
        return eq.is_ok() && eq.value();
    }
};

} // namespace

// L2-09: basic double integral ∫_0^1 ∫_0^1 1 dx dy = 1
TEST_F(MultipleIntegralTest, UnitSquareArea) {
    Symbol x("x"), y("y");
    ExprPtr one = make_int(1);
    ExprPtr zero = make_int(0);
    ExprPtr one_e = make_int(1);

    auto res = multiple_integral(one, {
        {x, zero, one_e},
        {y, zero, one_e}
    }, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    EXPECT_TRUE(equal_after_simplify(res.value(), one))
        << "∫_0^1 ∫_0^1 1 dx dy should equal 1";
}

// L2-09: ∫_0^2 ∫_0^3 x*y dx dy = (3²/2)*(2²/2) = (9/2)*2 = 9
TEST_F(MultipleIntegralTest, DoubleIntegralXTimesY) {
    Symbol x("x"), y("y");
    ExprPtr f = parse("x*y");
    ExprPtr zero = make_int(0);
    ExprPtr two = make_int(2);
    ExprPtr three = make_int(3);

    auto res = multiple_integral(f, {
        {x, zero, three},  // ∫_0^3 x*y dx = y*(9/2)
        {y, zero, two}     // ∫_0^2 y*(9/2) dy = (9/2)*(4/2) = 9
    }, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    EXPECT_TRUE(equal_after_simplify(res.value(), make_int(9)))
        << "∫_0^2 ∫_0^3 x*y dx dy should equal 9";
}

// L2-09: triple integral ∫_0^1 ∫_0^1 ∫_0^1 1 dx dy dz = 1
TEST_F(MultipleIntegralTest, TripleIntegralUnitCube) {
    Symbol x("x"), y("y"), z("z");
    ExprPtr one = make_int(1);
    ExprPtr zero = make_int(0);
    ExprPtr one_e = make_int(1);

    auto res = multiple_integral(one, {
        {x, zero, one_e},
        {y, zero, one_e},
        {z, zero, one_e}
    }, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    EXPECT_TRUE(equal_after_simplify(res.value(), one))
        << "∫_0^1 ∫_0^1 ∫_0^1 1 dx dy dz should equal 1";
}

// L2-09: inner integral result depends on outer variable.
// ∫_0^1 ∫_0^y x dx dy = ∫_0^1 y²/2 dy = 1/6
TEST_F(MultipleIntegralTest, InnerBoundDependsOnOuterVariable) {
    Symbol x("x"), y("y");
    ExprPtr f = parse("x");
    ExprPtr zero = make_int(0);
    ExprPtr one = make_int(1);
    ExprPtr y_expr = ctx.arena().make<Symbol>("y");

    // Inner: ∫_0^y x dx  (y appears as upper bound)
    auto inner = definite_integral(f, x, zero, y_expr, ctx);
    ASSERT_TRUE(inner.is_ok()) << "inner integral failed: " << inner.error().message;

    // Outer: ∫_0^1 (result) dy
    auto outer = definite_integral(inner.value(), y, zero, one, ctx);
    ASSERT_TRUE(outer.is_ok()) << "outer integral failed: " << outer.error().message;

    ExprPtr expected = ctx.arena().make<RationalLit>(BigInt(1), BigInt(6));
    EXPECT_TRUE(equal_after_simplify(outer.value(), expected))
        << "∫_0^1 ∫_0^y x dx dy should equal 1/6";
}

// L2-09: Fubini swap on rectangular domain.
// ∫_0^1 ∫_0^1 x*y dx dy = ∫_0^1 ∫_0^1 x*y dy dx = 1/4
TEST_F(MultipleIntegralTest, FubiniSwapRectangular) {
    Symbol x("x"), y("y");
    ExprPtr f = parse("x*y");
    ExprPtr zero = make_int(0);
    ExprPtr one = make_int(1);

    auto res = fubini_swap(f, x, zero, one, y, zero, one, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    ExprPtr expected = ctx.arena().make<RationalLit>(BigInt(1), BigInt(4));
    EXPECT_TRUE(equal_after_simplify(res.value(), expected))
        << "∫_0^1 ∫_0^1 x*y dx dy (via Fubini) should equal 1/4";
}

// L2-09: Fubini rejects non-rectangular domains (y-dependent x bound).
TEST_F(MultipleIntegralTest, FubiniRejectsNonRectangularDomain) {
    Symbol x("x"), y("y");
    ExprPtr f = parse("x*y");
    ExprPtr zero = make_int(0);
    ExprPtr one = make_int(1);
    ExprPtr y_sym = ctx.arena().make<Symbol>("y"); // ax depends on y

    auto res = fubini_swap(f, x, zero, y_sym, y, zero, one, ctx);
    EXPECT_TRUE(res.is_error())
        << "fubini_swap should fail when x-bounds depend on y";
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
}

// L2-09 anti-hardcode: ∫_0^a ∫_0^b 6*x*y^2 dx dy = a²·b³ (symbolic bounds)
TEST_F(MultipleIntegralTest, SymbolicBoundsDoubleIntegral) {
    Symbol x("x"), y("y");
    Symbol a("a"), b("b");
    ExprPtr f = parse("6*x*y^2");
    ExprPtr zero = make_int(0);
    ExprPtr a_sym = ctx.arena().make<Symbol>("a");
    ExprPtr b_sym = ctx.arena().make<Symbol>("b");

    auto res = multiple_integral(f, {
        {x, zero, a_sym},
        {y, zero, b_sym}
    }, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    // Expected: ∫_0^b ∫_0^a 6*x*y^2 dx dy = ∫_0^b 3*a²*y^2 dy = a²*b³
    ExprPtr expected = parse("a^2 * b^3");
    EXPECT_TRUE(equal_after_simplify(res.value(), expected))
        << "∫_0^b ∫_0^a 6xy² dx dy should equal a²b³";
}
