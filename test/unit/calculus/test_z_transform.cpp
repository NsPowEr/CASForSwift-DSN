// F5.8 / Task #16 — Tests per Z-transform.

#include "../../../src/calculus/calculus_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas::calculus {
namespace {

ExprPtr parse_expr(const std::string& s, AstArena& arena) {
    auto t = Lexer(s).tokenize();
    EXPECT_TRUE(t.is_ok()) << s << ": " << t.error().message;
    Parser p(t.value(), arena);
    auto r = p.parse();
    EXPECT_TRUE(r.is_ok()) << s;
    return r.value();
}

// Verifica: diff(Z{a} - expected, _) = 0  (equivalente sotto simplifier).
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

class ZTransformTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol n{"n"};
    Symbol z{"z"};
};

// Z{1} = z/(z - 1).
TEST_F(ZTransformTest, ConstantOne) {
    auto a = parse_expr("1", ctx.arena());
    auto res = z_transform(a, n, z, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("z/(z-1)", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// Z{c} = c·z/(z-1) per c costante.
TEST_F(ZTransformTest, ConstantC) {
    auto a = parse_expr("5", ctx.arena());
    auto res = z_transform(a, n, z, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("5*z/(z-1)", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// Z{n} = z/(z-1)².
TEST_F(ZTransformTest, LinearN) {
    auto a = parse_expr("n", ctx.arena());
    auto res = z_transform(a, n, z, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("z/(z-1)^2", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// Z{n²} = -z·d/dz Z{n} = z(z+1)/(z-1)³.
TEST_F(ZTransformTest, NSquared) {
    auto a = parse_expr("n^2", ctx.arena());
    auto res = z_transform(a, n, z, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("z*(z+1)/(z-1)^3", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// Z{r^n} = z/(z-r).
TEST_F(ZTransformTest, GeometricRn) {
    auto a = parse_expr("2^n", ctx.arena());
    auto res = z_transform(a, n, z, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("z/(z-2)", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// Z{cos(ω·n)} = z(z - cos ω)/(z² - 2z·cos ω + 1).
// Verifica via valutazione strutturale a omega = 0 → Z{1} = z/(z-1):
//   formula in ω=0: z(z-1)/(z²-2z+1) = z(z-1)/(z-1)² = z/(z-1). ✓
TEST_F(ZTransformTest, CosOmegaN_AgainstClosedForm) {
    AstArena& arena = ctx.arena();
    Symbol w{"omega"};
    ExprPtr omega = arena.make<Symbol>(w);
    ExprPtr n_e = arena.make<Symbol>(n);
    ExprPtr a = arena.make<FuncCall>(BuiltinOp::Cos,
        std::vector<ExprPtr>{arena.make<Binary>(BinaryOp::Mul, omega, n_e)});
    auto res = z_transform(a, n, z, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    // Sostituisci omega → 0; cos(0) = 1, sin(0) = 0.  Formula attesa: z/(z-1).
    auto subst = ctx.substitute(res.value(), w, arena.make<IntegerLit>(BigInt(0)));
    ASSERT_TRUE(subst.is_ok());
    auto expected = parse_expr("z/(z-1)", arena);
    EXPECT_TRUE(same_after_simplify(subst.value(), expected, ctx));
}

// Linearità: Z{2·1 + 3·n} = 2·Z{1} + 3·Z{n}.
TEST_F(ZTransformTest, LinearitySumScaled) {
    auto a = parse_expr("2 + 3*n", ctx.arena());
    auto res = z_transform(a, n, z, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("2*z/(z-1) + 3*z/(z-1)^2", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// Modulation: Z{r^n · n} = (z/r) / ((z/r) - 1)² = z·r / (z-r)².
TEST_F(ZTransformTest, Modulation_2nTimesN) {
    auto a = parse_expr("2^n * n", ctx.arena());
    auto res = z_transform(a, n, z, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    // Z{r^n·n}(z) = Z{n}(z/r) = (z/r)/((z/r)-1)² = z·r/(z-r)².
    auto expected = parse_expr("2*z/(z-2)^2", ctx.arena());
    EXPECT_TRUE(same_after_simplify(res.value(), expected, ctx));
}

// Pattern non riconosciuto: log(n) → Unimplemented diagnostico.
TEST_F(ZTransformTest, UnknownPatternUnimplemented) {
    auto a = parse_expr("ln(n)", ctx.arena());
    auto res = z_transform(a, n, z, ctx);
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
}

}  // namespace
}  // namespace cas::calculus
