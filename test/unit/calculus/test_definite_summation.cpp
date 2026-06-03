// F5.7 sub-block 0 — Definite hypergeometric summation via Gosper.
//
// The Newton-Leibniz finite-calculus identity
//   Σ_{k=a}^{b} t(k) = S(b+1) − S(a),  where S(k+1) − S(k) = t(k),
// turns Gosper antidifferences into definite-sum closed forms.

#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

#include <string>

namespace cas::test {
namespace {

ExprPtr parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    EXPECT_TRUE(tokens.is_ok()) << tokens.error().message;
    Parser parser(tokens.value(), arena);
    auto res = parser.parse();
    EXPECT_TRUE(res.is_ok()) << res.error().message;
    return res.value();
}

bool simplifies_to_zero(ExprPtr a, ExprPtr b, symbolic::CASContext& ctx) {
    auto diff = ctx.arena().make<Binary>(BinaryOp::Sub, a, b);
    auto expanded = algebra::expand(diff, ctx);
    if (expanded.is_error()) return false;
    auto s = ctx.simplify(expanded.value());
    if (s.is_error()) return false;
    if (const auto* il = expr_cast<IntegerLit>(s.value())) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(s.value())) return rl->numerator.is_zero();
    return false;
}

}  // namespace

class DefiniteSummationTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol k{"k"};
    Symbol n{"n"};
};

// Σ_{k=1}^{n} 1 = n.
TEST_F(DefiniteSummationTest, ConstantOne) {
    auto term = parse_expr("1", ctx.arena());
    auto lo = parse_expr("1", ctx.arena());
    auto hi = parse_expr("n", ctx.arena());
    auto res = calculus::sum(term, k, lo, hi, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("n", ctx.arena());
    EXPECT_TRUE(simplifies_to_zero(res.value(), ctx.simplify(expected).value(), ctx));
}

// Σ_{k=1}^{n} k = n·(n+1)/2.
TEST_F(DefiniteSummationTest, ArithmeticSeriesFirstN) {
    auto term = parse_expr("k", ctx.arena());
    auto lo = parse_expr("1", ctx.arena());
    auto hi = parse_expr("n", ctx.arena());
    auto res = calculus::sum(term, k, lo, hi, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("n*(n+1)/2", ctx.arena());
    EXPECT_TRUE(simplifies_to_zero(res.value(), ctx.simplify(expected).value(), ctx));
}

// Σ_{k=1}^{n} k² = n·(n+1)·(2n+1)/6  (closed via Gosper antidifference of k²
// after the F5.7-GOSPER-K2-NORMALIZATION rescaling).
TEST_F(DefiniteSummationTest, SumOfSquares) {
    auto term = parse_expr("k^2", ctx.arena());
    auto lo = parse_expr("1", ctx.arena());
    auto hi = parse_expr("n", ctx.arena());
    auto res = calculus::sum(term, k, lo, hi, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("n*(n+1)*(2*n+1)/6", ctx.arena());
    EXPECT_TRUE(simplifies_to_zero(res.value(), ctx.simplify(expected).value(), ctx));
}

// Σ_{k=1}^{n} 1/(k·(k+1)) = n/(n+1)  (telescoping; Gosper-summable rational).
TEST_F(DefiniteSummationTest, TelescopingRational) {
    auto term = parse_expr("1/(k*(k+1))", ctx.arena());
    auto lo = parse_expr("1", ctx.arena());
    auto hi = parse_expr("n", ctx.arena());
    auto res = calculus::sum(term, k, lo, hi, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto expected = parse_expr("n/(n+1)", ctx.arena());
    // Bring both through `together` first to collapse the partial-fraction
    // representation Gosper produces.
    auto res_tog = algebra::together(res.value(), ctx);
    ASSERT_TRUE(res_tog.is_ok());
    auto exp_tog = algebra::together(ctx.simplify(expected).value(), ctx);
    ASSERT_TRUE(exp_tog.is_ok());
    EXPECT_TRUE(simplifies_to_zero(res_tog.value(), exp_tog.value(), ctx));
}

// ── F5.7 sub-block 1 — Abramov polygamma rational summation ─────────────────

namespace {
bool is_digamma_call(ExprPtr e) {
    const auto* fc = expr_cast<FuncCall>(e);
    return fc && fc->func_id == BuiltinOp::Digamma;
}
bool is_polygamma_call(ExprPtr e) {
    const auto* fc = expr_cast<FuncCall>(e);
    return fc && fc->func_id == BuiltinOp::Polygamma;
}
}  // namespace

// Σ_{k=1}^{n} 1/k = ψ(n+1) − ψ(1) = ψ(n+1) + γ.  The polygamma path closes
// this where Gosper returns nullopt (no rational antidifference exists).
TEST_F(DefiniteSummationTest, HarmonicSum_ViaDigamma) {
    auto term = parse_expr("1/k", ctx.arena());
    auto lo = parse_expr("1", ctx.arena());
    auto hi = parse_expr("n", ctx.arena());
    auto res = calculus::sum(term, k, lo, hi, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    // Result must mention digamma (the only closed form).  Searching the
    // expression tree for a Digamma FuncCall confirms the polygamma path
    // produced the antidifference.
    bool found = false;
    auto walk = [&](auto self, ExprPtr e) -> void {
        if (!e) return;
        if (is_digamma_call(e)) { found = true; return; }
        if (const auto* bin = expr_cast<Binary>(e)) {
            self(self, bin->left); self(self, bin->right); return;
        }
        if (const auto* un = expr_cast<Unary>(e)) { self(self, un->operand); return; }
        if (const auto* sum = expr_cast<Sum>(e)) {
            for (ExprPtr t : sum->terms) self(self, t);
            return;
        }
        if (const auto* prod = expr_cast<Product>(e)) {
            for (ExprPtr t : prod->factors) self(self, t);
            return;
        }
    };
    walk(walk, res.value());
    EXPECT_TRUE(found) << "expected digamma antidifference in harmonic sum";
}

// F5.7 sub-block 2 — Abramov-Full: multi-atom rational summand.
// Σ_{k=1}^{n} 1/(k·(k+2)).  Partial fractions:  (1/2)/k − (1/2)/(k+2).
// Closed form  Σ_{k=1}^{n} 1/(k·(k+2)) = (1/2)·(1 + 1/2 − 1/(n+1) − 1/(n+2)) =
// (3/4) − (1/(2(n+1))) − (1/(2(n+2))).  Numeric check at n = 3:
// 1/3 + 1/8 + 1/15 = (40 + 15 + 8)/120 = 63/120 = 21/40 = 0.525.
TEST_F(DefiniteSummationTest, AbramovMultiAtom_RationalDecomposition) {
    auto term = parse_expr("1/(k*(k+2))", ctx.arena());
    auto lo = parse_expr("1", ctx.arena());
    auto hi = parse_expr("n", ctx.arena());
    auto res = calculus::sum(term, k, lo, hi, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    // Substitute n = 3 and verify numerical match with the hand-computed
    // value 21/40.
    auto at_three = ctx.substitute(res.value(), n,
        ctx.arena().make<IntegerLit>(BigInt(3)));
    ASSERT_TRUE(at_three.is_ok());
    auto simp = ctx.simplify(at_three.value());
    ASSERT_TRUE(simp.is_ok());
    auto expected = parse_expr("21/40", ctx.arena());
    auto exp_simp = ctx.simplify(expected);
    ASSERT_TRUE(exp_simp.is_ok());
    EXPECT_TRUE(simplifies_to_zero(simp.value(), exp_simp.value(), ctx))
        << "Abramov closed form does not evaluate to 21/40 at n = 3\n"
        << "got: " << debug_print(simp.value()) << "\n"
        << "raw S: " << debug_print(res.value());
}

// Σ_{k=1}^{n} 1/k² closes via polygamma ψ⁽¹⁾ (trigamma).
TEST_F(DefiniteSummationTest, BaselSumDefinite_ViaPolygamma) {
    auto term = parse_expr("1/k^2", ctx.arena());
    auto lo = parse_expr("1", ctx.arena());
    auto hi = parse_expr("n", ctx.arena());
    auto res = calculus::sum(term, k, lo, hi, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    bool found = false;
    auto walk = [&](auto self, ExprPtr e) -> void {
        if (!e) return;
        if (is_polygamma_call(e)) { found = true; return; }
        if (const auto* bin = expr_cast<Binary>(e)) {
            self(self, bin->left); self(self, bin->right); return;
        }
        if (const auto* un = expr_cast<Unary>(e)) { self(self, un->operand); return; }
        if (const auto* sum = expr_cast<Sum>(e)) {
            for (ExprPtr t : sum->terms) self(self, t);
            return;
        }
        if (const auto* prod = expr_cast<Product>(e)) {
            for (ExprPtr t : prod->factors) self(self, t);
            return;
        }
    };
    walk(walk, res.value());
    EXPECT_TRUE(found) << "expected polygamma antidifference in 1/k² sum";
}

}  // namespace cas::test
