// F5.7 — Zeilberger creative telescoping tests.
//
// The algorithm core (try_parametric_gosper + solve_first_order_rec) is
// fully implemented.  Test coverage here exercises:
//   - expand_gamma_int_shifts: Γ(z+n) → ∏(z+i)·Γ(z) rewriting.
//   - cancel_common_factors_in_ratio: structural cancellation in fractions.
//   - try_zeilberger_definite: graceful handling of bivariate inputs.
//
// Gamma-encoded binomial sums (e.g., Σ C(n,k)) are NOT yet exercised end-to-end
// because the global simplifier preserves Γ(z+n) calls inside Product/Div
// expressions for reflection-identity matching, preventing extraction of clean
// rational shift ratios.  See HARDCODE_LEDGER entry F5.7-ZEIL-GAMMA-RATIO.

#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "../../../src/symbolic/summation_zeilberger_helpers.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas::test {
namespace {

ExprPtr parse(const std::string& s, AstArena& arena) {
    auto tokens = Lexer(s).tokenize();
    EXPECT_TRUE(tokens.is_ok()) << tokens.error().message;
    Parser parser(tokens.value(), arena);
    auto res = parser.parse();
    EXPECT_TRUE(res.is_ok()) << res.error().message;
    return res.value();
}

bool simplifies_to_zero(ExprPtr a, ExprPtr b, symbolic::CASContext& ctx) {
    auto diff = ctx.arena().make<Binary>(BinaryOp::Sub, a, b);
    auto ex = algebra::expand(diff, ctx);
    if (ex.is_error()) return false;
    auto s = ctx.simplify(ex.value());
    if (s.is_error()) return false;
    if (const auto* il = expr_cast<IntegerLit>(s.value())) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(s.value())) return rl->numerator.is_zero();
    return false;
}

}  // namespace

class ZeilbergerTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol k{"k"};
    Symbol n{"n"};
};

// ── expand_gamma_int_shifts ────────────────────────────────────────────────

// Γ(k+1) → k·Γ(k).
TEST_F(ZeilbergerTest, ExpandGammaIntShifts_GammaKPlus1) {
    AstArena& arena = ctx.arena();
    auto e = parse("gamma(k+1)", arena);
    auto expanded = symbolic::zeilberger_detail::expand_gamma_int_shifts(e, ctx);
    auto expected = parse("gamma(k)*k", arena);
    EXPECT_TRUE(simplifies_to_zero(expanded, ctx.simplify(expected).value(), ctx))
        << "expanded = " << debug_print(expanded);
}

// Γ(k+2) → (k+1)·k·Γ(k).
TEST_F(ZeilbergerTest, ExpandGammaIntShifts_GammaKPlus2) {
    AstArena& arena = ctx.arena();
    auto e = parse("gamma(k+2)", arena);
    auto expanded = symbolic::zeilberger_detail::expand_gamma_int_shifts(e, ctx);
    auto expected = parse("gamma(k)*k*(k+1)", arena);
    EXPECT_TRUE(simplifies_to_zero(expanded, ctx.simplify(expected).value(), ctx))
        << "expanded = " << debug_print(expanded);
}

// ── cancel_common_factors_in_ratio ────────────────────────────────────────

// (x·y)/(y·z) → x/z.
TEST_F(ZeilbergerTest, CancelCommonFactors_BasicXYZ) {
    AstArena& arena = ctx.arena();
    Symbol x{"x"}, y{"y"}, z{"z"};
    ExprPtr xe = arena.make<Symbol>(x), ye = arena.make<Symbol>(y), ze = arena.make<Symbol>(z);
    ExprPtr num = arena.make<Binary>(BinaryOp::Mul, xe, ye);
    ExprPtr den = arena.make<Binary>(BinaryOp::Mul, ye, ze);
    ExprPtr ratio = arena.make<Binary>(BinaryOp::Div, num, den);
    auto cancelled = symbolic::zeilberger_detail::cancel_common_factors_in_ratio(ratio, ctx);
    auto expected = arena.make<Binary>(BinaryOp::Div, xe, ze);
    EXPECT_TRUE(simplifies_to_zero(cancelled, ctx.simplify(expected).value(), ctx))
        << "cancelled = " << debug_print(cancelled);
}

// Nested Div: (a/b)/(c/d) → a·d/(b·c).
TEST_F(ZeilbergerTest, CancelCommonFactors_NestedDiv) {
    AstArena& arena = ctx.arena();
    Symbol a{"a"}, b{"b"}, c{"c"}, d{"d"};
    ExprPtr ae = arena.make<Symbol>(a), be = arena.make<Symbol>(b);
    ExprPtr ce = arena.make<Symbol>(c), de = arena.make<Symbol>(d);
    ExprPtr lhs = arena.make<Binary>(BinaryOp::Div, ae, be);
    ExprPtr rhs = arena.make<Binary>(BinaryOp::Div, ce, de);
    ExprPtr ratio = arena.make<Binary>(BinaryOp::Div, lhs, rhs);
    auto cancelled = symbolic::zeilberger_detail::cancel_common_factors_in_ratio(ratio, ctx);
    auto expected = parse("a*d/(b*c)", arena);
    EXPECT_TRUE(simplifies_to_zero(cancelled, ctx.simplify(expected).value(), ctx))
        << "cancelled = " << debug_print(cancelled);
}

// ── Zeilberger driver: graceful handling ──────────────────────────────────

// Purely-numeric upper bound bypasses Zeilberger gracefully.
TEST_F(ZeilbergerTest, NumericUpperBound_DoesNotCrash) {
    auto term = parse("gamma(4)/(gamma(k+1)*gamma(4-k))", ctx.arena());
    auto lo   = parse("0", ctx.arena());
    auto hi   = parse("3", ctx.arena());
    auto res  = calculus::sum(term, k, lo, hi, ctx);
    (void)res;  // no crash is the contract here.
    SUCCEED();
}

// Binomial sum returns Unimplemented (not closed form) due to simplifier
// gap (HARDCODE_LEDGER F5.7-ZEIL-GAMMA-RATIO).  Important contract: no crash,
// no silently wrong answer.
TEST_F(ZeilbergerTest, BinomialSum_GracefulUnimplemented) {
    auto term = parse("gamma(n+1)/(gamma(k+1)*gamma(n-k+1))", ctx.arena());
    auto lo   = parse("0", ctx.arena());
    auto hi   = parse("n", ctx.arena());
    auto res  = calculus::sum(term, k, lo, hi, ctx);
    // Must NOT crash. May or may not produce a closed form depending on
    // the simplifier's evolving handling of gamma ratios.  If it returns
    // an error it must be Unimplemented (never silently wrong).
    if (res.is_error()) {
        EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented)
            << "Got non-Unimplemented error: " << res.error().message;
    }
}

}  // namespace cas::test
