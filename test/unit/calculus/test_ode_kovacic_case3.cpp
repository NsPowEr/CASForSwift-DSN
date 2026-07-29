// test_ode_kovacic_case3.cpp — Kovacic Case 3 acceptance tests.
//
// Reference: Kovacic 1986, §5, examples pp. 23-26.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Kovacic_Case3_SL2C.md
// Ledger: HC-KV-03.
//
// Test matrix (Kovacic_Case3_SL2C.md §"Test corpus minimo"):
//   C3-2  Paper Example 2 (p. 25):  r = −(5x + 27) / (36 (x − 1)²)
//                                    Case 3 with n = 4, octahedral group.
//   C3-3  Airy:  r = x  →  Case 3 cannot hold (ord(r,∞) = −1 < 2).
//   C3-4  Pole-of-order-3 input  →  necessary condition §2 violated.
//   C3-5  Algebraic certificate: returned minimal polynomial M(ω) is a non-
//         trivial polynomial of degree ∈ {4, 6, 12} in the fresh ω-symbol.

#include <gtest/gtest.h>
#include "cas/symbolic.hpp"
#include "cas/calculus.hpp"
#include "cas/ast_debug.hpp"
#include "calculus/ode_kovacic_internal.hpp"
#include <iostream>
#include <string>

using namespace cas;
using namespace cas::calculus;

namespace {

class OdeKovacicCase3Test : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

// ─── C3-3: Airy y'' = x·y — Case 3 cannot hold ───────────────────────────────
// ord(r=x, ∞) = −1 < 2.  Step 1 ∞ branch yields no family.
TEST_F(OdeKovacicCase3Test, Airy_Case3Unimplemented) {
    Symbol x("x");
    AstArena& a = ctx.arena();
    ExprPtr r = a.make<Symbol>("x");

    auto res = kovacic_impl::case3_omega(r, x, ctx);
    ASSERT_TRUE(res.is_error())
        << "Airy r=x violates ord(r,∞) ≥ 2; Case 3 must return Unimplemented.";
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
}

// ─── C3-4: pole of order 3 — Case 3 §2 necessary condition violated ─────────
// r = 1/x³.  ord at 0 = 3 > 2.  case3_omega rejects immediately.
TEST_F(OdeKovacicCase3Test, Order3Pole_Unimplemented) {
    Symbol x("x");
    AstArena& a = ctx.arena();
    ExprPtr r = a.make<Binary>(BinaryOp::Div,
        a.make<IntegerLit>(BigInt(1)),
        a.make<Binary>(BinaryOp::Pow, a.make<Symbol>("x"),
            a.make<IntegerLit>(BigInt(3))));

    auto res = kovacic_impl::case3_omega(r, x, ctx);
    ASSERT_TRUE(res.is_error())
        << "Pole of order 3 violates Case 3 §2; must return Unimplemented.";
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
    EXPECT_NE(res.error().message.find("order > 2"), std::string::npos)
        << "Diagnostic must mention pole order > 2: " << res.error().message;
}

// Build  r = −3/(16x²) − 2/(9(x−1)²) + 3/(16 x (x−1))   (Kovacic 1986 Ex.1, p.23)
static ExprPtr build_paper_example1_r(symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();
    ExprPtr xs = a.make<Symbol>("x");
    ExprPtr one = a.make<IntegerLit>(BigInt(1));
    ExprPtr xm1 = a.make<Binary>(BinaryOp::Sub, xs, one);
    ExprPtr xsq = a.make<Binary>(BinaryOp::Pow, xs, a.make<IntegerLit>(BigInt(2)));
    ExprPtr xm1sq = a.make<Binary>(BinaryOp::Pow, xm1, a.make<IntegerLit>(BigInt(2)));
    ExprPtr t1 = a.make<Unary>(UnaryOp::Neg,
        a.make<Binary>(BinaryOp::Div, a.make<IntegerLit>(BigInt(3)),
            a.make<Binary>(BinaryOp::Mul,
                a.make<IntegerLit>(BigInt(16)), xsq)));
    ExprPtr t2 = a.make<Unary>(UnaryOp::Neg,
        a.make<Binary>(BinaryOp::Div, a.make<IntegerLit>(BigInt(2)),
            a.make<Binary>(BinaryOp::Mul,
                a.make<IntegerLit>(BigInt(9)), xm1sq)));
    ExprPtr t3 = a.make<Binary>(BinaryOp::Div,
        a.make<IntegerLit>(BigInt(3)),
        a.make<Binary>(BinaryOp::Mul,
            a.make<IntegerLit>(BigInt(16)),
            a.make<Binary>(BinaryOp::Mul, xs, xm1)));
    return a.make<Binary>(BinaryOp::Add,
        a.make<Binary>(BinaryOp::Add, t1, t2), t3);
}

// ─── C3-2: Paper Example 1 (p. 23) — n=12 (icosahedral A₅) MUST succeed ──────
// Kovacic 1986 proves Case 3 with n=12 holds for this DE.  Since the
// HC-KV-06 closure (PolyExpr recurrence + general linear-factor pole
// extraction, 2026-07-08) the engine finds the degree-12 minimal polynomial
// of ω; a regression back to Unimplemented is a hard failure.
TEST_F(OdeKovacicCase3Test, Paper_Example1_p23_N12_Icosahedral) {
    Symbol x("x");
    ExprPtr r = build_paper_example1_r(ctx);

    auto res = kovacic_impl::case3_omega(r, x, ctx);
    ASSERT_TRUE(res.is_ok())
        << "Paper Example 1 (p. 23) is a proven n=12 Case 3 instance; "
           "case3_omega must succeed.  Got: " << res.error().message;
    auto* root = expr_cast<RootOf>(res.value().plus);
    ASSERT_NE(root, nullptr);
    EXPECT_FALSE(root->variable.name.empty());
}

// ─── C3-5: structural soundness of the returned minimal polynomial ───────────
// The wrapped minimal polynomial must reference the fresh ω-variable.
// Probe input: paper Example 1 (guaranteed-success since HC-KV-06 closure).
TEST_F(OdeKovacicCase3Test, MinPoly_Structural_Soundness) {
    Symbol x("x");
    ExprPtr r = build_paper_example1_r(ctx);

    auto res = kovacic_impl::case3_omega(r, x, ctx);
    ASSERT_TRUE(res.is_ok())
        << "case3_omega must succeed on paper Example 1: "
        << res.error().message;
    auto* root = expr_cast<RootOf>(res.value().plus);
    ASSERT_NE(root, nullptr);
    const std::string poly_str = debug_print(root->polynomial);
    EXPECT_NE(poly_str.find(root->variable.name), std::string::npos)
        << "Minimal polynomial does not reference ω fresh variable '"
        << root->variable.name << "': " << poly_str;
}

} // anonymous namespace
