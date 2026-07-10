// F5.7-B6BIS-QUADRATIC-M-GT-1 — polygamma antidifference for
// (A1·k + A0)/Q(k)^m with Q an irreducible quadratic and m > 1.
//
// The m>1 branch of try_quadratic_atom_antidiff was implemented but had no
// test pinning its mathematics (the ledger still listed it as OPEN). The
// load-bearing check is the partial-fraction identity behind it:
//
//   (A1·k + A0)/((k−α)^m (k−β)^m) = Σ_{j=1..m} [ C_j/(k−α)^j + D_j/(k−β)^j ]
//
// which is an identity in Q(α, β, k) — so it can be verified EXACTLY with
// rational α, β (no RootOf machinery involved), where simplify() can close
// the difference to a literal zero. If the Hermite-style coefficients in
// quadratic_pf_coeffs were wrong, these tests would fail with a nonzero
// rational function.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include "cas/symbolic.hpp"
#include "../../../src/calculus/summation_abramov_internal.hpp"

#include <string>
#include <vector>

using namespace cas;
using namespace cas::calculus::abramov_detail;

namespace {

class AbramovQuadraticMultiplicityTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol k{"k"};

    [[nodiscard]] ExprPtr integer(long long v) {
        return ctx.arena().make<IntegerLit>(BigInt(v));
    }

    // Σ_j [C_j/(k−α)^j + D_j/(k−β)^j]  −  (A1·k+A0)/((k−α)(k−β))^m  → must be 0.
    void expect_identity(long long a0, long long a1,
                         long long alpha_v, long long beta_v,
                         unsigned int m) {
        AstArena& arena = ctx.arena();
        ExprPtr A0 = integer(a0);
        ExprPtr A1 = integer(a1);
        ExprPtr alpha = integer(alpha_v);
        ExprPtr beta = integer(beta_v);
        ExprPtr k_e = arena.make<Symbol>(k);

        auto coeffs = quadratic_pf_coeffs(A0, A1, alpha, beta, m, arena);
        ASSERT_EQ(coeffs.size(), m);

        std::vector<ExprPtr> terms;
        for (unsigned int j = 1U; j <= m; ++j) {
            ExprPtr ka = arena.make<Binary>(BinaryOp::Sub, k_e, alpha);
            ExprPtr kb = arena.make<Binary>(BinaryOp::Sub, k_e, beta);
            ExprPtr ka_j = arena.make<Binary>(BinaryOp::Pow, ka,
                integer(static_cast<long long>(j)));
            ExprPtr kb_j = arena.make<Binary>(BinaryOp::Pow, kb,
                integer(static_cast<long long>(j)));
            terms.push_back(arena.make<Binary>(
                BinaryOp::Div, coeffs[j - 1U].first, ka_j));
            terms.push_back(arena.make<Binary>(
                BinaryOp::Div, coeffs[j - 1U].second, kb_j));
        }
        ExprPtr lhs = arena.make<Sum>(std::move(terms));

        ExprPtr num = arena.make<Binary>(BinaryOp::Add,
            arena.make<Binary>(BinaryOp::Mul, A1, k_e), A0);
        ExprPtr ka = arena.make<Binary>(BinaryOp::Sub, k_e, alpha);
        ExprPtr kb = arena.make<Binary>(BinaryOp::Sub, k_e, beta);
        ExprPtr den = arena.make<Binary>(BinaryOp::Pow,
            arena.make<Binary>(BinaryOp::Mul, ka, kb),
            integer(static_cast<long long>(m)));
        ExprPtr rhs = arena.make<Binary>(BinaryOp::Div, num, den);

        ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, lhs, rhs);
        auto together = algebra::together(delta, ctx);
        ASSERT_TRUE(together.is_ok()) << together.error().message;
        auto delta_s = ctx.simplify(together.value());
        ASSERT_TRUE(delta_s.is_ok()) << delta_s.error().message;
        const auto* il = expr_cast<IntegerLit>(delta_s.value());
        EXPECT_TRUE(il != nullptr && il->value.is_zero())
            << "partial-fraction identity violated for m=" << m
            << " A0=" << a0 << " A1=" << a1
            << " — residual: " << debug_print(delta_s.value());
    }
};

TEST_F(AbramovQuadraticMultiplicityTest, IdentityM1_Baseline) {
    expect_identity(/*a0=*/1, /*a1=*/0, /*alpha_v=*/2, /*beta_v=*/3, /*m=*/1U);
}

TEST_F(AbramovQuadraticMultiplicityTest, IdentityM2_ConstantNumerator) {
    expect_identity(1, 0, 2, 3, 2U);
}

TEST_F(AbramovQuadraticMultiplicityTest, IdentityM2_LinearNumerator) {
    expect_identity(5, 1, 2, 3, 2U);
}

TEST_F(AbramovQuadraticMultiplicityTest, IdentityM3_ConstantNumerator) {
    expect_identity(1, 0, -1, 4, 3U);
}

TEST_F(AbramovQuadraticMultiplicityTest, IdentityM3_LinearNumerator) {
    expect_identity(-7, 3, -1, 4, 3U);
}

// ── End-to-end: m=2 over a Q-irreducible quadratic fires the RootOf path ────

TEST_F(AbramovQuadraticMultiplicityTest, EndToEnd_M2_ProducesPolygammaAntidiff) {
    // 1/(k²+k+1)²: disc = −3 < 0 (irreducible over Q), multiplicity m=2.
    // The antidifference must contain ψ' (Polygamma order 1) at RootOf-shifted
    // arguments — the m>1 branch of try_quadratic_atom_antidiff.
    ExprPtr k_e = ctx.arena().make<Symbol>(k);
    ExprPtr quad = ctx.arena().make<Sum>(std::vector<ExprPtr>{
        ctx.arena().make<Binary>(BinaryOp::Pow, k_e, integer(2)),
        k_e,
        integer(1)});
    ExprPtr term = ctx.arena().make<Binary>(BinaryOp::Div, integer(1),
        ctx.arena().make<Binary>(BinaryOp::Pow, quad, integer(2)));

    auto anti = try_quadratic_atom_antidiff(term, k, ctx);
    ASSERT_TRUE(anti.has_value())
        << "m=2 quadratic atom must be handled (F5.7-B6BIS)";

    bool has_polygamma = false;
    auto walk = [&](auto self, ExprPtr e) -> void {
        if (!e || has_polygamma) return;
        if (const auto* fc = expr_cast<FuncCall>(e)) {
            if (fc->func_id == BuiltinOp::Polygamma) { has_polygamma = true; return; }
            for (ExprPtr a : fc->args) self(self, a);
            return;
        }
        if (const auto* bin = expr_cast<Binary>(e)) {
            self(self, bin->left); self(self, bin->right); return;
        }
        if (const auto* un = expr_cast<Unary>(e)) { self(self, un->operand); return; }
        if (const auto* s = expr_cast<Sum>(e)) {
            for (ExprPtr t2 : s->terms) self(self, t2);
            return;
        }
        if (const auto* p = expr_cast<Product>(e)) {
            for (ExprPtr t2 : p->factors) self(self, t2);
            return;
        }
    };
    walk(walk, anti.value());
    EXPECT_TRUE(has_polygamma)
        << "expected Polygamma(1, k - RootOf(...)) in the m=2 antidifference";
}

}  // namespace
