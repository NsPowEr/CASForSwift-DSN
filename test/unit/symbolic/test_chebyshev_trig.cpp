// F1.4b — Chebyshev trig generator tests.
//
// Covers six angle types:
//  1. cos(π/12) → constructible (Gauss: 12 = 4·3, Fermat primes {3})
//  2. cos(π/8)  → constructible (8 = 2³, pure power of 2)
//  3. cos(π/15) → constructible via angle subtraction (15 = 3·5, both Fermat)
//  4. cos(π/17) → constructible (Gauss heptadecagon, Fermat prime 17)
//  5. cos(2π/7) = cos(π·2/7) → non-constructible (7 not Fermat), emits RootOf
//  6. sin(π/30) → constructible (30 = 2·3·5)
//
// All tests use structural checks (not toString) as mandated by CLAUDE.md.
// RootOf test verifies the node kind, not a float approximation.

#include <gtest/gtest.h>
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::symbolic;

namespace {

// ── Helpers ──────────────────────────────────────────────────────────────────

[[nodiscard]] Result<ExprPtr> parse_and_simplify(const std::string& input, CASContext& ctx) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), ctx.arena());
    auto parsed = parser.parse();
    if (!parsed.is_ok()) return parsed;
    return ctx.simplify(parsed.value());
}

// Walk the AST to check whether any node is a RootOf.
[[nodiscard]] bool contains_rootof(ExprPtr e) {
    if (!e) return false;
    if (expr_kind(e) == ExprKind::RootOf) return true;
    bool found = false;
    visit_expr(e, [&](const auto& node) {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_base_of_v<ExprNode, T>) {
            (void)node; // handled by recursive walk below
        }
    });
    // Recursive manual walk for common node types.
    if (const auto* b = expr_cast<Binary>(e))
        return contains_rootof(b->left) || contains_rootof(b->right);
    if (const auto* u = expr_cast<Unary>(e))
        return contains_rootof(u->operand);
    if (const auto* s = expr_cast<Sum>(e)) {
        for (ExprPtr t : s->terms) if (contains_rootof(t)) return true;
    }
    if (const auto* p = expr_cast<Product>(e)) {
        for (ExprPtr f : p->factors) if (contains_rootof(f)) return true;
    }
    if (const auto* f = expr_cast<FuncCall>(e)) {
        for (ExprPtr a : f->args) if (contains_rootof(a)) return true;
    }
    if (expr_kind(e) == ExprKind::RootOf) return true;
    return found;
}

// Check that an expression does NOT contain RootOf (pure radical form).
[[nodiscard]] bool is_radical_form(ExprPtr e) {
    return !contains_rootof(e);
}

} // anonymous namespace

// ── Test 1: cos(π/12) — constructible (12 = 4·3) ────────────────────────────
// cos(π/12) = cos(15°) = (√6 + √2)/4
// Constructibility: 12 = 2² · 3; Fermat prime 3 ✓.
TEST(ChebyshevTrigTest, CosPiOver12_Constructible) {
    CASContext ctx;
    auto result = parse_and_simplify("cos(pi/12)", ctx);
    ASSERT_TRUE(result.is_ok()) << "Simplification failed for cos(pi/12)";
    EXPECT_TRUE(is_radical_form(result.value()))
        << "cos(pi/12) should simplify to a nested radical (constructible), not RootOf";
}

// ── Test 2: cos(π/8) — constructible (8 = 2³) ───────────────────────────────
// cos(π/8) = cos(22.5°) = √(2+√2)/2
// Constructibility: 8 = 2³, pure power of 2 ✓.
TEST(ChebyshevTrigTest, CosPiOver8_Constructible) {
    CASContext ctx;
    auto result = parse_and_simplify("cos(pi/8)", ctx);
    ASSERT_TRUE(result.is_ok()) << "Simplification failed for cos(pi/8)";
    EXPECT_TRUE(is_radical_form(result.value()))
        << "cos(pi/8) should be a nested radical, not RootOf";
    // Additional structural check: result must contain sqrt (not be trivial integer/rational).
    bool has_sqrt = false;
    std::function<void(ExprPtr)> find_sqrt = [&](ExprPtr e) {
        if (!e || has_sqrt) return;
        if (const auto* f = expr_cast<FuncCall>(e); f && f->func_id == BuiltinOp::Sqrt)
            has_sqrt = true;
        if (const auto* b = expr_cast<Binary>(e)) { find_sqrt(b->left); find_sqrt(b->right); }
        if (const auto* u = expr_cast<Unary>(e)) find_sqrt(u->operand);
        if (const auto* s = expr_cast<Sum>(e)) for (auto t : s->terms) find_sqrt(t);
        if (const auto* p = expr_cast<Product>(e)) for (auto t : p->factors) find_sqrt(t);
        if (const auto* f2 = expr_cast<FuncCall>(e)) for (auto a : f2->args) find_sqrt(a);
    };
    find_sqrt(result.value());
    EXPECT_TRUE(has_sqrt) << "cos(pi/8) should involve sqrt (nested radical)";
}

// ── Test 3: cos(π/15) — constructible via angle subtraction ─────────────────
// cos(π/15) = cos(12°): 15 = 3·5; both Fermat primes ✓.
// Result: constructible (angle subtraction from π/3 and π/5).
TEST(ChebyshevTrigTest, CosPiOver15_ConstructibleViaAngleCombination) {
    CASContext ctx;
    auto result = parse_and_simplify("cos(pi/15)", ctx);
    ASSERT_TRUE(result.is_ok()) << "Simplification failed for cos(pi/15)";
    EXPECT_TRUE(is_radical_form(result.value()))
        << "cos(pi/15) should be constructible (no RootOf)";
}

// ── Test 4: cos(π/17) — Gauss heptadecagon ────────────────────────────────────
// cos(π/17): 17 is the 3rd Fermat prime.
// With the F1.4c depth guard, try_angle_combination bails out at depth 3,
// preventing stack overflow. The engine falls through to the RootOf generator:
//   cos(π/17) → RootOf(Ψ_{34}(t), t, 0) / 2
// where Ψ_{34} is the minimal polynomial of 2cos(π/17) over Q (degree φ(34)/2=8).
//
// The full nested-radical expansion via Gauss period is a separate subproblem:
//   16cos(2π/17) = -1 + √17 + √(34-2√17) + 2√(17+3√17-√(34-2√17)-2√(34+2√17))
// (Gauss, Disquisitiones §VII art. 354, 1801).
// This closed form requires the Gauss period algorithm for Fermat primes,
// which is marked as Aperta permanente HPP-014c (see HARDCODE_LEDGER.md).
// Current canonical output: RootOf(Ψ_{34}, _tcc, 0) / 2 — structurally exact.
TEST(ChebyshevTrigTest, CosPiOver17_StackGuard_RootOf) {
    // F1.4c: depth guard prevents stack overflow.
    // Expected: ok result containing RootOf (no crash, no Unimplemented error).
    CASContext ctx;
    auto result = parse_and_simplify("cos(pi/17)", ctx);
    ASSERT_TRUE(result.is_ok())
        << "cos(pi/17) must not return an error after F1.4c depth guard; "
        << "got: " << (result.is_error() ? result.error().message : "(none)");
    EXPECT_TRUE(contains_rootof(result.value()))
        << "cos(pi/17): engine should fall through to RootOf(Psi_34, _tcc, 0)/2 "
        << "when try_angle_combination depth guard fires (Fermat prime, "
        << "closed-form Gauss period deferred to HPP-014c)";
}

// ── Test 5: cos(2π/7) = cos(π·2/7) — non-constructible ──────────────────────
// q = 7 is prime but NOT a Fermat prime (7 ≠ 2^(2^k)+1).
// The Galois group of Φ_{14}(x) has order φ(14)/2 = 3, which is odd and > 1
// → cos(2π/7) is NOT expressible in nested radicals over Q.
// The engine must emit RootOf(Ψ_{14}(t), t, 0)/2 — the canonical exact form.
//
// With the F1.4c depth guard, try_angle_combination returns nullptr at depth 3
// so cos(2π/7) = cos(π · 2/7) → Chebyshev T_2 applied to cos(π/7).
// cos(π/7) → RootOf(Ψ_{14}, _tcc, 0)/2; T_2 applied symbolically gives
// an expression in that RootOf — still contains_rootof = true.
// This is the mathematically correct representation: no nested radicals exist.
TEST(ChebyshevTrigTest, CosTwoPiOverSeven_StackGuard_RootOf) {
    // F1.4c: depth guard prevents stack overflow.
    // Expected: ok result containing RootOf (non-constructible, correct).
    CASContext ctx;
    auto result = parse_and_simplify("cos(2*pi/7)", ctx);
    ASSERT_TRUE(result.is_ok())
        << "cos(2*pi/7) must not error after F1.4c depth guard; "
        << (result.is_error() ? result.error().message : "");
    EXPECT_TRUE(contains_rootof(result.value()))
        << "cos(2pi/7) is non-constructible (7 not Fermat prime); "
        << "must emit RootOf(Psi_14, _tcc, 0)-based expression";
}

// ── Test 6: sin(π/30) — constructible (30 = 2·3·5) ──────────────────────────
// sin(π/30) = sin(6°): 30 = 2·3·5; Fermat primes {3,5} ✓.
// Result: constructible radical.
TEST(ChebyshevTrigTest, SinPiOver30_Constructible) {
    CASContext ctx;
    auto result = parse_and_simplify("sin(pi/30)", ctx);
    ASSERT_TRUE(result.is_ok()) << "Simplification failed for sin(pi/30)";
    EXPECT_TRUE(is_radical_form(result.value()))
        << "sin(pi/30) should be constructible (no RootOf); 30 = 2·3·5";
}

// ── Bonus: cos(π/5) is constructible ─────────────────────────────────────────
// Regression guard: cos(π/5) = (1+√5)/4 is a classic Gauss base angle.
TEST(ChebyshevTrigTest, CosPiOver5_Constructible_RegressionGuard) {
    CASContext ctx;
    auto result = parse_and_simplify("cos(pi/5)", ctx);
    ASSERT_TRUE(result.is_ok()) << "Simplification failed for cos(pi/5)";
    EXPECT_TRUE(is_radical_form(result.value()))
        << "cos(pi/5) must remain a nested radical (regression guard)";
}

// ── Bonus: cos(π/7) — same try_angle_combination issue ────────────────────────
// q = 7: Φ_{14} degree 6, Ψ_{14} degree 3. Same stack issue as 2π/7.
// With F1.4c depth guard, no crash. Engine emits RootOf(Ψ_{14}, _tcc, 0)/2.
TEST(ChebyshevTrigTest, CosPiOverSeven_StackGuard_RootOf) {
    // F1.4c: depth guard prevents stack overflow for cos(pi/7).
    // Expected: ok result with RootOf (7 is not a Fermat prime).
    CASContext ctx;
    auto result = parse_and_simplify("cos(pi/7)", ctx);
    ASSERT_TRUE(result.is_ok())
        << "cos(pi/7) must not error after F1.4c depth guard; "
        << (result.is_error() ? result.error().message : "");
    EXPECT_TRUE(contains_rootof(result.value()))
        << "cos(pi/7): q=7 is not a Fermat prime, must emit RootOf(Psi_14)";
}
