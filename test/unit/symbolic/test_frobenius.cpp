// L2-01: tests for solve_ode_frobenius_at_zero.
//
// Three scenarios:
//   1. Euler equation x^2 y'' - 6 y = 0     -> roots r = 3, r = -2.
//   2. 3 x^2 y'' - 4 x y' + 2 y = 0         -> roots r = 2, r = 1/3.
//   3. x^2 y'' + x y' - y = 0               -> roots r = 1, r = -1.

#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/ode.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

#include <memory>
#include <string>

namespace cas::test {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

// Structural search for x^r as a subexpression (looking for Binary Pow with x base
// and a target exponent expression matching `expected_exponent_src`).
bool contains_x_power(ExprPtr e, const std::string& var, ExprPtr expected_exponent,
                      symbolic::CASContext& ctx) {
    if (!e) return false;
    if (const auto* b = expr_cast<Binary>(e); b != nullptr && b->op == BinaryOp::Pow) {
        if (const auto* s = expr_cast<Symbol>(b->left); s != nullptr && s->name == var) {
            auto eq = symbolic::mathematically_equal(b->right, expected_exponent, ctx);
            if (eq.is_ok() && eq.value()) return true;
        }
    }
    // Also detect x itself when the exponent is 1.
    if (const auto* il = expr_cast<IntegerLit>(expected_exponent); il != nullptr && il->value == BigInt(1)) {
        if (const auto* s = expr_cast<Symbol>(e); s != nullptr && s->name == var) return true;
    }

    switch (e->kind) {
    case ExprKind::Unary:
        return contains_x_power(expr_ref<Unary>(e).operand, var, expected_exponent, ctx);
    case ExprKind::Binary: {
        const auto& b = expr_ref<Binary>(e);
        return contains_x_power(b.left, var, expected_exponent, ctx) ||
               contains_x_power(b.right, var, expected_exponent, ctx);
    }
    case ExprKind::Sum:
        for (auto t : expr_ref<Sum>(e).terms)
            if (contains_x_power(t, var, expected_exponent, ctx)) return true;
        return false;
    case ExprKind::Product:
        for (auto f : expr_ref<Product>(e).factors)
            if (contains_x_power(f, var, expected_exponent, ctx)) return true;
        return false;
    case ExprKind::FuncCall:
        for (auto a : expr_ref<FuncCall>(e).args)
            if (contains_x_power(a, var, expected_exponent, ctx)) return true;
        return false;
    default:
        return false;
    }
}

}  // namespace

class FrobeniusTest : public ::testing::Test {
protected:
    void SetUp() override { ctx = std::make_unique<symbolic::CASContext>(); }

    [[nodiscard]] ExprPtr E(const std::string& src) {
        auto r = parse_expr(src, ctx->arena());
        EXPECT_TRUE(r.is_ok()) << "parse: " << src
                               << " err=" << (r.is_error() ? r.error().message : std::string{});
        return r.is_ok() ? r.value() : ExprPtr{};
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

// 1. Euler:  x^2 y'' - 6 y = 0   ->   a_2 = x^2, a_1 = 0, a_0 = -6.
//    Indicial:  r(r-1) + 0 - 6 = r^2 - r - 6 = (r-3)(r+2),  roots r = 3, -2.
//    Recurrence: p_tilde = 0, q_tilde = -6 constants; all higher p_k,q_k = 0,
//    so c_n = 0 for n >= 1.  Solution = C1 * x^3 + C2 * x^(-2).
TEST_F(FrobeniusTest, EulerEquationProducesXCubedAndXMinusTwo) {
    ExprPtr a2 = E("x^2");
    ExprPtr a1 = E("0");
    ExprPtr a0 = E("-6");
    Symbol x("x");

    auto result = calculus::solve_ode_frobenius_at_zero(a2, a1, a0, x, 4U, *ctx);
    ASSERT_TRUE(result.is_ok()) << "err=" << (result.is_error() ? result.error().message : std::string{});

    ExprPtr y = result.value();
    ExprPtr exp_pos = ctx->arena().make<IntegerLit>(BigInt(3));
    ExprPtr exp_neg = ctx->arena().make<IntegerLit>(BigInt(-2));

    EXPECT_TRUE(contains_x_power(y, "x", exp_pos, *ctx))
        << "expected x^3 term in: " << debug_print(y);
    EXPECT_TRUE(contains_x_power(y, "x", exp_neg, *ctx))
        << "expected x^(-2) term in: " << debug_print(y);
}

// 2. 3 x^2 y'' - 4 x y' + 2 y = 0   ->   a_2 = 3 x^2,  a_1 = -4 x,  a_0 = 2.
//    p_tilde = x*p = x*(-4x)/(3x^2) = -4/3,  q_tilde = x^2*q = 2/3.  Both constant.
//    Indicial: r(r-1) - (4/3) r + 2/3 = r^2 - (7/3) r + 2/3 = 0.
//      Multiply by 3: 3 r^2 - 7 r + 2 = (3 r - 1)(r - 2),  roots r = 1/3, r = 2.
//    All p_k, q_k for k >= 1 are zero, so c_n = 0 for n >= 1.
//    Solution = C1 * x^2 + C2 * x^(1/3).
TEST_F(FrobeniusTest, RationalIndicialRoots) {
    ExprPtr a2 = E("3*x^2");
    ExprPtr a1 = E("-4*x");
    ExprPtr a0 = E("2");
    Symbol x("x");

    auto result = calculus::solve_ode_frobenius_at_zero(a2, a1, a0, x, 3U, *ctx);
    ASSERT_TRUE(result.is_ok()) << "err=" << (result.is_error() ? result.error().message : std::string{});

    ExprPtr y = result.value();
    ExprPtr exp_two = ctx->arena().make<IntegerLit>(BigInt(2));
    ExprPtr exp_third = E("1/3");

    EXPECT_TRUE(contains_x_power(y, "x", exp_two, *ctx))
        << "expected x^2 term in: " << debug_print(y);
    EXPECT_TRUE(contains_x_power(y, "x", exp_third, *ctx))
        << "expected x^(1/3) term in: " << debug_print(y);
}

// 3. Diagnostic / sanity:  x^2 y'' + x y' - y = 0.
//    p_tilde = x*p = x*(x)/(x^2) = 1,  q_tilde = x^2*q = x^2*(-1)/(x^2) = -1.
//    Indicial: r(r-1) + r - 1 = r^2 - 1 = 0,  roots r = 1, r = -1.
TEST_F(FrobeniusTest, IndicialPlusMinusOne) {
    ExprPtr a2 = E("x^2");
    ExprPtr a1 = E("x");
    ExprPtr a0 = E("-1");
    Symbol x("x");

    auto result = calculus::solve_ode_frobenius_at_zero(a2, a1, a0, x, 3U, *ctx);
    ASSERT_TRUE(result.is_ok()) << "err=" << (result.is_error() ? result.error().message : std::string{});

    ExprPtr y = result.value();
    ExprPtr exp_pos = ctx->arena().make<IntegerLit>(BigInt(1));
    ExprPtr exp_neg = ctx->arena().make<IntegerLit>(BigInt(-1));

    EXPECT_TRUE(contains_x_power(y, "x", exp_pos, *ctx))
        << "expected x^1 term in: " << debug_print(y);
    EXPECT_TRUE(contains_x_power(y, "x", exp_neg, *ctx))
        << "expected x^(-1) term in: " << debug_print(y);
}

// ── F5.3 / B2c — Logarithmic Frobenius branch (resonant integer gap) ────────
//
// Helper: detect a `Ln(x)` subexpression in the result, which is the
// hallmark of the second-solution log term  y_2 = c·ln(x)·y_1 + x^{r_2}·Σ b_n x^n.
bool contains_ln_x(ExprPtr e, const std::string& var) {
    if (!e) return false;
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        if (fc->name == "ln" || fc->func_id == BuiltinOp::Ln) {
            if (!fc->args.empty()) {
                if (const auto* s = expr_cast<Symbol>(fc->args[0]);
                    s && s->name == var) return true;
            }
        }
    }
    switch (e->kind) {
    case ExprKind::Unary:
        return contains_ln_x(expr_ref<Unary>(e).operand, var);
    case ExprKind::Binary: {
        const auto& b = expr_ref<Binary>(e);
        return contains_ln_x(b.left, var) || contains_ln_x(b.right, var);
    }
    case ExprKind::Sum:
        for (auto t : expr_ref<Sum>(e).terms) if (contains_ln_x(t, var)) return true;
        return false;
    case ExprKind::Product:
        for (auto f : expr_ref<Product>(e).factors) if (contains_ln_x(f, var)) return true;
        return false;
    case ExprKind::FuncCall:
        for (auto a : expr_ref<FuncCall>(e).args) if (contains_ln_x(a, var)) return true;
        return false;
    default:
        return false;
    }
}

// Bessel of order 1:  x² y'' + x y' + (x² − 1) y = 0.
// Indicial: r² − 1 = 0,  roots r = 1, r = −1,  gap N = 2.
// At r = −1 the recurrence hits I(1) = 0 with non-zero RHS (q_2 = 1 carries
// over from q_tilde = x² − 1).  The second solution must include a log term:
//     y_2 = c·ln(x)·y_1 + x^{−1}·Σ b_n x^n,    c = −1/2  (Coddington-Levinson §4.8).
TEST_F(FrobeniusTest, BesselOrder1_ResonantLogBranch) {
    ExprPtr a2 = E("x^2");
    ExprPtr a1 = E("x");
    ExprPtr a0 = E("x^2 - 1");
    Symbol x("x");

    auto result = calculus::solve_ode_frobenius_at_zero(a2, a1, a0, x, 4U, *ctx);
    ASSERT_TRUE(result.is_ok()) << "err=" << (result.is_error() ? result.error().message : std::string{});

    ExprPtr y = result.value();
    EXPECT_TRUE(contains_ln_x(y, "x"))
        << "expected ln(x) (log-branch) in: " << debug_print(y);

    ExprPtr exp_pos = ctx->arena().make<IntegerLit>(BigInt(1));
    ExprPtr exp_neg = ctx->arena().make<IntegerLit>(BigInt(-1));
    EXPECT_TRUE(contains_x_power(y, "x", exp_pos, *ctx))
        << "expected x^1 (y_1) in: " << debug_print(y);
    EXPECT_TRUE(contains_x_power(y, "x", exp_neg, *ctx))
        << "expected x^(-1) (y_2 leading) in: " << debug_print(y);
}

// ── A5 / F5.3 — Double indicial root (gap N = 0), logarithmic second solution ──
//
// Euler with a repeated indicial root:  x² y'' − x y' + y = 0.
//   p_tilde = x·(−x)/x² = −1  (p0 = −1),  q_tilde = x²·(1)/x² = 1  (q0 = 1).
//   Indicial: r² + (p0−1)r + q0 = r² − 2r + 1 = (r − 1)²,  double root r = 1.
//   All p_k, q_k (k ≥ 1) vanish ⇒ the power correction is zero and the general
//   solution is  C1·x + C2·x·ln(x).
TEST_F(FrobeniusTest, EulerDoubleRootProducesXLogX) {
    ExprPtr a2 = E("x^2");
    ExprPtr a1 = E("-x");
    ExprPtr a0 = E("1");
    Symbol x("x");

    auto result = calculus::solve_ode_frobenius_at_zero(a2, a1, a0, x, 4U, *ctx);
    ASSERT_TRUE(result.is_ok()) << "err=" << (result.is_error() ? result.error().message : std::string{});

    ExprPtr y = result.value();
    EXPECT_TRUE(contains_ln_x(y, "x"))
        << "expected x·ln(x) (double-root log term) in: " << debug_print(y);
    ExprPtr exp_one = ctx->arena().make<IntegerLit>(BigInt(1));
    EXPECT_TRUE(contains_x_power(y, "x", exp_one, *ctx))
        << "expected x^1 in: " << debug_print(y);
}

// Bessel of order 0:  x² y'' + x y' + x² y = 0.
//   p_tilde = 1 (p0 = 1),  q_tilde = x²  (q0 = 0, q_2 = 1).
//   Indicial: r² = 0,  double root r = 0.
//   y_1 = 1 − x²/4 + x⁴/64 − …  (J_0),
//   y_2 = J_0·ln(x) + x²/4 − 3x⁴/128 + …  (the logarithmic second solution).
TEST_F(FrobeniusTest, BesselOrder0_DoubleRootLogBranch) {
    ExprPtr a2 = E("x^2");
    ExprPtr a1 = E("x");
    ExprPtr a0 = E("x^2");
    Symbol x("x");

    auto result = calculus::solve_ode_frobenius_at_zero(a2, a1, a0, x, 4U, *ctx);
    ASSERT_TRUE(result.is_ok()) << "err=" << (result.is_error() ? result.error().message : std::string{});

    ExprPtr y = result.value();
    EXPECT_TRUE(contains_ln_x(y, "x"))
        << "expected ln(x) (double-root log term) in: " << debug_print(y);
    // y_1 carries the x² coefficient −1/4 (J_0 series); the x² power must appear.
    ExprPtr exp_two = ctx->arena().make<IntegerLit>(BigInt(2));
    EXPECT_TRUE(contains_x_power(y, "x", exp_two, *ctx))
        << "expected x^2 in: " << debug_print(y);
}

}  // namespace cas::test
