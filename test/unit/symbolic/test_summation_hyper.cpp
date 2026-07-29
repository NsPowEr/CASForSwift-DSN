// F5.7 — isolated tests for Petkovšek's Hyper (rational z) term-ratio solver.
#include <gtest/gtest.h>

#include "cas/symbolic.hpp"
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "../../../src/symbolic/summation_hyper.hpp"

using namespace cas;
using namespace cas::symbolic;

namespace {

// Numerically check that ρ is a hypergeometric solution of the recurrence
//   Σ_i p_i(n)·y(n+i)=0  by substituting  y(n+i) = ∏_{j<i} ρ(n+j)  at n=v.
// Returns true iff the recurrence collapses to 0 at every probe point.
bool ratio_solves_recurrence(const std::vector<ExprPtr>& p, ExprPtr rho,
                             const Symbol& n, CASContext& ctx) {
    AstArena& a = ctx.arena();
    for (long long v : {4LL, 5LL, 7LL}) {
        ExprPtr sum = a.make<IntegerLit>(BigInt(0));
        ExprPtr prod = a.make<IntegerLit>(BigInt(1));
        for (std::size_t i = 0; i < p.size(); ++i) {
            if (i > 0) {
                ExprPtr pt = a.make<IntegerLit>(BigInt(v + static_cast<long long>(i) - 1));
                auto rs = ctx.substitute(rho, n, pt);
                if (rs.is_error()) return false;
                prod = a.make<Binary>(BinaryOp::Mul, prod, rs.value());
            }
            ExprPtr vt = a.make<IntegerLit>(BigInt(v));
            auto ps = ctx.substitute(p[i], n, vt);
            if (ps.is_error()) return false;
            sum = a.make<Binary>(BinaryOp::Add, sum,
                a.make<Binary>(BinaryOp::Mul, ps.value(), prod));
        }
        auto s = ctx.simplify(sum);
        if (s.is_error()) return false;
        const auto* il = expr_cast<IntegerLit>(s.value());
        const auto* rl = expr_cast<RationalLit>(s.value());
        bool zero = (il && il->value.is_zero()) || (rl && rl->numerator.is_zero());
        if (!zero) return false;
    }
    return true;
}

// Build n + c as a polynomial coefficient.
ExprPtr lin(AstArena& a, long long coeff_n, long long c, const Symbol& n) {
    ExprPtr term = a.make<Binary>(BinaryOp::Mul, a.make<IntegerLit>(BigInt(coeff_n)), a.make<Symbol>(n));
    return a.make<Binary>(BinaryOp::Add, term, a.make<IntegerLit>(BigInt(c)));
}

}  // namespace

// y(n)=2^n:  y(n+1)−2y(n)=0  ⇒  p=[-2, 1],  ρ=2.
TEST(HyperTest, Geometric2n_FirstOrder) {
    CASContext ctx;
    Symbol n("n");
    AstArena& a = ctx.arena();
    std::vector<ExprPtr> p{a.make<IntegerLit>(BigInt(-2)), a.make<IntegerLit>(BigInt(1))};

    auto res = hyper_term_ratio(p, n, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    ASSERT_TRUE(res.value().ratio.has_value()) << "expected a hypergeometric ratio";
    EXPECT_TRUE(ratio_solves_recurrence(p, *res.value().ratio, n, ctx))
        << "ratio = " << debug_print(*res.value().ratio);
}

// y(n)=n!:  y(n+1)−(n+1)y(n)=0  ⇒  p=[-(n+1), 1],  ρ=n+1.
TEST(HyperTest, Factorial_FirstOrder) {
    CASContext ctx;
    Symbol n("n");
    AstArena& a = ctx.arena();
    ExprPtr neg_np1 = a.make<Unary>(UnaryOp::Neg, lin(a, 1, 1, n));
    std::vector<ExprPtr> p{neg_np1, a.make<IntegerLit>(BigInt(1))};

    auto res = hyper_term_ratio(p, n, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    ASSERT_TRUE(res.value().ratio.has_value());
    EXPECT_TRUE(ratio_solves_recurrence(p, *res.value().ratio, n, ctx))
        << "ratio = " << debug_print(*res.value().ratio);
}

// (E−3)(E−2):  y(n+2)−5y(n+1)+6y(n)=0  ⇒  p=[6,-5,1].  Solutions 2^n, 3^n.
TEST(HyperTest, Order2_TwoGeometricSolutions) {
    CASContext ctx;
    Symbol n("n");
    AstArena& a = ctx.arena();
    std::vector<ExprPtr> p{
        a.make<IntegerLit>(BigInt(6)),
        a.make<IntegerLit>(BigInt(-5)),
        a.make<IntegerLit>(BigInt(1))};

    auto res = hyper_term_ratio(p, n, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    ASSERT_TRUE(res.value().ratio.has_value()) << "order-2 geometric recurrence has rational solutions";
    EXPECT_TRUE(ratio_solves_recurrence(p, *res.value().ratio, n, ctx))
        << "ratio = " << debug_print(*res.value().ratio);
}

// (E−(n+1))²:  y(n+2)−(2n+3)y(n+1)+(n+1)²y(n)=0.  Solution n! (ratio n+1).
TEST(HyperTest, Order2_PolynomialCoefficients) {
    CASContext ctx;
    Symbol n("n");
    AstArena& a = ctx.arena();
    ExprPtr np1 = lin(a, 1, 1, n);
    ExprPtr p0 = a.make<Binary>(BinaryOp::Pow, np1, a.make<IntegerLit>(BigInt(2)));      // (n+1)²
    ExprPtr p1 = a.make<Unary>(UnaryOp::Neg, lin(a, 2, 3, n));                            // -(2n+3)
    std::vector<ExprPtr> p{p0, p1, a.make<IntegerLit>(BigInt(1))};

    auto res = hyper_term_ratio(p, n, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    ASSERT_TRUE(res.value().ratio.has_value());
    EXPECT_TRUE(ratio_solves_recurrence(p, *res.value().ratio, n, ctx))
        << "ratio = " << debug_print(*res.value().ratio);
}

// Fibonacci:  y(n+2)−y(n+1)−y(n)=0.  Char poly z²−z−1 has irrational roots ⇒
// no hypergeometric solution with rational z; needs_algebraic must be flagged.
TEST(HyperTest, Fibonacci_NoRationalHypergeometric) {
    CASContext ctx;
    Symbol n("n");
    AstArena& a = ctx.arena();
    std::vector<ExprPtr> p{
        a.make<IntegerLit>(BigInt(-1)),
        a.make<IntegerLit>(BigInt(-1)),
        a.make<IntegerLit>(BigInt(1))};

    auto res = hyper_term_ratio(p, n, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_FALSE(res.value().ratio.has_value());
    EXPECT_TRUE(res.value().needs_algebraic) << "irrational z should set needs_algebraic";
}

// Evaluate a closed form S(n) at an integer point as an exact rational.
namespace {
bool closed_form_equals_at(ExprPtr S, const Symbol& n, long long v,
                           long long num, long long den, CASContext& ctx) {
    AstArena& a = ctx.arena();
    auto sub = ctx.substitute(S, n, a.make<IntegerLit>(BigInt(v)));
    if (sub.is_error()) return false;
    auto s = ctx.simplify(sub.value());
    if (s.is_error()) return false;
    const auto* il = expr_cast<IntegerLit>(s.value());
    const auto* rl = expr_cast<RationalLit>(s.value());
    if (il) return il->value == BigInt(num) && den == 1;
    if (rl) return rl->numerator == BigInt(num) && rl->denominator == BigInt(den);
    return false;
}
}  // namespace

// Closed form of (E−3)(E−2): p=[6,-5,1], S(0)=2, S(1)=5 ⇒ S(n)=2^n+3^n.
TEST(HyperTest, ClosedForm_TwoGeometric_Combination) {
    CASContext ctx;
    Symbol n("n");
    AstArena& a = ctx.arena();
    std::vector<ExprPtr> p{
        a.make<IntegerLit>(BigInt(6)),
        a.make<IntegerLit>(BigInt(-5)),
        a.make<IntegerLit>(BigInt(1))};
    std::vector<ExprPtr> init{a.make<IntegerLit>(BigInt(2)), a.make<IntegerLit>(BigInt(5))};

    auto res = solve_recurrence_closed_form(p, init, n, 0, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    ASSERT_TRUE(res.value().has_value()) << "2^n+3^n is a hypergeometric combination";
    ExprPtr S = *res.value();
    // 2^n+3^n at n=2 → 13, n=3 → 35, n=4 → 97.
    EXPECT_TRUE(closed_form_equals_at(S, n, 2, 13, 1, ctx)) << debug_print(S);
    EXPECT_TRUE(closed_form_equals_at(S, n, 3, 35, 1, ctx)) << debug_print(S);
    EXPECT_TRUE(closed_form_equals_at(S, n, 4, 97, 1, ctx)) << debug_print(S);
}

// Closed form of (E−(n+1))²: p=[(n+1)²,-(2n+3),1], S(0)=1, S(1)=1 ⇒ S(n)=n!.
TEST(HyperTest, ClosedForm_Factorial_FromOrder2) {
    CASContext ctx;
    Symbol n("n");
    AstArena& a = ctx.arena();
    ExprPtr np1 = a.make<Binary>(BinaryOp::Add, a.make<Symbol>(n), a.make<IntegerLit>(BigInt(1)));
    ExprPtr p0 = a.make<Binary>(BinaryOp::Pow, np1, a.make<IntegerLit>(BigInt(2)));
    ExprPtr p1 = a.make<Unary>(UnaryOp::Neg, a.make<Binary>(BinaryOp::Add,
        a.make<Binary>(BinaryOp::Mul, a.make<IntegerLit>(BigInt(2)), a.make<Symbol>(n)),
        a.make<IntegerLit>(BigInt(3))));
    std::vector<ExprPtr> p{p0, p1, a.make<IntegerLit>(BigInt(1))};
    std::vector<ExprPtr> init{a.make<IntegerLit>(BigInt(1)), a.make<IntegerLit>(BigInt(1))};

    auto res = solve_recurrence_closed_form(p, init, n, 0, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    ASSERT_TRUE(res.value().has_value());
    ExprPtr S = *res.value();
    // n! at n=3 → 6, n=4 → 24, n=5 → 120.
    EXPECT_TRUE(closed_form_equals_at(S, n, 3, 6, 1, ctx)) << debug_print(S);
    EXPECT_TRUE(closed_form_equals_at(S, n, 4, 24, 1, ctx)) << debug_print(S);
    EXPECT_TRUE(closed_form_equals_at(S, n, 5, 120, 1, ctx)) << debug_print(S);
}

// End-to-end wrapper: Σ_{k=0}^{n} (2^k + 3^k) satisfies (E−1)(E−2)(E−3),
// p=[-6,11,-6,1].  Closed form = −3/2 + 2·2^n + (3/2)·3^n, cross-verified
// against the directly-computed sum.  Exercises sum_closed_form_from_recurrence,
// the 3-solution fit, and the direct-summation soundness gate.
TEST(HyperTest, SumClosedForm_FromRecurrence_ThreeSolutions) {
    CASContext ctx;
    Symbol n("n");
    Symbol k("k");
    AstArena& a = ctx.arena();
    // F(n,k) = 2^k + 3^k
    ExprPtr F = a.make<Binary>(BinaryOp::Add,
        a.make<Binary>(BinaryOp::Pow, a.make<IntegerLit>(BigInt(2)), a.make<Symbol>(k)),
        a.make<Binary>(BinaryOp::Pow, a.make<IntegerLit>(BigInt(3)), a.make<Symbol>(k)));
    std::vector<ExprPtr> p{
        a.make<IntegerLit>(BigInt(-6)), a.make<IntegerLit>(BigInt(11)),
        a.make<IntegerLit>(BigInt(-6)), a.make<IntegerLit>(BigInt(1))};
    ExprPtr lower = a.make<IntegerLit>(BigInt(0));

    auto res = sum_closed_form_from_recurrence(p, F, n, k, lower, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    ASSERT_TRUE(res.value().has_value()) << "Σ(2^k+3^k) has a hypergeometric closed form";
    ExprPtr S = *res.value();
    // Direct sums: S(3)=55, S(4)=152, S(5)=427.
    EXPECT_TRUE(closed_form_equals_at(S, n, 3, 55, 1, ctx)) << debug_print(S);
    EXPECT_TRUE(closed_form_equals_at(S, n, 4, 152, 1, ctx)) << debug_print(S);
    EXPECT_TRUE(closed_form_equals_at(S, n, 5, 427, 1, ctx)) << debug_print(S);
}

// Fibonacci has no rational hypergeometric closed form ⇒ ok(nullopt).
TEST(HyperTest, ClosedForm_Fibonacci_None) {
    CASContext ctx;
    Symbol n("n");
    AstArena& a = ctx.arena();
    std::vector<ExprPtr> p{
        a.make<IntegerLit>(BigInt(-1)),
        a.make<IntegerLit>(BigInt(-1)),
        a.make<IntegerLit>(BigInt(1))};
    std::vector<ExprPtr> init{a.make<IntegerLit>(BigInt(0)), a.make<IntegerLit>(BigInt(1))};

    auto res = solve_recurrence_closed_form(p, init, n, 0, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_FALSE(res.value().has_value()) << "Fibonacci is not ℚ-hypergeometric";
}
