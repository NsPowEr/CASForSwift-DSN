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
