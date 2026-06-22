// F3.3 — Cyclotomic polynomial via Möbius inversion: power-gain tests.
//
// Pre-fix: `compute_cyclotomic` used recursive identity with thread-
// unsafe static cache and recursion depth ≤ log₂(n).
// Post-fix: direct Möbius construction Φ_n(x) = ∏ (x^d − 1)^μ(n/d),
// no recursion, no cache. Verifies known closed forms and stress on
// large n where pre-fix would have hit the kCyclotomicMaxN cap.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "../../../src/algebra/polynomial_internal.hpp"

using namespace cas;
using namespace cas::algebra;
using cas::Symbol;

namespace {

[[nodiscard]] IntPoly poly(std::initializer_list<long long> coeffs) {
    std::vector<BigInt> bs;
    bs.reserve(coeffs.size());
    for (auto c : coeffs) bs.emplace_back(c);
    return IntPoly(std::move(bs));
}

}  // namespace

// Φ_1(x) = x - 1
TEST(CyclotomicMobiusTest, Phi1) {
    auto v = is_cyclotomic(poly({-1, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 1);
}

// Φ_2(x) = x + 1
TEST(CyclotomicMobiusTest, Phi2) {
    auto v = is_cyclotomic(poly({1, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 2);
}

// Φ_3(x) = x² + x + 1
TEST(CyclotomicMobiusTest, Phi3) {
    auto v = is_cyclotomic(poly({1, 1, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 3);
}

// Φ_4(x) = x² + 1
TEST(CyclotomicMobiusTest, Phi4) {
    auto v = is_cyclotomic(poly({1, 0, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 4);
}

// Φ_6(x) = x² - x + 1 (canonical "tricky" case: μ(6/d) involves
// μ(6) = +1, μ(2)/μ(3) = -1, μ(1) = +1).
TEST(CyclotomicMobiusTest, Phi6) {
    auto v = is_cyclotomic(poly({1, -1, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 6);
}

// Φ_8(x) = x^4 + 1
TEST(CyclotomicMobiusTest, Phi8) {
    auto v = is_cyclotomic(poly({1, 0, 0, 0, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 8);
}

// Φ_12(x) = x^4 - x² + 1
TEST(CyclotomicMobiusTest, Phi12) {
    auto v = is_cyclotomic(poly({1, 0, -1, 0, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 12);
}

// Φ_15(x) = x^8 − x^7 + x^5 − x^4 + x³ − x + 1
TEST(CyclotomicMobiusTest, Phi15) {
    auto v = is_cyclotomic(poly({1, -1, 0, 1, -1, 1, 0, -1, 1}));
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 15);
}

// Power-gain: Φ_97 (prime, degree 96). Pre-fix this would have spent
// ~97 recursive calls each allocating polynomial divisors; post-fix
// is a single Möbius pass. The polynomial is x^96 + x^95 + ... + 1.
TEST(CyclotomicMobiusTest, Phi97PrimeStillIdentified) {
    // Φ_p(x) = (x^p - 1)/(x - 1) = sum_{i=0}^{p-1} x^i  (for p prime)
    std::vector<BigInt> c97(97, BigInt(1));
    IntPoly phi_97(std::move(c97));
    auto v = is_cyclotomic(phi_97);
    ASSERT_TRUE(v.has_value());
    EXPECT_EQ(*v, 97);
}

// Power-gain: non-cyclotomic input must return nullopt without
// running forever or hitting the recursion cap.
TEST(CyclotomicMobiusTest, NonCyclotomicReturnsNullopt) {
    // (x+2) is not a cyclotomic polynomial.
    auto v = is_cyclotomic(poly({2, 1}));
    EXPECT_FALSE(v.has_value());
}

// ─── OOM cap guard (A5-LARGECYCLO) ───────────────────────────────────────────
//
// is_cyclotomic with an explicit small max_n that is less than the true
// cyclotomic order must return nullopt (cannot verify) rather than wrongly
// deciding the polynomial is or is not cyclotomic for those n > max_n.
//
// We use Φ_12(x) = x^4 - x^2 + 1 with max_n = 6: since n=12 > max_n, the
// scan cannot reach n=12 and must return nullopt.
TEST(CyclotomicMobiusTest, OomCapSmallMaxN_Phi12NotFound) {
    // Φ_12(x) = x^4 - x^2 + 1
    auto p = poly({1, 0, -1, 0, 1});
    // max_n = 6: scan only checks n = 1..6, none of which produce Φ_12.
    auto v = is_cyclotomic(p, /*max_n=*/6);
    EXPECT_FALSE(v.has_value())
        << "OomCapSmallMaxN: max_n=6 < 12, Phi_12 must not be found";
}

// is_cyclotomic correctly returns nullopt for degree-0 or empty input.
TEST(CyclotomicMobiusTest, DegenerateEmptyPoly) {
    IntPoly empty;
    auto v = is_cyclotomic(empty);
    EXPECT_FALSE(v.has_value());
}

// ─── A5-LARGECYCLO closed (T-022): deg > 724 detection ──────────────────────
//
// The former implementation materialised order-n intermediates and capped at
// n ≤ 2^20, so is_cyclotomic returned nullopt proactively for deg > 724. The
// squarefree-reduction compute_cyclotomic keeps intermediates at degree O(φ(n))
// and inverse-totient candidate enumeration replaces the O(d²) scan, so large
// degrees are now both computable and detectable.
//
// Φ_1009 (1009 prime) has degree 1008 > 724; Φ_2018 = Φ_2·1009 also degree 1008.
TEST(CyclotomicMobiusTest, RoundTripLargeDegreeAbove724) {
    for (int n : {105, 1155, 1009, 2018}) {
        IntPoly phi = compute_cyclotomic(n);
        ASSERT_FALSE(phi.empty()) << "compute_cyclotomic failed for n=" << n;
        auto v = is_cyclotomic(phi);
        ASSERT_TRUE(v.has_value())
            << "is_cyclotomic returned nullopt for Phi_" << n
            << " (deg " << phi.degree() << ")";
        EXPECT_EQ(*v, n) << "wrong order recovered for n=" << n;
    }
}

// Independent correctness check (not self-referential): Φ_105 = Φ_{3·5·7} is the
// smallest cyclotomic polynomial with a coefficient outside {-1,0,1}; it has two
// coefficients equal to -2. Validates the multi-prime squarefree fold.
TEST(CyclotomicMobiusTest, Phi105HasMinusTwoCoefficient) {
    IntPoly phi = compute_cyclotomic(105);
    ASSERT_EQ(phi.degree(), static_cast<std::size_t>(48));
    int minus_two_count = 0;
    for (std::size_t i = 0; i < phi.size(); ++i) {
        if (phi[i] == BigInt(-2)) ++minus_two_count;
    }
    EXPECT_EQ(minus_two_count, 2)
        << "Phi_105 must contain exactly two -2 coefficients";
}

// ─── cyclotomic_roots (production path via solve_polynomial wiring) ──────────
//
// cyclotomic_roots is called from solve_polynomial.cpp:493 when is_cyclotomic
// succeeds. These tests exercise it DIRECTLY (via polynomial_internal.hpp) and
// assert structural properties of the returned roots:
//
//   n=1: single root {1}, n=2: single root {-1}, n=k (prime): phi(k)=k-1 roots.
//
// We do NOT test exact symbolic form (implementation returns exp(2πik/n)),
// only COUNT of roots (phi(n) = Euler totient) and that the result vector
// is non-empty for n >= 1.
//
// The test also verifies that the integration path through solve_polynomial
// routes correctly: solve(x^4+1, x) should find that x^4+1 = Φ_8(x) and
// call cyclotomic_roots(8, ...), returning phi(8)=4 roots.

// Direct cyclotomic_roots: n=1 → 1 root {1}.
TEST(CyclotomicRootsTest, N1_SingleRootOne) {
    cas::symbolic::CASContext ctx;
    Symbol x("x");
    auto roots = cyclotomic_roots(1, x, ctx.arena());
    ASSERT_EQ(roots.size(), 1U) << "Phi_1 has phi(1)=1 primitive 1st root of unity";
    // Root is IntegerLit(1).
    ASSERT_TRUE(roots[0] != nullptr);
    const auto* lit = expr_cast<cas::IntegerLit>(roots[0]);
    ASSERT_NE(lit, nullptr) << "n=1: root must be integer literal 1";
    EXPECT_EQ(lit->value, cas::BigInt(1));
}

// Direct cyclotomic_roots: n=2 → 1 root {-1}.
TEST(CyclotomicRootsTest, N2_SingleRootMinusOne) {
    cas::symbolic::CASContext ctx;
    Symbol x("x");
    auto roots = cyclotomic_roots(2, x, ctx.arena());
    ASSERT_EQ(roots.size(), 1U) << "Phi_2 has phi(2)=1 primitive 2nd root of unity";
    ASSERT_TRUE(roots[0] != nullptr);
    const auto* lit = expr_cast<cas::IntegerLit>(roots[0]);
    ASSERT_NE(lit, nullptr) << "n=2: root must be integer literal -1";
    EXPECT_EQ(lit->value, cas::BigInt(-1));
}

// Direct cyclotomic_roots: n=3 → phi(3)=2 primitive cube roots.
TEST(CyclotomicRootsTest, N3_TwoPrimitiveCubeRoots) {
    cas::symbolic::CASContext ctx;
    Symbol x("x");
    auto roots = cyclotomic_roots(3, x, ctx.arena());
    // phi(3) = 2: k=1 and k=2 both have gcd(k,3)=1.
    EXPECT_EQ(roots.size(), 2U) << "Phi_3 has phi(3)=2 primitive 3rd roots";
    for (const auto& r : roots) {
        EXPECT_NE(r, nullptr) << "Each root must be non-null";
    }
}

// Direct cyclotomic_roots: n=4 → phi(4)=2 primitive 4th roots.
TEST(CyclotomicRootsTest, N4_TwoPrimitive4thRoots) {
    cas::symbolic::CASContext ctx;
    Symbol x("x");
    auto roots = cyclotomic_roots(4, x, ctx.arena());
    // phi(4) = 2: k=1 and k=3.
    EXPECT_EQ(roots.size(), 2U) << "Phi_4 has phi(4)=2 primitive 4th roots";
}

// Direct cyclotomic_roots: n=6 → phi(6)=2 primitive 6th roots.
TEST(CyclotomicRootsTest, N6_TwoPrimitive6thRoots) {
    cas::symbolic::CASContext ctx;
    Symbol x("x");
    auto roots = cyclotomic_roots(6, x, ctx.arena());
    // phi(6) = 2: k=1 and k=5 (gcd(1,6)=1, gcd(5,6)=1).
    EXPECT_EQ(roots.size(), 2U) << "Phi_6 has phi(6)=2 primitive 6th roots";
}

// Direct cyclotomic_roots: n=5 (prime) → phi(5)=4 primitive 5th roots.
TEST(CyclotomicRootsTest, N5_FourPrimitive5thRoots) {
    cas::symbolic::CASContext ctx;
    Symbol x("x");
    auto roots = cyclotomic_roots(5, x, ctx.arena());
    // phi(5) = 4: k=1,2,3,4 all coprime to 5.
    EXPECT_EQ(roots.size(), 4U) << "Phi_5 has phi(5)=4 primitive 5th roots";
    for (const auto& r : roots) {
        EXPECT_NE(r, nullptr);
        // Each root is exp(2*pi*i*k/5), a FuncCall node.
        const auto* fc = expr_cast<cas::FuncCall>(r);
        EXPECT_NE(fc, nullptr) << "n>2: each root must be a FuncCall(exp, ...)";
    }
}

// Wiring test: solve_polynomial on x^2 + x + 1 = Phi_3(x) routes through
// cyclotomic_roots and returns exactly phi(3)=2 roots.
TEST(CyclotomicRootsTest, SolvePolynomialWiring_Phi3) {
    cas::symbolic::CASContext ctx;
    Symbol x("x");

    // Parse x^2 + x + 1 = Phi_3(x).
    auto tokens = cas::Lexer("x^2 + x + 1").tokenize();
    ASSERT_TRUE(tokens.is_ok());
    cas::Parser parser(tokens.value(), ctx.arena());
    auto expr = parser.parse();
    ASSERT_TRUE(expr.is_ok());

    auto roots = cas::algebra::solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok()) << roots.error().message;
    // Must produce phi(3) = 2 roots via cyclotomic_roots path.
    EXPECT_EQ(roots.value().size(), 2U)
        << "solve(x^2+x+1, x) via cyclotomic_roots must return 2 roots";
    for (const auto& r : roots.value()) {
        EXPECT_NE(r, nullptr);
    }
}

// Wiring test: solve_polynomial on x^4 + 1 = Phi_8(x) → phi(8)=4 roots.
TEST(CyclotomicRootsTest, SolvePolynomialWiring_Phi8) {
    cas::symbolic::CASContext ctx;
    Symbol x("x");

    auto tokens = cas::Lexer("x^4 + 1").tokenize();
    ASSERT_TRUE(tokens.is_ok());
    cas::Parser parser(tokens.value(), ctx.arena());
    auto expr = parser.parse();
    ASSERT_TRUE(expr.is_ok());

    auto roots = cas::algebra::solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok()) << roots.error().message;
    // phi(8) = 4.
    EXPECT_EQ(roots.value().size(), 4U)
        << "solve(x^4+1, x) via cyclotomic_roots must return 4 roots";
}
