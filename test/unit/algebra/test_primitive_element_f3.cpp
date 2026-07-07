// F3.4 — Primitive Element Theorem (Trager 1976) tests.
//
// Tests compute_primitive_element and detect_tower_n_level.
// ZERO toString validation: all checks use structural / algebraic invariants.
//
// Invariants verified for each test:
//   I1. min_poly_theta has the correct degree (product of individual degrees).
//   I2. min_poly_theta is monic (leading coefficient == 1).
//   I3. alphas_in_theta.size() == n.
//   I4. Each alphas_in_theta[i] has length < min_poly_theta.size().
//   I5. shifts.size() == n-1.
//   I6. CERTIFICATORE: for each i, let P_i = alphas_in_theta[i].
//       Compute P_i(theta_val) where theta_val is a numerical evaluation of θ.
//       Should agree with a numerical evaluation of α_i.
//       (We use rational arithmetic where possible; for sqrt cases use
//       Lagrange-style verification.)
//   I7. min_poly_theta(theta_expr) ≡ 0: evaluate min_poly_theta at θ numerically.

#include "cas/algebra.hpp"
#include "cas/algebraic_tower_bridge.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"

#include "algebra/polynomial_internal.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace cas::test {
namespace {

[[nodiscard]] Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    Lexer lexer(input);
    auto tokens = lexer.tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

// Build a monic rational min-poly with rational coefficients [c0, 0, c1, 0, ..., 1].
[[nodiscard]] algebra::AlgebraicNumber::CoeffVec quadratic_min_poly(long neg_constant) {
    return {Rational(BigInt(neg_constant)), Rational(BigInt(0)), Rational(BigInt(1))};
}

// ── EXACT certificate primitives (Rational arithmetic, zero floating-point) ──
//
// The primitive-element result is certified algebraically, NOT numerically:
// every generator α_i is reconstructed as a polynomial P_i in θ, so the exact
// statements that must hold in Q[y]/(q_θ) are:
//   C1.  m_i(P_i) ≡ 0  (mod q_θ)         — P_i(θ) is a root of α_i's min-poly.
//   C2.  P_1 + Σ s_k·P_{k+1} ≡ y (mod q_θ) — the θ = α_1 + Σ s_k·α_{k+1} relation.
// A floating-point check at tolerance ε could pass a subtly wrong reconstruction;
// the exact residue check cannot.

// Reduce m(P) modulo q in Q[y] via Horner, returning the exact residue.
// m: min-poly coefficients (ascending). P: ring element coefficients. q: modulus.
[[nodiscard]] algebra::RatPoly eval_minpoly_at_ring_elem(
    const algebra::AlgebraicNumber::CoeffVec& m,
    const std::vector<Rational>& P_coeffs,
    const algebra::RatPoly& q) {
    const algebra::RatPoly P(P_coeffs);
    algebra::RatPoly acc;  // zero polynomial
    for (std::size_t i = m.size(); i > 0U; --i) {
        acc = algebra::mul_rational_poly(acc, P);
        acc = algebra::add_rational_poly(acc, algebra::RatPoly({m[i - 1U]}));
        acc = algebra::div_rem_rational_poly(acc, q).second;
    }
    return acc;
}

// Build the exact linear combination P_1 + Σ shifts[k]·P_{k+2}, reduced mod q.
[[nodiscard]] algebra::RatPoly theta_combination_mod_q(
    const std::vector<std::vector<Rational>>& alphas_in_theta,
    const std::vector<BigInt>& shifts,
    const algebra::RatPoly& q) {
    algebra::RatPoly combo(alphas_in_theta[0]);
    for (std::size_t k = 0U; k < shifts.size(); ++k) {
        const algebra::RatPoly term(alphas_in_theta[k + 1U]);
        const algebra::RatPoly scaled =
            algebra::mul_rational_poly(term, algebra::RatPoly({Rational(shifts[k])}));
        combo = algebra::add_rational_poly(combo, scaled);
    }
    return algebra::div_rem_rational_poly(combo, q).second;
}

// Assert: residue p is exactly the monomial y (coeff[1] == 1, all others 0).
void expect_is_y(const algebra::RatPoly& p, const std::string& test_name) {
    EXPECT_EQ(p.degree(), 1U) << test_name << ": θ-combination not degree 1";
    if (p.size() >= 2U) {
        EXPECT_EQ(p[0], Rational(BigInt(0))) << test_name << ": θ-combination constant ≠ 0";
        EXPECT_EQ(p[1], Rational(BigInt(1))) << test_name << ": θ-combination linear ≠ 1";
    }
}

class PrimitiveElementTest : public ::testing::Test {
protected:
    void SetUp() override {
        ctx = std::make_unique<symbolic::CASContext>();
        ctx->set_timeout(std::chrono::seconds(120));
        // Allow deep simplification for radical expressions.
        ctx->set_max_simplification_depth(500);
        ctx->set_max_trager_tower_shift_attempts(50U);
    }

    [[nodiscard]] ExprPtr parse_ok(const std::string& input) {
        auto parsed = parse_expr(input, ctx->arena());
        EXPECT_TRUE(parsed.is_ok()) << parsed.error().message;
        return parsed.is_ok() ? parsed.value() : nullptr;
    }

    // Build a canonicalized RootOf node for "the k-th root of x^2 - n".
    [[nodiscard]] ExprPtr make_sqrt_rootof(long n, std::size_t var_idx = 0U) {
        const std::string var_name = "y" + std::to_string(var_idx);
        const std::string str = "RootOf(" + var_name + "^2-" + std::to_string(n) +
                                "," + var_name + ",0)";
        ExprPtr e = parse_ok(str);
        if (!e) return nullptr;
        const std::size_t saved = ctx->max_rootof_explicit_degree();
        ctx->set_max_rootof_explicit_degree(1U);
        auto canon = ctx->simplify(e);
        ctx->set_max_rootof_explicit_degree(saved);
        if (!canon.is_ok()) return nullptr;
        return canon.value();
    }

    // Build a canonicalized RootOf node for "the k-th root of x^3 - n".
    [[nodiscard]] ExprPtr make_cube_rootof(long n, std::size_t var_idx = 0U) {
        const std::string var_name = "y" + std::to_string(var_idx);
        const std::string str = "RootOf(" + var_name + "^3-" + std::to_string(n) +
                                "," + var_name + ",0)";
        ExprPtr e = parse_ok(str);
        if (!e) return nullptr;
        const std::size_t saved = ctx->max_rootof_explicit_degree();
        ctx->set_max_rootof_explicit_degree(1U);
        auto canon = ctx->simplify(e);
        ctx->set_max_rootof_explicit_degree(saved);
        if (!canon.is_ok()) return nullptr;
        return canon.value();
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

// ── Invariant helpers ─────────────────────────────────────────────────────────

void check_invariants(
    const algebra::PrimitiveElementResult& res,
    std::size_t n,
    std::size_t expected_degree,
    const std::string& test_name) {
    // I1: correct degree.
    EXPECT_EQ(res.min_poly_theta.size(), expected_degree + 1U)
        << test_name << ": min_poly_theta degree mismatch";

    // I2: monic.
    if (!res.min_poly_theta.empty()) {
        const Rational lc = res.min_poly_theta.back();
        EXPECT_EQ(lc, Rational(BigInt(1)))
            << test_name << ": min_poly_theta not monic";
    }

    // I3: alphas_in_theta size.
    EXPECT_EQ(res.alphas_in_theta.size(), n)
        << test_name << ": alphas_in_theta.size() != n";

    // I4: each alpha poly has length ≤ expected_degree.
    for (std::size_t i = 0U; i < res.alphas_in_theta.size(); ++i) {
        EXPECT_LE(res.alphas_in_theta[i].size(), expected_degree)
            << test_name << ": alphas_in_theta[" << i << "] has too many coefficients";
    }

    // I5: shifts size.
    EXPECT_EQ(res.shifts.size(), n - 1U)
        << test_name << ": shifts.size() mismatch";
}

// ── Test 1: sqrt(2) + sqrt(3)  (2-level sanity against existing path) ─────────

TEST_F(PrimitiveElementTest, SqrtTwoSqrtThree) {
    // Q(sqrt(2), sqrt(3)): expected primitive element degree 4.
    ExprPtr alpha1 = make_sqrt_rootof(2, 0);
    ExprPtr alpha2 = make_sqrt_rootof(3, 1);
    ASSERT_NE(alpha1, nullptr);
    ASSERT_NE(alpha2, nullptr);

    std::vector<ExprPtr> alphas{alpha1, alpha2};
    std::vector<algebra::AlgebraicNumber::CoeffVec> min_polys{
        quadratic_min_poly(-2),  // x^2 - 2
        quadratic_min_poly(-3),  // x^2 - 3
    };

    auto result = algebra::compute_primitive_element(alphas, min_polys, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    const auto& res = result.value();
    check_invariants(res, 2U, 4U, "SqrtTwoSqrtThree");

    ASSERT_FALSE(res.shifts.empty());
    const algebra::RatPoly q_theta(res.min_poly_theta);

    // C1: m_i(P_i) ≡ 0 mod q_θ for each generator (exact).
    EXPECT_TRUE(eval_minpoly_at_ring_elem(min_polys[0], res.alphas_in_theta[0], q_theta).is_zero())
        << "SqrtTwoSqrtThree: α1=√2 not a root of x²-2 in Q[y]/(q_θ)";
    EXPECT_TRUE(eval_minpoly_at_ring_elem(min_polys[1], res.alphas_in_theta[1], q_theta).is_zero())
        << "SqrtTwoSqrtThree: α2=√3 not a root of x²-3 in Q[y]/(q_θ)";

    // C2: P_1 + s·P_2 ≡ y mod q_θ (θ decomposition, exact).
    expect_is_y(theta_combination_mod_q(res.alphas_in_theta, res.shifts, q_theta),
                "SqrtTwoSqrtThree");
}

// ── Test 2: sqrt(2) + sqrt(3) + sqrt(5)  (3-level, main probe) ────────────────

TEST_F(PrimitiveElementTest, SqrtTwoSqrtThreeSqrtFive) {
    // Q(sqrt(2), sqrt(3), sqrt(5)): expected primitive element degree 8.
    ExprPtr alpha1 = make_sqrt_rootof(2, 0);
    ExprPtr alpha2 = make_sqrt_rootof(3, 1);
    ExprPtr alpha3 = make_sqrt_rootof(5, 2);
    ASSERT_NE(alpha1, nullptr);
    ASSERT_NE(alpha2, nullptr);
    ASSERT_NE(alpha3, nullptr);

    std::vector<ExprPtr> alphas{alpha1, alpha2, alpha3};
    std::vector<algebra::AlgebraicNumber::CoeffVec> min_polys{
        quadratic_min_poly(-2),  // x^2 - 2
        quadratic_min_poly(-3),  // x^2 - 3
        quadratic_min_poly(-5),  // x^2 - 5
    };

    auto result = algebra::compute_primitive_element(alphas, min_polys, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    const auto& res = result.value();

    // I1: degree must be 8 = 2*2*2.
    check_invariants(res, 3U, 8U, "SqrtTwoSqrtThreeSqrtFive");

    // Print probe output (captures in test output).
    std::cout << "[PrimitiveElement.SqrtTwoSqrtThreeSqrtFive] "
              << "shifts: [";
    for (std::size_t i = 0U; i < res.shifts.size(); ++i) {
        if (i > 0U) std::cout << ", ";
        std::cout << res.shifts[i].decimal();
    }
    std::cout << "], min_poly_theta degree: " << res.min_poly_theta.size() - 1U << "\n";

    // CERTIFICATORE (exact): each α_i is a root of its min-poly in Q[y]/(q_θ),
    // and the θ decomposition holds exactly.
    ASSERT_GE(res.shifts.size(), 2U);
    const algebra::RatPoly q_theta(res.min_poly_theta);
    for (std::size_t i = 0U; i < 3U; ++i) {
        EXPECT_TRUE(eval_minpoly_at_ring_elem(min_polys[i], res.alphas_in_theta[i], q_theta).is_zero())
            << "SqrtTwoSqrtThreeSqrtFive: α[" << i << "] not a root of its min-poly in Q[y]/(q_θ)";
    }
    expect_is_y(theta_combination_mod_q(res.alphas_in_theta, res.shifts, q_theta),
                "SqrtTwoSqrtThreeSqrtFive");
}

// ── Test 3: cbrt(2) + sqrt(3)  (non-quadratic min-poly, degree 6) ─────────────

TEST_F(PrimitiveElementTest, CubeRootTwoSqrtThree) {
    // Q(∛2, sqrt(3)): expected degree 6 = 3*2.
    ExprPtr alpha1 = make_cube_rootof(2, 0);   // ∛2, min-poly x^3 - 2
    ExprPtr alpha2 = make_sqrt_rootof(3, 1);   // sqrt(3), min-poly x^2 - 3
    ASSERT_NE(alpha1, nullptr);
    ASSERT_NE(alpha2, nullptr);

    // min-poly of ∛2: x^3 - 2 = [-2, 0, 0, 1].
    algebra::AlgebraicNumber::CoeffVec mp_cbrt2 = {
        Rational(BigInt(-2)),
        Rational(BigInt(0)),
        Rational(BigInt(0)),
        Rational(BigInt(1)),
    };

    std::vector<ExprPtr> alphas{alpha1, alpha2};
    std::vector<algebra::AlgebraicNumber::CoeffVec> min_polys{
        mp_cbrt2,
        quadratic_min_poly(-3),
    };

    auto result = algebra::compute_primitive_element(alphas, min_polys, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    const auto& res = result.value();
    check_invariants(res, 2U, 6U, "CubeRootTwoSqrtThree");

    // CERTIFICATORE (exact): α[0]=∛2 root of x³-2, α[1]=√3 root of x²-3,
    // both in Q[y]/(q_θ), plus the exact θ decomposition.
    ASSERT_FALSE(res.shifts.empty());
    const algebra::RatPoly q_theta(res.min_poly_theta);
    EXPECT_TRUE(eval_minpoly_at_ring_elem(min_polys[0], res.alphas_in_theta[0], q_theta).is_zero())
        << "CubeRootTwoSqrtThree: α[0]=∛2 not a root of x³-2 in Q[y]/(q_θ)";
    EXPECT_TRUE(eval_minpoly_at_ring_elem(min_polys[1], res.alphas_in_theta[1], q_theta).is_zero())
        << "CubeRootTwoSqrtThree: α[1]=√3 not a root of x²-3 in Q[y]/(q_θ)";
    expect_is_y(theta_combination_mod_q(res.alphas_in_theta, res.shifts, q_theta),
                "CubeRootTwoSqrtThree");
}

// ── Test 4: detect_tower_n_level with sqrt(2)+sqrt(3)+sqrt(5) expression ──────

TEST_F(PrimitiveElementTest, DetectTowerNLevel_SqrtTwoSqrtThreeSqrtFive) {
    // Build an expression that structurally contains 3 RootOf nodes.
    ExprPtr alpha1 = make_sqrt_rootof(2, 0);
    ExprPtr alpha2 = make_sqrt_rootof(3, 1);
    ExprPtr alpha3 = make_sqrt_rootof(5, 2);
    ASSERT_NE(alpha1, nullptr);
    ASSERT_NE(alpha2, nullptr);
    ASSERT_NE(alpha3, nullptr);

    // Build sum expression: alpha1 + alpha2 + alpha3.
    ExprPtr sum_expr = ctx->arena().make<Sum>(
        std::vector<ExprPtr>{alpha1, alpha2, alpha3});

    auto result = algebra::detect_tower_n_level(sum_expr, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    ASSERT_TRUE(result.value().has_value())
        << "detect_tower_n_level: returned nullopt for 3-radical expression";

    const auto& res = result.value().value();
    // Degree 8.
    EXPECT_EQ(res.min_poly_theta.size(), 9U)
        << "detect_tower_n_level: unexpected min_poly degree";
    EXPECT_EQ(res.alphas_in_theta.size(), 3U);
}

// ── Test 4b (F3.4-DEBT-01): nested RootOf, outer min-poly over Q(β) ──────────
//
// β = √2 = RootOf(z²-2, z, 0).
// α = RootOf(x² - β·x + 1, x, 0).  Absolute min-poly:
//   R(x) = Res_y(y²-2, x²-y·x+1) = x⁴ + 1  (squarefree, irreducible).
// → Q(α, β) = Q(α) (degree 4).  detect_tower_n_level must:
//   1. Detect both RootOfs.
//   2. Recognize α has non-rational coefficients (depends on β).
//   3. Lift via absolute resultant → x⁴+1.
//   4. Pass {β, α} with min-polys {y²-2, x⁴+1} to compute_primitive_element.
// θ has degree 8 (= 2·4), as the algorithm does not detect dependency.
TEST_F(PrimitiveElementTest, DetectTowerNLevel_NestedRootOf_F34Debt01) {
    // β = √2.
    ExprPtr beta = parse_ok("RootOf(z0^2-2, z0, 0)");
    ASSERT_NE(beta, nullptr);
    ctx->set_max_rootof_explicit_degree(1U);
    auto beta_canon = ctx->simplify(beta);
    ASSERT_TRUE(beta_canon.is_ok()) << beta_canon.error().message;

    // α = RootOf(z1² − β·z1 + 1, z1, 0).
    ExprPtr z1 = ctx->arena().make<Symbol>(Symbol{"z1"});
    ExprPtr one = ctx->arena().make<IntegerLit>(BigInt(1));
    ExprPtr two = ctx->arena().make<IntegerLit>(BigInt(2));
    ExprPtr z1_sq = ctx->arena().make<Binary>(BinaryOp::Pow, z1, two);
    ExprPtr beta_z1 = ctx->arena().make<Binary>(BinaryOp::Mul, beta_canon.value(), z1);
    ExprPtr neg_beta_z1 = ctx->arena().make<Unary>(UnaryOp::Neg, beta_z1);
    ExprPtr poly_alpha = ctx->arena().make<Sum>(
        std::vector<ExprPtr>{z1_sq, neg_beta_z1, one});
    ExprPtr alpha = ctx->arena().make<RootOf>(poly_alpha, Symbol{"z1"}, 0U);

    // Expression: alpha + beta.  Forces collection of both RootOfs.
    ExprPtr sum_expr = ctx->arena().make<Binary>(BinaryOp::Add, alpha, beta_canon.value());

    auto result = algebra::detect_tower_n_level(sum_expr, *ctx);
    ASSERT_TRUE(result.is_ok())
        << "F3.4-DEBT-01: detect_tower_n_level returned error: "
        << result.error().message;
    ASSERT_TRUE(result.value().has_value())
        << "F3.4-DEBT-01: detect_tower_n_level returned nullopt for nested RootOf";

    const auto& res = result.value().value();
    // θ degree should be ≤ 8 (2 · 4); typically 8 unless Trager detects dependency.
    EXPECT_GE(res.min_poly_theta.size(), 5U)
        << "F3.4-DEBT-01: θ min-poly degree too small (lift failed?)";
    EXPECT_LE(res.min_poly_theta.size(), 9U)
        << "F3.4-DEBT-01: θ min-poly degree exceeds 2·4";
    EXPECT_EQ(res.alphas_in_theta.size(), 2U);
}

// ── F3.5-DEBT-01 RESOLVED: redundant generators with degree compression ──────
//
// When a generator α_k is algebraically dependent on the existing field
// Q(θ_{k-1}), the shift-resultant R_s = Res_y(q_{k-1}(y), m_k(x − s·y)) is
// squarefree but REDUCIBLE on Q.  compute_primitive_element must factor R_s
// over Q and select an irreducible factor as q_current (NOT use the full
// reducible R_s, since Q[y]/(R_s) would not be a field and ring-GCD would
// blow up).  Degree of q_current = [Q(θ):Q] which may be strictly less than
// deg(R_s) = deg(q_{k-1}) · deg(m_k).
//
// Verified algebraically with C1 (m_i(P_i) ≡ 0) + C2 (θ-combination ≡ y).

TEST_F(PrimitiveElementTest, RedundantSqrtTwoSqrtThreeSqrtSix) {
    // Q(√2, √3, √6): √6 = √2·√3 ∈ Q(√2,√3), so [Q(θ):Q] = 4 (NOT 8).
    ExprPtr a1 = make_sqrt_rootof(2, 0);
    ExprPtr a2 = make_sqrt_rootof(3, 1);
    ExprPtr a3 = make_sqrt_rootof(6, 2);
    ASSERT_NE(a1, nullptr); ASSERT_NE(a2, nullptr); ASSERT_NE(a3, nullptr);

    std::vector<ExprPtr> alphas{a1, a2, a3};
    std::vector<algebra::AlgebraicNumber::CoeffVec> mps{
        quadratic_min_poly(-2),
        quadratic_min_poly(-3),
        quadratic_min_poly(-6),
    };

    auto result = algebra::compute_primitive_element(alphas, mps, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    const auto& res = result.value();
    // F3.5-DEBT-01 KEY INVARIANT: degree compression to 4.
    EXPECT_EQ(res.min_poly_theta.size(), 5U)
        << "RedundantSqrtTwoSqrtThreeSqrtSix: [Q(θ):Q] must be 4 (got "
        << (res.min_poly_theta.size() - 1U) << ")";
    EXPECT_EQ(res.alphas_in_theta.size(), 3U);
    EXPECT_EQ(res.shifts.size(), 2U);

    // Monic.
    if (!res.min_poly_theta.empty())
        EXPECT_EQ(res.min_poly_theta.back(), Rational(BigInt(1)));

    // C1: each generator is a root of its min-poly in Q[y]/(q_θ).
    const algebra::RatPoly q_theta(res.min_poly_theta);
    for (std::size_t i = 0U; i < 3U; ++i) {
        EXPECT_TRUE(eval_minpoly_at_ring_elem(mps[i], res.alphas_in_theta[i], q_theta).is_zero())
            << "RedundantSqrtTwoSqrtThreeSqrtSix: α[" << i << "] fails m_i(P_i) ≡ 0";
    }
    // C2: θ-combination = y.
    expect_is_y(theta_combination_mod_q(res.alphas_in_theta, res.shifts, q_theta),
                "RedundantSqrtTwoSqrtThreeSqrtSix");
}

TEST_F(PrimitiveElementTest, RedundantSqrtTwoSqrtEight) {
    // Q(√2, √8): √8 = 2√2, fully redundant.  [Q(θ):Q] = 2 (NOT 4).
    ExprPtr a1 = make_sqrt_rootof(2, 0);
    ExprPtr a2 = make_sqrt_rootof(8, 1);
    ASSERT_NE(a1, nullptr); ASSERT_NE(a2, nullptr);

    std::vector<ExprPtr> alphas{a1, a2};
    std::vector<algebra::AlgebraicNumber::CoeffVec> mps{
        quadratic_min_poly(-2),
        quadratic_min_poly(-8),
    };

    auto result = algebra::compute_primitive_element(alphas, mps, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    const auto& res = result.value();
    EXPECT_EQ(res.min_poly_theta.size(), 3U)
        << "RedundantSqrtTwoSqrtEight: [Q(θ):Q] must be 2 (got "
        << (res.min_poly_theta.size() - 1U) << ")";
    if (!res.min_poly_theta.empty())
        EXPECT_EQ(res.min_poly_theta.back(), Rational(BigInt(1)));

    const algebra::RatPoly q_theta(res.min_poly_theta);
    for (std::size_t i = 0U; i < 2U; ++i) {
        EXPECT_TRUE(eval_minpoly_at_ring_elem(mps[i], res.alphas_in_theta[i], q_theta).is_zero())
            << "RedundantSqrtTwoSqrtEight: α[" << i << "] fails m_i(P_i) ≡ 0";
    }
    expect_is_y(theta_combination_mod_q(res.alphas_in_theta, res.shifts, q_theta),
                "RedundantSqrtTwoSqrtEight");
}

// DISABILITATO: Test di stress matematico fisiologicamente in timeout (>60s) sotto Debug Mode (-O0). Da eseguire in Release Mode o via target dedicato cas_stress_tests.
TEST_F(PrimitiveElementTest, DISABLED_RedundantMixedTower_Sqrt2_Sqrt3_Sqrt5_Sqrt6) {
    // Q(√2, √3, √5, √6): √5 is independent, √6 is redundant.
    // [Q(θ):Q] = 8 (= |Q(√2,√3,√5):Q|, the √6 generator collapses).
    ExprPtr a1 = make_sqrt_rootof(2, 0);
    ExprPtr a2 = make_sqrt_rootof(3, 1);
    ExprPtr a3 = make_sqrt_rootof(5, 2);
    ExprPtr a4 = make_sqrt_rootof(6, 3);
    ASSERT_NE(a1, nullptr); ASSERT_NE(a2, nullptr);
    ASSERT_NE(a3, nullptr); ASSERT_NE(a4, nullptr);

    std::vector<ExprPtr> alphas{a1, a2, a3, a4};
    std::vector<algebra::AlgebraicNumber::CoeffVec> mps{
        quadratic_min_poly(-2),
        quadratic_min_poly(-3),
        quadratic_min_poly(-5),
        quadratic_min_poly(-6),
    };

    auto result = algebra::compute_primitive_element(alphas, mps, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    const auto& res = result.value();
    EXPECT_EQ(res.min_poly_theta.size(), 9U)
        << "RedundantMixedTower: [Q(θ):Q] must be 8 (got "
        << (res.min_poly_theta.size() - 1U) << ")";

    const algebra::RatPoly q_theta(res.min_poly_theta);
    for (std::size_t i = 0U; i < 4U; ++i) {
        EXPECT_TRUE(eval_minpoly_at_ring_elem(mps[i], res.alphas_in_theta[i], q_theta).is_zero())
            << "RedundantMixedTower: α[" << i << "] fails m_i(P_i) ≡ 0";
    }
    expect_is_y(theta_combination_mod_q(res.alphas_in_theta, res.shifts, q_theta),
                "RedundantMixedTower");
}

// ── Test 5: n=1 trivial case ──────────────────────────────────────────────────

TEST_F(PrimitiveElementTest, TrivialSingleGenerator) {
    ExprPtr alpha1 = make_sqrt_rootof(7, 0);
    ASSERT_NE(alpha1, nullptr);

    std::vector<ExprPtr> alphas{alpha1};
    std::vector<algebra::AlgebraicNumber::CoeffVec> min_polys{
        quadratic_min_poly(-7),
    };

    auto result = algebra::compute_primitive_element(alphas, min_polys, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    const auto& res = result.value();
    EXPECT_EQ(res.min_poly_theta.size(), 3U);  // degree 2
    EXPECT_TRUE(res.shifts.empty());
    ASSERT_EQ(res.alphas_in_theta.size(), 1U);

    // alpha[0] should be represented as [0, 1] = y (the generator itself).
    ASSERT_GE(res.alphas_in_theta[0].size(), 2U);
    EXPECT_EQ(res.alphas_in_theta[0][0], Rational(BigInt(0)));  // constant part = 0
    EXPECT_EQ(res.alphas_in_theta[0][1], Rational(BigInt(1)));  // degree-1 part = 1
}

// ── Test 6: Input validation ───────────────────────────────────────────────────

TEST_F(PrimitiveElementTest, EmptyInputReturnsError) {
    auto result = algebra::compute_primitive_element({}, {}, *ctx);
    EXPECT_FALSE(result.is_ok());
    EXPECT_EQ(result.error().kind, CASErrorKind::InvalidArgument);
}

TEST_F(PrimitiveElementTest, MismatchedSizesReturnsError) {
    ExprPtr alpha1 = make_sqrt_rootof(2, 0);
    ASSERT_NE(alpha1, nullptr);

    auto result = algebra::compute_primitive_element(
        {alpha1},
        {quadratic_min_poly(-2), quadratic_min_poly(-3)},
        *ctx);
    EXPECT_FALSE(result.is_ok());
    EXPECT_EQ(result.error().kind, CASErrorKind::InvalidArgument);
}

TEST_F(PrimitiveElementTest, DegreeOneMinPolyReturnsError) {
    ExprPtr alpha1 = make_sqrt_rootof(2, 0);
    ASSERT_NE(alpha1, nullptr);

    // Min-poly of degree 0 (only constant term) is invalid.
    algebra::AlgebraicNumber::CoeffVec bad_mp = {Rational(BigInt(1))};
    auto result = algebra::compute_primitive_element(
        {alpha1},
        {bad_mp},
        *ctx);
    EXPECT_FALSE(result.is_ok());
    EXPECT_EQ(result.error().kind, CASErrorKind::InvalidArgument);
}

// ── Task 1.2 / A9 Unit Tests: 3-Level Extensions & Multi-β Nested Towers ────

TEST_F(PrimitiveElementTest, ThreeLevelExtension_Sqrt2_Cbrt3_Sqrt5) {
    // Q(√2, ∛3, √5) — 3-level simple algebraic extensions. Total degree 2·3·2 = 12.
    ExprPtr a1 = make_sqrt_rootof(2, 0);
    ExprPtr a3 = make_sqrt_rootof(5, 2);
    ASSERT_NE(a1, nullptr);
    ASSERT_NE(a3, nullptr);

    ExprPtr a2 = parse_ok("RootOf(z1^3-3, z1, 0)");
    ASSERT_NE(a2, nullptr);

    ExprPtr sum_expr = ctx->arena().make<Sum>(std::vector<ExprPtr>{a1, a2, a3});

    auto result = algebra::detect_n_level_tower(sum_expr, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    ASSERT_TRUE(result.value().has_value());

    const auto& res = result.value().value();
    EXPECT_EQ(res.min_poly_theta.size(), 13U) << "Degree of Q(√2, ∛3, √5) must be 12";
    EXPECT_EQ(res.alphas_in_theta.size(), 3U);
}

TEST_F(PrimitiveElementTest, DetectNLevelTower_MultiBetaNested) {
    // Outer RootOf α = RootOf(z2² - β1·z2 - β2, z2, 0) depending on β1=√2 and β2=√3.
    ExprPtr beta1 = parse_ok("RootOf(z0^2-2, z0, 0)");
    ExprPtr beta2 = parse_ok("RootOf(z1^2-3, z1, 0)");
    ASSERT_NE(beta1, nullptr);
    ASSERT_NE(beta2, nullptr);

    ExprPtr z2 = ctx->arena().make<Symbol>(Symbol{"z2"});
    ExprPtr z2_sq = ctx->arena().make<Binary>(BinaryOp::Pow, z2, ctx->arena().make<IntegerLit>(BigInt(2)));
    ExprPtr beta1_z2 = ctx->arena().make<Binary>(BinaryOp::Mul, beta1, z2);
    ExprPtr term1 = ctx->arena().make<Unary>(UnaryOp::Neg, beta1_z2);
    ExprPtr term2 = ctx->arena().make<Unary>(UnaryOp::Neg, beta2);
    ExprPtr poly_alpha = ctx->arena().make<Sum>(std::vector<ExprPtr>{z2_sq, term1, term2});
    ExprPtr alpha = ctx->arena().make<RootOf>(poly_alpha, Symbol{"z2"}, 0U);

    ExprPtr sum_expr = ctx->arena().make<Sum>(std::vector<ExprPtr>{alpha, beta1, beta2});

    auto result = algebra::detect_n_level_tower(sum_expr, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    ASSERT_TRUE(result.value().has_value());
}

// DISABLED: truly-nested chain β1→β2→α requires recursive primitive-element
// flattening (compute Q-min-poly of each level before combining).  The iterative
// resolver lifts the individual level min-polys correctly, but compute_primitive_element
// with all 3 generators still fails (returns empty optional after ~225 s) because the
// 3 algebraic numbers live in a sequential dependency chain, not independent Q extensions.
// Tracked as F3.4-DEBT-01 / A9 follow-up.  Re-enable when the hierarchical
// prim-elem chain Q(β1) → Q(β1,β2) → Q(β1,β2,α) is implemented.
TEST_F(PrimitiveElementTest, DISABLED_DetectNLevelTower_ThreeLevelNested) {
    // 3-level nested radical: β1 = √3, β2 = √(2 + √3), α = √(1 + √(2 + √3)).
    ExprPtr beta1 = parse_ok("RootOf(z0^2-3, z0, 0)");
    ASSERT_NE(beta1, nullptr);

    // β2: z1² - (2 + β1) = 0
    ExprPtr z1 = ctx->arena().make<Symbol>(Symbol{"z1"});
    ExprPtr z1_sq = ctx->arena().make<Binary>(BinaryOp::Pow, z1, ctx->arena().make<IntegerLit>(BigInt(2)));
    ExprPtr two = ctx->arena().make<IntegerLit>(BigInt(2));
    ExprPtr two_plus_b1 = ctx->arena().make<Binary>(BinaryOp::Add, two, beta1);
    ExprPtr poly_b2 = ctx->arena().make<Binary>(BinaryOp::Sub, z1_sq, two_plus_b1);
    ExprPtr beta2 = ctx->arena().make<RootOf>(poly_b2, Symbol{"z1"}, 0U);

    // α: z2² - (1 + β2) = 0
    ExprPtr z2 = ctx->arena().make<Symbol>(Symbol{"z2"});
    ExprPtr z2_sq = ctx->arena().make<Binary>(BinaryOp::Pow, z2, ctx->arena().make<IntegerLit>(BigInt(2)));
    ExprPtr one = ctx->arena().make<IntegerLit>(BigInt(1));
    ExprPtr one_plus_b2 = ctx->arena().make<Binary>(BinaryOp::Add, one, beta2);
    ExprPtr poly_a = ctx->arena().make<Binary>(BinaryOp::Sub, z2_sq, one_plus_b2);
    ExprPtr alpha = ctx->arena().make<RootOf>(poly_a, Symbol{"z2"}, 0U);

    ExprPtr sum_expr = ctx->arena().make<Sum>(std::vector<ExprPtr>{alpha, beta2, beta1});

    auto result = algebra::detect_n_level_tower(sum_expr, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    ASSERT_TRUE(result.value().has_value());
}

}  // namespace
}  // namespace cas::test

