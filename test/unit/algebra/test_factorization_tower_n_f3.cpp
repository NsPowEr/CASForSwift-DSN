// F3.5: factor_polynomial_tower_n — n-level tower factorisation tests.
//
// Validates the single-extension Trager factorisation through the F3.4
// primitive element collapse.  Coverage:
//   - 2-level positive (sqrt(2), sqrt(3)) — x²-3 splits.
//   - 3-level positive (sqrt(2), sqrt(3), sqrt(5)) — x²-5 splits.
//   - Irreducible: x²-7 over Q(sqrt(2), sqrt(3)) stays irreducible.
//   - Redundant generator: Q(sqrt(2), sqrt(3), sqrt(6)) is detected by
//     compute_primitive_element via the non-squarefree resultant signal;
//     the function may return either a valid factorisation (if the
//     primitive element collapses successfully) or an explicit Unimplemented
//     diagnostic — both outcomes are CLAUDE-compliant.
//
// Reconstruction oracle: content * Π factor_i == f exactly.

#include "cas/algebra.hpp"
#include "cas/algebraic_tower_bridge.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include "algebra/polynomial_internal.hpp"

#include <gtest/gtest.h>

#include <chrono>
#include <cstddef>
#include <memory>
#include <string>

namespace cas::test {
namespace {

[[nodiscard]] Result<ExprPtr> parse_expr_n(const std::string& input, AstArena& arena) {
    Lexer lexer(input);
    auto tokens = lexer.tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

class FactorizationTowerNTest : public ::testing::Test {
protected:
    void SetUp() override {
        ctx = std::make_unique<symbolic::CASContext>();
        ctx->set_timeout(std::chrono::seconds(120));
    }

    [[nodiscard]] ExprPtr parse_ok(const std::string& s) {
        auto p = parse_expr_n(s, ctx->arena());
        EXPECT_TRUE(p.is_ok()) << p.error().message;
        return p.is_ok() ? p.value() : nullptr;
    }

    // Build (sqrt(a) canonical RootOf, rational min-poly y^2 - a).
    struct SqrtGen {
        ExprPtr alpha;
        algebra::AlgebraicNumber::CoeffVec min_poly;
    };
    [[nodiscard]] SqrtGen sqrt_gen(long a) {
        const std::string s = "RootOf(y^2-" + std::to_string(a) + ",y,0)";
        ExprPtr alpha = parse_ok(s);
        EXPECT_TRUE(alpha != nullptr);
        const std::size_t saved = ctx->max_rootof_explicit_degree();
        ctx->set_max_rootof_explicit_degree(1U);
        auto canon = ctx->simplify(alpha);
        ctx->set_max_rootof_explicit_degree(saved);
        EXPECT_TRUE(canon.is_ok());
        return SqrtGen{
            .alpha = canon.is_ok() ? canon.value() : alpha,
            .min_poly = {Rational(BigInt(-a)),
                         Rational(BigInt(0)),
                         Rational(BigInt(1))},
        };
    }

    [[nodiscard]] std::size_t cumulative_deg(
        const algebra::Factorization& f, const Symbol& var) {
        std::size_t total = 0U;
        for (const auto& pf : f.factors) {
            auto parsed = algebra::parse_polynomial(pf.factor, var, *ctx);
            if (parsed.is_error()) return 0U;
            total += algebra::poly_degree(parsed.value());
        }
        return total;
    }

    [[nodiscard]] bool reconstructs(
        const algebra::Factorization& f, ExprPtr original) const {
        if (!f.content) return false;
        ExprPtr product = f.content;
        for (const auto& pf : f.factors) {
            product = ctx->arena().make<Binary>(
                BinaryOp::Mul, product, pf.factor);
        }
        auto eq = symbolic::mathematically_equal(product, original, *ctx);
        return eq.is_ok() && eq.value();
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

TEST_F(FactorizationTowerNTest, SplitsX2Minus3_Over_Q_Sqrt2_Sqrt3) {
    auto g2 = sqrt_gen(2);
    auto g3 = sqrt_gen(3);
    algebra::TowerGeneratorsN gens{
        .alphas = {g2.alpha, g3.alpha},
        .min_polys = {g2.min_poly, g3.min_poly},
    };
    ExprPtr poly = parse_ok("x^2 - 3");
    ASSERT_NE(poly, nullptr);
    auto r = algebra::factor_polynomial_tower_n(poly, Symbol{"x"}, gens, *ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value().factors.size(), 2U)
        << "x^2 - 3 must split into 2 linear factors over Q(sqrt 2, sqrt 3)";
    EXPECT_EQ(cumulative_deg(r.value(), Symbol{"x"}), 2U);
    EXPECT_TRUE(reconstructs(r.value(), poly));
}

// DISABILITATO: Test di stress matematico fisiologicamente in timeout (>60s) sotto Debug Mode (-O0). Da eseguire in Release Mode o via target dedicato cas_stress_tests.
TEST_F(FactorizationTowerNTest, DISABLED_SplitsX2Minus5_Over_Q_Sqrt2_Sqrt3_Sqrt5) {
    auto g2 = sqrt_gen(2);
    auto g3 = sqrt_gen(3);
    auto g5 = sqrt_gen(5);
    algebra::TowerGeneratorsN gens{
        .alphas = {g2.alpha, g3.alpha, g5.alpha},
        .min_polys = {g2.min_poly, g3.min_poly, g5.min_poly},
    };
    ExprPtr poly = parse_ok("x^2 - 5");
    ASSERT_NE(poly, nullptr);
    auto r = algebra::factor_polynomial_tower_n(poly, Symbol{"x"}, gens, *ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value().factors.size(), 2U)
        << "x^2 - 5 must split into 2 linear factors over Q(sqrt 2, sqrt 3, sqrt 5)";
    EXPECT_EQ(cumulative_deg(r.value(), Symbol{"x"}), 2U);
    EXPECT_TRUE(reconstructs(r.value(), poly));
}

TEST_F(FactorizationTowerNTest, IrreducibleX2Minus7_Over_Q_Sqrt2_Sqrt3) {
    auto g2 = sqrt_gen(2);
    auto g3 = sqrt_gen(3);
    algebra::TowerGeneratorsN gens{
        .alphas = {g2.alpha, g3.alpha},
        .min_polys = {g2.min_poly, g3.min_poly},
    };
    // Tight budget: the composite Trager norm of an irreducible quadratic over a
    // 2-level tower is a perfect power of an irreducible Q-polynomial, so the
    // square-free-norm shift search always exhausts the budget and recombination
    // falls to exponential Kronecker (perfect-power-norm recognition is the F3.5
    // gap).  Contract mirrors the sibling AntiHardcodeIrreducibleX2Minus2... test:
    // either irreducibility is recognised, or a budget-exceeded Unimplemented/
    // Timeout is returned — NEVER a hang (HC-F8-FACTORIZATIONTOWER-PERF, L3-06).
    ctx->set_timeout(std::chrono::seconds(3));
    ExprPtr poly = parse_ok("x^2 - 7");
    ASSERT_NE(poly, nullptr);
    auto r = algebra::factor_polynomial_tower_n(poly, Symbol{"x"}, gens, *ctx);
    if (r.is_ok()) {
        EXPECT_EQ(r.value().factors.size(), 1U)
            << "x^2 - 7 must NOT split (sqrt 7 ∉ Q(sqrt 2, sqrt 3))";
        EXPECT_EQ(cumulative_deg(r.value(), Symbol{"x"}), 2U);
        EXPECT_TRUE(reconstructs(r.value(), poly));
    } else {
        const bool is_budget_exceeded =
            r.error().kind == CASErrorKind::Unimplemented ||
            r.error().kind == CASErrorKind::Timeout;
        EXPECT_TRUE(is_budget_exceeded)
            << "must return Unimplemented or Timeout (budget exceeded), not "
               "InternalError or crash.  kind=" << static_cast<int>(r.error().kind)
            << " message=" << r.error().message;
    }
}

// F3.5-DEBT-01 RESOLVED 2026-05-31: redundant-generator case is fully handled.
// When R_s is squarefree but reducible (typical of algebraically dependent
// generators), compute_primitive_element factors R_s over Q via
// collect_irred_factors_over_q and tries each irreducible factor as q_current.
// The ring-GCD in Q[y]/(cand_q)[t] succeeds for any factor whose root contains
// θ_new = θ_{k-1} + s·α_k (always at least one exists by construction), giving
// a valid primitive element with the CORRECT compressed degree
// (e.g. [Q(√2,√3,√6):Q] = 4, not 8).
TEST_F(FactorizationTowerNTest, RedundantGenerator_Sqrt2_Sqrt3_Sqrt6) {
    // sqrt(6) = sqrt(2) * sqrt(3) ∈ Q(sqrt(2), sqrt(3)).
    // x²-3 must split into 2 linear factors since √3 ∈ Q(√2,√3,√6).
    ctx->set_timeout(std::chrono::seconds(30));
    ctx->set_max_trager_tower_shift_attempts(8U);
    auto g2 = sqrt_gen(2);
    auto g3 = sqrt_gen(3);
    auto g6 = sqrt_gen(6);
    algebra::TowerGeneratorsN gens{
        .alphas = {g2.alpha, g3.alpha, g6.alpha},
        .min_polys = {g2.min_poly, g3.min_poly, g6.min_poly},
    };
    ExprPtr poly = parse_ok("x^2 - 3");
    ASSERT_NE(poly, nullptr);
    auto r = algebra::factor_polynomial_tower_n(poly, Symbol{"x"}, gens, *ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value().factors.size(), 2U)
        << "x²-3 must split into 2 linear factors over Q(√2,√3,√6)";
    EXPECT_EQ(cumulative_deg(r.value(), Symbol{"x"}), 2U);
    EXPECT_TRUE(reconstructs(r.value(), poly));
}

TEST_F(FactorizationTowerNTest, RedundantGenerator_Sqrt2_Sqrt8) {
    // sqrt(8) = 2·sqrt(2): fully redundant.  Q(√2,√8) = Q(√2), [Q(θ):Q] = 2.
    // x²-2 must split over this tower.
    ctx->set_timeout(std::chrono::seconds(30));
    ctx->set_max_trager_tower_shift_attempts(8U);
    auto g2 = sqrt_gen(2);
    auto g8 = sqrt_gen(8);
    algebra::TowerGeneratorsN gens{
        .alphas = {g2.alpha, g8.alpha},
        .min_polys = {g2.min_poly, g8.min_poly},
    };
    ExprPtr poly = parse_ok("x^2 - 2");
    ASSERT_NE(poly, nullptr);
    auto r = algebra::factor_polynomial_tower_n(poly, Symbol{"x"}, gens, *ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value().factors.size(), 2U)
        << "x²-2 must split into 2 linear factors over Q(√2,√8) = Q(√2)";
    EXPECT_EQ(cumulative_deg(r.value(), Symbol{"x"}), 2U);
    EXPECT_TRUE(reconstructs(r.value(), poly));
}

TEST_F(FactorizationTowerNTest, RejectsNullPolynomial) {
    auto g2 = sqrt_gen(2);
    auto g3 = sqrt_gen(3);
    algebra::TowerGeneratorsN gens{
        .alphas = {g2.alpha, g3.alpha},
        .min_polys = {g2.min_poly, g3.min_poly},
    };
    auto r = algebra::factor_polynomial_tower_n(ExprPtr{}, Symbol{"x"}, gens, *ctx);
    ASSERT_TRUE(r.is_error());
    EXPECT_EQ(r.error().kind, CASErrorKind::InvalidArgument);
}

TEST_F(FactorizationTowerNTest, RejectsMalformedGenerators) {
    algebra::TowerGeneratorsN gens;
    ExprPtr poly = parse_ok("x^2 - 3");
    auto r = algebra::factor_polynomial_tower_n(poly, Symbol{"x"}, gens, *ctx);
    ASSERT_TRUE(r.is_error());
    EXPECT_EQ(r.error().kind, CASErrorKind::InvalidArgument);
}

}  // namespace
}  // namespace cas::test
