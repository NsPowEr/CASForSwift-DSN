#include "cas/algebraic_tower_bridge.hpp"

#include "cas/algebraic_number_bridge.hpp"
#include "cas/ast_debug.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "algebra/polynomial_internal.hpp"

#include <gtest/gtest.h>

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

class AlgebraicTowerBridgeTest : public ::testing::Test {
protected:
    void SetUp() override { ctx = std::make_unique<symbolic::CASContext>(); }

    [[nodiscard]] ExprPtr parse_ok(const std::string& input) {
        auto parsed = parse_expr(input, ctx->arena());
        EXPECT_TRUE(parsed.is_ok()) << parsed.error().message;
        return parsed.is_ok() ? parsed.value() : nullptr;
    }

    [[nodiscard]] algebra::AlgebraicNumber one_in_q_alpha(
        const algebra::AlgebraicNumber::CoeffVec& min_poly) const {
        return algebra::AlgebraicNumber({Rational(BigInt(1))}, min_poly);
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

TEST_F(AlgebraicTowerBridgeTest, DetectsIndependentSqrt2Sqrt3WithSortedGenerators) {
    ExprPtr expr = parse_ok("RootOf(y^2-2,y,0) + RootOf(z^2-3,z,0)");

    auto detected = algebra::detect_two_level_tower(expr, *ctx);
    ASSERT_TRUE(detected.is_ok()) << detected.error().message;
    ASSERT_TRUE(detected.value().has_value());

    const auto& gens = detected.value().value();
    const auto* root_1 = expr_cast<RootOf>(gens.alpha_1);
    const auto* root_2 = expr_cast<RootOf>(gens.alpha_2);
    ASSERT_NE(root_1, nullptr);
    ASSERT_NE(root_2, nullptr);

    auto min_poly_1 = algebra::rootof_min_poly(*root_1, *ctx);
    auto min_poly_2 = algebra::rootof_min_poly(*root_2, *ctx);
    ASSERT_TRUE(min_poly_1.is_ok()) << min_poly_1.error().message;
    ASSERT_TRUE(min_poly_2.is_ok()) << min_poly_2.error().message;
    EXPECT_EQ(min_poly_1.value(), algebra::AlgebraicNumber::CoeffVec({
        Rational(BigInt(-2)), Rational(BigInt(0)), Rational(BigInt(1))}));
    EXPECT_EQ(min_poly_2.value(), algebra::AlgebraicNumber::CoeffVec({
        Rational(BigInt(-3)), Rational(BigInt(0)), Rational(BigInt(1))}));
}

TEST_F(AlgebraicTowerBridgeTest, ExpressesSqrt2PlusSqrt3InTowerBasis) {
    ExprPtr expr = parse_ok("RootOf(y^2-2,y,0) + RootOf(z^2-3,z,0)");

    auto detected = algebra::detect_two_level_tower(expr, *ctx);
    ASSERT_TRUE(detected.is_ok()) << detected.error().message;
    ASSERT_TRUE(detected.value().has_value());

    auto tower = algebra::try_express_in_tower_two_level(expr, *detected.value(), *ctx);
    ASSERT_TRUE(tower.is_ok()) << tower.error().message;
    ASSERT_TRUE(tower.value().has_value());

    auto alpha_1 = algebra::alpha_from_rootof(expr_ref<RootOf>(detected.value()->alpha_1), *ctx);
    ASSERT_TRUE(alpha_1.is_ok()) << alpha_1.error().message;

    const auto& coeffs = tower.value()->value();
    ASSERT_EQ(coeffs.size(), 2U);
    EXPECT_EQ(coeffs[0], alpha_1.value());
    EXPECT_EQ(coeffs[1], one_in_q_alpha(detected.value()->min_poly_1));
}

TEST_F(AlgebraicTowerBridgeTest, RoundTripsTowerBackToExpression) {
    ExprPtr expr = parse_ok("(RootOf(y^2-2,y,0) + RootOf(z^2-3,z,0))^2");

    auto detected = algebra::detect_two_level_tower(expr, *ctx);
    ASSERT_TRUE(detected.is_ok()) << detected.error().message;
    ASSERT_TRUE(detected.value().has_value());

    auto tower = algebra::try_express_in_tower_two_level(expr, *detected.value(), *ctx);
    ASSERT_TRUE(tower.is_ok()) << tower.error().message;
    ASSERT_TRUE(tower.value().has_value());

    ExprPtr rebuilt_raw = algebra::tower_to_expr(*tower.value(), *detected.value(), ctx->arena());
    auto rebuilt = ctx->simplify(rebuilt_raw);
    ASSERT_TRUE(rebuilt.is_ok()) << rebuilt.error().message;

    auto redetected = algebra::detect_two_level_tower(rebuilt_raw, *ctx);
    ASSERT_TRUE(redetected.is_ok()) << redetected.error().message;
    ASSERT_TRUE(redetected.value().has_value());
    EXPECT_TRUE(structural_equal(redetected.value()->alpha_1, detected.value()->alpha_1));
    EXPECT_TRUE(structural_equal(redetected.value()->alpha_2, detected.value()->alpha_2));
    EXPECT_EQ(redetected.value()->min_poly_1, detected.value()->min_poly_1);
    EXPECT_EQ(redetected.value()->min_poly_2, detected.value()->min_poly_2);

    auto eq = symbolic::mathematically_equal(rebuilt.value(), expr, *ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value()) << debug_print(rebuilt.value());
}

TEST_F(AlgebraicTowerBridgeTest, RejectsThreeIndependentRootOfGenerators) {
    ExprPtr expr = parse_ok(
        "RootOf(a^2-2,a,0) + RootOf(b^2-3,b,0) + RootOf(c^2-5,c,0)");

    auto detected = algebra::detect_two_level_tower(expr, *ctx);
    ASSERT_TRUE(detected.is_ok()) << detected.error().message;
    EXPECT_FALSE(detected.value().has_value());
}

TEST_F(AlgebraicTowerBridgeTest, DetectsDependentOuterGeneratorOverQAlpha) {
    ExprPtr beta = parse_ok("RootOf(y^2-(3*RootOf(x^2-2,x,0)+1),y,0)");

    auto detected = algebra::detect_two_level_tower(beta, *ctx);
    ASSERT_TRUE(detected.is_ok()) << detected.error().message;
    ASSERT_TRUE(detected.value().has_value());

    auto alpha_min_poly = algebra::rootof_min_poly(expr_ref<RootOf>(detected.value()->alpha_1), *ctx);
    ASSERT_TRUE(alpha_min_poly.is_ok()) << alpha_min_poly.error().message;
    EXPECT_EQ(alpha_min_poly.value(), algebra::AlgebraicNumber::CoeffVec({
        Rational(BigInt(-2)), Rational(BigInt(0)), Rational(BigInt(1))}));

    ExprPtr expected_beta_poly = parse_ok("y^2-(3*RootOf(x^2-2,x,0)+1)");
    auto beta_eq = symbolic::mathematically_equal(
        expr_ref<RootOf>(detected.value()->alpha_2).polynomial,
        expected_beta_poly,
        *ctx);
    ASSERT_TRUE(beta_eq.is_ok()) << beta_eq.error().message;
    EXPECT_TRUE(beta_eq.value()) << debug_print(expr_ref<RootOf>(detected.value()->alpha_2).polynomial);

    auto tower = algebra::try_express_in_tower_two_level(beta, *detected.value(), *ctx);
    ASSERT_TRUE(tower.is_ok()) << tower.error().message;
    ASSERT_TRUE(tower.value().has_value());

    const auto& coeffs = tower.value()->value();
    ASSERT_EQ(coeffs.size(), 2U);
    EXPECT_TRUE(coeffs[0].is_zero());
    EXPECT_EQ(coeffs[1], one_in_q_alpha(detected.value()->min_poly_1));
}

TEST_F(AlgebraicTowerBridgeTest, DependentCoefficientIsExpressibleInQAlpha) {
    ExprPtr alpha = parse_ok("RootOf(x^2-2,x,0)");
    ExprPtr coeff = parse_ok("-(3*RootOf(x^2-2,x,0)+1)");

    auto min_poly = algebra::rootof_min_poly(expr_ref<RootOf>(alpha), *ctx);
    ASSERT_TRUE(min_poly.is_ok()) << min_poly.error().message;

    auto reduced = algebra::try_express_in_q_alpha(coeff, alpha, min_poly.value(), *ctx);
    ASSERT_TRUE(reduced.is_ok()) << reduced.error().message;
    ASSERT_TRUE(reduced.value().has_value());

    const auto& value = reduced.value()->value();
    ASSERT_EQ(value.size(), 2U);
    EXPECT_EQ(value[0], Rational(BigInt(-1)));
    EXPECT_EQ(value[1], Rational(BigInt(-3)));
}

TEST_F(AlgebraicTowerBridgeTest, SimplifyWithRootDegreeGuardKeepsDependentCoefficientInRootOfForm) {
    ExprPtr alpha = parse_ok("RootOf(x^2-2,x,0)");
    ExprPtr coeff = parse_ok("-(3*RootOf(x^2-2,x,0)+1)");

    ctx->set_max_rootof_explicit_degree(1U);
    auto simplified = ctx->simplify(coeff);
    ASSERT_TRUE(simplified.is_ok()) << simplified.error().message;

    auto min_poly = algebra::rootof_min_poly(expr_ref<RootOf>(alpha), *ctx);
    ASSERT_TRUE(min_poly.is_ok()) << min_poly.error().message;

    auto reduced = algebra::try_express_in_q_alpha(simplified.value(), alpha, min_poly.value(), *ctx);
    ASSERT_TRUE(reduced.is_ok()) << reduced.error().message;
    ASSERT_TRUE(reduced.value().has_value()) << debug_print(simplified.value());
}

TEST_F(AlgebraicTowerBridgeTest, ManualGeneratorsExpressDependentOuterRoot) {
    ExprPtr alpha = parse_ok("RootOf(x^2-2,x,0)");
    ExprPtr beta = parse_ok("RootOf(y^2-(3*RootOf(x^2-2,x,0)+1),y,0)");
    ctx->set_max_rootof_explicit_degree(1U);
    auto alpha_norm = ctx->simplify(alpha);
    auto beta_norm = ctx->simplify(beta);
    ASSERT_TRUE(alpha_norm.is_ok()) << alpha_norm.error().message;
    ASSERT_TRUE(beta_norm.is_ok()) << beta_norm.error().message;

    auto min_poly_1 = algebra::rootof_min_poly(expr_ref<RootOf>(alpha), *ctx);
    ASSERT_TRUE(min_poly_1.is_ok()) << min_poly_1.error().message;

    auto coeff0 = parse_ok("-(3*RootOf(x^2-2,x,0)+1)");
    auto coeff0_norm = ctx->simplify(coeff0);
    ASSERT_TRUE(coeff0_norm.is_ok()) << coeff0_norm.error().message;
    auto coeff0_an = algebra::try_express_in_q_alpha(coeff0_norm.value(), alpha_norm.value(), min_poly_1.value(), *ctx);
    ASSERT_TRUE(coeff0_an.is_ok()) << coeff0_an.error().message;
    ASSERT_TRUE(coeff0_an.value().has_value());

    algebra::TowerGenerators gens{
        .alpha_1 = alpha_norm.value(),
        .min_poly_1 = min_poly_1.value(),
        .alpha_2 = beta_norm.value(),
        .min_poly_2 = {
            coeff0_an.value().value(),
            algebra::AlgebraicNumber({Rational(BigInt(0))}, min_poly_1.value()),
            algebra::AlgebraicNumber({Rational(BigInt(1))}, min_poly_1.value()),
        },
    };

    auto tower = algebra::try_express_in_tower_two_level(beta_norm.value(), gens, *ctx);
    ASSERT_TRUE(tower.is_ok()) << tower.error().message;
    ASSERT_TRUE(tower.value().has_value());
}

TEST_F(AlgebraicTowerBridgeTest, DetectAndManualGeneratorsAgreeForDependentOuterRoot) {
    ExprPtr beta = parse_ok("RootOf(y^2-(3*RootOf(x^2-2,x,0)+1),y,0)");

    auto detected = algebra::detect_two_level_tower(beta, *ctx);
    ASSERT_TRUE(detected.is_ok()) << detected.error().message;
    ASSERT_TRUE(detected.value().has_value());

    ctx->set_max_rootof_explicit_degree(1U);
    ExprPtr alpha = parse_ok("RootOf(x^2-2,x,0)");
    auto alpha_norm = ctx->simplify(alpha);
    auto beta_norm = ctx->simplify(beta);
    ASSERT_TRUE(alpha_norm.is_ok()) << alpha_norm.error().message;
    ASSERT_TRUE(beta_norm.is_ok()) << beta_norm.error().message;

    auto min_poly_1 = algebra::rootof_min_poly(expr_ref<RootOf>(alpha), *ctx);
    ASSERT_TRUE(min_poly_1.is_ok()) << min_poly_1.error().message;

    ExprPtr coeff0 = parse_ok("-(3*RootOf(x^2-2,x,0)+1)");
    auto coeff0_norm = ctx->simplify(coeff0);
    ASSERT_TRUE(coeff0_norm.is_ok()) << coeff0_norm.error().message;
    auto coeff0_an = algebra::try_express_in_q_alpha(
        coeff0_norm.value(), alpha_norm.value(), min_poly_1.value(), *ctx);
    ASSERT_TRUE(coeff0_an.is_ok()) << coeff0_an.error().message;
    ASSERT_TRUE(coeff0_an.value().has_value());

    algebra::TowerGenerators manual{
        .alpha_1 = alpha_norm.value(),
        .min_poly_1 = min_poly_1.value(),
        .alpha_2 = beta_norm.value(),
        .min_poly_2 = {
            coeff0_an.value().value(),
            algebra::AlgebraicNumber({Rational(BigInt(0))}, min_poly_1.value()),
            algebra::AlgebraicNumber({Rational(BigInt(1))}, min_poly_1.value()),
        },
    };

    EXPECT_TRUE(structural_equal(detected.value()->alpha_1, manual.alpha_1));
    EXPECT_TRUE(structural_equal(detected.value()->alpha_2, manual.alpha_2));
    EXPECT_EQ(detected.value()->min_poly_1, manual.min_poly_1);
    EXPECT_EQ(detected.value()->min_poly_2, manual.min_poly_2);
}

}  // namespace
}  // namespace cas::test
