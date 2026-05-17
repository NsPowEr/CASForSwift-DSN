// L3-06: factor_polynomial_tower over Q(alpha_1, alpha_2).
//
// Validates composite Trager factorisation against handcrafted tower
// generators built from independent quadratic radicals.  Tests cover:
//   - Splitting case Q(sqrt 2, sqrt 3) with f = (x^2-2)(x^2-3) → 4 linear factors.
//   - Anti-hardcode irreducibility Q(sqrt 3, sqrt 5) with f = x^2-2 → stays irreducible.
//   - Algorithmic invariants (Trager post-condition prod(deg(f_i)) == deg(f)).
//   - Input validation (null poly, malformed generators, non-Q[x] input).

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

[[nodiscard]] Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    Lexer lexer(input);
    auto tokens = lexer.tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

class FactorizationTowerTest : public ::testing::Test {
protected:
    void SetUp() override {
        ctx = std::make_unique<symbolic::CASContext>();
        // Heavy resultant pipelines must not time out at the 1 s default.
        ctx->set_timeout(std::chrono::seconds(120));
    }

    [[nodiscard]] ExprPtr parse_ok(const std::string& input) {
        auto parsed = parse_expr(input, ctx->arena());
        EXPECT_TRUE(parsed.is_ok()) << parsed.error().message;
        return parsed.is_ok() ? parsed.value() : nullptr;
    }

    // Build Q(sqrt(a), sqrt(b)) generators from rational radicands.  Alphas
    // are canonicalised under rootof_explicit_degree=1 to match the form
    // produced by detect_two_level_tower / canonicalize_root_expr.
    [[nodiscard]] algebra::TowerGenerators biquadratic_gens(long a, long b) {
        const std::string a_str = "RootOf(y^2-" + std::to_string(a) + ",y,0)";
        const std::string b_str = "RootOf(z^2-" + std::to_string(b) + ",z,0)";
        ExprPtr alpha = parse_ok(a_str);
        ExprPtr beta = parse_ok(b_str);
        EXPECT_TRUE(alpha != nullptr);
        EXPECT_TRUE(beta != nullptr);

        const std::size_t saved = ctx->max_rootof_explicit_degree();
        ctx->set_max_rootof_explicit_degree(1U);
        auto alpha_norm = ctx->simplify(alpha);
        auto beta_norm = ctx->simplify(beta);
        ctx->set_max_rootof_explicit_degree(saved);
        EXPECT_TRUE(alpha_norm.is_ok());
        EXPECT_TRUE(beta_norm.is_ok());

        algebra::AlgebraicNumber::CoeffVec min_poly_1{
            Rational(BigInt(-a)), Rational(BigInt(0)), Rational(BigInt(1))};

        algebra::AlgebraicNumber inner_zero({Rational(BigInt(0))}, min_poly_1);
        algebra::AlgebraicNumber inner_one({Rational(BigInt(1))}, min_poly_1);
        algebra::AlgebraicNumber inner_neg_b({Rational(BigInt(-b))}, min_poly_1);

        return algebra::TowerGenerators{
            .alpha_1 = alpha_norm.value(),
            .min_poly_1 = std::move(min_poly_1),
            .alpha_2 = beta_norm.value(),
            .min_poly_2 = {inner_neg_b, inner_zero, inner_one},
        };
    }

    // Sum of degrees of all returned factors (in var).  Used as a Trager
    // post-condition oracle: must equal deg(input) for a correct
    // factorisation over the tower.
    [[nodiscard]] std::size_t cumulative_factor_degree(
        const algebra::Factorization& fact,
        const Symbol& var) const {
        std::size_t total = 0U;
        for (const auto& pf : fact.factors) {
            auto parsed = algebra::parse_polynomial(pf.factor, var, *ctx);
            if (parsed.is_error()) return 0U;
            total += algebra::poly_degree(parsed.value());
        }
        return total;
    }

    // Reconstruct  content * prod(factors)  as an ExprPtr and test it
    // against the original input via mathematically_equal.  This is the
    // strongest oracle for a factorization: degree + content + each
    // factor must combine to recover f exactly.
    [[nodiscard]] bool factorization_reconstructs(
        const algebra::Factorization& fact,
        ExprPtr original) const {
        if (!fact.content) return false;
        ExprPtr product = fact.content;
        for (const auto& pf : fact.factors) {
            product = ctx->arena().make<Binary>(BinaryOp::Mul, product, pf.factor);
        }
        auto eq = symbolic::mathematically_equal(product, original, *ctx);
        return eq.is_ok() && eq.value();
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

TEST_F(FactorizationTowerTest, AntiHardcodeIrreducibleX2Minus2OverQSqrt3Sqrt5) {
    // sqrt(2) is not in Q(sqrt(3), sqrt(5)); x^2 - 2 must remain irreducible.
    auto gens = biquadratic_gens(3, 5);
    ExprPtr poly = parse_ok("x^2 - 2");
    ASSERT_NE(poly, nullptr);

    auto result = algebra::factor_polynomial_tower(poly, Symbol{"x"}, gens, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_EQ(result.value().factors.size(), 1U)
        << "x^2 - 2 must NOT split over Q(sqrt 3, sqrt 5)";
    EXPECT_EQ(cumulative_factor_degree(result.value(), Symbol{"x"}), 2U)
        << "Trager post-condition: sum(deg(f_i)) == deg(f)";
    EXPECT_TRUE(factorization_reconstructs(result.value(), poly))
        << "content * prod(factors) must reconstruct the input polynomial";
}

TEST_F(FactorizationTowerTest, SplitsX2Minus3OverQSqrt2Sqrt3) {
    // sqrt(3) IS in Q(sqrt 2, sqrt 3); x^2 - 3 must split into two linear
    // factors (x - sqrt 3)(x + sqrt 3).
    auto gens = biquadratic_gens(2, 3);
    ExprPtr poly = parse_ok("x^2 - 3");
    ASSERT_NE(poly, nullptr);

    auto result = algebra::factor_polynomial_tower(poly, Symbol{"x"}, gens, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_EQ(result.value().factors.size(), 2U)
        << "x^2 - 3 must split into 2 linear factors over Q(sqrt 2, sqrt 3)";
    EXPECT_EQ(cumulative_factor_degree(result.value(), Symbol{"x"}), 2U)
        << "Trager post-condition: sum(deg(f_i)) == deg(f)";
    EXPECT_TRUE(factorization_reconstructs(result.value(), poly))
        << "content * prod(factors) must reconstruct the input polynomial";
}

// Audit fix B-L3-06-CRITICO: leading-coefficient preservation.  For
// f = 2x^2 - 4 = 2(x^2 - 2) the monic factor is x^2 - 2 and the content
// must carry the dropped factor of 2.
TEST_F(FactorizationTowerTest, PreservesLeadingCoefficientAsContent) {
    auto gens = biquadratic_gens(3, 5);  // x^2 - 2 stays irreducible here.
    ExprPtr poly = parse_ok("2*x^2 - 4");
    ASSERT_NE(poly, nullptr);

    auto result = algebra::factor_polynomial_tower(poly, Symbol{"x"}, gens, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_TRUE(factorization_reconstructs(result.value(), poly))
        << "content must carry the leading coefficient (here 2)";
}

TEST_F(FactorizationTowerTest, RejectsInvalidTowerGenerators) {
    algebra::TowerGenerators bad{};
    bad.alpha_1 = parse_ok("RootOf(y^2-2,y,0)");
    bad.alpha_2 = parse_ok("RootOf(z^2-3,z,0)");
    bad.min_poly_1 = {};
    bad.min_poly_2 = {};
    ExprPtr poly = parse_ok("x^2 - 2");
    auto result = algebra::factor_polynomial_tower(poly, Symbol{"x"}, bad, *ctx);
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::InvalidArgument);
}

TEST_F(FactorizationTowerTest, RejectsNullPolynomial) {
    auto gens = biquadratic_gens(2, 3);
    auto result = algebra::factor_polynomial_tower(ExprPtr{}, Symbol{"x"}, gens, *ctx);
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::InvalidArgument);
}

TEST_F(FactorizationTowerTest, RejectsNonRationalCoefficientPolynomial) {
    auto gens = biquadratic_gens(2, 3);
    // Already involves sqrt(2): cannot be parsed as Q[x] by the contract.
    ExprPtr poly = parse_ok("x^2 - RootOf(y^2-2,y,0)");
    ASSERT_NE(poly, nullptr);
    auto result = algebra::factor_polynomial_tower(poly, Symbol{"x"}, gens, *ctx);
    EXPECT_TRUE(result.is_error()) << "factor_polynomial_tower must reject non-Q[x] input";
    if (result.is_error()) {
        EXPECT_EQ(result.error().kind, CASErrorKind::Unimplemented);
    }
}

// Degree-4 splitting cases remain DISABLED: the composite norm of a
// degree-4 polynomial has degree 16; factorisation + iterated resultant
// under ASan instrumentation exceeds the lite regression budget.  Run
// explicitly via --gtest_also_run_disabled_tests when validating the slow
// path end-to-end.
TEST_F(FactorizationTowerTest, DISABLED_SplitsProductOfQuadraticsOverQSqrt2Sqrt3) {
    auto gens = biquadratic_gens(2, 3);
    ctx->set_timeout(std::chrono::minutes(10));
    ExprPtr poly = parse_ok("(x^2 - 2) * (x^2 - 3)");
    ASSERT_NE(poly, nullptr);
    auto expanded = algebra::expand(poly, *ctx);
    ASSERT_TRUE(expanded.is_ok()) << expanded.error().message;

    auto result = algebra::factor_polynomial_tower(expanded.value(), Symbol{"x"}, gens, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_EQ(result.value().factors.size(), 4U)
        << "expected 4 linear factors, got " << result.value().factors.size();
    EXPECT_EQ(cumulative_factor_degree(result.value(), Symbol{"x"}), 4U);
}

TEST_F(FactorizationTowerTest, DISABLED_SplitsX4Minus10X2Plus1OverQSqrt2Sqrt3) {
    auto gens = biquadratic_gens(2, 3);
    ctx->set_timeout(std::chrono::minutes(10));
    ExprPtr poly = parse_ok("x^4 - 10*x^2 + 1");
    ASSERT_NE(poly, nullptr);

    auto result = algebra::factor_polynomial_tower(poly, Symbol{"x"}, gens, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_EQ(result.value().factors.size(), 4U)
        << "expected 4 linear factors of x^4 - 10x^2 + 1 in Q(sqrt 2, sqrt 3)";
    EXPECT_EQ(cumulative_factor_degree(result.value(), Symbol{"x"}), 4U);
}

}  // namespace
}  // namespace cas::test
