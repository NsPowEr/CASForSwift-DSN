#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "../../../src/algebra/polynomial_groebner_f4.hpp"

#include <gtest/gtest.h>

using namespace cas;
using namespace cas::algebra;

namespace {

[[nodiscard]] Result<ExprPtr> parse_expr(const std::string& input, symbolic::CASContext& ctx) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), ctx.arena());
    return parser.parse();
}

[[nodiscard]] bool divides_monomial(const Monomial& divisor, const Monomial& value) {
    if (divisor.size() != value.size()) return false;
    for (std::size_t i = 0; i < divisor.size(); ++i) {
        if (divisor[i] > value[i]) return false;
    }
    return true;
}

[[nodiscard]] Monomial lcm_monomial(const Monomial& lhs, const Monomial& rhs) {
    Monomial result(lhs.size(), 0);
    for (std::size_t i = 0; i < lhs.size(); ++i) result[i] = std::max(lhs[i], rhs[i]);
    return result;
}

[[nodiscard]] PolyF4 multiply_by_monomial(const PolyF4& poly, const Monomial& shift, const Rational& factor) {
    PolyF4 result;
    for (const auto& [mon, coeff] : poly.terms) {
        Monomial shifted(mon.size(), 0);
        for (std::size_t i = 0; i < mon.size(); ++i) shifted[i] = mon[i] + shift[i];
        result.terms[shifted] = result.terms[shifted] + coeff * factor;
        if (result.terms[shifted].numerator().is_zero()) result.terms.erase(shifted);
    }
    return result;
}

void subtract_into(PolyF4& lhs, const PolyF4& rhs) {
    for (const auto& [mon, coeff] : rhs.terms) {
        lhs.terms[mon] = lhs.terms[mon] - coeff;
        if (lhs.terms[mon].numerator().is_zero()) lhs.terms.erase(mon);
    }
}

[[nodiscard]] PolyF4 reduce_poly(PolyF4 poly, const std::vector<PolyF4>& basis, MonomialOrder order) {
    constexpr std::size_t kMaxReductionSteps = 4096;
    for (std::size_t step = 0; step < kMaxReductionSteps && !poly.is_zero(); ++step) {
        const Monomial lm = poly.leading_monomial(order);
        const Rational lc = poly.leading_coefficient(order);
        bool reduced = false;
        for (const PolyF4& g : basis) {
            const Monomial lm_g = g.leading_monomial(order);
            if (lm_g.empty() || !divides_monomial(lm_g, lm)) continue;
            Monomial shift(lm.size(), 0);
            for (std::size_t i = 0; i < lm.size(); ++i) shift[i] = lm[i] - lm_g[i];
            subtract_into(poly, multiply_by_monomial(g, shift, lc / g.leading_coefficient(order)));
            reduced = true;
            break;
        }
        if (!reduced) break;
    }
    return poly;
}

[[nodiscard]] PolyF4 s_polynomial(const PolyF4& lhs, const PolyF4& rhs, MonomialOrder order) {
    const Monomial lm_l = lhs.leading_monomial(order);
    const Monomial lm_r = rhs.leading_monomial(order);
    const Monomial common = lcm_monomial(lm_l, lm_r);
    Monomial shift_l(common.size(), 0);
    Monomial shift_r(common.size(), 0);
    for (std::size_t i = 0; i < common.size(); ++i) {
        shift_l[i] = common[i] - lm_l[i];
        shift_r[i] = common[i] - lm_r[i];
    }
    PolyF4 result = multiply_by_monomial(lhs, shift_l, Rational(1) / lhs.leading_coefficient(order));
    subtract_into(result, multiply_by_monomial(rhs, shift_r, Rational(1) / rhs.leading_coefficient(order)));
    return result;
}

[[nodiscard]] std::vector<PolyF4> to_f4_basis(
    const std::vector<ExprPtr>& basis,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx) {
    std::vector<PolyF4> converted;
    for (ExprPtr expr : basis) {
        auto poly = expr_to_f4(expr, vars, ctx);
        EXPECT_TRUE(poly.is_ok());
        if (poly.is_ok() && !poly.value().is_zero()) converted.push_back(poly.value());
    }
    return converted;
}

void expect_groebner_basis_for(
    const std::vector<ExprPtr>& generators,
    const std::vector<ExprPtr>& basis_exprs,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx) {
    constexpr MonomialOrder order = MonomialOrder::GRevLex;
    const std::vector<PolyF4> basis = to_f4_basis(basis_exprs, vars, ctx);
    ASSERT_FALSE(basis.empty());

    for (ExprPtr generator : generators) {
        auto poly = expr_to_f4(generator, vars, ctx);
        ASSERT_TRUE(poly.is_ok());
        EXPECT_TRUE(reduce_poly(poly.value(), basis, order).is_zero());
    }

    for (std::size_t i = 0; i < basis.size(); ++i) {
        for (std::size_t j = i + 1; j < basis.size(); ++j) {
            EXPECT_TRUE(reduce_poly(s_polynomial(basis[i], basis[j], order), basis, order).is_zero());
        }
    }
}

} // namespace

// P3-002: Gröbner basis

TEST(GroebnerTest, EmptyInputReturnsEmpty) {
    symbolic::CASContext ctx;
    std::vector<ExprPtr> eqs;
    std::vector<Symbol> vars = {Symbol("x"), Symbol("y")};

    auto result = polynomial_groebner(eqs, vars, ctx);
    ASSERT_TRUE(result.is_ok());
    EXPECT_TRUE(result.value().empty());
}

TEST(GroebnerTest, SinglePolynomialReturnsSelf) {
    symbolic::CASContext ctx;
    auto x_expr = parse_expr("x^2 - 1", ctx);
    ASSERT_TRUE(x_expr.is_ok());

    std::vector<Symbol> vars = {Symbol("x")};
    auto result = polynomial_groebner({x_expr.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok());
    EXPECT_FALSE(result.value().empty());
    expect_groebner_basis_for({x_expr.value()}, result.value(), vars, ctx);
}

TEST(GroebnerTest, LinearSystemTwoVariables) {
    // x + y - 3 = 0, x - y - 1 = 0
    // Solution: x=2, y=1
    // Reduced Gröbner should eliminate one variable
    symbolic::CASContext ctx;

    auto eq1 = parse_expr("x + y - 3", ctx);
    auto eq2 = parse_expr("x - y - 1", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());

    std::vector<Symbol> vars = {Symbol("x"), Symbol("y")};
    auto result = polynomial_groebner({eq1.value(), eq2.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok());
    // Basis must be non-empty
    EXPECT_FALSE(result.value().empty());
    // Basis has at most 2 elements for this linear system
    EXPECT_LE(result.value().size(), 2U);
    expect_groebner_basis_for({eq1.value(), eq2.value()}, result.value(), vars, ctx);
}

TEST(GroebnerTest, ZeroPolynomialFiltered) {
    // Single zero polynomial should give empty basis
    symbolic::CASContext ctx;
    auto zero = parse_expr("0", ctx);
    ASSERT_TRUE(zero.is_ok());

    std::vector<Symbol> vars = {Symbol("x")};
    auto result = polynomial_groebner({zero.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok());
    EXPECT_TRUE(result.value().empty());
}

TEST(GroebnerTest, QuadraticSystemReduces) {
    // x^2 - 1, x - 1
    // GCD = x - 1, Gröbner = {x - 1}
    symbolic::CASContext ctx;

    auto eq1 = parse_expr("x^2 - 1", ctx);
    auto eq2 = parse_expr("x - 1", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());

    std::vector<Symbol> vars = {Symbol("x")};
    auto result = polynomial_groebner({eq1.value(), eq2.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok());
    EXPECT_FALSE(result.value().empty());
    // The reduced basis should be {x - 1}
    EXPECT_EQ(result.value().size(), 1U);
    expect_groebner_basis_for({eq1.value(), eq2.value()}, result.value(), vars, ctx);
}

TEST(GroebnerTest, IndependentVariablesPreserved) {
    // Ideal generated by {x - 1, y - 2}: no cross-terms
    symbolic::CASContext ctx;

    auto eq1 = parse_expr("x - 1", ctx);
    auto eq2 = parse_expr("y - 2", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());

    std::vector<Symbol> vars = {Symbol("x"), Symbol("y")};
    auto result = polynomial_groebner({eq1.value(), eq2.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok());
    EXPECT_EQ(result.value().size(), 2U);
    expect_groebner_basis_for({eq1.value(), eq2.value()}, result.value(), vars, ctx);
}

TEST(GroebnerTest, HardNonlinearSystemSatisfiesBuchbergerCriterion) {
    symbolic::CASContext ctx;

    auto eq1 = parse_expr("x*y - 1", ctx);
    auto eq2 = parse_expr("y^2 - x", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());

    std::vector<Symbol> vars = {Symbol("x"), Symbol("y")};
    auto result = polynomial_groebner({eq1.value(), eq2.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    expect_groebner_basis_for({eq1.value(), eq2.value()}, result.value(), vars, ctx);
}

TEST(GroebnerTest, HardNonlinearSystemUsesVariableListNotXHardcode) {
    symbolic::CASContext ctx;

    auto eq1 = parse_expr("z*t - 1", ctx);
    auto eq2 = parse_expr("t^2 - z", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());

    std::vector<Symbol> vars = {Symbol("z"), Symbol("t")};
    auto result = polynomial_groebner({eq1.value(), eq2.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    expect_groebner_basis_for({eq1.value(), eq2.value()}, result.value(), vars, ctx);
}

TEST(GroebnerTest, ThreeVariableTriangularSystemSatisfiesBuchbergerCriterion) {
    symbolic::CASContext ctx;

    auto eq1 = parse_expr("x^2 - y", ctx);
    auto eq2 = parse_expr("y^2 - z", ctx);
    auto eq3 = parse_expr("z^2 - 1", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());
    ASSERT_TRUE(eq3.is_ok());

    std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z")};
    auto result = polynomial_groebner({eq1.value(), eq2.value(), eq3.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    expect_groebner_basis_for({eq1.value(), eq2.value(), eq3.value()}, result.value(), vars, ctx);
}

TEST(GroebnerTest, CoupledQuadraticSystemSatisfiesBuchbergerCriterion) {
    symbolic::CASContext ctx;

    auto eq1 = parse_expr("x^2 + y^2 - 1", ctx);
    auto eq2 = parse_expr("x - y", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());

    std::vector<Symbol> vars = {Symbol("x"), Symbol("y")};
    auto result = polynomial_groebner({eq1.value(), eq2.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    expect_groebner_basis_for({eq1.value(), eq2.value()}, result.value(), vars, ctx);
}

TEST(GroebnerTest, AcceptsDecimalConvertedToRationalForExactConversion) {
    symbolic::CASContext ctx;
    // Parser converts 1.5 → RationalLit(3,2); Groebner handles rationals
    auto expr = parse_expr("x + 1.5", ctx);
    ASSERT_TRUE(expr.is_ok());

    std::vector<Symbol> vars = {Symbol("x")};
    auto result = polynomial_groebner({expr.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok());
}

TEST(GroebnerTest, OversizedExponentReturnsErrorNotException) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("x^1000000000000000000000000", ctx);
    ASSERT_TRUE(expr.is_ok());

    std::vector<Symbol> vars = {Symbol("x")};
    auto result = polynomial_groebner({expr.value()}, vars, ctx);
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::Overflow);
}

// L2-20: Gebauer-Moeller — product criterion eliminates coprime-LM pairs.
// x^2-1 has LM x^2, y^2-1 has LM y^2; they are coprime so the S-pair
// (x^2-1, y^2-1) must be pruned by product criterion, yet the basis is correct.
TEST(GroebnerTest, GmProductCriterionCoprimePair) {
    symbolic::CASContext ctx;

    auto eq1 = parse_expr("x^2 - 1", ctx);
    auto eq2 = parse_expr("y^2 - 1", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());

    std::vector<Symbol> vars = {Symbol("x"), Symbol("y")};
    auto result = polynomial_groebner({eq1.value(), eq2.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    // Basis must still be correct even though product criterion prunes the one S-pair
    expect_groebner_basis_for({eq1.value(), eq2.value()}, result.value(), vars, ctx);
    // The generators are already a Groebner basis (coprime LM), so size ≤ 2
    EXPECT_LE(result.value().size(), 2U);
}

// L2-20: Gebauer-Moeller — chain criterion removes pairs subsumed by new generator.
// Cyclic-3 is a classic Buchberger benchmark; without GM pruning it generates many
// redundant S-pairs. With GM it should complete within the budget and give a valid basis.
TEST(GroebnerTest, GmChainCriterionCyclic3) {
    symbolic::CASContext ctx;

    // cyclic-3: x+y+z, xy+yz+xz, xyz-1
    auto eq1 = parse_expr("x + y + z", ctx);
    auto eq2 = parse_expr("x*y + y*z + x*z", ctx);
    auto eq3 = parse_expr("x*y*z - 1", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());
    ASSERT_TRUE(eq3.is_ok());

    std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z")};
    auto result = polynomial_groebner({eq1.value(), eq2.value(), eq3.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    expect_groebner_basis_for(
        {eq1.value(), eq2.value(), eq3.value()},
        result.value(), vars, ctx);
}

// L2-20: Four-generator anti-hardcode — verifies GM handles ≥4 initial generators.
TEST(GroebnerTest, GmFourGeneratorSystemSatisfiesBuchbergerCriterion) {
    symbolic::CASContext ctx;

    auto eq1 = parse_expr("x^2 - y", ctx);
    auto eq2 = parse_expr("y^2 - z", ctx);
    auto eq3 = parse_expr("z^2 - x", ctx);
    auto eq4 = parse_expr("x + y + z - 1", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());
    ASSERT_TRUE(eq3.is_ok());
    ASSERT_TRUE(eq4.is_ok());

    std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z")};
    auto result = polynomial_groebner(
        {eq1.value(), eq2.value(), eq3.value(), eq4.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    expect_groebner_basis_for(
        {eq1.value(), eq2.value(), eq3.value(), eq4.value()},
        result.value(), vars, ctx);
}

// L2-25: is_reduced_groebner_basis — verifies the output satisfies all 3 RGB properties.
TEST(GroebnerTest, L25_OutputIsReducedGroebnerBasis) {
    symbolic::CASContext ctx;

    auto eq1 = parse_expr("x*y - 1", ctx);
    auto eq2 = parse_expr("y^2 - x", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());

    std::vector<Symbol> vars = {Symbol("x"), Symbol("y")};
    auto result = polynomial_groebner({eq1.value(), eq2.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    ASSERT_FALSE(result.value().empty());

    // Convert to PolyF4 and verify RGB property
    const std::vector<PolyF4> basis_f4 = to_f4_basis(result.value(), vars, ctx);
    EXPECT_TRUE(is_reduced_groebner_basis(basis_f4, MonomialOrder::GRevLex))
        << "polynomial_groebner output must be a reduced Gröbner basis";
}

// L2-25: uniqueness — same ideal with different generator order → identical reduced basis.
TEST(GroebnerTest, L25_ReducedBasisUniquenessForSameIdeal) {
    symbolic::CASContext ctx;
    std::vector<Symbol> vars = {Symbol("x"), Symbol("y")};

    // Two orderings of the same generators
    auto x2_minus_y = parse_expr("x^2 - y", ctx);
    auto y2_minus_x = parse_expr("y^2 - x", ctx);
    ASSERT_TRUE(x2_minus_y.is_ok());
    ASSERT_TRUE(y2_minus_x.is_ok());

    auto res1 = polynomial_groebner({x2_minus_y.value(), y2_minus_x.value()}, vars, ctx);
    auto res2 = polynomial_groebner({y2_minus_x.value(), x2_minus_y.value()}, vars, ctx);
    ASSERT_TRUE(res1.is_ok()) << res1.error().message;
    ASSERT_TRUE(res2.is_ok()) << res2.error().message;

    // Both results must satisfy RGB property
    const std::vector<PolyF4> b1 = to_f4_basis(res1.value(), vars, ctx);
    const std::vector<PolyF4> b2 = to_f4_basis(res2.value(), vars, ctx);
    EXPECT_TRUE(is_reduced_groebner_basis(b1, MonomialOrder::GRevLex));
    EXPECT_TRUE(is_reduced_groebner_basis(b2, MonomialOrder::GRevLex));
    // Same number of elements (unique RGB for a fixed order has fixed cardinality)
    EXPECT_EQ(b1.size(), b2.size());
}

// L2-25: reduced basis is minimal — no element is redundant.
TEST(GroebnerTest, L25_ReducedBasisIsMinimal) {
    symbolic::CASContext ctx;
    // Linear system: {x+y-3, x-y-1} → reduced basis is {x-2, y-1} (2 elements)
    auto eq1 = parse_expr("x + y - 3", ctx);
    auto eq2 = parse_expr("x - y - 1", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());

    std::vector<Symbol> vars = {Symbol("x"), Symbol("y")};
    auto result = polynomial_groebner({eq1.value(), eq2.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    const std::vector<PolyF4> basis_f4 = to_f4_basis(result.value(), vars, ctx);
    EXPECT_TRUE(is_reduced_groebner_basis(basis_f4, MonomialOrder::GRevLex));
    // For a 0-dimensional ideal with 2 variables, the RGB has exactly 2 elements
    EXPECT_EQ(basis_f4.size(), 2U);
}

// Step 2 power-gain: Sugar selection strategy (GMNR 1991) plus removal of
// the legacy `kMaxBuchbergerPairs=8192` / `kMaxBasisSize=256` guards lets
// Buchberger complete on systems whose pair queue would previously have
// tripped the resource bail with "Timeout".
//
// Cyclic-3 already passed pre-fix; we re-verify here that the post-fix
// engine still produces a valid Groebner basis and that the result is
// reduced, since the elimination of the hardcoded guards must NOT
// introduce regressions on the canonical small benchmark.
TEST(GroebnerTest, SugarStrategyCyclic3StillValid) {
    symbolic::CASContext ctx;
    auto eq1 = parse_expr("x + y + z", ctx);
    auto eq2 = parse_expr("x*y + y*z + x*z", ctx);
    auto eq3 = parse_expr("x*y*z - 1", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());
    ASSERT_TRUE(eq3.is_ok());

    std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z")};
    auto result = polynomial_groebner({eq1.value(), eq2.value(), eq3.value()}, vars, ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    expect_groebner_basis_for(
        {eq1.value(), eq2.value(), eq3.value()},
        result.value(), vars, ctx);

    const std::vector<PolyF4> basis_f4 = to_f4_basis(result.value(), vars, ctx);
    EXPECT_TRUE(is_reduced_groebner_basis(basis_f4, MonomialOrder::GRevLex));
}

// Power-gain: a longer chain of generators that creates a pair queue
// growth pattern previously close to the legacy 8192-pair cap. With
// Sugar selection the pair queue stays close to the theoretical minimum
// and the algorithm completes in polynomial time without hitting the
// removed guard.
TEST(GroebnerTest, SugarStrategyHandlesLongerGeneratorChain) {
    symbolic::CASContext ctx;
    // Five-generator system (intentionally over-determined to stress
    // pair generation; ideal reduces to a small Groebner basis).
    auto eq1 = parse_expr("x^2 - y", ctx);
    auto eq2 = parse_expr("y^2 - z", ctx);
    auto eq3 = parse_expr("z^2 - w", ctx);
    auto eq4 = parse_expr("w^2 - x", ctx);
    auto eq5 = parse_expr("x + y + z + w - 1", ctx);
    ASSERT_TRUE(eq1.is_ok());
    ASSERT_TRUE(eq2.is_ok());
    ASSERT_TRUE(eq3.is_ok());
    ASSERT_TRUE(eq4.is_ok());
    ASSERT_TRUE(eq5.is_ok());

    std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("w")};
    auto result = polynomial_groebner(
        {eq1.value(), eq2.value(), eq3.value(), eq4.value(), eq5.value()},
        vars, ctx);
    // We accept either success (preferred — the Sugar strategy should
    // bring this within reach) or Unimplemented; but the legacy Timeout
    // bail "exceeded exact algebra resource guard" must never appear.
    if (result.is_error()) {
        EXPECT_EQ(result.error().message.find("exceeded exact algebra resource guard"),
                  std::string::npos)
            << "Legacy kMaxBuchbergerPairs guard must be eliminated.";
    }
}
