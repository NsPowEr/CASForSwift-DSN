#include <gtest/gtest.h>

#include "cas/algebraic_number.hpp"
#include "cas/algebraic_tower.hpp"
#include "cas/rational.hpp"
#include "algebra/polynomial_internal.hpp"

using namespace cas;
using namespace cas::algebra;

namespace {

template <typename Coeff>
Coeff evaluate_poly(const std::vector<Coeff>& poly, const Coeff& x) {
    if (poly.empty()) return CoeffOps<Coeff>::zero_like(x);
    Coeff acc = CoeffOps<Coeff>::zero_like(x);
    for (std::size_t i = poly.size(); i > 0U; --i) {
        acc = acc * x + poly[i - 1U];
    }
    return acc;
}

[[nodiscard]] AlgebraicNumber make_q_alpha_zero(const AlgebraicNumber::CoeffVec& min_poly) {
    return AlgebraicNumber({Rational(BigInt(0))}, min_poly);
}

[[nodiscard]] AlgebraicNumber make_q_alpha_one(const AlgebraicNumber::CoeffVec& min_poly) {
    return AlgebraicNumber({Rational(BigInt(1))}, min_poly);
}

[[nodiscard]] AlgebraicNumber make_q_alpha_generator(const AlgebraicNumber::CoeffVec& min_poly) {
    return AlgebraicNumber({Rational(BigInt(0)), Rational(BigInt(1))}, min_poly);
}

}  // namespace

TEST(ResultantGenericTest, RationalMatchesKnownBaseline) {
    std::vector<Rational> f{
        Rational(BigInt(-2)),
        Rational(BigInt(0)),
        Rational(BigInt(1)),
    };
    std::vector<Rational> g{
        Rational(BigInt(-3)),
        Rational(BigInt(0)),
        Rational(BigInt(1)),
    };

    auto res = resultant_generic<Rational>(f, g);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_EQ(res.value(), Rational(BigInt(1)));
}

TEST(ResultantGenericTest, RejectsTwoZeroPolynomialsWithoutCoefficientWitness) {
    auto res = resultant_generic<Rational>({}, {});
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::InvalidArgument);
}

TEST(ResultantGenericTest, AlgebraicNumberResultantOverQAlphaIsNonZero) {
    const AlgebraicNumber::CoeffVec min_poly{
        Rational(BigInt(-2)),
        Rational(BigInt(0)),
        Rational(BigInt(1)),
    };
    const AlgebraicNumber alpha = make_q_alpha_generator(min_poly);

    std::vector<AlgebraicNumber> f{
        make_q_alpha_zero(min_poly),
        -alpha,
        make_q_alpha_one(min_poly),
    };
    std::vector<AlgebraicNumber> g{
        -(alpha + make_q_alpha_one(min_poly)),
        make_q_alpha_one(min_poly),
    };

    auto res = resultant_generic<AlgebraicNumber>(f, g);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    const AlgebraicNumber expected = evaluate_poly(f, alpha + make_q_alpha_one(min_poly));
    EXPECT_EQ(res.value(), expected);
    EXPECT_FALSE(res.value().is_zero());
}

TEST(ResultantGenericTest, TowerResultantMatchesSignedLinearEvaluationIdentity) {
    const AlgebraicNumber::CoeffVec min_poly_1{
        Rational(BigInt(-2)),
        Rational(BigInt(0)),
        Rational(BigInt(1)),
    };
    const AlgebraicNumber alpha = make_q_alpha_generator(min_poly_1);
    const AlgebraicNumber zero = make_q_alpha_zero(min_poly_1);
    const AlgebraicNumber one = make_q_alpha_one(min_poly_1);

    const std::vector<AlgebraicNumber> min_poly_2{
        -(alpha + one),
        zero,
        one,
    };
    const AlgebraicTowerTwoLevel beta({zero, one}, min_poly_2);
    const AlgebraicTowerTwoLevel tower_one({one}, min_poly_2);

    std::vector<AlgebraicTowerTwoLevel> f{
        -beta,
        tower_one,
    };
    std::vector<AlgebraicTowerTwoLevel> g{
        -(beta + tower_one),
        tower_one,
    };

    auto res = resultant_generic<AlgebraicTowerTwoLevel>(f, g);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    const AlgebraicTowerTwoLevel expected = -evaluate_poly(f, beta + tower_one);
    EXPECT_EQ(res.value(), expected);
    EXPECT_EQ(res.value(), -tower_one);
}

TEST(ResultantGenericTest, AntiHardcodeDegreeFiveInQAlphaMatchesSignedEvaluationAtGenerator) {
    const AlgebraicNumber::CoeffVec min_poly{
        Rational(BigInt(-2)),
        Rational(BigInt(0)),
        Rational(BigInt(1)),
    };
    const AlgebraicNumber alpha = make_q_alpha_generator(min_poly);
    const AlgebraicNumber one = make_q_alpha_one(min_poly);

    std::vector<AlgebraicNumber> f{
        one + alpha,
        one - alpha,
        AlgebraicNumber({Rational(BigInt(0))}, min_poly),
        one + alpha + alpha,
        AlgebraicNumber({Rational(BigInt(-5))}, min_poly),
        one,
    };
    std::vector<AlgebraicNumber> g{
        -alpha,
        one,
    };

    auto res = resultant_generic<AlgebraicNumber>(f, g);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    const AlgebraicNumber expected = -evaluate_poly(f, alpha);
    EXPECT_EQ(res.value(), expected);
}

TEST(ResultantGenericTest, ExprResultantRequiresContext) {
    AstArena arena;
    std::vector<ExprPtr> f{
        arena.make<IntegerLit>(-1),
        arena.make<IntegerLit>(1),
    };
    std::vector<ExprPtr> g{
        arena.make<IntegerLit>(-2),
        arena.make<IntegerLit>(1),
    };

    auto res = resultant_generic<ExprPtr>(f, g, nullptr);
    ASSERT_TRUE(res.is_error());
    EXPECT_EQ(res.error().kind, CASErrorKind::InvalidArgument);
}
