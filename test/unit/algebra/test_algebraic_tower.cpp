#include "cas/algebraic_number.hpp"
#include "cas/algebraic_tower.hpp"
#include "cas/rational.hpp"

#include <gtest/gtest.h>

#include <vector>

namespace cas {
namespace algebra {
namespace {

// Helpers for terse construction.
Rational R(int n) { return Rational(BigInt(static_cast<std::int64_t>(n))); }

// Build x^2 - n as a min_poly over Q.
std::vector<Rational> sqrt_min_poly(int n) {
    return {R(-n), R(0), R(1)};
}

// AlgebraicNumber convenience: alpha (i.e. value = x) in Q(sqrt(n)).
AlgebraicNumber make_sqrt_an(int n) {
    return AlgebraicNumber({R(0), R(1)}, sqrt_min_poly(n));
}

// AlgebraicNumber convenience: rational lift `q` into Q(sqrt(n)).
AlgebraicNumber lift_rational(Rational q, int n) {
    return AlgebraicNumber({q}, sqrt_min_poly(n));
}

// === Level-1 sanity (AlgebraicElement<Rational> = Q(α)) ====================

TEST(AlgebraicTowerTest, Level1_Sqrt2_DifferenceOfSquares) {
    // α = √2, m(x) = x² − 2.   (1 + α)(1 − α) = 1 − α² = 1 − 2 = −1.
    using E = AlgebraicElementQ;
    E one_plus_alpha ({R(1),  R(1)}, sqrt_min_poly(2));
    E one_minus_alpha({R(1), R(-1)}, sqrt_min_poly(2));
    E product = one_plus_alpha * one_minus_alpha;
    E expected({R(-1)}, sqrt_min_poly(2));
    EXPECT_EQ(product, expected);
}

TEST(AlgebraicTowerTest, Level1_Sqrt2_InverseOfAlpha) {
    // Inverse of √2 is √2 / 2 = (1/2)·α.
    using E = AlgebraicElementQ;
    E alpha({R(0), R(1)}, sqrt_min_poly(2));
    auto inv_res = alpha.inverse();
    ASSERT_TRUE(inv_res.is_ok()) << inv_res.error().message;
    E expected({R(0), Rational(BigInt(1), BigInt(2))}, sqrt_min_poly(2));
    EXPECT_EQ(inv_res.value(), expected);
    // α · α⁻¹ = 1.
    E product = alpha * inv_res.value();
    E one({R(1)}, sqrt_min_poly(2));
    EXPECT_EQ(product, one);
}

// === Level-2 tower (Q(√2, √3)) =============================================
//
// Inner extension: α₁ = √2 over Q, min_poly_1 = x² − 2.
// Outer extension: α₂ = √3 over Q(α₁), min_poly_2 = y² − 3
//   (coefficient `−3` lifted from Q into Q(α₁) as `−3 + 0·α₁`).
//
// The tower element type is AlgebraicElement<AlgebraicNumber>; its
// values and minimal polynomial use AlgebraicNumber (= Q(α₁)) as the
// coefficient ring.

namespace {

// Lift integer n into Q(α₁=√2).
AlgebraicNumber lift_int_sqrt2(int n) { return lift_rational(R(n), 2); }

// y² − 3  as a min_poly over Q(α₁=√2).
std::vector<AlgebraicNumber> outer_min_poly_y2_minus_3() {
    return {
        lift_int_sqrt2(-3),
        lift_int_sqrt2(0),
        lift_int_sqrt2(1),
    };
}

} // namespace

TEST(AlgebraicTowerTest, Level2_Sqrt2_Sqrt3_DifferenceOfSquares) {
    // E = Q(α₁=√2, α₂=√3); compute (α₁ + α₂)(α₁ − α₂) = α₁² − α₂² = 2 − 3 = −1.
    // α₁ embeds into the outer extension as the constant polynomial `α₁` (degree-0
    // in α₂); α₂ embeds as `0 + 1·y` (degree-1 in α₂).
    using Tower = AlgebraicTowerTwoLevel;

    AlgebraicNumber alpha1 = make_sqrt_an(2);                 // √2 ∈ Q(√2)
    AlgebraicNumber alpha1_neg = -alpha1;
    AlgebraicNumber zero_an = lift_int_sqrt2(0);
    AlgebraicNumber one_an  = lift_int_sqrt2(1);

    // α₁ + α₂   =   α₁ + 1·y   →   coeff vector [α₁, 1].
    Tower sum_e({alpha1, one_an}, outer_min_poly_y2_minus_3());
    // α₁ − α₂   =   α₁ − 1·y   →   coeff vector [α₁, −1].
    Tower diff_e({alpha1, -one_an}, outer_min_poly_y2_minus_3());

    Tower product = sum_e * diff_e;

    // Expected = −1 (rational) embedded in tower: [lift(−1)].
    Tower expected({lift_int_sqrt2(-1)}, outer_min_poly_y2_minus_3());
    EXPECT_EQ(product, expected);
    // Sanity: result is a "ground rational" — its value vector has size ≤ 1
    // and its single coefficient lies in Q (degree 0 in α₁).
    ASSERT_LE(product.value().size(), 1U);
    if (product.value().size() == 1U) {
        EXPECT_LE(product.value().front().value().size(), 1U);
    }
}

TEST(AlgebraicTowerTest, Level2_Sqrt2_Sqrt3_InverseExact) {
    // 1/(α₁ + α₂)  =  (α₁ − α₂) / ((α₁ + α₂)(α₁ − α₂))  =  (α₂ − α₁).
    //   ((α₁ + α₂)(α₁ − α₂) = −1.)
    using Tower = AlgebraicTowerTwoLevel;
    AlgebraicNumber alpha1 = make_sqrt_an(2);
    AlgebraicNumber one_an = lift_int_sqrt2(1);

    Tower sum_e({alpha1, one_an}, outer_min_poly_y2_minus_3());
    auto inv_res = sum_e.inverse();
    ASSERT_TRUE(inv_res.is_ok()) << inv_res.error().message;

    // Expected inverse: α₂ − α₁  =  [−α₁, 1].
    Tower expected({-alpha1, one_an}, outer_min_poly_y2_minus_3());
    EXPECT_EQ(inv_res.value(), expected);

    // Cross-check: sum_e · inv = 1.
    Tower one_tower({lift_int_sqrt2(1)}, outer_min_poly_y2_minus_3());
    EXPECT_EQ(sum_e * inv_res.value(), one_tower);
}

TEST(AlgebraicTowerTest, Level2_CoefficientsDependingOnInner) {
    // Outer min_poly with a *non-trivial* coefficient in Q(α₁):
    //   m_2(y) = y² − (3·α₁ + 1).
    // Construct β with β² = 3·√2 + 1, then verify β · β = 3·α₁ + 1.
    using Tower = AlgebraicTowerTwoLevel;
    AlgebraicNumber alpha1 = make_sqrt_an(2);
    AlgebraicNumber one_an = lift_int_sqrt2(1);

    // Build the AlgebraicNumber "3α₁ + 1" inside Q(α₁).
    AlgebraicNumber three_alpha1_plus_one(
        {R(1), R(3)}, sqrt_min_poly(2)); // value = 1 + 3·x.

    std::vector<AlgebraicNumber> mp_outer{
        -three_alpha1_plus_one,
        lift_int_sqrt2(0),
        one_an,
    };

    Tower beta({lift_int_sqrt2(0), one_an}, mp_outer);  // β = 0 + 1·y.
    Tower beta_squared = beta * beta;

    // β² mod m_2 = 3·α₁ + 1, embedded as the degree-0 element [3α₁+1].
    Tower expected({three_alpha1_plus_one}, mp_outer);
    EXPECT_EQ(beta_squared, expected);

    // Inversion test: 1/β exists (β ≠ 0) and  β · (1/β) = 1.
    auto inv_res = beta.inverse();
    ASSERT_TRUE(inv_res.is_ok()) << inv_res.error().message;
    Tower one_tower({one_an}, mp_outer);
    EXPECT_EQ(beta * inv_res.value(), one_tower);
}

// === Level-3 tower (anti-hardcode: recursion depth ≥ 3) ====================
//
// E = Q(α₁=√2, α₂=√3, α₃=√5); this exercises
// AlgebraicElement<AlgebraicElement<AlgebraicNumber>>.

TEST(AlgebraicTowerTest, Level3_Sqrt2_Sqrt3_Sqrt5_DifferenceOfSquares) {
    using Mid   = AlgebraicTowerTwoLevel;                       // Q(√2,√3)
    using Tower3 = AlgebraicElement<Mid>;                       // Q(√2,√3,√5)

    AlgebraicNumber alpha1 = make_sqrt_an(2);
    AlgebraicNumber one_an = lift_int_sqrt2(1);
    AlgebraicNumber neg_one_an = -one_an;

    // Mid-tower neutrals.
    Mid mid_one  ({one_an},      outer_min_poly_y2_minus_3());
    Mid mid_zero ({},            outer_min_poly_y2_minus_3());
    Mid mid_three({lift_int_sqrt2(3)}, outer_min_poly_y2_minus_3());
    Mid mid_neg_three({lift_int_sqrt2(-3)}, outer_min_poly_y2_minus_3());

    // α₃ = √5 in Q(√2,√3): min_poly_3(z) = z² − 5 (coefficient lifted from Q).
    Mid mid_neg_five({lift_int_sqrt2(-5)}, outer_min_poly_y2_minus_3());
    std::vector<Mid> mp3{mid_neg_five, mid_zero, mid_one};

    // α₂ ∈ Mid as a constant polynomial in z, i.e. the Mid element
    // representing √3:   value [0, 1] in the middle tower.
    Mid mid_alpha2({lift_int_sqrt2(0), one_an}, outer_min_poly_y2_minus_3());

    // Compute (α₂ + α₃)(α₂ − α₃) = α₂² − α₃² = 3 − 5 = −2.
    Tower3 sum_e ({mid_alpha2,  mid_one}, mp3);
    Tower3 diff_e({mid_alpha2, -mid_one}, mp3);
    Tower3 product = sum_e * diff_e;

    Mid mid_neg_two({lift_int_sqrt2(-2)}, outer_min_poly_y2_minus_3());
    Tower3 expected({mid_neg_two}, mp3);
    EXPECT_EQ(product, expected);
}

} // namespace
} // namespace algebra
} // namespace cas
