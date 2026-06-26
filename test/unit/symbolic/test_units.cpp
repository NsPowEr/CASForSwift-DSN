// CAS-L3-08 — Unit conversion tests.

#include <gtest/gtest.h>

#include "cas/units.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::units;

namespace {

class UnitsTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

TEST_F(UnitsTest, LookupBaseSI) {
    auto m = lookup_unit("m");
    ASSERT_TRUE(m.has_value());
    EXPECT_EQ(m->dimensions.m, 1);
    EXPECT_EQ(m->scale_to_si, Rational(BigInt(1), BigInt(1)));
}

TEST_F(UnitsTest, LookupCentimeter) {
    auto cm = lookup_unit("cm");
    ASSERT_TRUE(cm.has_value());
    EXPECT_EQ(cm->dimensions.m, 1);
    EXPECT_EQ(cm->scale_to_si, Rational(BigInt(1), BigInt(100)));
}

TEST_F(UnitsTest, LookupDerivedHz) {
    auto hz = lookup_unit("Hz");
    ASSERT_TRUE(hz.has_value());
    EXPECT_EQ(hz->dimensions.s, -1);
}

TEST_F(UnitsTest, LookupUnknown) {
    EXPECT_FALSE(lookup_unit("foo").has_value());
}

TEST_F(UnitsTest, MakeFromCentimeterScales) {
    // 5 cm → Quantity(5/100, {m=1}) = Quantity(1/20, ...).
    auto five = ctx.arena().make<IntegerLit>(BigInt(5));
    auto q = make_quantity_from_unit(five, "cm", ctx);
    ASSERT_TRUE(q.is_ok());
    auto* qty = expr_cast<Quantity>(q.value());
    ASSERT_NE(qty, nullptr);
    EXPECT_EQ(qty->dimensions.m, 1);
    auto* rl = expr_cast<RationalLit>(qty->value);
    ASSERT_NE(rl, nullptr);
    EXPECT_EQ(rl->numerator, BigInt(1));
    EXPECT_EQ(rl->denominator, BigInt(20));
}

TEST_F(UnitsTest, ConvertMeterToCentimeter) {
    // Build Quantity(2, {m=1}) → convert to cm → Quantity(200, {m=1}).
    SIDimensions m_dim;
    m_dim.m = 1;
    auto two = ctx.arena().make<IntegerLit>(BigInt(2));
    auto qty_2m = ctx.arena().make<Quantity>(static_cast<ExprPtr>(two), m_dim);
    auto cv = convert_quantity(qty_2m, "cm", ctx);
    ASSERT_TRUE(cv.is_ok());
    auto* qty = expr_cast<Quantity>(cv.value());
    ASSERT_NE(qty, nullptr);
    auto* il = expr_cast<IntegerLit>(qty->value);
    ASSERT_NE(il, nullptr);
    EXPECT_EQ(il->value, BigInt(200));
}

TEST_F(UnitsTest, ConvertFootToMeter) {
    // 1 ft → Quantity(3048/10000, {m=1}).
    auto one = ctx.arena().make<IntegerLit>(BigInt(1));
    auto q = make_quantity_from_unit(one, "ft", ctx);
    ASSERT_TRUE(q.is_ok());
    auto* qty = expr_cast<Quantity>(q.value());
    auto* rl = expr_cast<RationalLit>(qty->value);
    ASSERT_NE(rl, nullptr);
    // 3048/10000 = 381/1250 reduced.
    EXPECT_EQ(rl->numerator, BigInt(381));
    EXPECT_EQ(rl->denominator, BigInt(1250));
}

TEST_F(UnitsTest, ConvertDimensionalMismatchFails) {
    SIDimensions m_dim;
    m_dim.m = 1;
    auto two = ctx.arena().make<IntegerLit>(BigInt(2));
    auto qty_2m = ctx.arena().make<Quantity>(static_cast<ExprPtr>(two), m_dim);
    // Try to convert meter to second.
    auto cv = convert_quantity(qty_2m, "s", ctx);
    EXPECT_TRUE(cv.is_error())
        << "Expected error on dimensional mismatch (m → s)";
}

TEST_F(UnitsTest, AntiHardcodeMultipleUnits) {
    // Verify several units roundtrip correctly.
    struct Case { const char* unit; long long num; long long den; int dim_m; int dim_kg; };
    Case cases[] = {
        {"km", 1000, 1, 1, 0},
        {"mg", 1, 1000000, 0, 1},
        {"min", 60, 1, 0, 0},  // s dim, not m/kg
    };
    for (auto& c : cases) {
        auto info = lookup_unit(c.unit);
        ASSERT_TRUE(info.has_value()) << c.unit;
        EXPECT_EQ(info->scale_to_si, Rational(BigInt(c.num), BigInt(c.den)));
        EXPECT_EQ(info->dimensions.m, c.dim_m) << c.unit;
        EXPECT_EQ(info->dimensions.kg, c.dim_kg) << c.unit;
    }
}

// ── F6.6 — algorithmic SI prefixes (not in the explicit registry) ──
TEST_F(UnitsTest, SIPrefixGigahertz) {
    auto ghz = lookup_unit("GHz");  // giga (10^9) · Hz
    ASSERT_TRUE(ghz.has_value());
    EXPECT_EQ(ghz->dimensions.s, -1);
    EXPECT_EQ(ghz->scale_to_si, Rational(BigInt(1000000000), BigInt(1)));
}

TEST_F(UnitsTest, SIPrefixMicrometer) {
    auto um = lookup_unit("um");  // micro (10^-6) · m
    ASSERT_TRUE(um.has_value());
    EXPECT_EQ(um->dimensions.m, 1);
    EXPECT_EQ(um->scale_to_si, Rational(BigInt(1), BigInt(1000000)));
}

TEST_F(UnitsTest, SIPrefixMegagramIsThousandKg) {
    auto Mg = lookup_unit("Mg");  // mega (10^6) · gram = 10^6/1000 kg = 1000 kg
    ASSERT_TRUE(Mg.has_value());
    EXPECT_EQ(Mg->dimensions.kg, 1);
    EXPECT_EQ(Mg->scale_to_si, Rational(BigInt(1000), BigInt(1)));
}

TEST_F(UnitsTest, SIPrefixDecaVsDeci) {
    auto dam = lookup_unit("dam");  // deca (10^1) · m
    ASSERT_TRUE(dam.has_value());
    EXPECT_EQ(dam->scale_to_si, Rational(BigInt(10), BigInt(1)));
    auto dm = lookup_unit("dm");    // deci (10^-1) · m
    ASSERT_TRUE(dm.has_value());
    EXPECT_EQ(dm->scale_to_si, Rational(BigInt(1), BigInt(10)));
}

TEST_F(UnitsTest, ExactMatchWinsOverPrefix) {
    // "min" must stay minute (60 s), not milli-inch.
    auto mn = lookup_unit("min");
    ASSERT_TRUE(mn.has_value());
    EXPECT_EQ(mn->dimensions.s, 1);
    EXPECT_EQ(mn->scale_to_si, Rational(BigInt(60), BigInt(1)));
}

TEST_F(UnitsTest, PrefixOnlyAppliesToSISymbols) {
    // "kfoo" — foo is not a prefixable SI symbol → no decomposition.
    EXPECT_FALSE(lookup_unit("kfoo").has_value());
}

// ── F6.6 — exact physical constants (2019 SI) ──
TEST_F(UnitsTest, PhysicalConstantSpeedOfLight) {
    auto c = make_physical_constant("speed_of_light", ctx);
    ASSERT_TRUE(c.is_ok());
    auto* qty = expr_cast<Quantity>(c.value());
    ASSERT_NE(qty, nullptr);
    EXPECT_EQ(qty->dimensions.m, 1);
    EXPECT_EQ(qty->dimensions.s, -1);
    auto* il = expr_cast<IntegerLit>(qty->value);
    ASSERT_NE(il, nullptr);
    EXPECT_EQ(il->value, BigInt(299792458));
}

TEST_F(UnitsTest, PhysicalConstantElementaryChargeExactRational) {
    // e = 1.602176634e-19 C = 1602176634 / 10^28.
    auto e = make_physical_constant("elementary_charge", ctx);
    ASSERT_TRUE(e.is_ok());
    auto* qty = expr_cast<Quantity>(e.value());
    ASSERT_NE(qty, nullptr);
    EXPECT_EQ(qty->dimensions.A, 1);
    EXPECT_EQ(qty->dimensions.s, 1);  // Coulomb = A·s
    auto* rl = expr_cast<RationalLit>(qty->value);
    ASSERT_NE(rl, nullptr);
    // 1602176634 / 10^28, reduced (gcd with 10^28: 2 → 801088317 / 5·10^27).
    BigInt ten28(1);
    for (int i = 0; i < 28; ++i) ten28 = ten28 * BigInt(10);
    EXPECT_EQ(rl->numerator * ten28, BigInt(1602176634) * rl->denominator);
}

TEST_F(UnitsTest, PhysicalConstantUnknownFails) {
    EXPECT_TRUE(make_physical_constant("unobtainium", ctx).is_error());
}

}  // namespace
