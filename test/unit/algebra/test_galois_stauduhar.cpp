// A6 Brick 3c — first-layer Stauduhar descent, cross-validated against the
// INDEPENDENT deg-6/7 resolvent pipeline (galois_group, whose answers are
// themselves certified against the generated transitive lattice):
//
//     pipeline says S_n / A_n  ⟺  descent certifies G_f = ambient;
//     pipeline says a proper group ⟺ descent certifies a proper
//     first-layer containment (with the Frobenius-membership tripwire).
//
// Degree-8 smoke: Φ₁₆ (G ≅ (Z/16)^*, order 8, imprimitive — must leave the
// ambient) and the Osada trinomial x⁸−x−1 (G = S₈ — must certify the
// ambient by exhausting every maximal transitive candidate). These are the
// first certified Galois statements this engine makes at degree 8.

#include <gtest/gtest.h>

#include <cstddef>
#include <initializer_list>
#include <string>
#include <vector>

#include "../../../src/algebra/galois_stauduhar_internal.hpp"
#include "../../../src/algebra/polynomial_internal.hpp"
#include "cas/bigint.hpp"
#include "cas/error.hpp"
#include "cas/galois.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::algebra;
using namespace cas::algebra::galois_stauduhar;

namespace {

class GaloisStauduharTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};

    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }

    [[nodiscard]] static IntPoly make_poly(
        std::initializer_list<long long> coeffs) {
        std::vector<BigInt> c;
        c.reserve(coeffs.size());
        for (const long long v : coeffs) c.emplace_back(v);
        IntPoly p(std::move(c));
        p.normalize([](const BigInt& v) { return v.is_zero(); });
        return p;
    }

    // Cross-validation for degrees 6/7: the independent pipeline names the
    // group; the descent must agree on "ambient vs proper subgroup".
    void cross_check(const IntPoly& f, const std::string& expr) {
        auto truth = galois_group(parse(expr), x, ctx);
        ASSERT_TRUE(truth.is_ok()) << expr << ": " << truth.error().message;
        auto fl = stauduhar_first_layer(f, ctx);
        ASSERT_TRUE(fl.is_ok()) << expr << ": " << fl.error().message;
        const std::string& name = truth.value();
        const bool truth_ambient =
            name == fl.value().ambient_name;  // "S6"/"A6"/"S7"/"A7"
        EXPECT_EQ(fl.value().is_ambient, truth_ambient)
            << expr << " → pipeline " << name << ", descent ambient="
            << fl.value().is_ambient << " (" << fl.value().ambient_name
            << ")";
        if (!fl.value().is_ambient) {
            ASSERT_TRUE(fl.value().subgroup.has_value());
            EXPECT_LT(fl.value().subgroup->order(),
                      fl.value().disc_square
                          ? permgrp::factorial_u64(f.degree()) / 2U
                          : permgrp::factorial_u64(f.degree()));
        }
    }
};

TEST_F(GaloisStauduharTest, CrossCheckDegree6) {
    // Osada S6 (ambient), Φ9 → C6 and x^6−2 → D6 (proper, imprimitive).
    cross_check(make_poly({-1, -1, 0, 0, 0, 0, 1}), "x^6 - x - 1");
    cross_check(make_poly({1, 0, 0, 1, 0, 0, 1}), "x^6 + x^3 + 1");
    cross_check(make_poly({-2, 0, 0, 0, 0, 0, 1}), "x^6 - 2");
}

TEST_F(GaloisStauduharTest, CrossCheckDegree7) {
    // Osada S7 (ambient) and the Gaussian-period C7 poly for Q(ζ29)
    // (proper: C7 ⊆ F21 below the alternating ambient).
    cross_check(make_poly({-1, -1, 0, 0, 0, 0, 0, 1}), "x^7 - x - 1");
    cross_check(
        make_poly({1, -9, 14, 28, -7, -12, 1, 1}),
        "x^7 + x^6 - 12*x^5 - 7*x^4 + 28*x^3 + 14*x^2 - 9*x + 1");
}

TEST_F(GaloisStauduharTest, Degree8CyclotomicPhi16Descends) {
    // Φ₁₆ = x⁸ + 1: G ≅ (Z/16)^* has order 8 ⇒ far below the ambient; its
    // discriminant 2²⁴ is a square, so the walk starts at A₈ and must
    // certify a proper first-layer containment.
    const IntPoly f = make_poly({1, 0, 0, 0, 0, 0, 0, 0, 1});
    auto fl = stauduhar_first_layer(f, ctx);
    ASSERT_TRUE(fl.is_ok()) << fl.error().message;
    EXPECT_TRUE(fl.value().disc_square);
    EXPECT_EQ(fl.value().ambient_name, "A8");
    ASSERT_FALSE(fl.value().is_ambient);
    ASSERT_TRUE(fl.value().subgroup.has_value());
    EXPECT_LT(fl.value().subgroup->order(), 20160U);
    EXPECT_FALSE(fl.value().provenance.empty());
}

TEST_F(GaloisStauduharTest, Degree8OsadaTrinomialIsS8) {
    // x⁸ − x − 1: S₈ by Osada's theorem — every maximal transitive
    // candidate must be exhausted with a certified NotContained.
    const IntPoly f = make_poly({-1, -1, 0, 0, 0, 0, 0, 0, 1});
    auto fl = stauduhar_first_layer(f, ctx);
    ASSERT_TRUE(fl.is_ok()) << fl.error().message;
    EXPECT_FALSE(fl.value().disc_square);
    EXPECT_TRUE(fl.value().is_ambient);
    EXPECT_EQ(fl.value().ambient_name, "S8");
}

TEST_F(GaloisStauduharTest, StructuredFailures) {
    // Non-monic.
    auto bad = stauduhar_first_layer(make_poly({1, 0, 2}), ctx);
    ASSERT_TRUE(bad.is_error());
    EXPECT_EQ(bad.error().kind, CASErrorKind::InvalidArgument);
    // Degree out of the certified candidate range.
    auto deg3 = stauduhar_first_layer(make_poly({-2, 0, 0, 1}), ctx);
    ASSERT_TRUE(deg3.is_error());
    EXPECT_EQ(deg3.error().kind, CASErrorKind::Unimplemented);
    // Zero discriminant (not squarefree): (x−1)²·(x⁴+2) has degree 6.
    auto nsf = stauduhar_first_layer(
        make_poly({2, -4, 2, 0, 1, -2, 1}), ctx);
    ASSERT_TRUE(nsf.is_error());
    EXPECT_EQ(nsf.error().kind, CASErrorKind::InvalidArgument);
}

}  // namespace
