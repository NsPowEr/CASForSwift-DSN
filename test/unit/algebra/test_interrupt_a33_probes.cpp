// HC-F70-A33 — interrupt cancellation probes for pure combinatorial loops.
//
// These pin the contract added 2026-07-10: an interrupt inside the Zassenhaus
// subset recombination or the F5C S-pair loop surfaces as an ERROR, never as
// a "no factor found" / truncated basis. For the recombination this is a
// soundness contract, not just responsiveness: the wang_eez caller reads an
// empty result as a PROOF of irreducibility (exhaustive search), so an
// interrupt that decayed to nullopt would forge that proof.

#include <gtest/gtest.h>

#include "cas/symbolic.hpp"
#include "../../../src/algebra/polynomial_internal.hpp"
#include "../../../src/algebra/polynomial_groebner_f5.hpp"

#include <vector>

using namespace cas;
using namespace cas::algebra;

namespace {

// f = (x²+1)(x²+2) = x⁴ + 3x² + 2, with its two "lifted" factors given
// exactly. Any modulus above the Mignotte bound works; the subset search
// verifies true divisibility over Z, so the recombination must find x²+1.
struct RecombFixture {
    IntPoly f{std::vector<BigInt>{BigInt(2), BigInt(0), BigInt(3), BigInt(0), BigInt(1)}};
    std::vector<IntPoly> lifted{
        IntPoly{std::vector<BigInt>{BigInt(1), BigInt(0), BigInt(1)}},   // x²+1
        IntPoly{std::vector<BigInt>{BigInt(2), BigInt(0), BigInt(1)}}};  // x²+2
    BigInt modulus{BigInt(625)};
};

TEST(InterruptA33Probe, RecombinationFindsFactorWithoutInterrupt) {
    symbolic::CASContext ctx;
    RecombFixture fx;
    auto r = find_factor_by_hensel_recombination(fx.f, fx.lifted, fx.modulus, 2U, &ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    ASSERT_TRUE(r.value().has_value()) << "exhaustive search must find x^2+1";
    EXPECT_EQ(r.value()->degree(), 2U);
}

TEST(InterruptA33Probe, PreInterruptedRecombinationErrorsNeverNullopt) {
    symbolic::CASContext ctx;
    ctx.interrupt();
    RecombFixture fx;
    auto r = find_factor_by_hensel_recombination(fx.f, fx.lifted, fx.modulus, 2U, &ctx);
    ASSERT_FALSE(r.is_ok())
        << "interrupt must surface as an error: nullopt would forge an "
           "irreducibility proof in the wang_eez caller";
    EXPECT_EQ(r.error().kind, CASErrorKind::Timeout);
}

TEST(InterruptA33Probe, RecombinationNullCtxStillWorks) {
    // Legacy call shape (no ctx): behaviour unchanged, no poll.
    RecombFixture fx;
    auto r = find_factor_by_hensel_recombination(fx.f, fx.lifted, fx.modulus, 2U);
    ASSERT_TRUE(r.is_ok());
    ASSERT_TRUE(r.value().has_value());
}

// Two generators sharing variables → non-empty S-pair queue → the poll at the
// top of the F5C main loop fires before any reduction happens.
[[nodiscard]] std::vector<PolyF4> f5_probe_input() {
    // x² + y² - 1  and  x² - y  over GRevLex, 2 variables.
    PolyF4 p1, p2;
    p1.terms[{2U, 0U}] = Rational(1);
    p1.terms[{0U, 2U}] = Rational(1);
    p1.terms[{0U, 0U}] = Rational(-1);
    p2.terms[{2U, 0U}] = Rational(1);
    p2.terms[{0U, 1U}] = Rational(-1);
    return {p1, p2};
}

TEST(InterruptA33Probe, PreInterruptedF5GroebnerErrors) {
    symbolic::CASContext ctx;
    ctx.interrupt();
    auto r = f5c_groebner(f5_probe_input(), MonomialOrder::GRevLex, &ctx);
    ASSERT_FALSE(r.is_ok()) << "S-pair loop must observe the interrupt";
    EXPECT_EQ(r.error().kind, CASErrorKind::Timeout);
}

TEST(InterruptA33Probe, F5GroebnerNullCtxStillWorks) {
    auto r = f5c_groebner(f5_probe_input(), MonomialOrder::GRevLex);
    ASSERT_TRUE(r.is_ok());
    EXPECT_FALSE(r.value().basis.empty());
}

}  // namespace
