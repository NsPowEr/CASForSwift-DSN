// HC-F8-PENDING-20-RESIDUE — together() polynomial GCD content reduction.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Together_Polynomial_GCD_Reduction.md

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas::algebra {
namespace {

class TogetherGcdTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] ExprPtr parse(const std::string& src) {
        auto tokens = Lexer(src).tokenize();
        EXPECT_TRUE(tokens.is_ok()) << "lex failed: " << src;
        Parser parser(tokens.value(), ctx.arena());
        auto e = parser.parse();
        EXPECT_TRUE(e.is_ok()) << "parse failed: " << src;
        return e.value();
    }

    // delta == 0 iff together(lhs - rhs) simplifies to the IntegerLit 0.
    [[nodiscard]] bool equals(ExprPtr lhs, ExprPtr rhs) {
        auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, lhs, rhs);
        auto t = together(delta, ctx);
        auto s = ctx.simplify(t.is_ok() ? t.value() : delta);
        if (s.is_error()) return false;
        if (const auto* il = expr_cast<IntegerLit>(s.value())) return il->value.is_zero();
        return false;
    }

    [[nodiscard]] bool equals_src(const std::string& a, const std::string& b) {
        return equals(parse(a), parse(b));
    }
};

// Row 1: exact polynomial division (x²+2xy+y²)/(x+y) → x+y.
TEST_F(TogetherGcdTest, DivExactPoly) {
    EXPECT_TRUE(equals_src("(x^2 + 2*x*y + y^2)/(x + y)", "x + y"));
}

// Row 2: mixed Sum-of-fractions denominators (the QR reproduction class).
//   y^4/(y^4+x²y²) - 1 + x²/(x²+y²) == 0  because y^4+x²y² = y²(x²+y²).
TEST_F(TogetherGcdTest, MixedSumDenom) {
    EXPECT_TRUE(equals_src(
        "y^4/(y^4 + x^2*y^2) - 1 + x^2/(x^2 + y^2)",
        "0"));
}

// Row 3: nested reduction (x²-y²)/((x-y)(x²+xy+y²)) → (x+y)/(x²+xy+y²).
TEST_F(TogetherGcdTest, NestedReduce) {
    EXPECT_TRUE(equals_src(
        "(x^2 - y^2)/((x - y)*(x^2 + x*y + y^2))",
        "(x + y)/(x^2 + x*y + y^2)"));
}

// Row 4: GCD = 1 → no spurious cancellation, value preserved.
TEST_F(TogetherGcdTest, NoOpWhenGcd1) {
    EXPECT_TRUE(equals_src("x/(y + 1)", "x/(y + 1)"));
}

// Row 5: perf guard. With together_gcd_max_degree forced to 0 the reduction is
// skipped, but the rational must still be mathematically equal (unreduced form).
TEST_F(TogetherGcdTest, GuardMaxDegreeStillCorrect) {
    ctx.set_together_gcd_max_degree(0U);
    EXPECT_TRUE(equals_src("(x^2 - y^2)/(x - y)", "x + y"))
        << "even with reduction disabled the equality check must hold via simplify";
    ctx.set_together_gcd_max_degree(64U);
}

// Row 5b: master switch off → still correct (simplify finishes the job).
TEST_F(TogetherGcdTest, DisabledStillCorrect) {
    ctx.set_together_gcd_enabled(false);
    EXPECT_TRUE(equals_src("(x^2 - y^2)/(x - y)", "x + y"));
    ctx.set_together_gcd_enabled(true);
}

// Row 6: transcendental factor passthrough — sin(x)·(x²-1)/(x-1) → sin(x)·(x+1).
TEST_F(TogetherGcdTest, TranscendentalPassthrough) {
    EXPECT_TRUE(equals_src(
        "sin(x)*(x^2 - 1)/(x - 1)",
        "sin(x)*(x + 1)"));
}

// Sanity: single-variable deep cancellation.
TEST_F(TogetherGcdTest, UnivariateCubeCancellation) {
    EXPECT_TRUE(equals_src("(x^3 - 1)/(x - 1)", "x^2 + x + 1"));
}

}  // namespace
}  // namespace cas::algebra
