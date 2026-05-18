#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <cmath>

using namespace cas;
using namespace cas::algebra;

namespace {

class FsolveTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto toks = Lexer(s).tokenize();
        EXPECT_TRUE(toks.is_ok()) << s;
        Parser p(toks.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }

    // Extract roots from the (n×1) Matrix result
    [[nodiscard]] std::vector<double> extract_roots(ExprPtr result) {
        std::vector<double> roots;
        const auto* m = expr_cast<Matrix>(result);
        if (!m) return roots;
        for (ExprPtr e : m->elements) {
            if (const auto* d = expr_cast<DecimalLit>(e)) {
                roots.push_back(d->to_double());
            } else if (const auto* i = expr_cast<IntegerLit>(e)) {
                roots.push_back(i->value.to_double());
            } else if (const auto* r = expr_cast<RationalLit>(e)) {
                roots.push_back(r->numerator.to_double() / r->denominator.to_double());
            }
        }
        return roots;
    }
};

} // namespace

// L2-06 core test: sin(x) = x/2 has three roots in [-pi, pi]: -1.8955, 0, 1.8955
TEST_F(FsolveTest, SinXEqualsHalfXThreeRoots) {
    // f(x) = sin(x) - x/2
    ExprPtr f = parse("sin(x) - x/2");
    Symbol x("x");

    auto res = fsolve(f, x, ctx, -5.0, 5.0);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    auto roots = extract_roots(res.value());
    ASSERT_GE(roots.size(), 3U) << "Expected ≥3 roots for sin(x)=x/2";

    // Sort and check: roots ≈ -1.8955, 0, 1.8955
    std::sort(roots.begin(), roots.end());
    bool found_zero  = false;
    bool found_pos   = false;
    bool found_neg   = false;
    for (double r : roots) {
        if (std::abs(r) < 1e-4) found_zero = true;
        if (std::abs(r - 1.8955) < 1e-3) found_pos = true;
        if (std::abs(r + 1.8955) < 1e-3) found_neg = true;
    }
    EXPECT_TRUE(found_zero) << "Missing root at 0";
    EXPECT_TRUE(found_pos)  << "Missing root at ≈+1.8955";
    EXPECT_TRUE(found_neg)  << "Missing root at ≈-1.8955";
}

// L2-06: equation form f(x)=g(x) → internally converted to f-g=0
TEST_F(FsolveTest, EquationFormSinEqualsX) {
    // Parse as equation: this gives Binary(Equal, sin(x), x/2) — parser may not produce Equal,
    // so we build it manually via subtraction form
    ExprPtr f = parse("sin(x) - x");  // sin(x)=x has only x=0
    Symbol x("x");

    auto res = fsolve(f, x, ctx, -4.0, 4.0);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    auto roots = extract_roots(res.value());
    ASSERT_GE(roots.size(), 1U);
    bool found_zero = false;
    for (double r : roots) {
        if (std::abs(r) < 1e-4) found_zero = true;
    }
    EXPECT_TRUE(found_zero) << "sin(x)=x should have x=0 as root";
}

// L2-06: polynomial equation falls back to symbolic solve
TEST_F(FsolveTest, PolynomialEquationExactSolve) {
    // x^2 - 4 = 0 → symbolic solve gives ±2 (exact)
    ExprPtr f = parse("x^2 - 4");
    Symbol x("x");

    auto res = fsolve(f, x, ctx, -10.0, 10.0);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    auto roots = extract_roots(res.value());
    ASSERT_GE(roots.size(), 2U) << "x^2-4 should have 2 roots";
    bool has_pos2 = false, has_neg2 = false;
    for (double r : roots) {
        if (std::abs(r - 2.0) < 1e-6) has_pos2 = true;
        if (std::abs(r + 2.0) < 1e-6) has_neg2 = true;
    }
    EXPECT_TRUE(has_pos2) << "Missing root +2";
    EXPECT_TRUE(has_neg2) << "Missing root -2";
}

// L2-06: exp(x) = 2 has one root: x = ln(2) ≈ 0.6931
TEST_F(FsolveTest, ExpTranscendentalOneRoot) {
    ExprPtr f = parse("exp(x) - 2");
    Symbol x("x");

    auto res = fsolve(f, x, ctx, -5.0, 5.0);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    auto roots = extract_roots(res.value());
    ASSERT_GE(roots.size(), 1U);
    bool found = false;
    for (double r : roots) {
        if (std::abs(r - std::log(2.0)) < 1e-6) found = true;
    }
    EXPECT_TRUE(found) << "Expected root at ln(2)≈0.6931";
}

// L2-06: no roots → returns empty matrix (no exception, no error)
TEST_F(FsolveTest, NoRootsReturnsEmptyMatrix) {
    // exp(x) + 1 > 0 always, so exp(x) + 1 = 0 has no real root
    ExprPtr f = parse("exp(x) + 1");
    Symbol x("x");

    auto res = fsolve(f, x, ctx, -20.0, 20.0);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    auto roots = extract_roots(res.value());
    EXPECT_TRUE(roots.empty()) << "exp(x)+1 should have no real roots";
}

// L2-06 anti-hardcode: x*sin(x) - 1 = 0 has roots at ~1.1142, ~2.7726 in [0,4]
TEST_F(FsolveTest, XTimesSinXMinusOneRoots) {
    ExprPtr f = parse("x*sin(x) - 1");
    Symbol x("x");

    auto res = fsolve(f, x, ctx, 0.0, 5.0);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    auto roots = extract_roots(res.value());
    ASSERT_GE(roots.size(), 1U) << "x*sin(x)=1 should have ≥1 root in [0,5]";
    // The smallest positive root is ≈ 1.1142
    bool found = false;
    for (double r : roots) {
        if (std::abs(r - 1.1142) < 0.01) found = true;
    }
    EXPECT_TRUE(found) << "Missing root near 1.1142 for x*sin(x)=1";
}
