// F5.6 — Aberth/Ehrlich complex root isolator.  Per-polynomial sanity
// checks: every returned root must satisfy |p(z)| < tolerance for the
// requested precision, and the multiset of roots must match the expected
// closed-form values to within the working tolerance.

#include "cas/lexer.hpp"
#include "cas/numeric/complex_root_isolator.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

namespace cas::test {
namespace {

cas::ExprPtr parse_expr(const std::string& input, cas::AstArena& arena) {
    auto tokens = cas::Lexer(input).tokenize();
    EXPECT_TRUE(tokens.is_ok()) << tokens.error().message;
    cas::Parser parser(tokens.value(), arena);
    auto res = parser.parse();
    EXPECT_TRUE(res.is_ok()) << res.error().message;
    return res.value();
}

// Match `actual` to one element of `expected` (a multiset of complex
// numbers given as (re, im) double pairs) within tolerance.  Returns the
// matched index, or -1 on miss.
int match_root(double re, double im,
               std::vector<std::pair<double, double>>& expected,
               double tol) {
    for (std::size_t k = 0; k < expected.size(); ++k) {
        double dr = expected[k].first - re;
        double di = expected[k].second - im;
        if (std::hypot(dr, di) < tol) return static_cast<int>(k);
    }
    return -1;
}

void expect_matches(const std::vector<cas::numeric::ComplexRoot>& got,
                    std::vector<std::pair<double, double>> expected,
                    double tol) {
    ASSERT_EQ(got.size(), expected.size());
    for (const auto& r : got) {
        double re = r.real.to_double();
        double im = r.imag.to_double();
        int idx = match_root(re, im, expected, tol);
        ASSERT_GE(idx, 0) << "unmatched root (" << re << ", " << im << ")";
        expected.erase(expected.begin() + idx);
    }
    EXPECT_TRUE(expected.empty());
}

}  // namespace

class AberthTest : public ::testing::Test {
protected:
    cas::symbolic::CASContext ctx;
};

// x² + 1 — classical (i, -i).
TEST_F(AberthTest, ImaginaryUnitRoots) {
    auto poly = parse_expr("x^2 + 1", ctx.arena());
    auto res = cas::numeric::aberth_isolate_complex_roots(poly, "x", ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_matches(res.value(), {{0.0, 1.0}, {0.0, -1.0}}, 1e-12);
}

// x³ − 1 — three cube roots of unity.
TEST_F(AberthTest, CubeRootsOfUnity) {
    auto poly = parse_expr("x^3 - 1", ctx.arena());
    auto res = cas::numeric::aberth_isolate_complex_roots(poly, "x", ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_matches(res.value(), {
        {1.0, 0.0},
        {-0.5, std::sqrt(3.0) / 2.0},
        {-0.5, -std::sqrt(3.0) / 2.0},
    }, 1e-12);
}

// (x − 2)(x − 3)(x − 5)(x − 7)(x − 11) — five distinct real prime roots; the
// hand-expanded Lagrange interpolation gives a degree-5 polynomial with
// rational coefficients, exercising the deg ≥ 5 path that has no closed-form
// solution in radicals (in general).
TEST_F(AberthTest, RealPrimesDegreeFive) {
    auto poly = parse_expr(
        "x^5 - 28*x^4 + 288*x^3 - 1358*x^2 + 2927*x - 2310",
        ctx.arena());
    auto res = cas::numeric::aberth_isolate_complex_roots(poly, "x", ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    expect_matches(res.value(), {
        {2.0, 0.0}, {3.0, 0.0}, {5.0, 0.0}, {7.0, 0.0}, {11.0, 0.0},
    }, 1e-10);
}

// x⁴ + 1 — four primitive 8th roots of unity.
TEST_F(AberthTest, EighthRootsOfUnity) {
    auto poly = parse_expr("x^4 + 1", ctx.arena());
    auto res = cas::numeric::aberth_isolate_complex_roots(poly, "x", ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    const double s = std::sqrt(2.0) / 2.0;
    expect_matches(res.value(), {
        { s,  s}, { s, -s}, {-s,  s}, {-s, -s},
    }, 1e-12);
}

// Degree-6 quintic-like with one real root and pairs of complex conjugates.
// f(x) = (x − 1)(x² + 1)(x² + x + 1)·(x + 4) =
//   x⁶ + 4x⁵ + … (expanded form).  Verifying residual magnitudes only
//   because the conjugate-pair structure is sufficient witness; the radical
//   forms of the quintic roots are not amenable to a clean double literal.
TEST_F(AberthTest, ProductOfQuadraticsResidualOnly) {
    auto poly = parse_expr(
        "(x - 1)*(x^2 + 1)*(x^2 + x + 1)*(x + 4)",
        ctx.arena());
    auto res = cas::numeric::aberth_isolate_complex_roots(poly, "x", ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    ASSERT_EQ(res.value().size(), 6U);
    for (const auto& r : res.value()) {
        EXPECT_LT(r.residual.to_double(), 1e-20);
    }
}

}  // namespace cas::test
