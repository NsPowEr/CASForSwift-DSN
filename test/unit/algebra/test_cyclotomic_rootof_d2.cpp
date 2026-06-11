// F7.5.A1 / HC-F75-CYCLOTOMIC-ROOTOF — coverage for the geometric
// RootOf expansion that lets `mathematically_equal` recognise the
// canonical roots of (x^n - c^n)/(x - c) as equal to
// `c · exp(2πi·m/n)`.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include "algebra/algebra_internal.hpp"
#include "cas/ast_debug.hpp"

using namespace cas;

namespace {

class CyclotomicRootOfTest : public ::testing::Test {
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
};

// CAS produces RootOf(x^4 + 2x^3 + 4x^2 + 8x + 16, x, k) for k=0..3 from
// solve(x^5-32). Each of those roots equals 2·exp(2πi·m/5) for some
// m=1..4. The geometric-RootOf expansion should make
// `mathematically_equal` accept ANY such pairing (set-level matching;
// the bipartite assignment in compare_solve_sets picks a consistent
// bijection).
TEST_F(CyclotomicRootOfTest, RootOfXFifthMinusThirtyTwoQuarticEqualsExpForm) {
    auto poly = parse("x^4 + 2*x^3 + 4*x^2 + 8*x + 16");
    for (std::size_t k = 0; k < 4; ++k) {
        ExprPtr root = ctx.arena().make<RootOf>(
            poly, x, std::optional<std::size_t>{k});
        bool any_match = false;
        for (std::size_t m = 1; m <= 4; ++m) {
            // Build 2·exp(2*pi*i*m/5) in the same canonical form Maxima
            // produces after normalise_maxima_output.
            std::string s = "2 * exp(2*pi*i*" + std::to_string(m) + "/5)";
            // Substitute 'pi' and 'i' with the canonical MathConstant
            // symbols would require a parser hook; instead build by
            // hand so we depend only on the AST.
            auto pi = ctx.arena().make<Constant>(MathConstant::Pi);
            auto img = ctx.arena().make<Constant>(MathConstant::I);
            auto two = ctx.arena().make<IntegerLit>(BigInt(2));
            auto m_lit = ctx.arena().make<IntegerLit>(
                BigInt(static_cast<long long>(m)));
            auto five = ctx.arena().make<IntegerLit>(BigInt(5));
            auto num = ctx.arena().make<Product>(
                std::vector<ExprPtr>{two, m_lit, pi, img});
            auto arg = ctx.arena().make<Binary>(BinaryOp::Div, num, five);
            auto exp_term = ctx.arena().make<FuncCall>(
                BuiltinOp::Exp, std::vector<ExprPtr>{arg});
            auto rhs = ctx.arena().make<Product>(
                std::vector<ExprPtr>{two, exp_term});
            auto eq = cas::symbolic::mathematically_equal(root, rhs, ctx);
            ASSERT_TRUE(eq.is_ok()) << "k=" << k << " m=" << m;
            if (eq.value()) { any_match = true; break; }
        }
        EXPECT_TRUE(any_match) << "RootOf index k=" << k
                               << " did not match any 5th-root form";
    }
}

TEST_F(CyclotomicRootOfTest, NonGeometricRootOfFallsThrough) {
    // x^3 + x + 1 is not (x^n - c^n)/(x - c) for any rational c. The
    // helper must return nullopt and the standard mathematically_equal
    // comparison (which would return false for a structurally different
    // RHS) must be reached. We verify nullopt indirectly: a known
    // non-equal RHS returns false (no spurious true).
    auto poly = parse("x^3 + x + 1");
    ExprPtr root = ctx.arena().make<RootOf>(
        poly, x, std::optional<std::size_t>{0});
    auto rhs = parse("0");
    auto eq = cas::symbolic::mathematically_equal(root, rhs, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_FALSE(eq.value());
}

TEST_F(CyclotomicRootOfTest, GeometricRootOfRejectsUnrelatedRhs) {
    // RootOf of the geometric quartic above is in
    // {2·ω, 2·ω^2, 2·ω^3, 2·ω^4} (ω = exp(2πi/5)), none of which equals 7.
    auto poly = parse("x^4 + 2*x^3 + 4*x^2 + 8*x + 16");
    ExprPtr root = ctx.arena().make<RootOf>(
        poly, x, std::optional<std::size_t>{0});
    auto rhs = parse("7");
    auto eq = cas::symbolic::mathematically_equal(root, rhs, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_FALSE(eq.value());
}

// Direct coverage for the geometric enumerator helper: a geometric
// RootOf must expand to exactly `degree` closed-form roots. A
// non-geometric RootOf must return nullopt. This pins the helper
// surface without depending on the rest of mathematically_equal.
TEST_F(CyclotomicRootOfTest, EnumeratorDirectGeometric) {
    auto poly = parse("x^4 + 2*x^3 + 4*x^2 + 8*x + 16");
    RootOf node(poly, x, std::optional<std::size_t>{0});
    auto enumerated = cas::algebra::enumerate_geometric_rootof(node, ctx);
    ASSERT_TRUE(enumerated.has_value());
    // The enumerator emits both the positive-angle representative
    // (m * 2π/n) and the symmetric negative-angle representative
    // (m - n) * 2π/n for each of the d cyclotomic roots, so the
    // candidate set holds 2 * d entries. The duplication is harmless
    // for set-level matching and necessary to align with Maxima's
    // canonical angle range (-π, π].
    const std::size_t d = 4U;
    EXPECT_EQ(enumerated->size(), 2U * d);
}

TEST_F(CyclotomicRootOfTest, EnumeratorDirectNonGeometric) {
    auto poly = parse("x^3 + x + 1");
    RootOf node(poly, x, std::optional<std::size_t>{0});
    auto enumerated = cas::algebra::enumerate_geometric_rootof(node, ctx);
    EXPECT_FALSE(enumerated.has_value());
}

TEST_F(CyclotomicRootOfTest, DistinctRootOfIndicesCompareUnequal) {
    // Two RootOf nodes that share the same minimal polynomial but
    // carry distinct concrete indices represent distinct roots; the
    // RootOf-specific dispatch must short-circuit the comparison to
    // false before polynomial_normal_form treats both as the same
    // opaque generator.
    auto poly = parse("x^4 + 2*x^3 + 4*x^2 + 8*x + 16");
    ExprPtr r0 = ctx.arena().make<RootOf>(
        poly, x, std::optional<std::size_t>{0});
    ExprPtr r1 = ctx.arena().make<RootOf>(
        poly, x, std::optional<std::size_t>{1});
    auto eq = cas::symbolic::mathematically_equal(r0, r1, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_FALSE(eq.value());
}

TEST_F(CyclotomicRootOfTest, EnumeratorRejectsNonMonic) {
    // 2x^2 + 4x + 8 — the leading coefficient is not 1.
    auto poly = parse("2*x^2 + 4*x + 8");
    RootOf node(poly, x, std::optional<std::size_t>{0});
    auto enumerated = cas::algebra::enumerate_geometric_rootof(node, ctx);
    EXPECT_FALSE(enumerated.has_value());
}

}  // namespace
