// L2-08 Complex Polar/Log: basic identities on a + b*i form.
//   - abs(a + b*i) = sqrt(a^2 + b^2)
//   - abs(i) = 1, abs(-i) = 1
//   - ln(i) = i*pi/2
//   - ln(-1) = i*pi (principal branch)

#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <string>

namespace cas::test {

namespace {
Result<ExprPtr> parse_expr(const std::string& s, AstArena& arena) {
    auto t = Lexer(s).tokenize();
    if (t.is_error()) return fail<ExprPtr>(t.error());
    Parser p(t.value(), arena);
    return p.parse();
}
}  // namespace

class ComplexPolarTest : public ::testing::Test {
protected:
    void SetUp() override { ctx = std::make_unique<symbolic::CASContext>(); }
    void expect_equiv(const std::string& src, const std::string& expected_src) {
        auto e = parse_expr(src, ctx->arena());
        ASSERT_TRUE(e.is_ok()) << src;
        auto exp = parse_expr(expected_src, ctx->arena());
        ASSERT_TRUE(exp.is_ok());
        auto sa = ctx->simplify(e.value());
        ASSERT_TRUE(sa.is_ok());
        auto sb = ctx->simplify(exp.value());
        ASSERT_TRUE(sb.is_ok());
        auto eq = symbolic::mathematically_equal(sa.value(), sb.value(), *ctx);
        ASSERT_TRUE(eq.is_ok());
        EXPECT_TRUE(eq.value())
            << src << " ≠ " << expected_src
            << "  (got: " << debug_print(sa.value()) << ")";
    }
    std::unique_ptr<symbolic::CASContext> ctx;
};

TEST_F(ComplexPolarTest, AbsimaginaryUnit) {
    expect_equiv("abs(i)", "1");
}

TEST_F(ComplexPolarTest, AbsNegativeimaginaryUnit) {
    expect_equiv("abs(-i)", "1");
}

TEST_F(ComplexPolarTest, AbsOnePlusi) {
    expect_equiv("abs(1 + i)", "sqrt(2)");
}

TEST_F(ComplexPolarTest, AbsThreePlusFouri) {
    expect_equiv("abs(3 + 4*i)", "5");
}

TEST_F(ComplexPolarTest, LnOfiisHalfiPi) {
    expect_equiv("ln(i)", "i*pi/2");
}

TEST_F(ComplexPolarTest, LnOfMinusOneisiPi) {
    expect_equiv("ln(-1)", "i*pi");
}

TEST_F(ComplexPolarTest, ConjOfComplexFlipsimag) {
    // already working per existing extract_complex; regression check.
    expect_equiv("conj(3 + 4*i)", "3 - 4*i");
}

// arg(z) — principal argument
TEST_F(ComplexPolarTest, ArgOfPositiveRealIsZero) {
    expect_equiv("arg(5)", "0");
}
TEST_F(ComplexPolarTest, ArgOfNegativeRealIsPi) {
    expect_equiv("arg(-3)", "pi");
}
TEST_F(ComplexPolarTest, ArgOfImaginaryUnit) {
    expect_equiv("arg(i)", "pi/2");
}
TEST_F(ComplexPolarTest, ArgOfNegativeImaginaryUnit) {
    expect_equiv("arg(-i)", "-pi/2");
}
TEST_F(ComplexPolarTest, ArgOfOnePlusiIsPiOverFour) {
    expect_equiv("arg(1 + i)", "pi/4");
}
TEST_F(ComplexPolarTest, ArgOfMinusOnePlusiIsThreePiOverFour) {
    // a < 0, b > 0:  atan(1/-1) + pi = atan(-1) + pi = -pi/4 + pi = 3*pi/4
    expect_equiv("arg(-1 + i)", "3*pi/4");
}

}  // namespace cas::test
