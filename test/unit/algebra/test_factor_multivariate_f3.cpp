// test_factor_multivariate_f3.cpp — GTest coverage for F3.2 Wang multivariate
// factorization (factor_multivariate).
//
// INVARIANTS:
//   - ZERO toString() validation, ZERO floating-point certificates.
//   - Certify EXACTLY: reconstruct product(content, factors^mult) as a
//     MultivariatePolynomial and assert structural equality with the input.

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "cas/formatter.hpp"
#include "algebra/algebra_internal.hpp"

#include <gtest/gtest.h>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

namespace cas::algebra {
namespace {

using Monomial  = std::vector<unsigned int>;
using SparseMap = std::map<Monomial, BigInt>;

[[nodiscard]] std::vector<Symbol> union_vars(const std::vector<Symbol>& a,
                                             const std::vector<Symbol>& b) {
    std::vector<Symbol> r = a;
    for (const auto& s : b) {
        bool present = false;
        for (const auto& t : r) if (t.name == s.name) { present = true; break; }
        if (!present) r.push_back(s);
    }
    std::sort(r.begin(), r.end(), [](const Symbol& x, const Symbol& y){ return x.name < y.name; });
    return r;
}

[[nodiscard]] SparseMap to_canonical(const MultivariatePolynomial& p,
                                     const std::vector<Symbol>& vars) {
    SparseMap m;
    for (const auto& term : p.terms()) {
        Monomial mono(vars.size(), 0U);
        for (const auto& [sym, exp] : term.factors)
            for (std::size_t i = 0; i < vars.size(); ++i)
                if (vars[i].name == sym.name) mono[i] += exp;
        m[mono] += term.coefficient;
        if (m[mono].is_zero()) m.erase(mono);
    }
    return m;
}

[[nodiscard]] Result<MultivariatePolynomial> parse_to_mv(
        const std::string& s, symbolic::CASContext& ctx) {
    auto tokens = Lexer(s).tokenize();
    if (tokens.is_error()) return fail<MultivariatePolynomial>(tokens.error());
    Parser parser(tokens.value(), ctx.arena());
    auto expr = parser.parse();
    if (expr.is_error()) return fail<MultivariatePolynomial>(expr.error());
    auto expanded = expand(expr.value(), ctx);
    if (expanded.is_error()) return fail<MultivariatePolynomial>(expanded.error());
    return parse_multivariate_polynomial(expanded.value(), ctx);
}

[[nodiscard]] Result<MultivariatePolynomial> parse_to_mv_expr(
        ExprPtr e, symbolic::CASContext& ctx) {
    auto expanded = expand(e, ctx);
    if (expanded.is_error()) return fail<MultivariatePolynomial>(expanded.error());
    return parse_multivariate_polynomial(expanded.value(), ctx);
}

// Reconstruct the product of all returned factors (with multiplicities) times the
// content, then compare structurally to the (expanded) input.
[[nodiscard]] ::testing::AssertionResult certify_factorization(
        const std::string& input, symbolic::CASContext& ctx) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return ::testing::AssertionFailure() << "lex";
    Parser parser(tokens.value(), ctx.arena());
    auto expr = parser.parse();
    if (expr.is_error()) return ::testing::AssertionFailure() << "parse";

    auto fac = factor_multivariate(expr.value(), ctx);
    if (fac.is_error())
        return ::testing::AssertionFailure() << "factor error: " << fac.error().message;

    // Build product expression: content * prod(factor_i ^ mult_i).
    ExprPtr prod = fac.value().content;
    for (const auto& pf : fac.value().factors) {
        for (unsigned int m = 0; m < pf.multiplicity; ++m) {
            prod = ctx.arena().make<Binary>(BinaryOp::Mul, prod, pf.factor);
        }
    }
    auto prod_mv = parse_to_mv_expr(prod, ctx);
    if (prod_mv.is_error())
        return ::testing::AssertionFailure() << "product parse: " << prod_mv.error().message;

    auto in_mv = parse_to_mv(input, ctx);
    if (in_mv.is_error())
        return ::testing::AssertionFailure() << "input parse";

    auto vars = union_vars(prod_mv.value().variables(), in_mv.value().variables());
    if (to_canonical(prod_mv.value(), vars) == to_canonical(in_mv.value(), vars))
        return ::testing::AssertionSuccess();
    return ::testing::AssertionFailure()
        << "product of factors != input (factors=" << fac.value().factors.size() << ")";
}

class FactorMultivariateTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

// (a) bivariate reducible: x^2 - y^2 = (x - y)(x + y)
TEST_F(FactorMultivariateTest, DifferenceOfSquares) {
    EXPECT_TRUE(certify_factorization("x^2 - y^2", ctx));
    auto fac = factor_multivariate(
        [&]{ auto t = Lexer("x^2 - y^2").tokenize(); Parser p(t.value(), ctx.arena());
             return p.parse().value(); }(), ctx);
    ASSERT_TRUE(fac.is_ok());
    EXPECT_EQ(fac.value().factors.size(), 2U);  // (x-y),(x+y)
}

// (a) bivariate reducible: (x*y + 1)(x + y)
TEST_F(FactorMultivariateTest, BivariateProduct) {
    EXPECT_TRUE(certify_factorization("(x*y + 1)*(x + y)", ctx));
}

// (b) trivariate reducible: (x + y + z)(x - y + z)
TEST_F(FactorMultivariateTest, Trivariate) {
    EXPECT_TRUE(certify_factorization("(x + y + z)*(x - y + z)", ctx));
}

// (b) trivariate with repeated structure: (x + y)(x + z)
TEST_F(FactorMultivariateTest, TrivariateTwoLinear) {
    EXPECT_TRUE(certify_factorization("(x + y)*(x + z)", ctx));
}

// (c) non-trivial (non-constant) leading coefficient in main var x:
//     (x*y + 1)*(x + y) has lc_x(first) = y → exercises the Wang single-factor
//     non-constant LC distribution + multivariate Hensel lift.
TEST_F(FactorMultivariateTest, NonTrivialLeadingCoeff) {
    EXPECT_TRUE(certify_factorization("(x*y + 1)*(x + y)", ctx));
}

// (c') Documented limitations.  When the leading-coefficient distribution across
//      MULTIPLE factors with non-constant lc_x and interacting integer content is
//      required (Wang's leading-coefficient correction, GCL §6.6 — OPEN here),
//      the algorithm MUST NOT return a silently-wrong answer: it returns an
//      explicit Unimplemented diagnostic, or certifies if it happens to succeed.
TEST_F(FactorMultivariateTest, HardLeadingCoeffNeverSilentlyWrong) {
    for (const std::string& s : {std::string("(x*y + 1)*(x*y + 2)"),
                                 std::string("(y*x + 1)*(2*x + y)")}) {
        auto t = Lexer(s).tokenize();
        Parser p(t.value(), ctx.arena());
        auto fac = factor_multivariate(p.parse().value(), ctx);
        if (fac.is_ok()) {
            EXPECT_TRUE(certify_factorization(s, ctx)) << s;
        } else {
            EXPECT_EQ(fac.error().kind, CASErrorKind::Unimplemented) << s;
        }
    }
}

// F3.2-WANG-LC-CORRECTION closure: Wang's leading-coefficient correction
// (GCL §6.6 Algorithm 6.4) now handles points where lc-factors share images.
// Each case below was previously Unimplemented; all are now certified.
TEST_F(FactorMultivariateTest, WangLcCorrectionXyPlus1XyPlus2) {
    EXPECT_TRUE(certify_factorization("(x*y + 1)*(x*y + 2)", ctx));
    auto t = Lexer("(x*y + 1)*(x*y + 2)").tokenize();
    Parser p(t.value(), ctx.arena());
    auto fac = factor_multivariate(p.parse().value(), ctx);
    ASSERT_TRUE(fac.is_ok());
    EXPECT_EQ(fac.value().factors.size(), 2U);
}
TEST_F(FactorMultivariateTest, WangLcCorrection2Xy1_3Xy1) {
    EXPECT_TRUE(certify_factorization("(2*x*y + 1)*(3*x*y + 1)", ctx));
    auto t = Lexer("(2*x*y + 1)*(3*x*y + 1)").tokenize();
    Parser p(t.value(), ctx.arena());
    auto fac = factor_multivariate(p.parse().value(), ctx);
    ASSERT_TRUE(fac.is_ok());
    EXPECT_EQ(fac.value().factors.size(), 2U);
}
TEST_F(FactorMultivariateTest, WangLcCorrectionThreeFactors) {
    EXPECT_TRUE(certify_factorization("(x*y + 1)*(x*y + 2)*(x*y + 3)", ctx));
    auto t = Lexer("(x*y + 1)*(x*y + 2)*(x*y + 3)").tokenize();
    Parser p(t.value(), ctx.arena());
    auto fac = factor_multivariate(p.parse().value(), ctx);
    ASSERT_TRUE(fac.is_ok());
    EXPECT_EQ(fac.value().factors.size(), 3U);
}
TEST_F(FactorMultivariateTest, WangLcCorrectionMixedLc) {
    // (y*x + 1)(2*x + y) — non-trivial lc_x interaction
    EXPECT_TRUE(certify_factorization("(y*x + 1)*(2*x + y)", ctx));
}

// (d) irreducible multivariate must return itself (not crash).
TEST_F(FactorMultivariateTest, IrreducibleReturnsItself) {
    auto t = Lexer("x^2 + y^2 + 1").tokenize();
    Parser p(t.value(), ctx.arena());
    auto fac = factor_multivariate(p.parse().value(), ctx);
    ASSERT_TRUE(fac.is_ok());
    EXPECT_EQ(fac.value().factors.size(), 1U);
    EXPECT_TRUE(certify_factorization("x^2 + y^2 + 1", ctx));
}

// (a) squared factor with multiplicity: (x + y)^2
TEST_F(FactorMultivariateTest, RepeatedFactor) {
    EXPECT_TRUE(certify_factorization("(x + y)^2", ctx));
}

// content extraction: 6*x^2 - 6*y^2 = 6*(x-y)(x+y)
TEST_F(FactorMultivariateTest, IntegerContent) {
    EXPECT_TRUE(certify_factorization("6*x^2 - 6*y^2", ctx));
}

// (e) input validation: zero polynomial → error, null → error.
TEST_F(FactorMultivariateTest, ZeroPolynomialRejected) {
    auto fac = factor_multivariate(ctx.arena().make<IntegerLit>(BigInt(0)), ctx);
    EXPECT_TRUE(fac.is_error());
    auto facn = factor_multivariate(nullptr, ctx);
    EXPECT_TRUE(facn.is_error());
}

}  // namespace
}  // namespace cas::algebra
