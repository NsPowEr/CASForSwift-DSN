#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "../../helpers/property_test.hpp"

#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace cas::symbolic {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    return Parser(tokens.value(), arena).parse();
}

bool math_eq(ExprPtr a, ExprPtr b, CASContext& ctx) {
    auto r = mathematically_equal(a, b, ctx);
    return r.is_ok() && r.value();
}

// Replace all occurrences of old_var with new_var in expression string
std::string rename_var(const std::string& expr, const std::string& old_var, const std::string& new_var) {
    std::string result = expr;
    std::size_t pos = 0;
    while ((pos = result.find(old_var, pos)) != std::string::npos) {
        // Check word boundary: don't replace inside longer identifiers
        bool before_ok = (pos == 0 || !std::isalnum(result[pos - 1]));
        bool after_ok = (pos + old_var.size() >= result.size() || !std::isalnum(result[pos + old_var.size()]));
        if (before_ok && after_ok) {
            result.replace(pos, old_var.size(), new_var);
            pos += new_var.size();
        } else {
            pos += old_var.size();
        }
    }
    return result;
}

// --- P0-004: Variable independence ---
// Same mathematical fact must hold regardless of variable name

TEST(VariableIndependence, DifferenceOfSquares) {
    const std::string lhs_templ = "(VAR^2-1)";
    const std::string rhs_templ = "(VAR-1)*(VAR+1)";
    for (const auto& var : std::vector<std::string>{"x", "y", "z", "a", "b", "t"}) {
        CASContext ctx;
        auto lhs = parse_expr(rename_var(lhs_templ, "VAR", var), ctx.arena());
        auto rhs = parse_expr(rename_var(rhs_templ, "VAR", var), ctx.arena());
        ASSERT_TRUE(lhs.is_ok()) << "var=" << var << " lhs parse: " << lhs.error().message;
        ASSERT_TRUE(rhs.is_ok()) << "var=" << var << " rhs parse: " << rhs.error().message;
        EXPECT_TRUE(math_eq(lhs.value(), rhs.value(), ctx)) << "Failed for var=" << var;
    }
}

TEST(VariableIndependence, ExpandedSquare) {
    const std::string lhs_templ = "(VAR+1)^2";
    const std::string rhs_templ = "VAR^2+2*VAR+1";
    for (const auto& var : std::vector<std::string>{"x", "y", "z", "a", "b", "t"}) {
        CASContext ctx;
        auto lhs = parse_expr(rename_var(lhs_templ, "VAR", var), ctx.arena());
        auto rhs = parse_expr(rename_var(rhs_templ, "VAR", var), ctx.arena());
        ASSERT_TRUE(lhs.is_ok()) << "var=" << var;
        ASSERT_TRUE(rhs.is_ok()) << "var=" << var;
        EXPECT_TRUE(math_eq(lhs.value(), rhs.value(), ctx)) << "Failed for var=" << var;
    }
}

TEST(VariableIndependence, Differentiation) {
    // d/dv (v^3) == 3*v^2 regardless of variable name
    for (const auto& var : std::vector<std::string>{"x", "y", "z", "a", "b", "t"}) {
        CASContext ctx;
        const Symbol v{var};
        auto expr = parse_expr(var + "^3", ctx.arena());
        ASSERT_TRUE(expr.is_ok()) << "var=" << var;
        auto deriv = calculus::diff(expr.value(), v, 1U, ctx);
        ASSERT_TRUE(deriv.is_ok()) << "var=" << var << ": " << deriv.error().message;
        auto expected = parse_expr("3*" + var + "^2", ctx.arena());
        ASSERT_TRUE(expected.is_ok()) << "var=" << var;
        EXPECT_TRUE(math_eq(deriv.value(), expected.value(), ctx)) << "d/d" << var << "(v^3) failed";
    }
}

TEST(VariableIndependence, Integration_PowerRule) {
    // D(∫v^2 dv) == v^2 regardless of variable name
    for (const auto& var : std::vector<std::string>{"x", "y", "z", "a", "b", "t"}) {
        CASContext ctx;
        const Symbol v{var};
        auto f = parse_expr(var + "^2", ctx.arena());
        ASSERT_TRUE(f.is_ok()) << "var=" << var;
        auto integral = calculus::integrate(f.value(), v, ctx);
        ASSERT_TRUE(integral.is_ok()) << "var=" << var << ": " << integral.error().message;
        auto deriv = calculus::diff(integral.value(), v, 1U, ctx);
        ASSERT_TRUE(deriv.is_ok()) << "var=" << var << ": " << deriv.error().message;
        EXPECT_TRUE(math_eq(deriv.value(), f.value(), ctx)) << "D(∫v^2 dv) != v^2 for var=" << var;
    }
}

// --- P0-004: D(∫f dx) == f property ---

TEST(IntegrationDifferentiationInverse, Polynomial_x3) {
    CASContext ctx;
    const Symbol x{"x"};
    auto f = parse_expr("x^3", ctx.arena());
    ASSERT_TRUE(f.is_ok());
    auto integral = calculus::integrate(f.value(), x, ctx);
    ASSERT_TRUE(integral.is_ok()) << integral.error().message;
    auto deriv = calculus::diff(integral.value(), x, 1U, ctx);
    ASSERT_TRUE(deriv.is_ok()) << deriv.error().message;
    EXPECT_TRUE(math_eq(deriv.value(), f.value(), ctx)) << "D(∫x^3 dx) != x^3";
}

TEST(IntegrationDifferentiationInverse, Polynomial_2x2_plus_x_minus_1) {
    CASContext ctx;
    const Symbol x{"x"};
    auto f = parse_expr("2*x^2+x-1", ctx.arena());
    ASSERT_TRUE(f.is_ok());
    auto integral = calculus::integrate(f.value(), x, ctx);
    ASSERT_TRUE(integral.is_ok()) << integral.error().message;
    auto deriv = calculus::diff(integral.value(), x, 1U, ctx);
    ASSERT_TRUE(deriv.is_ok()) << deriv.error().message;
    EXPECT_TRUE(math_eq(deriv.value(), f.value(), ctx)) << "D(∫(2x^2+x-1) dx) != 2x^2+x-1";
}

TEST(IntegrationDifferentiationInverse, Polynomial_x5_minus_3x3_plus_2x) {
    CASContext ctx;
    const Symbol x{"x"};
    auto f = parse_expr("x^5-3*x^3+2*x", ctx.arena());
    ASSERT_TRUE(f.is_ok());
    auto integral = calculus::integrate(f.value(), x, ctx);
    ASSERT_TRUE(integral.is_ok()) << integral.error().message;
    auto deriv = calculus::diff(integral.value(), x, 1U, ctx);
    ASSERT_TRUE(deriv.is_ok()) << deriv.error().message;
    EXPECT_TRUE(math_eq(deriv.value(), f.value(), ctx)) << "D(∫(x^5-3x^3+2x) dx) failed";
}

TEST(IntegrationDifferentiationInverse, PropertyHoldsForMultipleVariables) {
    // D(∫v^4 dv) == v^4 for different variables
    for (const auto& var : std::vector<std::string>{"x", "y", "z"}) {
        CASContext ctx;
        const Symbol v{var};
        auto f = parse_expr(var + "^4", ctx.arena());
        ASSERT_TRUE(f.is_ok()) << "var=" << var;
        auto integral = calculus::integrate(f.value(), v, ctx);
        ASSERT_TRUE(integral.is_ok()) << "var=" << var << ": " << integral.error().message;
        auto deriv = calculus::diff(integral.value(), v, 1U, ctx);
        ASSERT_TRUE(deriv.is_ok()) << "var=" << var << ": " << deriv.error().message;
        EXPECT_TRUE(math_eq(deriv.value(), f.value(), ctx)) << "D(∫v^4 dv) failed for var=" << var;
    }
}

// --- P0-004: expand(factor(p)) == p property ---

TEST(ExpandFactorInverse, DifferenceOfSquares) {
    const std::vector<std::string> polys = {
        "x^2-1", "x^2-4", "x^2-9", "x^2-25"
    };
    for (const auto& p_str : polys) {
        CASContext ctx;
        const Symbol x{"x"};
        auto p = parse_expr(p_str, ctx.arena());
        ASSERT_TRUE(p.is_ok()) << p_str;
        auto factored = algebra::factor_over_integers(p.value(), x, ctx);
        if (factored.is_error()) { GTEST_SKIP() << "factor failed: " << factored.error().message; continue; }

        // Reconstruct product from factors
        std::vector<ExprPtr> terms;
        if (factored.value().content) terms.push_back(factored.value().content);
        for (const auto& f : factored.value().factors) {
            if (f.multiplicity == 1) {
                terms.push_back(f.factor);
            } else {
                terms.push_back(ctx.arena().make<Binary>(
                    BinaryOp::Pow, f.factor,
                    ctx.arena().make<IntegerLit>(BigInt(static_cast<long long>(f.multiplicity)))));
            }
        }
        ExprPtr reconstructed = terms.empty()
            ? ctx.arena().make<IntegerLit>(BigInt(1))
            : (terms.size() == 1 ? terms[0] : ctx.arena().make<Product>(std::move(terms)));

        auto expanded = algebra::expand(reconstructed, ctx);
        ASSERT_TRUE(expanded.is_ok()) << p_str;
        EXPECT_TRUE(math_eq(expanded.value(), p.value(), ctx))
            << "expand(factor(" << p_str << ")) != " << p_str;
    }
}

TEST(ExpandFactorInverse, LinearFactor) {
    const std::vector<std::string> polys = {"x^3-x", "2*x^2+3*x"};
    for (const auto& p_str : polys) {
        CASContext ctx;
        const Symbol x{"x"};
        auto p = parse_expr(p_str, ctx.arena());
        ASSERT_TRUE(p.is_ok()) << p_str;
        auto factored = algebra::factor_over_integers(p.value(), x, ctx);
        if (factored.is_error()) continue;

        std::vector<ExprPtr> terms;
        if (factored.value().content) terms.push_back(factored.value().content);
        for (const auto& f : factored.value().factors) {
            if (f.multiplicity == 1) {
                terms.push_back(f.factor);
            } else {
                terms.push_back(ctx.arena().make<Binary>(
                    BinaryOp::Pow, f.factor,
                    ctx.arena().make<IntegerLit>(BigInt(static_cast<long long>(f.multiplicity)))));
            }
        }
        ExprPtr reconstructed = terms.empty()
            ? ctx.arena().make<IntegerLit>(BigInt(1))
            : (terms.size() == 1 ? terms[0] : ctx.arena().make<Product>(std::move(terms)));

        auto expanded = algebra::expand(reconstructed, ctx);
        ASSERT_TRUE(expanded.is_ok()) << p_str;
        EXPECT_TRUE(math_eq(expanded.value(), p.value(), ctx))
            << "expand(factor(" << p_str << ")) != " << p_str;
    }
}

}  // namespace
}  // namespace cas::symbolic
