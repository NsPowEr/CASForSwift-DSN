#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace cas::algebra {
namespace {

[[nodiscard]] Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

[[nodiscard]] bool is_expr_zero(ExprPtr e) {
    if (!e) return false;
    if (const auto* il = expr_cast<IntegerLit>(e)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(e)) return rl->numerator.is_zero();
    return false;
}

[[nodiscard]] bool verify_is_root(
    ExprPtr poly, const Symbol& var, ExprPtr root, symbolic::CASContext& ctx) {
    auto subs = substitute(poly, var, root, ctx);
    if (subs.is_error()) return false;
    auto simp = ctx.simplify(subs.value());
    if (simp.is_error()) return false;
    return is_expr_zero(simp.value());
}

[[nodiscard]] bool contains_symbol_named(ExprPtr expr, const std::string& name) {
    if (!expr) return false;
    if (const auto* symbol = expr_cast<Symbol>(expr)) return symbol->name == name;
    if (const auto* unary = expr_cast<Unary>(expr)) return contains_symbol_named(unary->operand, name);
    if (const auto* binary = expr_cast<Binary>(expr)) {
        return contains_symbol_named(binary->left, name) || contains_symbol_named(binary->right, name);
    }
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        for (ExprPtr arg : call->args) {
            if (contains_symbol_named(arg, name)) return true;
        }
        return false;
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr term : sum->terms) {
            if (contains_symbol_named(term, name)) return true;
        }
        return false;
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr factor : product->factors) {
            if (contains_symbol_named(factor, name)) return true;
        }
        return false;
    }
    if (const auto* root = expr_cast<RootOf>(expr)) return contains_symbol_named(root->polynomial, name);
    if (const auto* matrix = expr_cast<Matrix>(expr)) {
        for (ExprPtr element : matrix->elements) {
            if (contains_symbol_named(element, name)) return true;
        }
    }
    return false;
}

[[nodiscard]] bool contains_imaginary_constant(ExprPtr expr) {
    if (!expr) return false;
    if (const auto* constant = expr_cast<Constant>(expr)) return constant->value == MathConstant::I;
    if (const auto* unary = expr_cast<Unary>(expr)) return contains_imaginary_constant(unary->operand);
    if (const auto* binary = expr_cast<Binary>(expr)) {
        return contains_imaginary_constant(binary->left) || contains_imaginary_constant(binary->right);
    }
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        for (ExprPtr arg : call->args) {
            if (contains_imaginary_constant(arg)) return true;
        }
        return false;
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr term : sum->terms) {
            if (contains_imaginary_constant(term)) return true;
        }
        return false;
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr factor : product->factors) {
            if (contains_imaginary_constant(factor)) return true;
        }
        return false;
    }
    if (const auto* root = expr_cast<RootOf>(expr)) return contains_imaginary_constant(root->polynomial);
    if (const auto* matrix = expr_cast<Matrix>(expr)) {
        for (ExprPtr element : matrix->elements) {
            if (contains_imaginary_constant(element)) return true;
        }
    }
    return false;
}

[[nodiscard]] bool is_imaginary_unit(ExprPtr expr) {
    const auto* constant = expr_cast<Constant>(expr);
    return constant != nullptr && constant->value == MathConstant::I;
}

[[nodiscard]] bool is_negative_imaginary_unit(ExprPtr expr) {
    const auto* unary = expr_cast<Unary>(expr);
    return unary != nullptr &&
        unary->op == UnaryOp::Neg &&
        is_imaginary_unit(unary->operand);
}

// ─── solve_polynomial ───────────────────────────────────────────────────

TEST(SolvePolynomial, ConstantNonZeroNoRoots) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("5", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto roots = solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok());
    EXPECT_EQ(roots.value().size(), 0U);
}

TEST(SolvePolynomial, LinearExact) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("2*x - 7", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto roots = solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok());
    ASSERT_EQ(roots.value().size(), 1U);
    EXPECT_TRUE(verify_is_root(expr.value(), x, roots.value()[0], ctx));
}

TEST(SolvePolynomial, LinearNegativeRoot) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("x + 3", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto roots = solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok());
    ASSERT_EQ(roots.value().size(), 1U);
    EXPECT_TRUE(verify_is_root(expr.value(), x, roots.value()[0], ctx));
}

TEST(SolvePolynomial, QuadraticTwoRationalRoots) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("x^2 - 5*x + 6", ctx.arena()); // (x-2)(x-3)
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto roots = solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok());
    ASSERT_EQ(roots.value().size(), 2U);
    for (ExprPtr root : roots.value()) {
        EXPECT_TRUE(verify_is_root(expr.value(), x, root, ctx));
    }
}

TEST(SolvePolynomial, QuadraticDoubleRoot) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("x^2 - 2*x + 1", ctx.arena()); // (x-1)^2
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto roots = solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok());
    ASSERT_GE(roots.value().size(), 1U);
    EXPECT_TRUE(verify_is_root(expr.value(), x, roots.value()[0], ctx));
}

TEST(SolvePolynomial, QuadraticIrrationalRoots) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("x^2 - 2", ctx.arena()); // roots ±sqrt(2), D=8 not perfect square
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto roots = solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok());
    ASSERT_EQ(roots.value().size(), 2U);
    for (ExprPtr root : roots.value()) {
        EXPECT_FALSE(expr_is<RootOf>(root));
    }
}

TEST(SolvePolynomial, QuadraticNegativeDiscriminant) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("x^2 + 1", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto roots = solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok());
    ASSERT_EQ(roots.value().size(), 2U);
    bool found_i = false;
    bool found_neg_i = false;
    for (ExprPtr root : roots.value()) {
        EXPECT_TRUE(verify_is_root(expr.value(), x, root, ctx));
        found_i = found_i || is_imaginary_unit(root);
        found_neg_i = found_neg_i || is_negative_imaginary_unit(root);
    }
    EXPECT_TRUE(found_i);
    EXPECT_TRUE(found_neg_i);
}

TEST(SolvePolynomial, CubicThreeRationalRoots) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("x^3 - x", ctx.arena()); // x(x-1)(x+1)
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto roots = solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok());
    ASSERT_EQ(roots.value().size(), 3U);
    for (ExprPtr root : roots.value()) {
        EXPECT_TRUE(verify_is_root(expr.value(), x, root, ctx));
    }
}

TEST(SolvePolynomial, CubicOneRationalOneIrreducibleQuadratic) {
    symbolic::CASContext ctx;
    // x^3 - 1 = (x-1)(x^2+x+1), quadratic factor has complex roots
    auto expr = parse_expr("x^3 - 1", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto roots = solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok());
    ASSERT_EQ(roots.value().size(), 3U);
    // First root should be rational (x=1)
    EXPECT_TRUE(verify_is_root(expr.value(), x, roots.value()[0], ctx));
    EXPECT_TRUE(contains_imaginary_constant(roots.value()[1]) || contains_imaginary_constant(roots.value()[2]));
    EXPECT_FALSE(contains_symbol_named(roots.value()[1], "i"));
    EXPECT_FALSE(contains_symbol_named(roots.value()[2], "i"));
}

TEST(SolvePolynomial, CubicFormulaUsesCanonicalImaginaryConstant) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("x^3 - 2", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto roots = solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok());
    ASSERT_EQ(roots.value().size(), 3U);

    bool found_complex_branch = false;
    for (ExprPtr root : roots.value()) {
        EXPECT_FALSE(contains_symbol_named(root, "i"));
        found_complex_branch = found_complex_branch || contains_imaginary_constant(root);
    }
    EXPECT_TRUE(found_complex_branch);
}

TEST(SolvePolynomial, QuarticFourRationalRoots) {
    symbolic::CASContext ctx;
    // (x-1)(x+1)(x-2)(x+2) = x^4 - 5x^2 + 4
    auto expr = parse_expr("x^4 - 5*x^2 + 4", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto roots = solve_polynomial(expr.value(), x, ctx);
    ASSERT_TRUE(roots.is_ok());
    ASSERT_EQ(roots.value().size(), 4U);
    for (ExprPtr root : roots.value()) {
        EXPECT_TRUE(verify_is_root(expr.value(), x, root, ctx));
    }
}

// ─── factor_over_integers (Kronecker) ───────────────────────────────────

TEST(FactorOverIntegers, KroneckerSophieGermain) {
    symbolic::CASContext ctx;
    // x^4 + 4 = (x^2+2x+2)(x^2-2x+2) — Sophie Germain identity
    auto expr = parse_expr("x^4 + 4", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto factorization = factor_over_integers(expr.value(), x, ctx);
    ASSERT_TRUE(factorization.is_ok());
    // Expect 2 quadratic factors (not 1 quartic residual)
    EXPECT_EQ(factorization.value().factors.size(), 2U);
    for (const auto& pf : factorization.value().factors) {
        EXPECT_EQ(pf.multiplicity, 1U);
    }
}

TEST(FactorOverIntegers, IrreducibleQuarticStaysIrreducible) {
    symbolic::CASContext ctx;
    // x^4 + 1 is irreducible over Q
    auto expr = parse_expr("x^4 + 1", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto factorization = factor_over_integers(expr.value(), x, ctx);
    ASSERT_TRUE(factorization.is_ok());
    // Should remain as one irreducible factor
    EXPECT_EQ(factorization.value().factors.size(), 1U);
}

// ─── partial_fractions (linear + quadratic) ─────────────────────────────

TEST(PartialFractions, SimpleLinear) {
    symbolic::CASContext ctx;
    auto expr = parse_expr("1/((x-1)*(x+1))", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto terms = partial_fractions(expr.value(), x, ctx);
    ASSERT_TRUE(terms.is_ok());
    EXPECT_EQ(terms.value().size(), 2U);
}

TEST(PartialFractions, LinearPlusMixedQuadratic) {
    symbolic::CASContext ctx;
    // 1/((x-1)*(x^2+1)): one linear + one irreducible quadratic factor
    auto expr = parse_expr("1/((x-1)*(x^2+1))", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto terms = partial_fractions(expr.value(), x, ctx);
    ASSERT_TRUE(terms.is_ok());
    // Should produce 2 terms: A/(x-1) and (Bx+C)/(x^2+1)
    ASSERT_EQ(terms.value().size(), 2U);

    // Verify equality: reconstruct from terms via together and compare together(original)
    std::vector<ExprPtr> term_vec = terms.value();
    ExprPtr sum_expr = ctx.arena().make<Sum>(std::move(term_vec));
    auto together_sum = together(sum_expr, ctx);
    auto together_orig = together(expr.value(), ctx);
    ASSERT_TRUE(together_sum.is_ok() && together_orig.is_ok());

    // Substitute x=2 (non-singular: 2-1=1, 4+1=5) into both and compare
    ExprPtr x_two = ctx.arena().make<IntegerLit>(BigInt(2LL));
    auto s1 = ctx.simplify(substitute(together_sum.value(), x, x_two, ctx).value());
    auto s2 = ctx.simplify(substitute(together_orig.value(), x, x_two, ctx).value());
    ASSERT_TRUE(s1.is_ok() && s2.is_ok());
    EXPECT_TRUE(structural_equal(s1.value(), s2.value()));
}

TEST(PartialFractions, PureQuadraticDenominator) {
    symbolic::CASContext ctx;
    // x/(x^2+1): already in partial fraction form (single quadratic factor)
    auto expr = parse_expr("x/(x^2+1)", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto terms = partial_fractions(expr.value(), x, ctx);
    ASSERT_TRUE(terms.is_ok());
    ASSERT_EQ(terms.value().size(), 1U);
}

TEST(PartialFractions, RepeatedLinear) {
    symbolic::CASContext ctx;
    // 1/(x*(x-1)^2)
    auto expr = parse_expr("1/(x*(x-1)^2)", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto terms = partial_fractions(expr.value(), x, ctx);
    ASSERT_TRUE(terms.is_ok());
    // Expect 3 terms: A/x + B/(x-1) + C/(x-1)^2
    ASSERT_EQ(terms.value().size(), 3U);
}

} // namespace
} // namespace cas::algebra
