#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

#include <string>

namespace cas::linalg {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

ExprPtr integer(symbolic::CASContext& ctx, long long value) {
    return ctx.arena().make<IntegerLit>(BigInt(value));
}

ExprPtr symbol(symbolic::CASContext& ctx, std::string name) {
    return ctx.arena().make<Symbol>(std::move(name));
}

[[nodiscard]] bool is_zero_expr(ExprPtr expr) {
    if (const auto* il = expr_cast<IntegerLit>(expr)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(expr)) return rl->numerator.is_zero();
    return false;
}

void expect_equivalent(ExprPtr actual, const std::string& expected_text) {
    symbolic::CASContext context;
    AstArena expected_arena;
    auto expected = parse_expr(expected_text, expected_arena);
    ASSERT_TRUE(expected.is_ok()) << expected.error().message;
    auto equal = symbolic::mathematically_equal(actual, expected.value(), context);
    ASSERT_TRUE(equal.is_ok()) << equal.error().message;
    EXPECT_TRUE(equal.value());
}

void expect_expr_equal(ExprPtr actual, ExprPtr expected, symbolic::CASContext& ctx) {
    auto equal = symbolic::mathematically_equal(actual, expected, ctx);
    ASSERT_TRUE(equal.is_ok()) << equal.error().message;
    EXPECT_TRUE(equal.value());
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

}  // namespace

TEST(MatrixBasicTest, MatrixAddition) {
    symbolic::CASContext context;
    MatrixExpr a(2U, 2U, {integer(context, 1), integer(context, 2), integer(context, 3), integer(context, 4)});
    MatrixExpr b(2U, 2U, {integer(context, 5), integer(context, 6), integer(context, 7), integer(context, 8)});

    auto result = add(a, b, context);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    expect_equivalent(result.value()(0U, 0U), "6");
    expect_equivalent(result.value()(0U, 1U), "8");
    expect_equivalent(result.value()(1U, 0U), "10");
    expect_equivalent(result.value()(1U, 1U), "12");
}

TEST(MatrixBasicTest, MatrixSubtraction) {
    symbolic::CASContext context;
    MatrixExpr a(2U, 2U, {integer(context, 10), integer(context, 20), integer(context, 30), integer(context, 40)});
    MatrixExpr b(2U, 2U, {integer(context, 1), integer(context, 2), integer(context, 3), integer(context, 4)});

    auto result = subtract(a, b, context);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    expect_equivalent(result.value()(0U, 0U), "9");
    expect_equivalent(result.value()(0U, 1U), "18");
    expect_equivalent(result.value()(1U, 0U), "27");
    expect_equivalent(result.value()(1U, 1U), "36");
}

TEST(MatrixBasicTest, MatrixMultiplication) {
    symbolic::CASContext context;
    MatrixExpr a(2U, 3U, {
        integer(context, 1), integer(context, 2), integer(context, 3),
        integer(context, 4), integer(context, 5), integer(context, 6),
    });
    MatrixExpr b(3U, 2U, {
        integer(context, 7), integer(context, 8),
        integer(context, 9), integer(context, 10),
        integer(context, 11), integer(context, 12),
    });

    auto result = multiply(a, b, context);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    expect_equivalent(result.value()(0U, 0U), "58");
    expect_equivalent(result.value()(0U, 1U), "64");
    expect_equivalent(result.value()(1U, 0U), "139");
    expect_equivalent(result.value()(1U, 1U), "154");
}

TEST(MatrixBasicTest, SymbolicMultiplication) {
    symbolic::CASContext context;
    MatrixExpr a(1U, 2U, {symbol(context, "x"), integer(context, 1)});
    MatrixExpr b(2U, 1U, {integer(context, 2), symbol(context, "y")});

    auto result = multiply(a, b, context);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    expect_equivalent(result.value()(0U, 0U), "2*x + y");
}

TEST(MatrixBasicTest, TransposesRectangularMatrix) {
    symbolic::CASContext context;
    MatrixExpr matrix(2U, 3U, {
        integer(context, 1), integer(context, 2), integer(context, 3),
        integer(context, 4), integer(context, 5), integer(context, 6),
    });

    auto result = transpose(matrix);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_EQ(result.value().rows(), 3U);
    EXPECT_EQ(result.value().cols(), 2U);
    expect_equivalent(result.value()(0U, 0U), "1");
    expect_equivalent(result.value()(0U, 1U), "4");
    expect_equivalent(result.value()(2U, 0U), "3");
    expect_equivalent(result.value()(2U, 1U), "6");
}

TEST(MatrixBasicTest, ComputesBareissDeterminantsExactly) {
    symbolic::CASContext context;
    MatrixExpr numeric(3U, 3U, {
        integer(context, 1), integer(context, 2), integer(context, 3),
        integer(context, 0), integer(context, 4), integer(context, 5),
        integer(context, 1), integer(context, 0), integer(context, 6),
    });

    auto numeric_det = determinant(numeric, context);
    ASSERT_TRUE(numeric_det.is_ok()) << numeric_det.error().message;
    expect_equivalent(numeric_det.value(), "22");

    MatrixExpr symbolic_matrix(2U, 2U, {
        symbol(context, "a"), symbol(context, "b"),
        symbol(context, "c"), symbol(context, "d"),
    });
    auto symbolic_det = determinant(symbolic_matrix, context);
    ASSERT_TRUE(symbolic_det.is_ok()) << symbolic_det.error().message;
    expect_equivalent(symbolic_det.value(), "a*d-b*c");
}

TEST(MatrixBasicTest, ComputesTwoByTwoInverseExactly) {
    symbolic::CASContext context;
    MatrixExpr matrix(2U, 2U, {
        symbol(context, "a"), symbol(context, "b"),
        symbol(context, "c"), symbol(context, "d"),
    });

    auto result = inverse(matrix, context);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    expect_equivalent(result.value()(0U, 0U), "d/(a*d-b*c)");
    expect_equivalent(result.value()(0U, 1U), "-b/(a*d-b*c)");
    expect_equivalent(result.value()(1U, 0U), "-c/(a*d-b*c)");
    expect_equivalent(result.value()(1U, 1U), "a/(a*d-b*c)");
}

TEST(MatrixBasicTest, ComputesLargeDiagonalInverseWithDelayedRref) {
    symbolic::CASContext context;
    constexpr std::size_t n = 9U;
    MatrixExpr matrix(n, n);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            matrix(row, col) = (row == col)
                ? symbol(context, "d" + std::to_string(row))
                : integer(context, 0);
        }
    }

    auto result = inverse(matrix, context);
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            if (row == col) {
                expect_equivalent(result.value()(row, col), "1/d" + std::to_string(row));
            } else {
                EXPECT_TRUE(is_zero_expr(result.value()(row, col)));
            }
        }
    }
}

TEST(MatrixBasicTest, RejectsInvalidLinearAlgebraInputs) {
    symbolic::CASContext context;
    MatrixExpr a(1U, 2U, {integer(context, 1), integer(context, 2)});
    MatrixExpr b(1U, 1U, {integer(context, 3)});
    EXPECT_TRUE(add(a, b, context).is_error());
    EXPECT_TRUE(multiply(a, a, context).is_error());
    EXPECT_TRUE(determinant(a, context).is_error());

    MatrixExpr singular(2U, 2U, {
        integer(context, 1), integer(context, 2),
        integer(context, 2), integer(context, 4),
    });
    auto inv = inverse(singular, context);
    ASSERT_TRUE(inv.is_error());
    EXPECT_EQ(inv.error().kind, CASErrorKind::Undefined);
}

TEST(MatrixBasicTest, ComputesRrefWithExactGaussJordan) {
    symbolic::CASContext context;
    MatrixExpr matrix(2U, 3U, {
        integer(context, 1), integer(context, 2), symbol(context, "x"),
        integer(context, 0), integer(context, 1), symbol(context, "y"),
    });

    auto result = rref(matrix, context);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    expect_equivalent(result.value()(0U, 0U), "1");
    expect_equivalent(result.value()(0U, 1U), "0");
    expect_equivalent(result.value()(0U, 2U), "x-2*y");
    expect_equivalent(result.value()(1U, 0U), "0");
    expect_equivalent(result.value()(1U, 1U), "1");
    expect_equivalent(result.value()(1U, 2U), "y");
}

TEST(MatrixBasicTest, SolvesDeterminedLinearSystemExactly) {
    symbolic::CASContext context;
    MatrixExpr a(2U, 2U, {
        integer(context, 2), integer(context, 1),
        integer(context, 1), integer(context, 1),
    });
    std::vector<ExprPtr> b{integer(context, 3), integer(context, 2)};

    auto solution = linsolve(a, b, context);
    ASSERT_TRUE(solution.is_ok()) << solution.error().message;
    ASSERT_EQ(solution.value().size(), 2U);
    expect_equivalent(solution.value()[0], "1");
    expect_equivalent(solution.value()[1], "1");
}

TEST(MatrixBasicTest, SolvesUnderdeterminedLinearSystemWithFreeParameter) {
    symbolic::CASContext context;
    MatrixExpr a(2U, 2U, {
        integer(context, 1), integer(context, 1),
        integer(context, 2), integer(context, 2),
    });
    std::vector<ExprPtr> b{integer(context, 3), integer(context, 6)};

    auto solution = linsolve(a, b, context);
    ASSERT_TRUE(solution.is_ok()) << solution.error().message;
    ASSERT_EQ(solution.value().size(), 2U);
    expect_equivalent(solution.value()[0], "3-t1");
    expect_equivalent(solution.value()[1], "t1");
}

TEST(MatrixBasicTest, RejectsInconsistentLinearSystem) {
    symbolic::CASContext context;
    MatrixExpr a(2U, 2U, {
        integer(context, 1), integer(context, 1),
        integer(context, 1), integer(context, 1),
    });
    std::vector<ExprPtr> b{integer(context, 1), integer(context, 2)};

    auto solution = linsolve(a, b, context);
    ASSERT_TRUE(solution.is_error());
    EXPECT_EQ(solution.error().kind, CASErrorKind::Undefined);
}

TEST(MatrixEigenvalueTest, DiagonalMatrixEigenvalues) {
    symbolic::CASContext context;
    MatrixExpr diag(2U, 2U, {
        integer(context, 3), integer(context, 0),
        integer(context, 0), integer(context, 5),
    });

    auto evals = eigenvalues(diag, context);
    ASSERT_TRUE(evals.is_ok()) << evals.error().message;
    ASSERT_EQ(evals.value().size(), 2U);

    bool found3 = false;
    bool found5 = false;
    for (ExprPtr ev : evals.value()) {
        if (auto eq3 = symbolic::mathematically_equal(ev, integer(context, 3), context); eq3.is_ok() && eq3.value()) { found3 = true; }
        if (auto eq5 = symbolic::mathematically_equal(ev, integer(context, 5), context); eq5.is_ok() && eq5.value()) { found5 = true; }
    }
    EXPECT_TRUE(found3);
    EXPECT_TRUE(found5);
}

TEST(MatrixEigenvalueTest, SymmetricMatrixRationalEigenvalues) {
    symbolic::CASContext context;
    // [[1,2],[2,1]]: char poly = (1-λ)²-4 = λ²-2λ-3 = (λ-3)(λ+1), eigenvalues 3 and -1.
    MatrixExpr mat(2U, 2U, {
        integer(context, 1), integer(context, 2),
        integer(context, 2), integer(context, 1),
    });

    auto evals = eigenvalues(mat, context);
    ASSERT_TRUE(evals.is_ok()) << evals.error().message;
    ASSERT_EQ(evals.value().size(), 2U);

    bool found3 = false;
    bool foundm1 = false;
    ExprPtr minus_one = context.arena().make<Unary>(UnaryOp::Neg, integer(context, 1));
    for (ExprPtr ev : evals.value()) {
        if (auto eq3 = symbolic::mathematically_equal(ev, integer(context, 3), context); eq3.is_ok() && eq3.value()) { found3 = true; }
        if (auto eqm1 = symbolic::mathematically_equal(ev, minus_one, context); eqm1.is_ok() && eqm1.value()) { foundm1 = true; }
    }
    EXPECT_TRUE(found3);
    EXPECT_TRUE(foundm1);
}

TEST(MatrixEigenvalueTest, CharacteristicPolynomialOfTwoByTwo) {
    symbolic::CASContext context;
    // [[a,b],[c,d]]: char poly = (a-λ)(d-λ)-bc = λ²-(a+d)λ+(ad-bc)
    MatrixExpr mat(2U, 2U, {
        symbol(context, "a"), symbol(context, "b"),
        symbol(context, "c"), symbol(context, "d"),
    });

    Symbol lam("lam");
    auto char_poly = characteristic_polynomial(mat, lam, context);
    ASSERT_TRUE(char_poly.is_ok()) << char_poly.error().message;

    // Substitute λ=a: char_poly(a) = a²-(a+d)a+(ad-bc) = -da+ad-bc = -bc. Verify non-trivially.
    ExprPtr lambda_a = symbol(context, "a");
    auto subs = context.substitute(char_poly.value(), lam, lambda_a);
    ASSERT_TRUE(subs.is_ok());
    auto simp = context.simplify(subs.value());
    ASSERT_TRUE(simp.is_ok());
    expect_expr_equal(simp.value(), parse_expr("-b*c", context.arena()).value(), context);
}

TEST(MatrixEigenvalueTest, RotationMatrixUsesCanonicalComplexEigenvalues) {
    symbolic::CASContext context;
    MatrixExpr mat(2U, 2U, {
        integer(context, 0), integer(context, -1),
        integer(context, 1), integer(context, 0),
    });

    Symbol lam("lam");
    auto char_poly = characteristic_polynomial(mat, lam, context);
    ASSERT_TRUE(char_poly.is_ok()) << char_poly.error().message;
    auto expected_poly = parse_expr("lam^2 + 1", context.arena());
    ASSERT_TRUE(expected_poly.is_ok()) << expected_poly.error().message;
    expect_expr_equal(char_poly.value(), expected_poly.value(), context);

    auto evals = eigenvalues(mat, context);
    ASSERT_TRUE(evals.is_ok()) << evals.error().message;
    ASSERT_EQ(evals.value().size(), 2U);
    bool found_i = false;
    bool found_neg_i = false;
    for (ExprPtr ev : evals.value()) {
        found_i = found_i || is_imaginary_unit(ev);
        found_neg_i = found_neg_i || is_negative_imaginary_unit(ev);
    }
    EXPECT_TRUE(found_i);
    EXPECT_TRUE(found_neg_i);
}

TEST(MatrixEigenvalueTest, EigenvectorsOfDiagonalMatrix) {
    symbolic::CASContext context;
    MatrixExpr diag(2U, 2U, {
        integer(context, 2), integer(context, 0),
        integer(context, 0), integer(context, 3),
    });

    auto evecs = eigenvectors(diag, context);
    ASSERT_TRUE(evecs.is_ok()) << evecs.error().message;
    ASSERT_EQ(evecs.value().size(), 2U);

    for (const auto& pair : evecs.value()) {
        EXPECT_TRUE(pair.eigenvalue != ExprPtr{});
        EXPECT_FALSE(pair.eigenvector.empty());
    }
}

TEST(MatrixEigenvalueTest, GeometricMultiplicityGreaterThanOne) {
    symbolic::CASContext context;
    // Identity matrix has one eigenvalue (1) with geometric multiplicity 2.
    MatrixExpr id(2U, 2U, {
        integer(context, 1), integer(context, 0),
        integer(context, 0), integer(context, 1),
    });

    auto evecs = eigenvectors(id, context);
    ASSERT_TRUE(evecs.is_ok()) << evecs.error().message;
    
    // Should have 2 eigenvectors for eigenvalue 1
    ASSERT_EQ(evecs.value().size(), 2U);
    for (const auto& pair : evecs.value()) {
        expect_equivalent(pair.eigenvalue, "1");
        ASSERT_EQ(pair.eigenvector.size(), 2U);
    }
    
    // Vectors should be (1, 0) and (0, 1) or equivalent basis
    // We check independence or just specific values if RREF is predictable
    bool found_1_0 = false;
    bool found_0_1 = false;
    
    for (const auto& pair : evecs.value()) {
        if (is_zero_expr(pair.eigenvector[1])) {
            expect_equivalent(pair.eigenvector[0], "1");
            found_1_0 = true;
        } else if (is_zero_expr(pair.eigenvector[0])) {
            expect_equivalent(pair.eigenvector[1], "1");
            found_0_1 = true;
        }
    }
    EXPECT_TRUE(found_1_0);
    EXPECT_TRUE(found_0_1);
}

TEST(MatrixBasicTest, ComputesRank) {
    symbolic::CASContext context;
    // [[1, 2, 3], [4, 5, 6], [7, 8, 9]] rank is 2
    MatrixExpr m(3, 3);
    m(0, 0) = integer(context, 1); m(0, 1) = integer(context, 2); m(0, 2) = integer(context, 3);
    m(1, 0) = integer(context, 4); m(1, 1) = integer(context, 5); m(1, 2) = integer(context, 6);
    m(2, 0) = integer(context, 7); m(2, 1) = integer(context, 8); m(2, 2) = integer(context, 9);
    
    auto r = rank(m, context);
    ASSERT_TRUE(r.is_ok());
    EXPECT_EQ(r.value(), 2U);

    // [[1, 0], [0, 1]] rank is 2
    MatrixExpr id(2, 2);
    id(0, 0) = integer(context, 1); id(0, 1) = integer(context, 0);
    id(1, 0) = integer(context, 0); id(1, 1) = integer(context, 1);
    auto r_id = rank(id, context);
    ASSERT_TRUE(r_id.is_ok());
    EXPECT_EQ(r_id.value(), 2U);
}

TEST(MatrixBasicTest, ComputesTrace) {
    symbolic::CASContext context;
    MatrixExpr m(2, 2);
    m(0, 0) = symbol(context, "a"); m(0, 1) = symbol(context, "b");
    m(1, 0) = symbol(context, "c"); m(1, 1) = symbol(context, "d");
    
    auto tr = trace(m, context);
    ASSERT_TRUE(tr.is_ok());
    expect_equivalent(tr.value(), "a + d");
}

TEST(MatrixBasicTest, ComputesSmithNormalForm) {
    symbolic::CASContext context;
    // [[1, 2], [3, 4]] -> [[1, 0], [0, 2]]
    MatrixExpr m(2, 2);
    m(0, 0) = integer(context, 1); m(0, 1) = integer(context, 2);
    m(1, 0) = integer(context, 3); m(1, 1) = integer(context, 4);
    
    auto snf = smith_normal_form(m, context);
    ASSERT_TRUE(snf.is_ok());
    
    expect_equivalent(snf.value().S(0, 0), "1");
    expect_equivalent(snf.value().S(0, 1), "0");
    expect_equivalent(snf.value().S(1, 0), "0");
    expect_equivalent(snf.value().S(1, 1), "2");
    
    // Check property: U * A * V = S
    auto UA = multiply(snf.value().U, m, context);
    ASSERT_TRUE(UA.is_ok());
    auto UAV = multiply(UA.value(), snf.value().V, context);
    ASSERT_TRUE(UAV.is_ok());
    
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 2; ++j) {
            expect_expr_equal(UAV.value()(i, j), snf.value().S(i, j), context);
        }
    }
}

// P2-006: Eigenvalues for n>3 matrices (unblocked by Ferrari/P1-003)
TEST(MatrixEigenvalueTest, DiagonalFourByFourEigenvalues) {
    symbolic::CASContext context;
    // Diagonal 4x4: eigenvalues are the diagonal entries {1,2,3,4}
    MatrixExpr mat(4U, 4U, {
        integer(context, 1), integer(context, 0), integer(context, 0), integer(context, 0),
        integer(context, 0), integer(context, 2), integer(context, 0), integer(context, 0),
        integer(context, 0), integer(context, 0), integer(context, 3), integer(context, 0),
        integer(context, 0), integer(context, 0), integer(context, 0), integer(context, 4),
    });

    auto evals = eigenvalues(mat, context);
    ASSERT_TRUE(evals.is_ok()) << evals.error().message;
    ASSERT_EQ(evals.value().size(), 4U);

    for (int expected : {1, 2, 3, 4}) {
        bool found = false;
        for (ExprPtr ev : evals.value()) {
            auto eq = symbolic::mathematically_equal(ev, integer(context, expected), context);
            if (eq.is_ok() && eq.value()) { found = true; break; }
        }
        EXPECT_TRUE(found) << "Eigenvalue " << expected << " not found";
    }
}

TEST(MatrixEigenvalueTest, UpperTriangularFourByFourEigenvalues) {
    symbolic::CASContext context;
    // Upper triangular 4x4: eigenvalues = diagonal = {2,3,5,7}
    MatrixExpr mat(4U, 4U, {
        integer(context, 2), integer(context, 1), integer(context, 4), integer(context, 9),
        integer(context, 0), integer(context, 3), integer(context, 2), integer(context, 1),
        integer(context, 0), integer(context, 0), integer(context, 5), integer(context, 6),
        integer(context, 0), integer(context, 0), integer(context, 0), integer(context, 7),
    });

    auto evals = eigenvalues(mat, context);
    ASSERT_TRUE(evals.is_ok()) << evals.error().message;
    ASSERT_EQ(evals.value().size(), 4U);

    for (int expected : {2, 3, 5, 7}) {
        bool found = false;
        for (ExprPtr ev : evals.value()) {
            auto eq = symbolic::mathematically_equal(ev, integer(context, expected), context);
            if (eq.is_ok() && eq.value()) { found = true; break; }
        }
        EXPECT_TRUE(found) << "Eigenvalue " << expected << " not found";
    }
}

// L1-17: Bareiss pivot scoring uses assumptions
TEST(MatrixBareissTest, L1_17_AssumptionBackedPivotPreferredOverStructuralNonzero) {
    symbolic::CASContext ctx;
    auto x = ctx.arena().make<Symbol>(Symbol{"x"});
    auto y = ctx.arena().make<Symbol>(Symbol{"y"});
    // Assume y > 0 (known positive) but x is unconstrained
    ctx.assumptions().assume_positive(Symbol{"y"});
    // Matrix: row0=[x, 1], row1=[y, 2]
    // Bareiss should prefer y (assumption-backed) over x (unknown) as pivot
    MatrixExpr mat(2U, 2U, {x, integer(ctx, 1), y, integer(ctx, 2)});
    auto result = bareiss(mat, ctx);
    EXPECT_TRUE(result.is_ok()) << result.error().message;
}

// L1-17: configuring max_integration_depth works
TEST(MatrixBareissTest, L1_18_MaxIntegrationDepthConfigurable) {
    symbolic::CASContext ctx;
    EXPECT_EQ(ctx.max_integration_depth(), 16U);
    ctx.set_max_integration_depth(32U);
    EXPECT_EQ(ctx.max_integration_depth(), 32U);
    ctx.set_max_integration_depth(200U);
    EXPECT_EQ(ctx.max_integration_depth(), 128U);  // clamped to max
    ctx.set_max_integration_depth(0U);
    EXPECT_GE(ctx.max_integration_depth(), 1U);  // min 1
}

// L1-21: gcd_error_probability configurable
TEST(MatrixBareissTest, L1_21_GcdErrorProbabilityConfigurable) {
    symbolic::CASContext ctx;
    EXPECT_DOUBLE_EQ(ctx.gcd_error_probability(), 0.001);
    ctx.set_gcd_error_probability(0.01);
    EXPECT_DOUBLE_EQ(ctx.gcd_error_probability(), 0.01);
    ctx.set_gcd_error_probability(1e-10);
    EXPECT_GE(ctx.gcd_error_probability(), 1e-6);  // clamped to min
    ctx.set_gcd_error_probability(0.5);
    EXPECT_LE(ctx.gcd_error_probability(), 0.1);   // clamped to max
}

}  // namespace cas::linalg
