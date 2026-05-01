#include <gtest/gtest.h>
#include "cas/calculus.hpp"
#include "cas/algebra.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/formatter.hpp"
#include "cas/ast_debug.hpp"
#include "cas/numtheory.hpp"
#include "algebra/algebra_internal.hpp"
#include <chrono>
#include <iostream>

using namespace cas;
using namespace cas::calculus;
using namespace cas::algebra;
using namespace cas::linalg;
using namespace cas::symbolic;

class AcidTest : public ::testing::Test {
protected:
    CASContext ctx;
};

// TEST 1: GRUNTZ LIMIT TRAP
TEST_F(AcidTest, Test1_GruntzLimit) {
    Symbol acid_x("acid_x");
    ExprPtr exp_x = ctx.arena().make<FuncCall>("exp", std::vector<ExprPtr>{ctx.arena().make<Symbol>("acid_x")});
    ExprPtr one = ctx.arena().make<IntegerLit>(BigInt(1));
    ExprPtr x_ptr = ctx.arena().make<Symbol>("acid_x");
    ExprPtr x2_2 = ctx.arena().make<Binary>(BinaryOp::Div, 
        ctx.arena().make<Binary>(BinaryOp::Pow, x_ptr, ctx.arena().make<IntegerLit>(BigInt(2))),
        ctx.arena().make<IntegerLit>(BigInt(2)));
    
    ExprPtr num = ctx.arena().make<Binary>(BinaryOp::Sub, 
        ctx.arena().make<Binary>(BinaryOp::Sub, 
            ctx.arena().make<Binary>(BinaryOp::Sub, exp_x, one), 
            x_ptr), 
        x2_2);
    
    ExprPtr den = ctx.arena().make<Binary>(BinaryOp::Pow, x_ptr, ctx.arena().make<IntegerLit>(BigInt(3)));
    ExprPtr expr = ctx.arena().make<Binary>(BinaryOp::Div, num, den);

    auto res = limit(expr, acid_x, ctx.arena().make<IntegerLit>(BigInt(0)), LimitDirection::Both, ctx);
    ASSERT_TRUE(res.is_ok()) << (res.is_error() ? res.error().message : "");

    const auto* rat = expr_cast<RationalLit>(res.value());
    ASSERT_NE(rat, nullptr) << "Result is not a rational: " << debug_print(res.value());
    EXPECT_EQ(rat->numerator, BigInt(1));
    EXPECT_EQ(rat->denominator, BigInt(6));
}

// TEST 2: RISCH ALGORITHM
TEST_F(AcidTest, Test2_RischIntegration) {
    Symbol x("x");
    ExprPtr lnx = ctx.arena().make<FuncCall>("ln", std::vector<ExprPtr>{ctx.arena().make<Symbol>("x")});
    ExprPtr expr = ctx.arena().make<Binary>(BinaryOp::Pow, lnx, ctx.arena().make<IntegerLit>(BigInt(2)));

    auto res = integrate(expr, x, ctx);
    ASSERT_TRUE(res.is_ok()) << (res.is_error() ? res.error().message : "");
    
    std::cout << "[ ACID DEBUG ] Test 2 Result: " << debug_print(res.value()) << std::endl;

    auto d_res = diff(res.value(), x, 1U, ctx);
    ASSERT_TRUE(d_res.is_ok());
    auto simplified = ctx.simplify(d_res.value());
    
    auto mat_eq = mathematically_equal(simplified.value(), expr, ctx);
    EXPECT_TRUE(mat_eq.is_ok() && mat_eq.value());
}

// TEST 3: EULER ALIASING & TRIG IDENTITY
TEST_F(AcidTest, Test3_TrigSimplification) {
    Symbol x("x");
    ExprPtr sinx = ctx.arena().make<FuncCall>("sin", std::vector<ExprPtr>{ctx.arena().make<Symbol>("x")});
    ExprPtr cosx = ctx.arena().make<FuncCall>("cos", std::vector<ExprPtr>{ctx.arena().make<Symbol>("x")});
    
    ExprPtr s4 = ctx.arena().make<Binary>(BinaryOp::Pow, sinx, ctx.arena().make<IntegerLit>(BigInt(4)));
    ExprPtr c4 = ctx.arena().make<Binary>(BinaryOp::Pow, cosx, ctx.arena().make<IntegerLit>(BigInt(4)));
    ExprPtr s2 = ctx.arena().make<Binary>(BinaryOp::Pow, sinx, ctx.arena().make<IntegerLit>(BigInt(2)));
    ExprPtr c2 = ctx.arena().make<Binary>(BinaryOp::Pow, cosx, ctx.arena().make<IntegerLit>(BigInt(2)));

    ExprPtr expr = ctx.arena().make<Sum>(std::vector<ExprPtr>{
        s4, 
        ctx.arena().make<Unary>(UnaryOp::Neg, c4),
        ctx.arena().make<Unary>(UnaryOp::Neg, s2),
        c2
    });

    auto res = ctx.simplify(expr);
    ASSERT_TRUE(res.is_ok());
    std::cout << "[ ACID DEBUG ] Test 3 Simplified: " << debug_print(res.value()) << std::endl;
    EXPECT_TRUE(poly_is_zero_expr(res.value()));
}

// TEST 4: SYMBOLIC MATRIX EIGEN-DECOMPOSITION
TEST_F(AcidTest, Test4_SymbolicEigenvals) {
    MatrixExpr m(2, 2);
    m(0, 0) = ctx.arena().make<IntegerLit>(BigInt(0));
    m(0, 1) = ctx.arena().make<IntegerLit>(BigInt(1));
    m(1, 0) = ctx.arena().make<IntegerLit>(BigInt(-1));
    m(1, 1) = ctx.arena().make<IntegerLit>(BigInt(0));

    auto res = eigenvalues(m, ctx);
    ASSERT_TRUE(res.is_ok()) << (res.is_error() ? res.error().message : "");
    std::cout << "[ ACID DEBUG ] Test 4 Result size: " << res.value().size() << std::endl;
    for (auto v : res.value()) std::cout << "  - " << debug_print(v) << std::endl;

    ASSERT_EQ(res.value().size(), 2U);

    bool found_i = false;
    bool found_neg_i = false;

    for (auto val : res.value()) {
        if (const auto* c = expr_cast<Constant>(val)) {
            if (c->value == MathConstant::I) found_i = true;
        } else if (const auto* u = expr_cast<Unary>(val)) {
            if (u->op == UnaryOp::Neg && expr_is<Constant>(u->operand) && expr_cast<Constant>(u->operand)->value == MathConstant::I) {
                found_neg_i = true;
            }
        }
    }

    EXPECT_TRUE(found_i);
    EXPECT_TRUE(found_neg_i);
}

// TEST 5: STRESS TEST / MEMORY ARENA
TEST_F(AcidTest, Test5_ExpansionStress) {
    ctx.set_timeout(std::chrono::milliseconds(60000));
    Symbol x("x");
    Symbol y("y");
    ExprPtr sum = ctx.arena().make<Binary>(BinaryOp::Add, 
        ctx.arena().make<Symbol>("x"), 
        ctx.arena().make<Symbol>("y"));
    ExprPtr expr = ctx.arena().make<Binary>(BinaryOp::Pow, sum, ctx.arena().make<IntegerLit>(BigInt(100)));

    auto start = std::chrono::steady_clock::now();
    auto res = expand(expr, ctx);
    auto end = std::chrono::steady_clock::now();

    ASSERT_TRUE(res.is_ok()) << (res.is_error() ? res.error().message : "");
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "[ PERF ] (x+y)^100 expansion took " << duration << "ms" << std::endl;
    EXPECT_LT(duration, 500);
}

// --- TEST 6: TRIGONOMETRIC LINEARIZATION ---
// sin(x)^3 → 3/4*sin(x) - 1/4*sin(3*x)
// Weak CAS leaves sin(x)^3 unchanged. Strong CAS applies multiple-angle identity.
// NOTE: no dedicated tcollect() API; tests ctx.simplify() pipeline.
// The math_equal check is always valid (both forms are identical functions);
// structural_equal check is the real discriminator — it FAILS if simplifier is passive.
TEST_F(AcidTest, Test6_TrigLinearization) {
    ExprPtr x_ptr = ctx.arena().make<Symbol>("x");
    ExprPtr sinx   = ctx.arena().make<FuncCall>("sin", std::vector<ExprPtr>{x_ptr});
    ExprPtr expr   = ctx.arena().make<Binary>(BinaryOp::Pow, sinx, ctx.arena().make<IntegerLit>(BigInt(3)));

    auto res = ctx.simplify(expr);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    // Expected: 3/4*sin(x) - 1/4*sin(3*x)
    ExprPtr three   = ctx.arena().make<IntegerLit>(BigInt(3));
    ExprPtr four    = ctx.arena().make<IntegerLit>(BigInt(4));
    ExprPtr one     = ctx.arena().make<IntegerLit>(BigInt(1));
    ExprPtr sin3x   = ctx.arena().make<FuncCall>("sin",
        std::vector<ExprPtr>{ctx.arena().make<Binary>(BinaryOp::Mul, three, x_ptr)});
    ExprPtr c_3_4   = ctx.arena().make<Binary>(BinaryOp::Div, three, four);
    ExprPtr c_1_4   = ctx.arena().make<Binary>(BinaryOp::Div, one, four);
    ExprPtr expected = ctx.arena().make<Binary>(BinaryOp::Sub,
        ctx.arena().make<Binary>(BinaryOp::Mul, c_3_4, sinx),
        ctx.arena().make<Binary>(BinaryOp::Mul, c_1_4, sin3x));

    // Must be mathematically correct (always true for any valid simplification)
    auto eq = mathematically_equal(res.value(), expected, ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value()) << "Result not mathematically equal to 3/4*sin(x) - 1/4*sin(3*x)";

    // Must differ structurally from the original — otherwise the simplifier is passive
    // on trig linearization (this is the real discriminator for weak vs strong CAS).
    EXPECT_FALSE(structural_equal(res.value(), expr))
        << "WEAK CAS: simplifier returned sin(x)^3 unchanged — trig linearization not implemented";
}

// --- TEST 7: GAUSSIAN INTEGRAL ---
// int(exp(-x^2), x, -inf, +inf) = sqrt(pi)
// Tests recognition of non-elementary primitives over infinite domains.
// Engine must NOT try to evaluate erf(±∞) symbolically — it should pattern-match
// the Gaussian and return sqrt(pi) directly.
TEST_F(AcidTest, Test7_GaussianIntegral) {
    ExprPtr x_ptr  = ctx.arena().make<Symbol>("x");
    ExprPtr neg_x2 = ctx.arena().make<Unary>(UnaryOp::Neg,
        ctx.arena().make<Binary>(BinaryOp::Pow, x_ptr, ctx.arena().make<IntegerLit>(BigInt(2))));
    ExprPtr integrand = ctx.arena().make<FuncCall>("exp", std::vector<ExprPtr>{neg_x2});

    ExprPtr pos_inf = ctx.arena().make<Constant>(MathConstant::Infinity);
    ExprPtr neg_inf = ctx.arena().make<Unary>(UnaryOp::Neg, pos_inf);

    Symbol x("x");
    auto res = definite_integral(integrand, x, neg_inf, pos_inf, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    std::cout << "[ ACID DEBUG ] Test 7 Result: " << debug_print(res.value()) << std::endl;

    // Expected: sqrt(pi) = pi^(1/2)
    ExprPtr pi   = ctx.arena().make<Constant>(MathConstant::Pi);
    ExprPtr half = ctx.arena().make<Binary>(BinaryOp::Div,
        ctx.arena().make<IntegerLit>(BigInt(1)), ctx.arena().make<IntegerLit>(BigInt(2)));
    ExprPtr sqrt_pi = ctx.arena().make<Binary>(BinaryOp::Pow, pi, half);

    auto eq = mathematically_equal(res.value(), sqrt_pi, ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value()) << "Result not equal to sqrt(pi)";
}

// --- TEST 8: PARTIAL FRACTIONS (Repeated Linear Factor) ---
// partfrac((x^2+1) / (x^3-x^2-x+1))
// Denominator = (x-1)^2*(x+1) → 3-term decomposition required.
// NOTE: User spec "1/(x+1) + 1/(x-1)^2" is mathematically INCORRECT.
// Correct: (1/2)/(x+1) + (1/2)/(x-1) + 1/(x-1)^2.  Verified below via oracle.
TEST_F(AcidTest, Test8_PartialFractionsRepeatedFactor) {
    ExprPtr x_ptr = ctx.arena().make<Symbol>("x");
    Symbol x("x");

    // Numerator: x^2 + 1
    ExprPtr num = ctx.arena().make<Binary>(BinaryOp::Add,
        ctx.arena().make<Binary>(BinaryOp::Pow, x_ptr, ctx.arena().make<IntegerLit>(BigInt(2))),
        ctx.arena().make<IntegerLit>(BigInt(1)));

    // Denominator: ((x^3 - x^2) - x) + 1  i.e. x^3 - x^2 - x + 1
    ExprPtr den = ctx.arena().make<Binary>(BinaryOp::Add,
        ctx.arena().make<Binary>(BinaryOp::Sub,
            ctx.arena().make<Binary>(BinaryOp::Sub,
                ctx.arena().make<Binary>(BinaryOp::Pow, x_ptr, ctx.arena().make<IntegerLit>(BigInt(3))),
                ctx.arena().make<Binary>(BinaryOp::Pow, x_ptr, ctx.arena().make<IntegerLit>(BigInt(2)))),
            x_ptr),
        ctx.arena().make<IntegerLit>(BigInt(1)));

    ExprPtr expr = ctx.arena().make<Binary>(BinaryOp::Div, num, den);

    auto terms = partial_fractions(expr, x, ctx);
    ASSERT_TRUE(terms.is_ok()) << terms.error().message;

    std::cout << "[ ACID DEBUG ] Test 8: " << terms.value().size() << " PF terms" << std::endl;
    for (auto& t : terms.value()) {
        std::cout << "  " << debug_print(t) << std::endl;
    }

    // Denominator = (x-1)^2*(x+1) → 3 terms: A/(x+1), B/(x-1), C/(x-1)^2
    ASSERT_EQ(terms.value().size(), 3U) << "Expected 3 partial fraction terms for repeated root";

    // Oracle: sum of terms must equal original rational expression
    ExprPtr rebuilt = terms.value()[0];
    for (std::size_t i = 1; i < terms.value().size(); ++i) {
        rebuilt = ctx.arena().make<Binary>(BinaryOp::Add, rebuilt, terms.value()[i]);
    }
    auto eq = mathematically_equal(rebuilt, expr, ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value()) << "Sum of PF terms != original (x^2+1)/(x^3-x^2-x+1)";
}

// --- TEST 9: CYCLOTOMIC ROOTS ---
// solve(x^6 - 1 = 0) → 6 roots: ±1, ±1/2 ± i*sqrt(3)/2
// Weak CAS returns decimals or crashes on degree-6.
// Strong CAS returns 2 integer roots + 4 algebraic/exact complex roots.
TEST_F(AcidTest, Test9_CyclotomicRoots) {
    ExprPtr x_ptr = ctx.arena().make<Symbol>("x");
    Symbol x("x");

    // x^6 - 1
    ExprPtr expr = ctx.arena().make<Binary>(BinaryOp::Sub,
        ctx.arena().make<Binary>(BinaryOp::Pow, x_ptr, ctx.arena().make<IntegerLit>(BigInt(6))),
        ctx.arena().make<IntegerLit>(BigInt(1)));

    auto roots = solve_polynomial(expr, x, ctx);
    ASSERT_TRUE(roots.is_ok()) << roots.error().message;

    std::cout << "[ ACID DEBUG ] Test 9: " << roots.value().size() << " roots" << std::endl;
    for (auto& r : roots.value()) {
        std::cout << "  " << debug_print(r) << std::endl;
    }

    ASSERT_EQ(roots.value().size(), 6U) << "x^6-1 must have exactly 6 roots";

    // Oracle: each root r must satisfy r^6 - 1 = 0
    ExprPtr zero = ctx.arena().make<IntegerLit>(BigInt(0));
    for (ExprPtr root : roots.value()) {
        auto subst = ctx.substitute(expr, x, root);
        ASSERT_TRUE(subst.is_ok()) << "substitute failed for root: " << debug_print(root);
        auto simp = ctx.simplify(subst.value());
        ASSERT_TRUE(simp.is_ok()) << "simplify failed for root: " << debug_print(root);
        auto is_zero = mathematically_equal(simp.value(), zero, ctx);
        ASSERT_TRUE(is_zero.is_ok()) << is_zero.error().message;
        EXPECT_TRUE(is_zero.value())
            << "Root " << debug_print(root) << " does not satisfy x^6-1=0 — not exact";
    }
}

// --- TEST 10: EXTENDED EUCLIDEAN ALGORITHM ---
// egcd(240, 46) → [gcd=2, a=-9, b=47] such that 240*(-9) + 46*47 = 2
// Fundamental for RSA/Bézout. Must return exact integer coefficients, not floats.
TEST_F(AcidTest, Test10_ExtendedGCD) {
    auto result = numtheory::extended_gcd(numtheory::Integer(240), numtheory::Integer(46));
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    const auto& [g, a, b] = result.value();
    std::cout << "[ ACID DEBUG ] Test 10: gcd=" << g.decimal()
              << " a=" << a.decimal() << " b=" << b.decimal() << std::endl;

    EXPECT_EQ(g.decimal(), "2") << "GCD(240,46) must be 2";

    // Bézout identity: 240*a + 46*b = 2
    numtheory::Integer check = numtheory::Integer(240) * a + numtheory::Integer(46) * b;
    EXPECT_EQ(check.decimal(), "2")
        << "Bezout identity failed: 240*" << a.decimal() << " + 46*" << b.decimal()
        << " = " << check.decimal();

    // Canonical coefficients (HP Prime spec): a=-9, b=47
    EXPECT_EQ(a.decimal(), "-9") << "Expected canonical Bezout coefficient a=-9";
    EXPECT_EQ(b.decimal(), "47") << "Expected canonical Bezout coefficient b=47";
}

// --- TEST 11: MACLAURIN SERIES WITH ERROR TERM ---
// taylor(sin(x)/x, x, 0, 5) → 1 - x^2/6 + x^4/120 + O(x^6)
// sin(x)/x has a removable singularity at 0; the engine must handle it via
// series division (sin series ÷ x) rather than direct differentiation at 0.
// If O(x^6) remainder is absent, the result is numeric, not symbolic.
TEST_F(AcidTest, Test11_MaclaurinSeriesWithRemainder) {
    ExprPtr x_ptr = ctx.arena().make<Symbol>("x");
    Symbol x("x");
    ExprPtr sinx  = ctx.arena().make<FuncCall>("sin", std::vector<ExprPtr>{x_ptr});
    ExprPtr expr  = ctx.arena().make<Binary>(BinaryOp::Div, sinx, x_ptr);
    ExprPtr zero  = ctx.arena().make<IntegerLit>(BigInt(0));

    auto expansion = taylor_series(expr, x, zero, 5U, ctx);
    ASSERT_TRUE(expansion.is_ok()) << expansion.error().message;
    ASSERT_EQ(expansion.value().computed_order, 5U);

    std::cout << "[ ACID DEBUG ] Test 11 poly: " << debug_print(expansion.value().polynomial) << std::endl;

    // Expected polynomial: 1 - x^2/6 + x^4/120
    ExprPtr one = ctx.arena().make<IntegerLit>(BigInt(1));
    ExprPtr expected = ctx.arena().make<Sum>(std::vector<ExprPtr>{
        one,
        ctx.arena().make<Unary>(UnaryOp::Neg,
            ctx.arena().make<Binary>(BinaryOp::Div,
                ctx.arena().make<Binary>(BinaryOp::Pow, x_ptr, ctx.arena().make<IntegerLit>(BigInt(2))),
                ctx.arena().make<IntegerLit>(BigInt(6)))),
        ctx.arena().make<Binary>(BinaryOp::Div,
            ctx.arena().make<Binary>(BinaryOp::Pow, x_ptr, ctx.arena().make<IntegerLit>(BigInt(4))),
            ctx.arena().make<IntegerLit>(BigInt(120)))
    });

    auto eq = mathematically_equal(expansion.value().polynomial, expected, ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value()) << "Polynomial != 1 - x^2/6 + x^4/120";

    // Remainder must be present: O(x^6) is a symbolic claim, not a numeric approximation
    EXPECT_NE(expansion.value().remainder, nullptr)
        << "WEAK CAS: O(x^6) remainder missing — output is numeric, not symbolic";
}

// --- TEST 12: BRANCH CUT TRAP ---
// simplify(ln(-e) - ln(e)) = I*pi
// Principal value: ln(-e) = 1 + I*pi, ln(e) = 1 → difference = I*pi.
// A broken CAS returns 1 (applies real log incorrectly) or crashes.
TEST_F(AcidTest, Test12_BranchCutTrap) {
    ExprPtr e_const = ctx.arena().make<Constant>(MathConstant::E);
    ExprPtr neg_e   = ctx.arena().make<Unary>(UnaryOp::Neg, e_const);

    ExprPtr ln_neg_e = ctx.arena().make<FuncCall>("ln", std::vector<ExprPtr>{neg_e});
    ExprPtr ln_e     = ctx.arena().make<FuncCall>("ln", std::vector<ExprPtr>{e_const});
    ExprPtr expr     = ctx.arena().make<Binary>(BinaryOp::Sub, ln_neg_e, ln_e);

    auto res = ctx.simplify(expr);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    std::cout << "[ ACID DEBUG ] Test 12 Result: " << debug_print(res.value()) << std::endl;

    // Expected: I*pi
    ExprPtr I_const  = ctx.arena().make<Constant>(MathConstant::I);
    ExprPtr pi_const = ctx.arena().make<Constant>(MathConstant::Pi);
    ExprPtr expected = ctx.arena().make<Binary>(BinaryOp::Mul, I_const, pi_const);

    auto eq = mathematically_equal(res.value(), expected, ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value()) << "Result not equal to I*pi — branch cut not handled correctly";
}

// --- TEST 16: BASEL PROBLEM (Symbolic Infinite Summation) ---
// sum(1/k^2, k, 1, infinity) = pi^2/6
TEST_F(AcidTest, Test16_BaselProblem) {
    Symbol k("k");
    ExprPtr k_ptr = ctx.arena().make<Symbol>("k");
    ExprPtr term = ctx.arena().make<Binary>(BinaryOp::Div,
        ctx.arena().make<IntegerLit>(BigInt(1)),
        ctx.arena().make<Binary>(BinaryOp::Pow, k_ptr, ctx.arena().make<IntegerLit>(BigInt(2))));
    
    ExprPtr lower = ctx.arena().make<IntegerLit>(BigInt(1));
    ExprPtr upper = ctx.arena().make<Constant>(MathConstant::Infinity);
    
    auto res = sum(term, k, lower, upper, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    
    std::cout << "[ ACID DEBUG ] Test 16 Result: " << debug_print(res.value()) << std::endl;
    
    ExprPtr pi = ctx.arena().make<Constant>(MathConstant::Pi);
    ExprPtr expected = ctx.arena().make<Binary>(BinaryOp::Div,
        ctx.arena().make<Binary>(BinaryOp::Pow, pi, ctx.arena().make<IntegerLit>(BigInt(2))),
        ctx.arena().make<IntegerLit>(BigInt(6)));
        
    auto eq = mathematically_equal(res.value(), expected, ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value()) << "Result != pi^2 / 6";
}

// --- TEST 17: SQUEEZE THEOREM LIMIT (Calculus Blindspot) ---
// limit(x * sin(1/x), x, 0) = 0
TEST_F(AcidTest, Test17_SqueezeTheoremLimit) {
    Symbol x("x");
    ExprPtr x_ptr = ctx.arena().make<Symbol>("x");
    ExprPtr sin_inv_x = ctx.arena().make<FuncCall>("sin", std::vector<ExprPtr>{
        ctx.arena().make<Binary>(BinaryOp::Div, 
            ctx.arena().make<IntegerLit>(BigInt(1)),
            x_ptr)
    });
    ExprPtr expr = ctx.arena().make<Binary>(BinaryOp::Mul, x_ptr, sin_inv_x);
    
    auto res = limit(expr, x, ctx.arena().make<IntegerLit>(BigInt(0)), LimitDirection::Both, ctx);
    ASSERT_TRUE(res.is_ok()) << (res.is_error() ? res.error().message : "");
    
    std::cout << "[ ACID DEBUG ] Test 17 Result: " << debug_print(res.value()) << std::endl;
    
    EXPECT_TRUE(poly_is_zero_expr(res.value()));
}

// --- TEST 18: SYMBOLIC COMPLEX RESIDUES (Advanced Calculus) ---
// residue(1 / (z^2 + 1)^2, z, I) = -I / 4
TEST_F(AcidTest, Test18_SymbolicComplexResidues) {
    Symbol z("z");
    ExprPtr z_ptr = ctx.arena().make<Symbol>("z");
    ExprPtr den = ctx.arena().make<Binary>(BinaryOp::Pow,
        ctx.arena().make<Binary>(BinaryOp::Add,
            ctx.arena().make<Binary>(BinaryOp::Pow, z_ptr, ctx.arena().make<IntegerLit>(BigInt(2))),
            ctx.arena().make<IntegerLit>(BigInt(1))),
        ctx.arena().make<IntegerLit>(BigInt(2)));
    ExprPtr expr = ctx.arena().make<Binary>(BinaryOp::Div, ctx.arena().make<IntegerLit>(BigInt(1)), den);
    
    ExprPtr pole = ctx.arena().make<Constant>(MathConstant::I);
    
    auto res = residue(expr, z, pole, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    
    std::cout << "[ ACID DEBUG ] Test 18 Result: " << debug_print(res.value()) << std::endl;
    
    // Expected: -I / 4
    ExprPtr I_const = ctx.arena().make<Constant>(MathConstant::I);
    ExprPtr expected = ctx.arena().make<Binary>(BinaryOp::Div,
        ctx.arena().make<Unary>(UnaryOp::Neg, I_const),
        ctx.arena().make<IntegerLit>(BigInt(4)));
        
    auto eq = mathematically_equal(res.value(), expected, ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value()) << "Result != -I / 4";
}

// --- TEST 19: PATHOLOGICAL GRÖBNER SYSTEM (Algebraic Geometry) ---
// csolve([x^2 + y^2 = 1, x^2 - y = 0], [x, y])
TEST_F(AcidTest, Test19_GroebnerSystem) {
    Symbol x("x");
    Symbol y("y");
    ExprPtr x_ptr = ctx.arena().make<Symbol>("x");
    ExprPtr y_ptr = ctx.arena().make<Symbol>("y");
    ExprPtr one = ctx.arena().make<IntegerLit>(BigInt(1));
    ExprPtr zero = ctx.arena().make<IntegerLit>(BigInt(0));

    // eqs = [x^2 + y^2 = 1, x^2 - y = 0]
    ExprPtr eq1 = ctx.arena().make<Binary>(BinaryOp::Equal,
        ctx.arena().make<Sum>(std::vector<ExprPtr>{
            ctx.arena().make<Binary>(BinaryOp::Pow, x_ptr, ctx.arena().make<IntegerLit>(BigInt(2))),
            ctx.arena().make<Binary>(BinaryOp::Pow, y_ptr, ctx.arena().make<IntegerLit>(BigInt(2)))
        }),
        one);
    ExprPtr eq2 = ctx.arena().make<Binary>(BinaryOp::Equal,
        ctx.arena().make<Binary>(BinaryOp::Sub,
            ctx.arena().make<Binary>(BinaryOp::Pow, x_ptr, ctx.arena().make<IntegerLit>(BigInt(2))),
            y_ptr),
        zero);
    
    ExprPtr eqs = ctx.arena().make<cas::Matrix>(1, 2, std::vector<ExprPtr>{eq1, eq2});
    ExprPtr vars = ctx.arena().make<cas::Matrix>(1, 2, std::vector<ExprPtr>{x_ptr, y_ptr});

    auto res = csolve(eqs, vars, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    std::cout << "[ ACID DEBUG ] Test 19 Result: " << debug_print(res.value()) << std::endl;

    const auto* m = expr_cast<cas::Matrix>(res.value());    ASSERT_NE(m, nullptr);
    // x^2 + y^2 = 1, x^2 = y => y^2 + y - 1 = 0 => y = (-1 ± sqrt(5))/2
    // For each y, x = ±sqrt(y). 
    // Total 4 solutions.
    EXPECT_GE(m->rows, 2U); // Almeno le soluzioni reali
    
    // Verifica una soluzione: x^2 = y
    for (std::size_t i = 0; i < m->rows; ++i) {
        ExprPtr sol_x = m->elements[i * 2];
        ExprPtr sol_y = m->elements[i * 2 + 1];
        
        // x^2 - y should be 0
        auto x2 = ctx.arena().make<Binary>(BinaryOp::Pow, sol_x, ctx.arena().make<IntegerLit>(BigInt(2)));
        auto diff = ctx.arena().make<Binary>(BinaryOp::Sub, x2, sol_y);
        auto simp = ctx.simplify(diff);
        ASSERT_TRUE(simp.is_ok());
        
        auto is_zero = mathematically_equal(simp.value(), zero, ctx);
        ASSERT_TRUE(is_zero.is_ok());
        EXPECT_TRUE(is_zero.value()) << "Solution " << i << " does not satisfy x^2 = y: "
                                     << debug_print(sol_x) << ", " << debug_print(sol_y);
    }
}

// ==============================================================================
// BATCH 5: THE NUCLEAR OPTION (IDENTITY & CONSISTENCY)
// ==============================================================================

// --- TEST 25: THE IDENTITY TRAP (Zero-Recognition) ---
// DIAGNOSI: ln(sqrt(e)) = ln(e^(1/2)) = 1/2  ⟹  ln(sqrt(e)) - 1/2 = 0
// Se il risultato non è esattamente 0, il simplifier ha un leak
// nelle identità logaritmiche/esponenziali.
TEST_F(AcidTest, Test25_IdentityTrapLnSqrtE) {
    // Costruisce: ln(sqrt(e)) - 1/2
    ExprPtr e_const  = ctx.arena().make<Constant>(MathConstant::E);
    ExprPtr sqrt_e   = ctx.arena().make<FuncCall>("sqrt", std::vector<ExprPtr>{e_const});
    ExprPtr ln_sqrt_e = ctx.arena().make<FuncCall>("ln", std::vector<ExprPtr>{sqrt_e});
    ExprPtr half = ctx.arena().make<Binary>(BinaryOp::Div,
        ctx.arena().make<IntegerLit>(BigInt(1)),
        ctx.arena().make<IntegerLit>(BigInt(2)));
    ExprPtr expr = ctx.arena().make<Binary>(BinaryOp::Sub, ln_sqrt_e, half);

    auto res = ctx.simplify(expr);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    std::cout << "[ ACID DEBUG ] Test 25 Result: " << debug_print(res.value()) << std::endl;
    EXPECT_TRUE(poly_is_zero_expr(res.value()))
        << "FAIL: ln(sqrt(e)) - 1/2 non si semplifica a 0. Got: " << debug_print(res.value());
}

// --- TEST 26: SCHANUEL'S CONJECTURE LIMITS ---
// DIAGNOSI: e^(pi + ln(x)) / x = e^pi * x / x = e^pi  ⟹  espressione - e^pi = 0
// Richiede la regola e^(ln(x)) = x (per variabili simboliche).
TEST_F(AcidTest, Test26_SchanuelLimits) {
    ExprPtr x_ptr    = ctx.arena().make<Symbol>("x");
    ExprPtr pi_const = ctx.arena().make<Constant>(MathConstant::Pi);
    ExprPtr e_const  = ctx.arena().make<Constant>(MathConstant::E);

    // e^(pi + ln(x))
    ExprPtr ln_x        = ctx.arena().make<FuncCall>("ln", std::vector<ExprPtr>{x_ptr});
    ExprPtr pi_plus_lnx = ctx.arena().make<Binary>(BinaryOp::Add, pi_const, ln_x);
    ExprPtr e_pow       = ctx.arena().make<Binary>(BinaryOp::Pow, e_const, pi_plus_lnx);

    // e^(pi + ln(x)) / x
    ExprPtr frac = ctx.arena().make<Binary>(BinaryOp::Div, e_pow, x_ptr);

    // e^pi
    ExprPtr e_pi = ctx.arena().make<Binary>(BinaryOp::Pow, e_const, pi_const);

    // espressione completa: frac - e^pi
    ExprPtr expr = ctx.arena().make<Binary>(BinaryOp::Sub, frac, e_pi);

    auto res = ctx.simplify(expr);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    std::cout << "[ ACID DEBUG ] Test 26 Result: " << debug_print(res.value()) << std::endl;
    EXPECT_TRUE(poly_is_zero_expr(res.value()))
        << "FAIL: e^(pi+ln(x))/x - e^pi non si semplifica a 0. Got: " << debug_print(res.value());
}

// --- TEST 27: MULTIVARIATE POLYNOMIAL GCD ---
// DIAGNOSI: gcd(x^2 - y^2, x^2 + 2*x*y + y^2) = x + y
// Richiede algoritmo EZ-GCD o Heuristic GCD multivariato.
TEST_F(AcidTest, Test27_MultivariatePolynomialGCD) {
    ExprPtr x_ptr = ctx.arena().make<Symbol>("x");
    ExprPtr y_ptr = ctx.arena().make<Symbol>("y");
    Symbol x("x"), y("y");

    // p = x^2 - y^2
    ExprPtr x2 = ctx.arena().make<Binary>(BinaryOp::Pow, x_ptr, ctx.arena().make<IntegerLit>(BigInt(2)));
    ExprPtr y2 = ctx.arena().make<Binary>(BinaryOp::Pow, y_ptr, ctx.arena().make<IntegerLit>(BigInt(2)));
    ExprPtr p = ctx.arena().make<Binary>(BinaryOp::Sub, x2, y2);

    // q = x^2 + 2*x*y + y^2 = (x+y)^2
    ExprPtr two_xy = ctx.arena().make<Binary>(BinaryOp::Mul,
        ctx.arena().make<IntegerLit>(BigInt(2)),
        ctx.arena().make<Binary>(BinaryOp::Mul, x_ptr, y_ptr));
    ExprPtr q = ctx.arena().make<Binary>(BinaryOp::Add,
        ctx.arena().make<Binary>(BinaryOp::Add, x2, two_xy), y2);

    auto res = polynomial_gcd_multivariate(p, q, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    std::cout << "[ ACID DEBUG ] Test 27 GCD: " << debug_print(res.value()) << std::endl;

    // Expected: x + y  (o equivalente, es. y + x)
    ExprPtr expected = ctx.arena().make<Binary>(BinaryOp::Add, x_ptr, y_ptr);
    auto eq = mathematically_equal(res.value(), expected, ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value())
        << "FAIL: gcd(x^2-y^2, x^2+2xy+y^2) != x+y. Got: " << debug_print(res.value());
}

// --- TEST 28: LIMITS OF NESTED LOGARITHMS (Hard Gruntz) ---
// DIAGNOSI: limit(ln(ln(x + e)) / ln(ln(x)), x -> inf) = 1
// Se restituisce 1 o infinito (ERRATO) il calcolo gerarchie di crescita è fallato.
TEST_F(AcidTest, Test28_NestedLogLimit) {
    ctx.set_timeout(std::chrono::milliseconds(30000));
    Symbol x("x");
    ExprPtr x_ptr   = ctx.arena().make<Symbol>("x");
    ExprPtr e_const = ctx.arena().make<Constant>(MathConstant::E);

    // ln(ln(x + e))
    ExprPtr x_plus_e = ctx.arena().make<Binary>(BinaryOp::Add, x_ptr, e_const);
    ExprPtr ln_xpe   = ctx.arena().make<FuncCall>("ln", std::vector<ExprPtr>{x_plus_e});
    ExprPtr ln_ln_xpe = ctx.arena().make<FuncCall>("ln", std::vector<ExprPtr>{ln_xpe});

    // ln(ln(x))
    ExprPtr ln_x   = ctx.arena().make<FuncCall>("ln", std::vector<ExprPtr>{x_ptr});
    ExprPtr ln_ln_x = ctx.arena().make<FuncCall>("ln", std::vector<ExprPtr>{ln_x});

    ExprPtr expr = ctx.arena().make<Binary>(BinaryOp::Div, ln_ln_xpe, ln_ln_x);

    ExprPtr inf = ctx.arena().make<Constant>(MathConstant::Infinity);
    auto res = limit(expr, x, inf, LimitDirection::Both, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    std::cout << "[ ACID DEBUG ] Test 28 Result: " << debug_print(res.value()) << std::endl;

    // Expected: 1
    ExprPtr one = ctx.arena().make<IntegerLit>(BigInt(1));
    auto eq = mathematically_equal(res.value(), one, ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value())
        << "FAIL: limit(ln(ln(x+e))/ln(ln(x)), inf) != 1. Got: " << debug_print(res.value());
}

// --- TEST 29: TENSOR-LIKE PRODUCT CONSISTENCY ---
// DIAGNOSI: A*B - B*A ≠ 0 per matrici simboliche 2x2 generiche.
// Verifica la non-commutatività del prodotto matriciale.
TEST_F(AcidTest, Test29_MatrixNonCommutativity) {
    // Costruisce A = [[a, b], [c, d]]
    ExprPtr a = ctx.arena().make<Symbol>("a");
    ExprPtr b = ctx.arena().make<Symbol>("b");
    ExprPtr c = ctx.arena().make<Symbol>("c");
    ExprPtr d = ctx.arena().make<Symbol>("d");
    ExprPtr e_sym = ctx.arena().make<Symbol>("e_m");
    ExprPtr f = ctx.arena().make<Symbol>("f");
    ExprPtr g = ctx.arena().make<Symbol>("g");
    ExprPtr h = ctx.arena().make<Symbol>("h");

    linalg::MatrixExpr A(2, 2);
    A(0,0) = a; A(0,1) = b;
    A(1,0) = c; A(1,1) = d;

    linalg::MatrixExpr B(2, 2);
    B(0,0) = e_sym; B(0,1) = f;
    B(1,0) = g;     B(1,1) = h;

    // AB
    auto AB_res = linalg::multiply(A, B, ctx);
    ASSERT_TRUE(AB_res.is_ok()) << AB_res.error().message;

    // BA
    auto BA_res = linalg::multiply(B, A, ctx);
    ASSERT_TRUE(BA_res.is_ok()) << BA_res.error().message;

    // AB - BA
    auto diff_res = linalg::subtract(AB_res.value(), BA_res.value(), ctx);
    ASSERT_TRUE(diff_res.is_ok()) << diff_res.error().message;

    std::cout << "[ ACID DEBUG ] Test 29 AB-BA[0,0]: " << debug_print(diff_res.value()(0,0)) << std::endl;
    std::cout << "[ ACID DEBUG ] Test 29 AB-BA[0,1]: " << debug_print(diff_res.value()(0,1)) << std::endl;

    // Semplifica ogni elemento della matrice differenza
    bool all_zero = true;
    for (std::size_t r = 0; r < 2; ++r) {
        for (std::size_t col = 0; col < 2; ++col) {
            auto simp = ctx.simplify(diff_res.value()(r, col));
            ASSERT_TRUE(simp.is_ok());
            if (!poly_is_zero_expr(simp.value())) all_zero = false;
        }
    }

    // Il risultato NON deve essere la matrice zero: A e B non commutano in generale
    EXPECT_FALSE(all_zero)
        << "FAIL: A*B - B*A == 0 per matrici simboliche generiche. "
        << "Il motore semplifica erroneamente a zero (falsa commutatività).";

    // Verifica specifica: [0,0] = b*g - c*f  (deve essere non-zero simbolicamente)
    ExprPtr expected_00 = ctx.arena().make<Binary>(BinaryOp::Sub,
        ctx.arena().make<Binary>(BinaryOp::Mul, b, g),
        ctx.arena().make<Binary>(BinaryOp::Mul, c, f));
    auto simp_00 = ctx.simplify(diff_res.value()(0, 0));
    ASSERT_TRUE(simp_00.is_ok());

    auto eq_00 = mathematically_equal(simp_00.value(), expected_00, ctx);
    ASSERT_TRUE(eq_00.is_ok()) << eq_00.error().message;
    EXPECT_TRUE(eq_00.value())
        << "FAIL: (AB-BA)[0,0] != b*g - c*f. Got: " << debug_print(simp_00.value());
}
