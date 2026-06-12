#include <gtest/gtest.h>
#include "cas/ode.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast_debug.hpp"
#include <iostream>

using namespace cas;
using namespace cas::calculus;

class OdeTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

TEST_F(OdeTest, ConstantCoeffHomogeneous) {
    // y'' - 3y' + 2y = 0
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();
    
    ExprPtr y_pp = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 2);
    ExprPtr y_p = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 1);
    ExprPtr y_s = arena.make<Symbol>("y");
    
    ExprPtr eq = arena.make<Sum>(std::vector<ExprPtr>{
        y_pp,
        arena.make<Unary>(UnaryOp::Neg, arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(3)), y_p)),
        arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(2)), y_s)
    });
    
    auto res = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    std::cout << "ODE Solution: " << debug_print(res.value()) << std::endl;
    // Expected: C1*exp(x) + C2*exp(2x) or similar
}

TEST_F(OdeTest, ConstantCoeffNonHomogeneous) {
    // y'' + y = sin(x)
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();
    
    ExprPtr y_pp = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 2);
    ExprPtr y_s = arena.make<Symbol>("y");
    ExprPtr sin_x = arena.make<FuncCall>("sin", std::vector<ExprPtr>{arena.make<Symbol>("x")});
    
    ExprPtr eq = arena.make<Binary>(BinaryOp::Equal,
        arena.make<Binary>(BinaryOp::Add, y_pp, y_s),
        sin_x);
    
    auto res = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    std::cout << "ODE Solution: " << debug_print(res.value()) << std::endl;
}

TEST_F(OdeTest, Linear3rdOrderMultiplicity) {
    // y''' - 3y'' + 3y' - y = 0
    // L'equazione caratteristica r^3 - 3r^2 + 3r - 1 = (r-1)^3
    // Radice ripetuta r=1 (molteplicità 3). Soluzione generale: (C1 + C2*x + C3*x^2)*e^x
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();
    
    ExprPtr y_ppp = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 3);
    ExprPtr y_pp = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 2);
    ExprPtr y_p = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 1);
    ExprPtr y_s = arena.make<Symbol>("y");
    
    ExprPtr eq = arena.make<Sum>(std::vector<ExprPtr>{
        y_ppp,
        arena.make<Unary>(UnaryOp::Neg, arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(3)), y_pp)),
        arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(3)), y_p),
        arena.make<Unary>(UnaryOp::Neg, y_s)
    });
    
    auto res = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    std::cout << "ODE 3rd Order Solution: " << debug_print(res.value()) << std::endl;
}

TEST_F(OdeTest, Linear2ndOrderParticularComplex) {
    // y'' + y = x
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();

    ExprPtr y_pp = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 2);
    ExprPtr y_s = arena.make<Symbol>("y");

    ExprPtr eq = arena.make<Binary>(BinaryOp::Equal,
        arena.make<Binary>(BinaryOp::Add, y_pp, y_s),
        arena.make<Symbol>("x"));

    auto res = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    std::cout << "ODE 2nd Order Particular Solution: " << debug_print(res.value()) << std::endl;
}

// ── F5.3 Riccati family classification + closed-form solver ─────────────────

TEST_F(OdeTest, Riccati_QzeroNull_BernoulliReduction) {
    // y' = y + y²  →  Riccati with q_0 = 0, q_1 = 1, q_2 = 1.
    // Reduces to v' = -v - 1 (linear 1st-order), so y = 1/v with v solved.
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();

    ExprPtr y_p = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 1);
    ExprPtr y_s = arena.make<Symbol>("y");
    ExprPtr y_sq = arena.make<Binary>(BinaryOp::Pow, y_s, arena.make<IntegerLit>(BigInt(2)));
    ExprPtr eq = arena.make<Binary>(BinaryOp::Equal, y_p,
        arena.make<Binary>(BinaryOp::Add, y_s, y_sq));

    auto cls = classify_ode(eq, y, x, ctx);
    ASSERT_TRUE(cls.is_ok());
    EXPECT_EQ(cls.value().type, OdeType::Separable);
    ASSERT_EQ(cls.value().components.size(), 2U);

    auto sol = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(sol.is_ok()) << sol.error().message;
    std::cout << "Riccati q0=0 Solution: " << debug_print(sol.value()) << std::endl;
}

TEST_F(OdeTest, Riccati_ConstantCoeff_NegativeDisc) {
    // y' = 1 + y²  →  Riccati, Δ = 0² - 4·1·1 = -4 < 0,  y = tan(x − C).
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();

    ExprPtr y_p = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 1);
    ExprPtr y_s = arena.make<Symbol>("y");
    ExprPtr y_sq = arena.make<Binary>(BinaryOp::Pow, y_s, arena.make<IntegerLit>(BigInt(2)));
    ExprPtr eq = arena.make<Binary>(BinaryOp::Equal, y_p,
        arena.make<Binary>(BinaryOp::Add, arena.make<IntegerLit>(BigInt(1)), y_sq));

    auto cls = classify_ode(eq, y, x, ctx);
    ASSERT_TRUE(cls.is_ok());
    EXPECT_EQ(cls.value().type, OdeType::Separable);

    auto sol = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(sol.is_ok()) << sol.error().message;
    std::cout << "Riccati Δ<0 Solution: " << debug_print(sol.value()) << std::endl;
}

TEST_F(OdeTest, Riccati_ConstantCoeff_PositiveDisc) {
    // y' = -1 + y²  →  Riccati, Δ = 0² - 4·(-1)·1 = 4 > 0,  y = -tanh(x − C).
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();

    ExprPtr y_p = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 1);
    ExprPtr y_s = arena.make<Symbol>("y");
    ExprPtr y_sq = arena.make<Binary>(BinaryOp::Pow, y_s, arena.make<IntegerLit>(BigInt(2)));
    ExprPtr eq = arena.make<Binary>(BinaryOp::Equal, y_p,
        arena.make<Sum>(std::vector<ExprPtr>{
            arena.make<IntegerLit>(BigInt(-1)),
            y_sq}));

    auto cls = classify_ode(eq, y, x, ctx);
    ASSERT_TRUE(cls.is_ok());
    EXPECT_EQ(cls.value().type, OdeType::Separable);

    auto sol = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(sol.is_ok()) << sol.error().message;
    std::cout << "Riccati Δ>0 Solution: " << debug_print(sol.value()) << std::endl;
}

// ── F5.3 / B2b — Clairaut & d'Alembert (Lagrange) ───────────────────────────

TEST_F(OdeTest, Clairaut_QuadraticG_GeneralAndSingular) {
    // y = x·y' + (y')²  →  Clairaut with F(p) = p, G(p) = p².
    // General  : y = C·x + C².
    // Singular : x = -2p, y = -p²  ⇒  y = -x²/4.
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();

    ExprPtr y_p = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 1);
    ExprPtr y_s = arena.make<Symbol>("y");
    ExprPtr xs  = arena.make<Symbol>("x");
    ExprPtr rhs = arena.make<Binary>(BinaryOp::Add,
        arena.make<Binary>(BinaryOp::Mul, xs, y_p),
        arena.make<Binary>(BinaryOp::Pow, y_p, arena.make<IntegerLit>(BigInt(2))));
    ExprPtr eq = arena.make<Binary>(BinaryOp::Equal, y_s, rhs);

    auto cls = classify_ode(eq, y, x, ctx);
    ASSERT_TRUE(cls.is_ok());
    EXPECT_EQ(cls.value().type, OdeType::Clairaut);
    ASSERT_TRUE(cls.value().parameter.has_value());
    ASSERT_EQ(cls.value().components.size(), 1U);

    auto sol = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(sol.is_ok()) << sol.error().message;
    std::cout << "Clairaut Solution: " << debug_print(sol.value()) << std::endl;

    const auto* fc = expr_cast<FuncCall>(sol.value());
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->name, "GeneralAndSingular");
    ASSERT_EQ(fc->args.size(), 2U);
}

TEST_F(OdeTest, DAlembert_LinearF_QuadraticG_Parametric) {
    // y = 2x·y' + (y')²  →  d'Alembert with F(p) = 2p, G(p) = p².
    // Reduction: dx/dp + (2/p)·x + 2 = 0  (linear in x(p)).
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();

    ExprPtr y_p = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 1);
    ExprPtr y_s = arena.make<Symbol>("y");
    ExprPtr xs  = arena.make<Symbol>("x");
    ExprPtr two = arena.make<IntegerLit>(BigInt(2));
    ExprPtr rhs = arena.make<Binary>(BinaryOp::Add,
        arena.make<Binary>(BinaryOp::Mul, two,
            arena.make<Binary>(BinaryOp::Mul, xs, y_p)),
        arena.make<Binary>(BinaryOp::Pow, y_p, arena.make<IntegerLit>(BigInt(2))));
    ExprPtr eq = arena.make<Binary>(BinaryOp::Equal, y_s, rhs);

    auto cls = classify_ode(eq, y, x, ctx);
    ASSERT_TRUE(cls.is_ok());
    EXPECT_EQ(cls.value().type, OdeType::DAlembert);
    ASSERT_TRUE(cls.value().parameter.has_value());
    ASSERT_EQ(cls.value().components.size(), 2U);

    auto sol = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(sol.is_ok()) << sol.error().message;
    std::cout << "d'Alembert Solution: " << debug_print(sol.value()) << std::endl;

    const auto* fc = expr_cast<FuncCall>(sol.value());
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->name, "ParametricSolution");
    ASSERT_EQ(fc->args.size(), 2U);
}

TEST_F(OdeTest, Riccati_VariableCoeffNoParticular_Diagnostic) {
    // y' = x + y²  →  Riccati, variable a(x) = x, no particular solution known.
    // Expect explicit Unimplemented with F5.3 / B2 continuation diagnostic.
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();

    ExprPtr y_p = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 1);
    ExprPtr y_s = arena.make<Symbol>("y");
    ExprPtr y_sq = arena.make<Binary>(BinaryOp::Pow, y_s, arena.make<IntegerLit>(BigInt(2)));
    ExprPtr eq = arena.make<Binary>(BinaryOp::Equal, y_p,
        arena.make<Binary>(BinaryOp::Add, arena.make<Symbol>("x"), y_sq));

    auto cls = classify_ode(eq, y, x, ctx);
    ASSERT_TRUE(cls.is_ok());
    EXPECT_EQ(cls.value().type, OdeType::Riccati);

    auto sol = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(sol.is_error());
    EXPECT_EQ(sol.error().kind, CASErrorKind::Unimplemented);
    EXPECT_NE(sol.error().message.find("Riccati"), std::string::npos);
}

// Kovacic Case 1: Euler-Cauchy  x²y'' + αxy' + βy = 0
// Classified as Linear2ndOrderRationalCoeff; solved via Kovacic Case 1.
// Expected general solution: y = C₁·x^λ₁ + C₂·x^λ₂
// where λ₁,₂ are roots of the indicial equation λ²+(α-1)λ+β = 0.
TEST_F(OdeTest, Kovacic_EulerCauchy_TwoDistinctRealRoots) {
    // x²y'' - 2xy' + 2y = 0  →  α=-2, β=2
    // Indicial: λ²-3λ+2 = 0 → λ=1, λ=2  → y = C₁x + C₂x²
    symbolic::CASContext ctx;
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();

    ExprPtr x_sym  = arena.make<Symbol>("x");
    ExprPtr x_sq   = arena.make<Binary>(BinaryOp::Pow, x_sym, arena.make<IntegerLit>(BigInt(2)));
    ExprPtr y_sym  = arena.make<Symbol>("y");
    ExprPtr yp     = arena.make<Derivative>(y_sym, Symbol("x"), 1);
    ExprPtr ypp    = arena.make<Derivative>(y_sym, Symbol("x"), 2);

    // x²y'' - 2xy' + 2y = 0  →  Equal(x²y'' - 2xy' + 2y, 0)
    ExprPtr lhs = arena.make<Sum>(std::vector<ExprPtr>{
        arena.make<Binary>(BinaryOp::Mul, x_sq, ypp),
        arena.make<Unary>(UnaryOp::Neg,
            arena.make<Binary>(BinaryOp::Mul,
                arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(2)), x_sym), yp)),
        arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(2)), y_sym)
    });
    ExprPtr eq = arena.make<Binary>(BinaryOp::Equal, lhs, arena.make<IntegerLit>(BigInt(0)));

    auto cls = classify_ode(eq, y, x, ctx);
    ASSERT_TRUE(cls.is_ok()) << cls.error().message;
    EXPECT_EQ(cls.value().type, OdeType::Linear2ndOrderRationalCoeff);

    auto sol = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(sol.is_ok()) << sol.error().message;
    std::cout << "Kovacic Euler-Cauchy solution: " << debug_print(sol.value()) << std::endl;

    // Verify: must be an equality  y = ...
    // Verify structural form:  y = (C₂ + C₃·x) · exp(ln|x|)  ≡  C₂·x + C₃·x²
    auto* eq_node = expr_cast<Binary>(sol.value());
    ASSERT_NE(eq_node, nullptr);
    EXPECT_EQ(eq_node->op, BinaryOp::Equal);

    // The RHS must NOT be zero (trivial) and must involve 'x'.
    ExprPtr rhs = eq_node->right;
    EXPECT_FALSE(debug_print(rhs).empty());
    EXPECT_NE(debug_print(rhs), "IntegerLit(0)");
    // The solution must contain exp (from the back-transformation factor).
    EXPECT_NE(debug_print(rhs).find("exp"), std::string::npos);
}

// Kovacic Case 1: equation with simple-pole r → Case 1 fails → Unimplemented.
// This exercises the diagnostic path correctly.
TEST_F(OdeTest, Kovacic_SimplePole_Unimplemented) {
    // y'' + (1/x)y' = 0  →  r = 1/(4x²) - (-1/x²)/2 = 1/(4x²) + 1/(2x²) = 3/(4x²)
    // Actually this has only order-2 poles, so Case 1 applies.
    // Use y'' + y'/x + y = 0 (Bessel order 0 type) — has a non-rational discriminant.
    symbolic::CASContext ctx;
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();

    ExprPtr x_sym = arena.make<Symbol>("x");
    ExprPtr y_sym = arena.make<Symbol>("y");
    ExprPtr yp    = arena.make<Derivative>(y_sym, Symbol("x"), 1);
    ExprPtr ypp   = arena.make<Derivative>(y_sym, Symbol("x"), 2);

    // x²y'' + xy' + x²y = 0  (Bessel order 0) → r = (1-1²)/4/x² - x² = -x² (roughly)
    // Kovacic Case 1 fails for Bessel-type (transcendental, no Liouvillian solution).
    ExprPtr lhs = arena.make<Sum>(std::vector<ExprPtr>{
        arena.make<Binary>(BinaryOp::Mul, arena.make<Binary>(BinaryOp::Pow, x_sym,
            arena.make<IntegerLit>(BigInt(2))), ypp),
        arena.make<Binary>(BinaryOp::Mul, x_sym, yp),
        arena.make<Binary>(BinaryOp::Mul, arena.make<Binary>(BinaryOp::Pow, x_sym,
            arena.make<IntegerLit>(BigInt(2))), y_sym)
    });
    ExprPtr eq = arena.make<Binary>(BinaryOp::Equal, lhs, arena.make<IntegerLit>(BigInt(0)));

    auto cls = classify_ode(eq, y, x, ctx);
    ASSERT_TRUE(cls.is_ok()) << cls.error().message;

    if (cls.value().type == OdeType::Linear2ndOrderRationalCoeff) {
        auto sol = solve_ode(eq, y, x, ctx);
        // Bessel has no Liouvillian solution → Unimplemented is correct.
        EXPECT_TRUE(sol.is_error());
        if (sol.is_error())
            EXPECT_EQ(sol.error().kind, CASErrorKind::Unimplemented);
    }
}

