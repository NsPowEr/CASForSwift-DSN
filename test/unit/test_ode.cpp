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
    EXPECT_EQ(cls.value().type, OdeType::Riccati);
    ASSERT_EQ(cls.value().components.size(), 3U);

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
    EXPECT_EQ(cls.value().type, OdeType::Riccati);

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
    EXPECT_EQ(cls.value().type, OdeType::Riccati);

    auto sol = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(sol.is_ok()) << sol.error().message;
    std::cout << "Riccati Δ>0 Solution: " << debug_print(sol.value()) << std::endl;
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

