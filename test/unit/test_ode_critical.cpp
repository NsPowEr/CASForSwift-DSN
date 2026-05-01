#include <gtest/gtest.h>
#include "cas/ode.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast_debug.hpp"
#include "cas/algebra.hpp"
#include <iostream>

using namespace cas;
using namespace cas::calculus;

class OdeCriticalTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

TEST_F(OdeCriticalTest, NegativeSignHandling) {
    std::cout << "TEST_START: NegativeSignHandling" << std::endl;
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();
    
    ExprPtr y_pp = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 2);
    ExprPtr y_s = arena.make<Symbol>("y");
    
    ExprPtr eq_sub = arena.make<Binary>(BinaryOp::Sub, y_pp, y_s);
    auto res_sub = classify_ode(eq_sub, y, x, ctx);
    ASSERT_TRUE(res_sub.is_ok());
    EXPECT_EQ(res_sub.value().type, OdeType::Linear2ndOrderConstantCoeff);
    std::cout << "a0 (Sub): " << debug_print(res_sub.value().components[2]) << std::endl;
}

TEST_F(OdeCriticalTest, ComplexCoefficients) {
    std::cout << "TEST_START: ComplexCoefficients_NEW_VERSION" << std::endl;
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();
    
    // 7*y'' + 5*y' + 11*y = 0
    ExprPtr y_pp = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 2);
    ExprPtr y_p = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 1);
    ExprPtr y_s = arena.make<Symbol>("y");
    
    ExprPtr term2 = arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(7)), y_pp);
    ExprPtr term1 = arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(5)), y_p);
    ExprPtr term0 = arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(11)), y_s);
    
    ExprPtr eq = arena.make<Sum>(std::vector<ExprPtr>{term2, term1, term0});
    
    auto res = classify_ode(eq, y, x, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_EQ(res.value().type, OdeType::Linear2ndOrderConstantCoeff);
    std::cout << "a2: " << debug_print(res.value().components[0]) << std::endl;
    std::cout << "a1: " << debug_print(res.value().components[1]) << std::endl;
    std::cout << "a0: " << debug_print(res.value().components[2]) << std::endl;
}

TEST_F(OdeCriticalTest, SolveHomogeneous) {
    std::cout << "TEST_START: SolveHomogeneous" << std::endl;
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();
    ExprPtr y_pp = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 2);
    ExprPtr y_s = arena.make<Symbol>("y");
    // y'' - y = 0
    ExprPtr eq = arena.make<Binary>(BinaryOp::Equal, 
        arena.make<Binary>(BinaryOp::Sub, y_pp, y_s),
        arena.make<IntegerLit>(BigInt(0)));
    
    auto sol = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(sol.is_ok());
    std::cout << "Solution: " << debug_print(sol.value()) << std::endl;
    // Expected something like y = C1*exp(x) + C2*exp(-x)
}

TEST_F(OdeCriticalTest, SolveNonHomogeneous) {
    std::cout << "TEST_START: SolveNonHomogeneous" << std::endl;
    Symbol y("y");
    Symbol x("x");
    AstArena& arena = ctx.arena();
    
    // y'' + 3*y' + 2*y = exp(x)
    ExprPtr y_pp = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 2);
    ExprPtr y_p = arena.make<Derivative>(arena.make<Symbol>("y"), Symbol("x"), 1);
    ExprPtr y_s = arena.make<Symbol>("y");
    
    ExprPtr term2 = y_pp;
    ExprPtr term1 = arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(3)), y_p);
    ExprPtr term0 = arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(2)), y_s);
    ExprPtr f = arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{arena.make<Symbol>("x")});
    
    ExprPtr eq = arena.make<Binary>(BinaryOp::Equal, 
        arena.make<Sum>(std::vector<ExprPtr>{term2, term1, term0}),
        f);
    
    auto sol = solve_ode(eq, y, x, ctx);
    ASSERT_TRUE(sol.is_ok());
    std::cout << "Solution: " << debug_print(sol.value()) << std::endl;
    // Expected y = C1*exp(-x) + C2*exp(-2*x) + (1/6)*exp(x)
}
