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
