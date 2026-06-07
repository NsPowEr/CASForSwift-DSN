#include <gtest/gtest.h>
#include "cas/ode.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast_debug.hpp"
#include "cas/algebra.hpp"
using namespace cas;
using namespace cas::calculus;

class OdeClassifierStrictTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    ExprPtr parse(const std::string& /*expr*/) {
        // Dummy parse, we'll build AST manually.
        return nullptr;
    }
};

TEST_F(OdeClassifierStrictTest, Separable) {
    auto y = Symbol("y");
    auto x = Symbol("x");
    AstArena& arena = ctx.arena();
    auto y_sym = arena.make<Symbol>(y.name);
    auto x_sym = arena.make<Symbol>(x.name);
    auto y_prime = arena.make<Derivative>(y_sym, x, 1);
    
    // y * y' - x = 0  => N(y)=y, M(x)=x
    auto eq = arena.make<Binary>(BinaryOp::Sub,
        arena.make<Binary>(BinaryOp::Mul, y_sym, y_prime),
        x_sym);
        
    auto res = classify_ode(eq, y, x, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_EQ(res.value().type, OdeType::Separable);
    EXPECT_EQ(res.value().components.size(), 2);
}

TEST_F(OdeClassifierStrictTest, Exact) {
    auto y = Symbol("y");
    auto x = Symbol("x");
    AstArena& arena = ctx.arena();
    auto y_sym = arena.make<Symbol>(y.name);
    auto x_sym = arena.make<Symbol>(x.name);
    auto y_prime = arena.make<Derivative>(y_sym, x, 1);
    
    // (2xy + x^2) + (x^2 - y)y' = 0
    auto M = arena.make<Binary>(BinaryOp::Add,
        arena.make<Binary>(BinaryOp::Mul, 
            arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(2)), x_sym), y_sym),
        arena.make<Binary>(BinaryOp::Pow, x_sym, arena.make<IntegerLit>(BigInt(2))));
    auto N = arena.make<Binary>(BinaryOp::Sub,
        arena.make<Binary>(BinaryOp::Pow, x_sym, arena.make<IntegerLit>(BigInt(2))),
        y_sym);
    
    auto eq = arena.make<Binary>(BinaryOp::Add, M, arena.make<Binary>(BinaryOp::Mul, N, y_prime));
        
    auto res = classify_ode(eq, y, x, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_EQ(res.value().type, OdeType::Exact);
    EXPECT_EQ(res.value().components.size(), 2);
}

TEST_F(OdeClassifierStrictTest, Bernoulli) {
    auto y = Symbol("y");
    auto x = Symbol("x");
    AstArena& arena = ctx.arena();
    auto y_sym = arena.make<Symbol>(y.name);
    auto x_sym = arena.make<Symbol>(x.name);
    auto y_prime = arena.make<Derivative>(y_sym, x, 1);
    
    // y' + x y - x^2 y^3 = 0  => P = x, Q = x^2, n = 3
    auto Py = arena.make<Binary>(BinaryOp::Mul, x_sym, y_sym);
    auto Qy3 = arena.make<Binary>(BinaryOp::Mul,
        arena.make<Binary>(BinaryOp::Pow, x_sym, arena.make<IntegerLit>(BigInt(2))),
        arena.make<Binary>(BinaryOp::Pow, y_sym, arena.make<IntegerLit>(BigInt(3))));
    
    auto eq = arena.make<Binary>(BinaryOp::Sub,
        arena.make<Binary>(BinaryOp::Add, y_prime, Py),
        Qy3);
        
    auto res = classify_ode(eq, y, x, ctx);
    ASSERT_TRUE(res.is_ok());
    EXPECT_EQ(res.value().type, OdeType::Bernoulli);
    EXPECT_EQ(res.value().components.size(), 3);
}

TEST_F(OdeClassifierStrictTest, Homogeneous) {
    auto y = Symbol("y");
    auto x = Symbol("x");
    AstArena& arena = ctx.arena();
    auto y_sym = arena.make<Symbol>(y.name);
    auto x_sym = arena.make<Symbol>(x.name);
    auto y_prime = arena.make<Derivative>(y_sym, x, 1);
    
    // x y' - y = 0
    auto eq = arena.make<Binary>(BinaryOp::Sub,
        arena.make<Binary>(BinaryOp::Mul, x_sym, y_prime),
        y_sym);
        
    auto res = classify_ode(eq, y, x, ctx);
    ASSERT_TRUE(res.is_ok());
    // In this specific case, it might also be separable, but we check homogenous logic
    // Depending on order, separable is checked first. Let's make it not separable.
    // Wait, x y' - y = 0 is separable.
    // xy y' - (y^2 + x^2) = 0
    auto A = arena.make<Binary>(BinaryOp::Mul, x_sym, y_sym);
    auto y2 = arena.make<Binary>(BinaryOp::Pow, y_sym, arena.make<IntegerLit>(BigInt(2)));
    auto x2 = arena.make<Binary>(BinaryOp::Pow, x_sym, arena.make<IntegerLit>(BigInt(2)));
    auto B = arena.make<Unary>(UnaryOp::Neg, arena.make<Binary>(BinaryOp::Add, y2, x2));
    
    auto eq2 = arena.make<Binary>(BinaryOp::Add,
        arena.make<Binary>(BinaryOp::Mul, A, y_prime),
        B);
        
    auto res2 = classify_ode(eq2, y, x, ctx);
    ASSERT_TRUE(res2.is_ok());
    EXPECT_EQ(res2.value().type, OdeType::Homogeneous);
    EXPECT_EQ(res2.value().components.size(), 1);
}
