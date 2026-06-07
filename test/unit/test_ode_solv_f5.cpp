#include <gtest/gtest.h>
#include "cas/ode.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast_debug.hpp"
using namespace cas;
using namespace cas::calculus;

class OdeSolverStrictTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

TEST_F(OdeSolverStrictTest, Separable) {
    auto& arena = ctx.arena();
    Symbol x("x");
    Symbol y("y");

    ExprPtr Mx = arena.make<Symbol>("x"); // M(x) = x
    ExprPtr Ny = arena.make<Symbol>("y"); // N(y) = y

    OdeClassification cls(OdeType::Separable, nullptr, y, x);
    cls.components = {Ny, Mx};

    auto res = solve_ode_1st_order(cls, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    std::string s = debug_print(res.value());
    std::cout << "Separable output: " << s << std::endl;
    EXPECT_TRUE(s.find("equal") != std::string::npos);
}

TEST_F(OdeSolverStrictTest, Exact) {
    auto& arena = ctx.arena();
    Symbol x("x");
    Symbol y("y");

    // M(x,y) = 2xy
    ExprPtr two = arena.make<IntegerLit>(BigInt(2));
    ExprPtr sx = arena.make<Symbol>("x");
    ExprPtr sy = arena.make<Symbol>("y");
    ExprPtr M = arena.make<Binary>(BinaryOp::Mul, arena.make<Binary>(BinaryOp::Mul, two, sx), sy);

    // N(x,y) = x^2 + 3y^2
    ExprPtr x2 = arena.make<Binary>(BinaryOp::Pow, sx, arena.make<IntegerLit>(BigInt(2)));
    ExprPtr three = arena.make<IntegerLit>(BigInt(3));
    ExprPtr y2 = arena.make<Binary>(BinaryOp::Pow, sy, arena.make<IntegerLit>(BigInt(2)));
    ExprPtr N = arena.make<Binary>(BinaryOp::Add, x2, arena.make<Binary>(BinaryOp::Mul, three, y2));

    OdeClassification cls(OdeType::Exact, nullptr, y, x);
    cls.components = {M, N};

    auto res = solve_ode_1st_order(cls, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    std::string s = debug_print(res.value());
    std::cout << "Exact output: " << s << std::endl;
    EXPECT_TRUE(s.find("equal") != std::string::npos);
}

TEST_F(OdeSolverStrictTest, Bernoulli) {
    auto& arena = ctx.arena();
    Symbol x("x");
    Symbol y("y");
    
    ctx.assumptions().assume_positive(x);

    // P(x) = 1/x
    ExprPtr one = arena.make<IntegerLit>(BigInt(1));
    ExprPtr sx = arena.make<Symbol>("x");
    ExprPtr P = arena.make<Binary>(BinaryOp::Div, one, sx);

    // Q(x) = x
    ExprPtr Q = sx;

    // n = 2
    ExprPtr n_expr = arena.make<IntegerLit>(BigInt(2));

    OdeClassification cls(OdeType::Bernoulli, nullptr, y, x);
    cls.components = {P, Q, n_expr};

    auto res = solve_ode_1st_order(cls, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    std::string s = debug_print(res.value());
    std::cout << "Bernoulli output: " << s << std::endl;
    EXPECT_TRUE(s.find("equal") != std::string::npos);
}

TEST_F(OdeSolverStrictTest, Homogeneous) {
    auto& arena = ctx.arena();
    Symbol x("x");
    Symbol y("y");
    
    ctx.assumptions().assume_positive(x);

    // F(v) = 1/v + v
    Symbol v("v");
    ctx.assumptions().assume_positive(v);
    ExprPtr sv = arena.make<Symbol>("v");
    ExprPtr one = arena.make<IntegerLit>(BigInt(1));
    ExprPtr inv_v = arena.make<Binary>(BinaryOp::Div, one, sv);
    ExprPtr F_v = arena.make<Binary>(BinaryOp::Add, inv_v, sv);

    OdeClassification cls(OdeType::Homogeneous, nullptr, y, x);
    cls.components = {F_v};
    cls.parameter = v;

    auto res = solve_ode_1st_order(cls, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;

    std::string s = debug_print(res.value());
    std::cout << "Homogeneous output: " << s << std::endl;
    EXPECT_TRUE(s.find("equal") != std::string::npos);
}
