#include <gtest/gtest.h>
#include "cas/calculus.hpp"
#include "cas/algebra.hpp"
#include "cas/formatter.hpp"
#include "cas/ast_debug.hpp"
#include <iostream>

using namespace cas;
using namespace cas::calculus;
using namespace cas::algebra;
using namespace cas::symbolic;

namespace {
// Parser mock removed as unused
}

TEST(SupremeStressTest, Test6_Residues) {
    CASContext ctx;
    Symbol x("x");
    // x^2 / (x^4 + 1)
    ExprPtr num = ctx.arena().make<Binary>(BinaryOp::Pow, ctx.arena().make<Symbol>("x"), ctx.arena().make<IntegerLit>(BigInt(2)));
    ExprPtr den = ctx.arena().make<Binary>(BinaryOp::Add, 
        ctx.arena().make<Binary>(BinaryOp::Pow, ctx.arena().make<Symbol>("x"), ctx.arena().make<IntegerLit>(BigInt(4))),
        ctx.arena().make<IntegerLit>(BigInt(1)));
    ExprPtr f = ctx.arena().make<Binary>(BinaryOp::Div, num, den);
    
    ExprPtr inf = ctx.arena().make<Constant>(MathConstant::Infinity);
    ExprPtr n_inf = ctx.arena().make<Unary>(UnaryOp::Neg, inf);
    
    auto res = definite_integral(f, x, n_inf, inf, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    std::cout << "Result Test 6: " << formatter::TextFormatter{}.format(res.value()) << std::endl;
}

TEST(SupremeStressTest, Test7_AsymptoticLimit) {
    CASContext ctx;
    Symbol x("x");
    // limit(x * ((1 + 1/x)^x - e), x, infinity)
    ExprPtr one = ctx.arena().make<IntegerLit>(BigInt(1));
    ExprPtr inv_x = ctx.arena().make<Binary>(BinaryOp::Div, one, ctx.arena().make<Symbol>("x"));
    ExprPtr base = ctx.arena().make<Binary>(BinaryOp::Add, one, inv_x);
    ExprPtr term = ctx.arena().make<Binary>(BinaryOp::Pow, base, ctx.arena().make<Symbol>("x"));
    ExprPtr e_const = ctx.arena().make<Constant>(MathConstant::E);
    ExprPtr diff_expr = ctx.arena().make<Binary>(BinaryOp::Sub, term, e_const);
    ExprPtr f = ctx.arena().make<Binary>(BinaryOp::Mul, ctx.arena().make<Symbol>("x"), diff_expr);
    
    ExprPtr inf = ctx.arena().make<Constant>(MathConstant::Infinity);
    
    auto res = limit(f, x, inf, LimitDirection::Both, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    std::cout << "Result Test 7: " << formatter::TextFormatter{}.format(res.value()) << std::endl;
}

TEST(SupremeStressTest, Test10_CSE_Derivative) {
    CASContext ctx;
    Symbol x("x");
    // diff(x^x, x, 4)
    ExprPtr x_sym = ctx.arena().make<Symbol>("x");
    ExprPtr f = ctx.arena().make<Binary>(BinaryOp::Pow, x_sym, x_sym);
    
    auto res = diff(f, x, 4, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    
    std::string output = formatter::TextFormatter{}.format(res.value());
    std::cout << "Result Test 10 (CSE Output): " << output << std::endl;
    
    // Verifica che l'output contenga 'where' o variabili v1, v2 indicando successo CSE
    EXPECT_TRUE(output.find("where") != std::string::npos || output.find("v1") != std::string::npos);
}
