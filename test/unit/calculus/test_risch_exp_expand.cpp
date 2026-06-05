// F5.1 / Task #23 — Tests cap.9 sub-case exp-decomposition.
//
// Verifica `expand_exp_args_via_decomposition` come pre-processore per il
// wiring cap.8 esponenziale: l'integrale di exp(2x) viene riconosciuto
// come polinomio in t = exp(x) dopo decomposizione.

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include "../../../src/calculus/calculus_internal.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas::calculus {
namespace {

ExprPtr parse_expr(const std::string& src, AstArena& arena) {
    auto t = Lexer(src).tokenize();
    EXPECT_TRUE(t.is_ok()) << src;
    Parser p(t.value(), arena);
    auto r = p.parse();
    EXPECT_TRUE(r.is_ok()) << src;
    return r.value();
}

class RischExpExpandTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};
};

// exp(2x) → exp(x)^2 dopo decomposizione.
TEST_F(RischExpExpandTest, ExpOfTwoXDecomposes) {
    auto e = parse_expr("exp(2*x)", ctx.arena());
    ExprPtr decomp = expand_exp_args_via_decomposition(e, x, ctx);

    // Atteso: Pow(exp(x), 2).
    const auto* bin = expr_cast<Binary>(decomp);
    ASSERT_NE(bin, nullptr) << "decomposed form must be Binary(Pow,...)";
    EXPECT_EQ(bin->op, BinaryOp::Pow);
    const auto* base = expr_cast<FuncCall>(bin->left);
    ASSERT_NE(base, nullptr);
    EXPECT_EQ(base->func_id, BuiltinOp::Exp);
    const auto* il = expr_cast<IntegerLit>(bin->right);
    ASSERT_NE(il, nullptr);
    EXPECT_EQ(il->value, BigInt(2));
}

// exp(x + y) → exp(x) * exp(y).
TEST_F(RischExpExpandTest, ExpOfSumDecomposes) {
    auto e = parse_expr("exp(x + y)", ctx.arena());
    ExprPtr decomp = expand_exp_args_via_decomposition(e, x, ctx);
    const auto* prod = expr_cast<Product>(decomp);
    ASSERT_NE(prod, nullptr) << "expected Product after Sum decomposition";
    EXPECT_EQ(prod->factors.size(), 2U);
    int exp_count = 0;
    for (ExprPtr f : prod->factors) {
        if (const auto* fc = expr_cast<FuncCall>(f); fc && fc->func_id == BuiltinOp::Exp) ++exp_count;
    }
    EXPECT_EQ(exp_count, 2);
}

// exp(x) inalterato (coeff ±1, no scomposizione).
TEST_F(RischExpExpandTest, ExpOfXUnchanged) {
    auto e = parse_expr("exp(x)", ctx.arena());
    ExprPtr decomp = expand_exp_args_via_decomposition(e, x, ctx);
    const auto* fc = expr_cast<FuncCall>(decomp);
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->func_id, BuiltinOp::Exp);
}

// Walk ricorsivo: somma con exp(2x) interno viene espansa.
TEST_F(RischExpExpandTest, RecursiveWalkInsideSum) {
    auto e = parse_expr("exp(x) + exp(3*x)", ctx.arena());
    ExprPtr decomp = expand_exp_args_via_decomposition(e, x, ctx);
    const auto* sum = expr_cast<Sum>(decomp);
    ASSERT_NE(sum, nullptr);
    bool found_pow = false;
    for (ExprPtr t : sum->terms) {
        if (const auto* bin = expr_cast<Binary>(t); bin && bin->op == BinaryOp::Pow) {
            if (const auto* base = expr_cast<FuncCall>(bin->left);
                base && base->func_id == BuiltinOp::Exp) {
                found_pow = true;
            }
        }
    }
    EXPECT_TRUE(found_pow) << "exp(3*x) should decompose to exp(x)^3";
}

// Integrazione end-to-end: ∫ exp(2x) dx = (1/2)·exp(2x).
// Validation: la derivata della primitiva deve uguagliare exp(2x).
TEST_F(RischExpExpandTest, IntegrateExp2xEndToEnd) {
    auto e = parse_expr("exp(2*x)", ctx.arena());
    auto integ = integrate(e, x, ctx);
    ASSERT_TRUE(integ.is_ok()) << integ.error().message;

    auto d = diff(integ.value(), x, 1U, ctx);
    ASSERT_TRUE(d.is_ok());
    ExprPtr delta = ctx.arena().make<Binary>(BinaryOp::Sub, d.value(), e);
    auto tog = algebra::together(delta, ctx);
    ASSERT_TRUE(tog.is_ok());
    auto simp = ctx.simplify(tog.value());
    ASSERT_TRUE(simp.is_ok());
    bool is_zero = false;
    if (const auto* il = expr_cast<IntegerLit>(simp.value())) is_zero = il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(simp.value())) is_zero = rl->numerator.is_zero();
    EXPECT_TRUE(is_zero) << "d/dx(integral) should equal exp(2x)";
}

}  // namespace
}  // namespace cas::calculus
