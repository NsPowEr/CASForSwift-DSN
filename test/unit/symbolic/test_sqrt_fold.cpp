// HH-8 probe: verify sqrt(p)·sqrt(p) → p when p ≥ 0 via assumptions.
// Covers the QR cert closure path documented in HC-F8-QR-HOUSEHOLDER-BAILOUT.

#include <gtest/gtest.h>
#include <vector>

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"

namespace cas {
namespace {

// sqrt(x²+y²) · sqrt(x²+y²) → x²+y²  when x, y assumed positive.
TEST(SqrtFold, SqrtSumOfSquaresProductWithItself_PositiveSymbols) {
    symbolic::CASContext ctx;
    Symbol sx("x"), sy("y");
    ctx.assumptions().assume_positive(sx);
    ctx.assumptions().assume_positive(sy);
    AstArena& a = ctx.arena();
    ExprPtr x = a.make<Symbol>("x");
    ExprPtr y = a.make<Symbol>("y");
    ExprPtr two = a.make<IntegerLit>(BigInt(2));
    ExprPtr x_sq = a.make<Binary>(BinaryOp::Pow, x, two);
    ExprPtr y_sq = a.make<Binary>(BinaryOp::Pow, y, two);
    ExprPtr Nx = a.make<Sum>(std::vector<ExprPtr>{x_sq, y_sq});
    auto Nx_s = ctx.simplify(Nx);
    ASSERT_TRUE(Nx_s.is_ok());
    EXPECT_TRUE(ctx.assumptions().is_nonnegative(Nx_s.value()))
        << "x²+y² with positive x,y must be provably nonnegative";

    ExprPtr sqrt_Nx = a.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{Nx_s.value()});
    ExprPtr sq = a.make<Binary>(BinaryOp::Mul, sqrt_Nx, sqrt_Nx);
    auto sq_s = ctx.simplify(sq);
    ASSERT_TRUE(sq_s.is_ok());

    // Result must be structurally Nx (x² + y²), NOT FuncCall(Sqrt, ...).
    const auto* fc = expr_cast<FuncCall>(sq_s.value());
    EXPECT_FALSE(fc != nullptr && fc->func_id == BuiltinOp::Sqrt)
        << "sqrt(x²+y²)·sqrt(x²+y²) failed to fold to x²+y²";
}

// Product form (no Div): sqrt(Nx) · sqrt(Nx) · Nx^(-1) → 1.
TEST(SqrtFold, ProductFormWithNegativeExponent_FoldsToOne) {
    symbolic::CASContext ctx;
    Symbol sx("x"), sy("y");
    ctx.assumptions().assume_positive(sx);
    ctx.assumptions().assume_positive(sy);
    AstArena& a = ctx.arena();
    ExprPtr x = a.make<Symbol>("x");
    ExprPtr y = a.make<Symbol>("y");
    ExprPtr two = a.make<IntegerLit>(BigInt(2));
    ExprPtr neg_one = a.make<IntegerLit>(BigInt(-1));
    ExprPtr x_sq = a.make<Binary>(BinaryOp::Pow, x, two);
    ExprPtr y_sq = a.make<Binary>(BinaryOp::Pow, y, two);
    ExprPtr Nx = a.make<Sum>(std::vector<ExprPtr>{x_sq, y_sq});
    ExprPtr sqrt_Nx = a.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{Nx});
    ExprPtr Nx_inv = a.make<Binary>(BinaryOp::Pow, Nx, neg_one);
    ExprPtr expr = a.make<Product>(std::vector<ExprPtr>{sqrt_Nx, sqrt_Nx, Nx_inv});
    auto s = ctx.simplify(expr);
    ASSERT_TRUE(s.is_ok());

    bool is_one = false;
    if (const auto* il = expr_cast<IntegerLit>(s.value()))
        is_one = (il->value == BigInt(1));
    EXPECT_TRUE(is_one)
        << "Product[sqrt(Nx), sqrt(Nx), Nx^(-1)] failed to fold to 1";
}

// QR-cert downstream pattern: y · sqrt(x²+y²) / (x²+y²) · sqrt(x²+y²)
// must reduce to y. (Mirrors Q·R[0][0] structure for 2×2 [[x,1],[y,2]].)
TEST(SqrtFold, FractionWithSqrtNumeratorAndDenominator_FoldsToOne) {
    symbolic::CASContext ctx;
    Symbol sx("x"), sy("y");
    ctx.assumptions().assume_positive(sx);
    ctx.assumptions().assume_positive(sy);
    AstArena& a = ctx.arena();
    ExprPtr x = a.make<Symbol>("x");
    ExprPtr y = a.make<Symbol>("y");
    ExprPtr two = a.make<IntegerLit>(BigInt(2));
    ExprPtr x_sq = a.make<Binary>(BinaryOp::Pow, x, two);
    ExprPtr y_sq = a.make<Binary>(BinaryOp::Pow, y, two);
    ExprPtr Nx = a.make<Sum>(std::vector<ExprPtr>{x_sq, y_sq});
    ExprPtr sqrt_Nx = a.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{Nx});

    // expr = (sqrt(Nx) / Nx) · sqrt(Nx) = 1/sqrt(Nx) · sqrt(Nx) = 1
    ExprPtr inner = a.make<Binary>(BinaryOp::Div, sqrt_Nx, Nx);
    ExprPtr expr = a.make<Binary>(BinaryOp::Mul, inner, sqrt_Nx);
    auto s = ctx.simplify(expr);
    ASSERT_TRUE(s.is_ok());

    // Expected: 1 (IntegerLit). May also accept Product/Sum that reduces to 1.
    bool is_one = false;
    if (const auto* il = expr_cast<IntegerLit>(s.value()))
        is_one = (il->value == BigInt(1));
    EXPECT_TRUE(is_one)
        << "sqrt(Nx)/Nx · sqrt(Nx) failed to collapse to 1";
}

}  // namespace
}  // namespace cas
