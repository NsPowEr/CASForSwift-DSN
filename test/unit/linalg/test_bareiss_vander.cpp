#include "cas/linalg/Matrix.hpp"
#include "cas/symbolic.hpp"
#include <gtest/gtest.h>
#include <string>

namespace cas::linalg {
namespace {
ExprPtr integer(symbolic::CASContext& ctx, long long value) {
    return ctx.arena().make<IntegerLit>(BigInt(value));
}

ExprPtr symbol(symbolic::CASContext& ctx, std::string name) {
    return ctx.arena().make<Symbol>(std::move(name));
}

ExprPtr pow_expr(symbolic::CASContext& ctx, ExprPtr base, long long exp) {
    return ctx.arena().make<Binary>(BinaryOp::Pow, base, integer(ctx, exp));
}

[[nodiscard]] bool is_zero_expr(ExprPtr expr) {
    if (const auto* il = expr_cast<IntegerLit>(expr)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(expr)) return rl->numerator.is_zero();
    return false;
}
}

TEST(MatrixBareissTest, Vandermonde4x4) {
    symbolic::CASContext context;
    // Vandermonde 4x4:
    // [1, x0, x0^2, x0^3]
    // [1, x1, x1^2, x1^3]
    // [1, x2, x2^2, x2^3]
    // [1, x3, x3^2, x3^3]
    // Det = (x1-x0)(x2-x0)(x3-x0)(x2-x1)(x3-x1)(x3-x2)
    
    MatrixExpr m(4, 4);
    for (int i = 0; i < 4; ++i) {
        ExprPtr xi = symbol(context, "x" + std::to_string(i));
        for (int j = 0; j < 4; ++j) {
            if (j == 0) m(i, j) = integer(context, 1);
            else if (j == 1) m(i, j) = xi;
            else m(i, j) = pow_expr(context, xi, j);
        }
    }
    
    auto det = determinant(m, context);
    ASSERT_TRUE(det.is_ok()) << det.error().message;
    
    // Check it's not zero
    EXPECT_FALSE(is_zero_expr(det.value()));
    
    // Evaluate at specific points to verify
    // x0=0, x1=1, x2=2, x3=3 -> Det=12
    for (int i = 0; i < 4; ++i) {
        context.define(Symbol("x" + std::to_string(i)), integer(context, i));
    }
    
    // Wait, context.define defines symbols for lookup during simplification or substitution.
    // I should use substitute.
    ExprPtr result = det.value();
    for (int i = 0; i < 4; ++i) {
        auto sub_res = context.substitute(result, Symbol("x" + std::to_string(i)), integer(context, i));
        ASSERT_TRUE(sub_res.is_ok());
        result = sub_res.value();
    }
    
    auto final_simp = context.simplify(result);
    ASSERT_TRUE(final_simp.is_ok());
    
    const auto* res_int = expr_cast<IntegerLit>(final_simp.value());
    ASSERT_NE(res_int, nullptr);
    EXPECT_EQ(res_int->value, BigInt(12));
}
}
