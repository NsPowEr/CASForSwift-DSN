// CAS-F4.2d — Test companion matrix.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::linalg;

namespace {

class CompanionTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    [[nodiscard]] ExprPtr lit(long long v) {
        return ctx.arena().make<IntegerLit>(BigInt(v));
    }
    [[nodiscard]] Symbol sym(const std::string& name) { return Symbol(name); }
    [[nodiscard]] bool entries_equal(ExprPtr a, ExprPtr b) {
        auto eq = symbolic::mathematically_equal(a, b, ctx);
        return eq.is_ok() && eq.value();
    }
};

TEST_F(CompanionTest, MonicDeg2_x2_minus_3x_plus_2) {
    // p(x) = x^2 − 3x + 2.  Coeffs [2, -3, 1].
    // Companion: [[0, -2],[1, 3]]
    Symbol x = sym("xc");
    ExprPtr x_sym = ctx.arena().make<Symbol>(x.name);
    ExprPtr x2 = ctx.arena().make<Binary>(BinaryOp::Pow, x_sym, lit(2));
    ExprPtr term1 = ctx.arena().make<Binary>(BinaryOp::Mul, lit(3), x_sym);
    ExprPtr p = ctx.arena().make<Sum>(std::vector<ExprPtr>{
        x2,
        ctx.arena().make<Unary>(UnaryOp::Neg, term1),
        lit(2)
    });
    auto C = companion_matrix(p, x, ctx);
    ASSERT_TRUE(C.is_ok()) << C.error().message;
    EXPECT_EQ(C.value().rows(), 2U);
    EXPECT_EQ(C.value().cols(), 2U);
    EXPECT_TRUE(entries_equal(C.value()(0, 0), lit(0)));
    EXPECT_TRUE(entries_equal(C.value()(0, 1), lit(-2)));
    EXPECT_TRUE(entries_equal(C.value()(1, 0), lit(1)));
    EXPECT_TRUE(entries_equal(C.value()(1, 1), lit(3)));
}

TEST_F(CompanionTest, Certificator_det_lambdaI_minus_C_equals_p) {
    // Per p(x) = x^3 − 6x^2 + 11x − 6 (radici 1,2,3), det(λI − C) ≡ p(λ).
    Symbol x = sym("xc");
    ExprPtr x_sym = ctx.arena().make<Symbol>(x.name);
    ExprPtr x3 = ctx.arena().make<Binary>(BinaryOp::Pow, x_sym, lit(3));
    ExprPtr x2 = ctx.arena().make<Binary>(BinaryOp::Pow, x_sym, lit(2));
    ExprPtr term2 = ctx.arena().make<Binary>(BinaryOp::Mul, lit(6), x2);
    ExprPtr term1 = ctx.arena().make<Binary>(BinaryOp::Mul, lit(11), x_sym);
    ExprPtr p = ctx.arena().make<Sum>(std::vector<ExprPtr>{
        x3,
        ctx.arena().make<Unary>(UnaryOp::Neg, term2),
        term1,
        ctx.arena().make<Unary>(UnaryOp::Neg, lit(6))
    });
    auto C = companion_matrix(p, x, ctx);
    ASSERT_TRUE(C.is_ok()) << C.error().message;

    // characteristic_polynomial in CAS computes det(A - λI). Standard companion
    // matrix soddisfa det(λI - C) = p(λ), quindi det(C - λI) = (-1)^n · p(λ).
    // Per n=3 dispari, cp atteso = -p.
    auto cp = characteristic_polynomial(C.value(), x, ctx);
    ASSERT_TRUE(cp.is_ok());

    auto cp_coeffs = algebra::univariate_coefficients(cp.value(), x, ctx);
    auto p_coeffs = algebra::univariate_coefficients(p, x, ctx);
    ASSERT_TRUE(cp_coeffs.is_ok());
    ASSERT_TRUE(p_coeffs.is_ok());
    ASSERT_EQ(cp_coeffs.value().size(), p_coeffs.value().size());
    const std::size_t n = p_coeffs.value().size() - 1U;
    const int sign = (n % 2U == 0U) ? 1 : -1;
    for (std::size_t k = 0; k < p_coeffs.value().size(); ++k) {
        ExprPtr expected = (sign == 1) ? p_coeffs.value()[k]
            : static_cast<ExprPtr>(ctx.arena().make<Unary>(UnaryOp::Neg, p_coeffs.value()[k]));
        EXPECT_TRUE(entries_equal(cp_coeffs.value()[k], expected))
            << "coefficient mismatch at x^" << k;
    }
}

TEST_F(CompanionTest, NonMonicNormalized_2x2_minus_8x_plus_6) {
    // p(x) = 2x^2 − 8x + 6 → normalizzato: x^2 − 4x + 3. Companion [[0,-3],[1,4]].
    Symbol x = sym("xc");
    ExprPtr x_sym = ctx.arena().make<Symbol>(x.name);
    ExprPtr x2 = ctx.arena().make<Binary>(BinaryOp::Pow, x_sym, lit(2));
    ExprPtr t1 = ctx.arena().make<Binary>(BinaryOp::Mul, lit(2), x2);
    ExprPtr t2 = ctx.arena().make<Binary>(BinaryOp::Mul, lit(8), x_sym);
    ExprPtr p = ctx.arena().make<Sum>(std::vector<ExprPtr>{
        t1,
        ctx.arena().make<Unary>(UnaryOp::Neg, t2),
        lit(6)
    });
    auto C = companion_matrix(p, x, ctx);
    ASSERT_TRUE(C.is_ok());
    EXPECT_TRUE(entries_equal(C.value()(0, 1), lit(-3)));
    EXPECT_TRUE(entries_equal(C.value()(1, 1), lit(4)));
}

TEST_F(CompanionTest, Deg0Rejected) {
    auto C = companion_matrix(ctx.arena().make<IntegerLit>(BigInt(5)),
                              sym("xc"), ctx);
    EXPECT_TRUE(C.is_error());
}

}  // namespace
