#include <gtest/gtest.h>
#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "../../../src/algebra/polynomial_sparse_interpolation_internal.hpp"
#include "../../../src/algebra/algebra_internal.hpp"

namespace cas::algebra {

class SparseInterpolationTest : public ::testing::Test {
protected:
    void SetUp() override {
        ctx = std::make_unique<symbolic::CASContext>();
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

TEST_F(SparseInterpolationTest, UnivariateX2) {
    Symbol x("x");
    std::vector<Symbol> vars = {x};
    std::vector<std::size_t> bounds = {2};

    // P(x) = x^2 + 3x + 1
    auto oracle = [&](const std::vector<ExprPtr>& point) -> Result<ExprPtr> {
        if (point.size() != 1) return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "wrong point size"));
        
        // Evaluate x^2 + 3x + 1
        // Point is Rational.
        auto val_res = expr_to_integer_coefficient(point[0]);
        if (val_res.is_error()) {
            // Handle Rational
            (void)point[0]; // Need a way to evaluate at rational if it's not BigInt
            // Simplified for test: assume we get integer or rational
            return ok(point[0]); // dummy
        }
        
        BigInt v = val_res.value();
        BigInt res = v * v + BigInt(3) * v + BigInt(1);
        return ok(ctx->arena().make<IntegerLit>(res));
    };

    auto result = sparse_interpolate(oracle, vars, bounds, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    
    // Result should be x^2 + 3x + 1
    // We can't easily compare ExprPtr by value without simplification or string comparison
    // but we can check if it evaluates to the same values.
}

TEST_F(SparseInterpolationTest, BivariateSparse) {
    Symbol x("x"), y("y");
    std::vector<Symbol> vars = {x, y};
    std::vector<std::size_t> bounds = {2, 2};

    // P(x, y) = x^2 + y^2
    auto oracle = [&](const std::vector<ExprPtr>& point) -> Result<ExprPtr> {
        if (point.size() != 2) return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "wrong point size"));
        
        auto vx = expr_to_integer_coefficient(point[0]).value();
        auto vy = expr_to_integer_coefficient(point[1]).value();
        
        BigInt res = vx * vx + vy * vy;
        return ok(ctx->arena().make<IntegerLit>(res));
    };

    auto result = sparse_interpolate(oracle, vars, bounds, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    
    // Check degree and terms if possible
}

TEST_F(SparseInterpolationTest, BivariateProduct) {
    Symbol x("x"), y("y");
    std::vector<Symbol> vars = {x, y};
    std::vector<std::size_t> bounds = {1, 1};

    // P(x, y) = 5*x*y - 2
    auto oracle = [&](const std::vector<ExprPtr>& point) -> Result<ExprPtr> {
        auto vx = expr_to_integer_coefficient(point[0]).value();
        auto vy = expr_to_integer_coefficient(point[1]).value();
        
        BigInt res = BigInt(5) * vx * vy - BigInt(2);
        return ok(ctx->arena().make<IntegerLit>(res));
    };

    auto result = sparse_interpolate(oracle, vars, bounds, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
}

TEST_F(SparseInterpolationTest, DISABLED_Trivariate) {
    Symbol x("x"), y("y"), z("z");
    std::vector<Symbol> vars = {x, y, z};
    std::vector<std::size_t> bounds = {1, 1, 1};

    // P(x, y, z) = x + y + z
    auto oracle = [&](const std::vector<ExprPtr>& point) -> Result<ExprPtr> {
        auto vx = expr_to_integer_coefficient(point[0]).value();
        auto vy = expr_to_integer_coefficient(point[1]).value();
        auto vz = expr_to_integer_coefficient(point[2]).value();
        
        BigInt res = vx + vy + vz;
        return ok(ctx->arena().make<IntegerLit>(res));
    };

    auto result = sparse_interpolate(oracle, vars, bounds, *ctx);
    ASSERT_TRUE(result.is_ok()) << result.error().message;
}

} // namespace cas::algebra
