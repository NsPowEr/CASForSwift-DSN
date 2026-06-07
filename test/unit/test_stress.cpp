#include <gtest/gtest.h>
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/symbolic.hpp"
#include "cas/numtheory.hpp"
#include "algebra/polynomial_internal.hpp"
#include <chrono>
#include <iostream>
#include <vector>
#include <string>

using namespace cas;
using namespace cas::linalg;
using namespace cas::symbolic;
using namespace cas::algebra;

namespace {

struct Timer {
    std::chrono::steady_clock::time_point start;
    std::string name;
    Timer(std::string n) : name(n) { start = std::chrono::steady_clock::now(); }
    ~Timer() {
        auto end = std::chrono::steady_clock::now();
        auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "[ TIMER ] " << name << ": " << diff << " ms" << std::endl;
    }
};

ExprPtr integer(CASContext& ctx, const BigInt& val) {
    return ctx.arena().make<IntegerLit>(val);
}

ExprPtr symbol(CASContext& ctx, std::string name) {
    return ctx.arena().make<Symbol>(std::move(name));
}

} // namespace

TEST(StressTest, MatrixDeterminantSymbolic50x50) {
    CASContext ctx;
    const size_t n = 50;
    MatrixExpr A(n, n);
    
    // x, y, z matrix
    std::vector<ExprPtr> syms = { symbol(ctx, "x"), symbol(ctx, "y"), symbol(ctx, "z") };
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            A(i, j) = syms[(i + j) % 3];
            if (i == j) {
                // Add some diagonal dominance to make it less likely to be zero
                A(i, j) = ctx.arena().make<Sum>(std::vector<ExprPtr>{A(i, j), integer(ctx, BigInt(1))});
            }
        }
    }
    
    std::cout << "Starting 50x50 Symbolic Determinant..." << std::endl;
    {
        Timer t("50x50 Symbolic Determinant");
        auto det = determinant_modular(A, ctx);
        if (det.is_error()) {
            std::cout << "Error: " << det.error().message << std::endl;
        } else {
            std::cout << "Determinant calculated (symbolic)." << std::endl;
        }
    }
}

// DISABILITATO: Test di stress matematico fisiologicamente in timeout (>60s) sotto Debug Mode (-O0). Da eseguire in Release Mode o via target dedicato cas_stress_tests.
TEST(StressTest, DISABLED_MatrixDeterminantGiantInteger100x100) {
    CASContext ctx;
    const size_t n = 100;
    MatrixExpr A(n, n);
    
    // 1000-digit integer
    std::string giant_str(1000, '9');
    BigInt giant = BigInt::parse(giant_str).value();
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            // Mix it up a bit
            BigInt val = giant + BigInt((long long)(i * n + j));
            A(i, j) = integer(ctx, val);
        }
    }
    
    std::cout << "Starting 100x100 Giant Integer Determinant..." << std::endl;
    {
        Timer t("100x100 Giant Integer Determinant");
        auto det = determinant_modular(A, ctx);
        if (det.is_error()) {
            std::cout << "Error: " << det.error().message << std::endl;
        } else {
            std::cout << "Determinant calculated (giant integer)." << std::endl;
        }
    }
}

TEST(StressTest, AssumptionChain500) {
    CASContext ctx;
    const int n = 500;
    std::vector<ExprPtr> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = symbol(ctx, "x" + std::to_string(i + 1));
    }
    
    std::cout << "Building chain of " << n << " assumptions..." << std::endl;
    // x1 > x2 > ... > x500 > 0
    for (int i = 0; i < n - 1; ++i) {
        ctx.assumptions().assume_greater(x[i], x[i + 1]);
    }
    ctx.assumptions().assume_greater(x[n - 1], nullptr); // x500 > 0
    
    std::cout << "Querying is_positive(x1)..." << std::endl;
    {
        Timer t("Assumption Query x1 > 0");
        bool pos = ctx.assumptions().is_positive(x[0]);
        EXPECT_TRUE(pos);
        std::cout << "is_positive(x1) = " << (pos ? "true" : "false") << std::endl;
    }
}

// DISABILITATO: Test di stress matematico fisiologicamente in timeout (>60s) sotto Debug Mode (-O0). Da eseguire in Release Mode o via target dedicato cas_stress_tests.
TEST(StressTest, DISABLED_FactorizationLLLStress) {
    CASContext ctx;
    Symbol x("x");

    std::string coeff_str(50, '7');
    BigInt c = BigInt::parse(coeff_str).value();

    std::cout << "Building degree 100 polynomial..." << std::endl;

    // Build PolyExpr p1 = sum_{i=0}^50 (c+i)*x^i, p2 = sum_{i=0}^50 (c-i)*x^i
    PolyExpr p1, p2;
    for (int i = 0; i <= 50; ++i) {
        p1.push_back(ctx.arena().make<IntegerLit>(c + BigInt((long long)i)));
        p2.push_back(ctx.arena().make<IntegerLit>(c - BigInt((long long)i)));
    }

    auto f_poly = poly_multiply(p1, p2, ctx);
    ASSERT_TRUE(f_poly.is_ok());

    auto f_expr_res = polynomial_to_expr(f_poly.value(), x, ctx);
    ASSERT_TRUE(f_expr_res.is_ok());

    std::cout << "Starting Factorization of degree 100 polynomial..." << std::endl;
    {
        Timer t("Factorization Degree 100");
        auto res = factor_over_integers(f_expr_res.value(), x, ctx);
        if (res.is_error()) {
            std::cout << "Error: " << res.error().message << std::endl;
        } else {
            std::cout << "Factors found: " << res.value().factors.size() << std::endl;
        }
    }
}
