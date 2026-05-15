// Validates eigenvalue + eigenvector computation for matrices whose
// characteristic polynomial has irrational roots represented as RootOf.
// Updated from the previous tautological version (which accepted 0
// eigenvectors): we now verify A*v = lambda*v explicitly using the
// algebraic-extension reduction path.

#include "cas/algebraic_number_bridge.hpp"
#include "cas/ast_debug.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

using namespace cas;
using namespace cas::linalg;

TEST(EigenTest, EigenvaluesDimension4) {
    symbolic::CASContext ctx;
    AstArena& arena = ctx.arena();

    // Companion matrix for lambda^4 - 2 = 0
    // [0 0 0 2]
    // [1 0 0 0]
    // [0 1 0 0]
    // [0 0 1 0]
    MatrixExpr M(4, 4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            M(i, j) = arena.make<IntegerLit>(BigInt(0));
        }
    }
    M(0, 3) = arena.make<IntegerLit>(BigInt(2));
    M(1, 0) = arena.make<IntegerLit>(BigInt(1));
    M(2, 1) = arena.make<IntegerLit>(BigInt(1));
    M(3, 2) = arena.make<IntegerLit>(BigInt(1));

    // Eigenvalues are computed as RootOf instances since n > 3.
    auto ev_res = eigenvalues(M, ctx);
    ASSERT_TRUE(ev_res.is_ok()) << "Eigenvalues failed: " << ev_res.error().message;
    auto ev = ev_res.value();
    EXPECT_EQ(ev.size(), 4U);

    int root_of_count = 0;
    for (auto val : ev) {
        if (expr_is<RootOf>(val)) ++root_of_count;
    }
    EXPECT_EQ(root_of_count, 4) << "Expected all 4 eigenvalues to be represented as RootOf.";

    // Now compute eigenvectors and verify A*v = lambda*v over Q(lambda).
    auto evec_res = eigenvectors(M, ctx);
    ASSERT_TRUE(evec_res.is_ok()) << "Eigenvectors failed: " << evec_res.error().message;
    ASSERT_GT(evec_res.value().size(), 0U)
        << "Expected at least one eigenpair via algebraic extension path";

    auto residual_is_zero = [&](ExprPtr lambda, const std::vector<ExprPtr>& v) -> bool {
        const std::size_t n = M.rows();
        for (std::size_t i = 0; i < n; ++i) {
            ExprPtr accum = arena.make<IntegerLit>(BigInt(0));
            for (std::size_t j = 0; j < n; ++j) {
                ExprPtr prod = arena.make<Binary>(BinaryOp::Mul, M(i, j), v[j]);
                accum = arena.make<Binary>(BinaryOp::Add, accum, prod);
            }
            ExprPtr lam_vi = arena.make<Binary>(BinaryOp::Mul, lambda, v[i]);
            ExprPtr residual = arena.make<Binary>(BinaryOp::Sub, accum, lam_vi);
            auto reduced = cas::algebra::simplify_in_q_alpha(residual, ctx);
            if (reduced.is_error()) return false;
            const auto* il = expr_cast<IntegerLit>(reduced.value());
            if (!il || !il->value.is_zero()) {
                ADD_FAILURE() << "Component " << i << " residual non-zero: "
                              << debug_print(reduced.value());
                return false;
            }
        }
        return true;
    };

    for (const auto& pair : evec_res.value()) {
        EXPECT_EQ(pair.eigenvector.size(), 4U);
        EXPECT_TRUE(residual_is_zero(pair.eigenvalue, pair.eigenvector));
    }
}
