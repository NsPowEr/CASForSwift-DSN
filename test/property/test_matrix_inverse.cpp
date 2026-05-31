// F0.4 — Property: A * inv(A) ≡ I  for random rational 3×3 matrices.
//
// Seed corpus: 5 rational 3×3 matrices with known inverses.
// Verification: multiply A by inv(A) and check each entry of the
// product against the identity matrix via simplify.
//
// TODO: scale to 5×5 (plan: use rapidcheck generator for rational entries
//       once matrix_inverse robustness on symbolic rationals is confirmed).

#include "cas/linalg/Matrix.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <rapidcheck/gtest.h>
#include <string>
#include <vector>

namespace cas::property {
namespace {

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

[[nodiscard]] ExprPtr must_parse(const std::string& s, symbolic::CASContext& ctx) {
    auto t = Lexer(s).tokenize();
    if (!t.is_ok()) throw std::runtime_error("lex: " + t.error().message);
    auto e = Parser(t.value(), ctx.arena()).parse();
    if (!e.is_ok()) throw std::runtime_error("parse: " + e.error().message);
    return e.value();
}

[[nodiscard]] bool is_zero(ExprPtr e, symbolic::CASContext& ctx) {
    auto s = ctx.simplify(e);
    if (!s.is_ok()) return false;
    ExprPtr sv = s.value();
    if (const auto* il = expr_cast<IntegerLit>(sv)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(sv)) return rl->numerator.is_zero();
    return false;
}

[[nodiscard]] bool is_one(ExprPtr e, symbolic::CASContext& ctx) {
    auto s = ctx.simplify(e);
    if (!s.is_ok()) return false;
    ExprPtr sv = s.value();
    if (const auto* il = expr_cast<IntegerLit>(sv)) return il->value == BigInt(1);
    if (const auto* rl = expr_cast<RationalLit>(sv))
        return rl->numerator == rl->denominator;
    return false;
}

// Build a 3×3 MatrixExpr from row-major string entries.
[[nodiscard]] linalg::MatrixExpr make_matrix3(
    const std::vector<std::string>& entries,
    symbolic::CASContext& ctx)
{
    linalg::MatrixExpr m(3, 3);
    for (std::size_t i = 0; i < 9; ++i) {
        m(i / 3, i % 3) = must_parse(entries[i], ctx);
    }
    return m;
}

// Check A * inv(A) == I (all diagonal == 1, off-diagonal == 0).
[[nodiscard]] bool verify_inverse(
    const linalg::MatrixExpr& A, symbolic::CASContext& ctx)
{
    auto inv_res = linalg::inverse(A, ctx);
    if (!inv_res.is_ok()) return false;  // singular or unimplemented
    auto prod_res = linalg::multiply(A, inv_res.value(), ctx);
    if (!prod_res.is_ok()) return false;
    const auto& P = prod_res.value();
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            if (r == c) {
                if (!is_one(P(r, c), ctx)) return false;
            } else {
                if (!is_zero(P(r, c), ctx)) return false;
            }
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Seed matrices (rational 3×3, all invertible)
// ---------------------------------------------------------------------------

using MatrixSpec = std::vector<std::string>;

static const std::vector<MatrixSpec> kSeedMatrices = {
    // Identity
    {"1","0","0", "0","1","0", "0","0","1"},
    // Simple invertible
    {"2","1","0", "1","2","1", "0","1","2"},
    // Rational entries
    {"1","1/2","1/3", "0","1","1/2", "0","0","1"},
    // Permutation-like
    {"0","1","0", "0","0","1", "1","0","0"},
    // Symmetric positive-definite small integers (det=8)
    {"3","1","2", "1","2","1", "2","1","3"},
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

class MatrixInverseTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

TEST_F(MatrixInverseTest, IdentityMatrix) {
    auto A = make_matrix3(kSeedMatrices[0], ctx);
    EXPECT_TRUE(verify_inverse(A, ctx));
}

TEST_F(MatrixInverseTest, TridiagonalLike) {
    auto A = make_matrix3(kSeedMatrices[1], ctx);
    EXPECT_TRUE(verify_inverse(A, ctx));
}

TEST_F(MatrixInverseTest, UpperTriangularRational) {
    auto A = make_matrix3(kSeedMatrices[2], ctx);
    EXPECT_TRUE(verify_inverse(A, ctx));
}

TEST_F(MatrixInverseTest, PermutationMatrix) {
    auto A = make_matrix3(kSeedMatrices[3], ctx);
    EXPECT_TRUE(verify_inverse(A, ctx));
}

TEST_F(MatrixInverseTest, SmallIntegers) {
    auto A = make_matrix3(kSeedMatrices[4], ctx);
    EXPECT_TRUE(verify_inverse(A, ctx));
}

// rapidcheck: all seeds pass.
RC_GTEST_FIXTURE_PROP(MatrixInverseTest, AllSeedsPass, ()) {
    for (const auto& spec : kSeedMatrices) {
        symbolic::CASContext lctx;
        auto A = make_matrix3(spec, lctx);
        RC_ASSERT(verify_inverse(A, lctx));
    }
}

}  // namespace
}  // namespace cas::property
