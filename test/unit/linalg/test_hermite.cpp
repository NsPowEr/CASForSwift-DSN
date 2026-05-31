// CAS-F4.2c — Hermite Normal Form test.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::linalg;

namespace {

class HermiteTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    [[nodiscard]] ExprPtr lit(long long v) {
        return ctx.arena().make<IntegerLit>(BigInt(v));
    }
    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        if (t.is_error()) return nullptr;
        Parser p(t.value(), ctx.arena());
        auto e = p.parse();
        if (e.is_error()) return nullptr;
        auto simp = ctx.simplify(e.value());
        return simp.is_ok() ? simp.value() : e.value();
    }
    [[nodiscard]] bool entries_equal(ExprPtr a, ExprPtr b) {
        auto eq = symbolic::mathematically_equal(a, b, ctx);
        return eq.is_ok() && eq.value();
    }
};

TEST_F(HermiteTest, Z_2x2_simple) {
    // A = [[2, 4], [6, 10]] → HNF: upper triangular, 0 ≤ H[0][1] < H[0][0].
    MatrixExpr A(2, 2, {lit(2), lit(4), lit(6), lit(10)});
    auto h = hermite_normal_form(A, ctx);
    ASSERT_TRUE(h.is_ok()) << h.error().message;
    // Cert: U·A == H
    auto UA = multiply(h.value().U, A, ctx);
    ASSERT_TRUE(UA.is_ok());
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            EXPECT_TRUE(entries_equal(UA.value()(i, j), h.value().H(i, j)))
                << "U·A != H at (" << i << "," << j << ")";
    // H upper-triangular: H[1][0] = 0
    EXPECT_TRUE(entries_equal(h.value().H(1, 0), lit(0)));
}

TEST_F(HermiteTest, Z_3x3_random) {
    MatrixExpr A(3, 3, {
        lit(2), lit(3), lit(5),
        lit(1), lit(2), lit(4),
        lit(3), lit(7), lit(11),
    });
    auto h = hermite_normal_form(A, ctx);
    ASSERT_TRUE(h.is_ok());
    auto UA = multiply(h.value().U, A, ctx);
    ASSERT_TRUE(UA.is_ok());
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            EXPECT_TRUE(entries_equal(UA.value()(i, j), h.value().H(i, j)));
    // Upper-triangular
    EXPECT_TRUE(entries_equal(h.value().H(1, 0), lit(0)));
    EXPECT_TRUE(entries_equal(h.value().H(2, 0), lit(0)));
    EXPECT_TRUE(entries_equal(h.value().H(2, 1), lit(0)));
}

TEST_F(HermiteTest, Qx_diagonal) {
    auto x = parse("x");
    auto x2 = parse("x^2");
    MatrixExpr A(2, 2, {x, lit(0), lit(0), x2});
    auto h = hermite_normal_form(A, ctx);
    ASSERT_TRUE(h.is_ok()) << h.error().message;
    // U·A == H
    auto UA = multiply(h.value().U, A, ctx);
    ASSERT_TRUE(UA.is_ok());
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            EXPECT_TRUE(entries_equal(UA.value()(i, j), h.value().H(i, j)));
}

}  // namespace
