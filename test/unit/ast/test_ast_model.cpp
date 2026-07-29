#include "cas/ast.hpp"

#include <gtest/gtest.h>

namespace cas {
namespace {

TEST(AstArenaTest, ArenaAllocatesStableExprNodes) {
    AstArena arena;
    const auto one = arena.make<IntegerLit>(BigInt(1));
    const auto two = arena.make<IntegerLit>(BigInt(2));

    ASSERT_TRUE(one);
    ASSERT_TRUE(two);
    EXPECT_NE(one.get(), two.get());
    EXPECT_EQ(arena.size(), 2U);
}

TEST(AstArenaTest, HashConsingInternsIdenticalNodes) {
    AstArena arena;
    const auto x1 = arena.make<Symbol>(std::string("x"));
    const auto x2 = arena.make<Symbol>(std::string("x"));
    
    EXPECT_EQ(x1.get(), x2.get());
    EXPECT_TRUE(is_equal(x1, x2));
    
    const auto one1 = arena.make<IntegerLit>(BigInt(1));
    const auto one2 = arena.make<IntegerLit>(BigInt(1));
    EXPECT_EQ(one1.get(), one2.get());
    
    const auto sum1 = arena.make<Binary>(BinaryOp::Add, x1, one1);
    const auto sum2 = arena.make<Binary>(BinaryOp::Add, x2, one2);
    EXPECT_EQ(sum1.get(), sum2.get());
    EXPECT_TRUE(is_equal(sum1, sum2));
}

TEST(AstArenaTest, DynamicShardingConfigurable) {
    AstArena arena(32);
    EXPECT_EQ(arena.num_shards(), 32U);

    const auto x1 = arena.make<Symbol>(std::string("x"));
    const auto x2 = arena.make<Symbol>(std::string("x"));
    EXPECT_EQ(x1.get(), x2.get());

    AstArena default_arena;
    EXPECT_GT(default_arena.num_shards(), 0U);
    // Should be power of 2
    EXPECT_EQ(default_arena.num_shards() & (default_arena.num_shards() - 1U), 0U);
}

TEST(AstArenaTest, ReconfigureShardsEmpty) {
    AstArena arena(8);
    EXPECT_EQ(arena.num_shards(), 8U);

    arena.configure_shards(64);
    EXPECT_EQ(arena.num_shards(), 64U);

    // After allocation, configure_shards should do nothing
    const auto x = arena.make<Symbol>(std::string("x"));
    EXPECT_TRUE(x);
    arena.configure_shards(128);
    EXPECT_EQ(arena.num_shards(), 64U);
}

TEST(AstModelTest, DistinguishesExactAndDecimalLiterals) {
    AstArena arena;
    const auto exact = arena.make<IntegerLit>(BigInt(314));
    const auto decimal = arena.make<DecimalLit>(std::string("3.14"));

    EXPECT_EQ(expr_kind(exact), ExprKind::IntegerLit);
    EXPECT_EQ(expr_kind(decimal), ExprKind::DecimalLit);
    EXPECT_FALSE(structural_equal(exact, decimal));
}

TEST(AstModelTest, NullExprPtrHasDedicatedKind) {
    const ExprPtr null_expr;

    EXPECT_EQ(expr_kind(null_expr), ExprKind::Null);
    EXPECT_EQ(expr_kind_name(expr_kind(null_expr)), "Null");
}

TEST(AstModelTest, PreservesDedicatedMathConstants) {
    AstArena arena;
    const auto pi = arena.make<Constant>(MathConstant::Pi);
    const auto symbol = arena.make<Symbol>(std::string("pi"));

    EXPECT_EQ(expr_kind(pi), ExprKind::Constant);
    EXPECT_EQ(expr_kind(symbol), ExprKind::Symbol);
    EXPECT_FALSE(structural_equal(pi, symbol));
}

TEST(AstModelTest, SupportsDedicatedRootOfNode) {
    AstArena arena;
    const auto x = arena.make<Symbol>(std::string("x"));
    const auto exponent = arena.make<IntegerLit>(BigInt(5));
    const auto polynomial = arena.make<Binary>(BinaryOp::Pow, x, exponent);
    const auto root = arena.make<RootOf>(polynomial, Symbol{"x"}, std::optional<std::size_t>{2U});

    EXPECT_EQ(expr_kind(root), ExprKind::RootOf);
    EXPECT_EQ(expr_kind_name(expr_kind(root)), "RootOf");
}

TEST(AstModelTest, StructuralEqualityRecursesAcrossChildNodes) {
    AstArena arena;
    const auto x1 = arena.make<Symbol>(std::string("x"));
    const auto x2 = arena.make<Symbol>(std::string("x"));
    const auto one1 = arena.make<IntegerLit>(BigInt(1));
    const auto one2 = arena.make<IntegerLit>(BigInt(1));
    const auto left = arena.make<Binary>(BinaryOp::Add, x1, one1);
    const auto right = arena.make<Binary>(BinaryOp::Add, x2, one2);

    EXPECT_TRUE(structural_equal(left, right));
}

TEST(AstModelTest, StoresInertCalculusNodesWithoutEvaluation) {
    AstArena arena;
    const auto x = arena.make<Symbol>(std::string("x"));
    const auto expr = arena.make<FuncCall>(std::string("sin"), std::vector<ExprPtr>{x});
    const auto derivative = arena.make<Derivative>(expr, Symbol{"x"}, 2U);
    const auto integral = arena.make<Integral>(expr, Symbol{"x"}, std::nullopt, std::nullopt);
    const auto infinity = arena.make<Constant>(MathConstant::Infinity);
    const auto limit = arena.make<Limit>(expr, Symbol{"x"}, infinity, LimitDirection::Left);

    EXPECT_EQ(expr_kind(derivative), ExprKind::Derivative);
    EXPECT_EQ(expr_kind(integral), ExprKind::Integral);
    EXPECT_EQ(expr_kind(limit), ExprKind::Limit);
}

TEST(AstModelTest, SupportsNArySumAndProductNodes) {
    AstArena arena;
    const auto x = arena.make<Symbol>(std::string("x"));
    const auto y = arena.make<Symbol>(std::string("y"));
    const auto z = arena.make<Symbol>(std::string("z"));

    const auto sum = arena.make<Sum>(std::vector<ExprPtr>{x, y, z});
    const auto product = arena.make<Product>(std::vector<ExprPtr>{x, y, z});

    EXPECT_EQ(expr_kind(sum), ExprKind::Sum);
    EXPECT_EQ(expr_kind(product), ExprKind::Product);
}

TEST(AstModelTest, MatrixCarriesShapeAndElements) {
    AstArena arena;
    const auto one = arena.make<IntegerLit>(BigInt(1));
    const auto zero = arena.make<IntegerLit>(BigInt(0));
    const auto matrix = arena.make<Matrix>(2U, 2U, std::vector<ExprPtr>{one, zero, zero, one});

    ASSERT_EQ(expr_kind(matrix), ExprKind::Matrix);
    const auto* matrix_value = expr_cast<Matrix>(matrix);
    ASSERT_NE(matrix_value, nullptr);
    EXPECT_EQ(matrix_value->rows, 2U);
    EXPECT_EQ(matrix_value->cols, 2U);
    EXPECT_EQ(matrix_value->elements.size(), 4U);
}

}  // namespace
}  // namespace cas
