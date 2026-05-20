#include <gtest/gtest.h>
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "cas/formatter.hpp"

namespace cas {

class QuantityTest : public ::testing::Test {
protected:
    void SetUp() override {
        arena = std::make_unique<AstArena>();
    }

    std::unique_ptr<AstArena> arena;
};

TEST_F(QuantityTest, StructuralEqual) {
    SIDimensions m; m.m = 1;
    SIDimensions s; s.s = 1;
    SIDimensions ms2; ms2.m = 1; ms2.s = -2;

    auto val5 = arena->make<IntegerLit>(BigInt(5));
    auto val10 = arena->make<IntegerLit>(BigInt(10));

    auto q1 = arena->make<Quantity>(val5, m);
    auto q2 = arena->make<Quantity>(val5, m);
    auto q3 = arena->make<Quantity>(val10, m);
    auto q4 = arena->make<Quantity>(val5, s);
    auto q5 = arena->make<Quantity>(val5, ms2);

    EXPECT_TRUE(structural_equal(q1, q2));
    EXPECT_FALSE(structural_equal(q1, q3));
    EXPECT_FALSE(structural_equal(q1, q4));
    EXPECT_FALSE(structural_equal(q1, q5));
}

TEST_F(QuantityTest, Hashing) {
    SIDimensions m; m.m = 1;
    auto val5 = arena->make<IntegerLit>(BigInt(5));

    auto q1 = arena->make<Quantity>(val5, m);
    auto q2 = arena->make<Quantity>(val5, m);

    EXPECT_EQ(expr_hash(q1), expr_hash(q2));
}

TEST_F(QuantityTest, Formatting) {
    formatter::TextFormatter fmt;
    
    SIDimensions m; m.m = 1;
    SIDimensions ms2; ms2.m = 1; ms2.s = -2;
    SIDimensions kgm2s2; kgm2s2.kg = 1; kgm2s2.m = 2; kgm2s2.s = -2;

    auto val5 = arena->make<IntegerLit>(BigInt(5));
    auto val981 = arena->make<DecimalLit>("9.81");

    EXPECT_EQ(fmt.format(arena->make<Quantity>(val5, m)), "5[m]");
    EXPECT_EQ(fmt.format(arena->make<Quantity>(val981, ms2)), "9.81[m*s^-2]");
    EXPECT_EQ(fmt.format(arena->make<Quantity>(val5, kgm2s2)), "5[m^2*kg*s^-2]");
}

TEST_F(QuantityTest, Cloning) {
    AstArena target_arena;
    SIDimensions ms2; ms2.m = 1; ms2.s = -2;
    auto val981 = arena->make<DecimalLit>("9.81");
    auto q_source = arena->make<Quantity>(val981, ms2);

    std::unordered_map<ExprPtr, ExprPtr> cache;
    auto q_cloned = clone_into_arena(q_source, target_arena, cache);

    EXPECT_TRUE(structural_equal(q_source, q_cloned));
    EXPECT_NE(q_source, q_cloned); // Should be different pointers (different arenas)
    EXPECT_EQ(expr_kind(q_cloned), ExprKind::Quantity);
}

TEST_F(QuantityTest, Multiplication) {
    SIDimensions m; m.m = 1;
    SIDimensions s_inv; s_inv.s = -1;
    
    auto val5 = arena->make<IntegerLit>(BigInt(5));
    auto val2 = arena->make<IntegerLit>(BigInt(2));
    
    auto q_m = arena->make<Quantity>(val5, m);
    auto q_s_inv = arena->make<Quantity>(val2, s_inv);
    
    // Test: 5[m] * 2[s^-1] -> 10[m*s^-1]
    auto prod = arena->make<Product>(std::vector<ExprPtr>{q_m, q_s_inv});
    auto simplified = symbolic::simplify(prod, *arena);
    
    ASSERT_TRUE(simplified.is_ok());
    EXPECT_EQ(expr_kind(simplified.value()), ExprKind::Quantity);
    
    const auto& res = expr_ref<Quantity>(simplified.value());
    EXPECT_EQ(res.dimensions.m, 1);
    EXPECT_EQ(res.dimensions.s, -1);
    
    formatter::TextFormatter fmt;
    EXPECT_EQ(fmt.format(simplified.value()), "10[m*s^-1]");
}

TEST_F(QuantityTest, Division) {
    SIDimensions m; m.m = 1;
    SIDimensions s; s.s = 1;
    
    auto val10 = arena->make<IntegerLit>(BigInt(10));
    auto val2 = arena->make<IntegerLit>(BigInt(2));
    
    auto q_m = arena->make<Quantity>(val10, m);
    auto q_s = arena->make<Quantity>(val2, s);
    
    // Test: 10[m] / 2[s] -> 5[m*s^-1]
    auto div = arena->make<Binary>(BinaryOp::Div, q_m, q_s);
    auto simplified = symbolic::simplify(div, *arena);
    
    ASSERT_TRUE(simplified.is_ok());
    formatter::TextFormatter fmt;
    EXPECT_EQ(fmt.format(simplified.value()), "5[m*s^-1]");
}

TEST_F(QuantityTest, MultiplicationSameUnit) {
    SIDimensions m; m.m = 1;
    auto val5 = arena->make<IntegerLit>(BigInt(5));
    auto val2 = arena->make<IntegerLit>(BigInt(2));
    
    auto q1 = arena->make<Quantity>(val5, m);
    auto q2 = arena->make<Quantity>(val2, m);
    
    // Test: 5[m] * 2[m] -> 10[m^2]
    auto prod = arena->make<Product>(std::vector<ExprPtr>{q1, q2});
    auto simplified = symbolic::simplify(prod, *arena);
    
    ASSERT_TRUE(simplified.is_ok());
    formatter::TextFormatter fmt;
    EXPECT_EQ(fmt.format(simplified.value()), "10[m^2]");
}

TEST_F(QuantityTest, AdditionSameUnit) {
    SIDimensions m; m.m = 1;
    auto val5 = arena->make<IntegerLit>(BigInt(5));
    auto val3 = arena->make<IntegerLit>(BigInt(3));
    
    auto q1 = arena->make<Quantity>(val5, m);
    auto q2 = arena->make<Quantity>(val3, m);
    
    // Test: 5[m] + 3[m] -> 8[m]
    auto sum = arena->make<Sum>(std::vector<ExprPtr>{q1, q2});
    auto simplified = symbolic::simplify(sum, *arena);
    
    ASSERT_TRUE(simplified.is_ok());
    formatter::TextFormatter fmt;
    EXPECT_EQ(fmt.format(simplified.value()), "8[m]");
}

TEST_F(QuantityTest, SubtractionSameUnit) {
    SIDimensions m; m.m = 1;
    auto valx = arena->make<Symbol>("x");
    
    auto q1 = arena->make<Quantity>(valx, m);
    auto neg_q1 = arena->make<Unary>(UnaryOp::Neg, q1);
    
    // Test: x[m] - x[m] -> 0
    auto sum = arena->make<Sum>(std::vector<ExprPtr>{q1, neg_q1});
    auto simplified = symbolic::simplify(sum, *arena);
    
    ASSERT_TRUE(simplified.is_ok());
    // Check the simplified result is structurally zero (IntegerLit 0) or simplifies to 0
    if (const auto* lit = expr_cast<IntegerLit>(simplified.value())) {
        EXPECT_TRUE(lit->value.is_zero());
    } else {
        symbolic::CASContext ctx;
        auto zero = ctx.arena().make<IntegerLit>(BigInt(0));
        auto eq = symbolic::mathematically_equal(simplified.value(), zero, ctx);
        EXPECT_TRUE(eq.is_ok() && eq.value());
    }
}

TEST_F(QuantityTest, AdditionDifferentUnits) {
    SIDimensions m; m.m = 1;
    SIDimensions s; s.s = 1;
    auto val5 = arena->make<IntegerLit>(BigInt(5));
    auto val2 = arena->make<IntegerLit>(BigInt(2));
    
    auto q1 = arena->make<Quantity>(val5, m);
    auto q2 = arena->make<Quantity>(val2, s);
    
    // Test: 5[m] + 2[s] -> 5[m] + 2[s] (stay separate)
    auto sum = arena->make<Sum>(std::vector<ExprPtr>{q1, q2});
    auto simplified = symbolic::simplify(sum, *arena);
    
    ASSERT_TRUE(simplified.is_ok());
    formatter::TextFormatter fmt;
    // Note: TextFormatter sorts Sum terms. Quantity precedence is 100.
    // canonical_compare might put m before s or vice versa.
    std::string result = fmt.format(simplified.value());
    EXPECT_TRUE(result == "5[m] + 2[s]" || result == "2[s] + 5[m]");
}

} // namespace cas
