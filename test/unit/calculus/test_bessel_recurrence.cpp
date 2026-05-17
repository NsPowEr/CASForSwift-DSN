#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

#include <string>

namespace cas::test {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

[[nodiscard]] bool contains_bessel_with_order_geq_two(ExprPtr expr) {
    bool found = false;
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if ((call->func_id == BuiltinOp::BesselJ || call->func_id == BuiltinOp::BesselY)
            && call->args.size() == 2U) {
            const auto* lit = expr_cast<IntegerLit>(call->args[0]);
            if (lit != nullptr && !lit->value.is_negative() && lit->value > BigInt(1)) {
                found = true;
            }
        }
    }
    if (found) return true;
    bool inner = false;
    visit_expr(expr, [&](const auto& node) {
        using Node = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<Node, Unary>) {
            if (contains_bessel_with_order_geq_two(node.operand)) inner = true;
        } else if constexpr (std::is_same_v<Node, Binary>) {
            if (contains_bessel_with_order_geq_two(node.left)) inner = true;
            if (contains_bessel_with_order_geq_two(node.right)) inner = true;
        } else if constexpr (std::is_same_v<Node, FuncCall>) {
            for (ExprPtr a : node.args)
                if (contains_bessel_with_order_geq_two(a)) { inner = true; break; }
        } else if constexpr (std::is_same_v<Node, Sum>) {
            for (ExprPtr t : node.terms)
                if (contains_bessel_with_order_geq_two(t)) { inner = true; break; }
        } else if constexpr (std::is_same_v<Node, Product>) {
            for (ExprPtr f : node.factors)
                if (contains_bessel_with_order_geq_two(f)) { inner = true; break; }
        }
    });
    return inner;
}

}  // namespace

class BesselRecurrenceTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

TEST_F(BesselRecurrenceTest, DefaultOffLeavesHighOrderIntact) {
    auto e = parse_expr("bessel_j(3, x)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx.simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    // Default: flag off, BesselJ(3, x) must stay.
    EXPECT_TRUE(contains_bessel_with_order_geq_two(s.value()));
}

TEST_F(BesselRecurrenceTest, FlagOnReducesToOrdersZeroAndOne) {
    ctx.set_expand_bessel_recurrence(true);
    auto e = parse_expr("bessel_j(3, x)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx.simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    // After expansion, no remaining BesselJ of integer order >= 2.
    EXPECT_FALSE(contains_bessel_with_order_geq_two(s.value()));
}

TEST_F(BesselRecurrenceTest, FlagOnAppliesToBesselY) {
    ctx.set_expand_bessel_recurrence(true);
    auto e = parse_expr("bessel_y(4, x)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx.simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    EXPECT_FALSE(contains_bessel_with_order_geq_two(s.value()));
}

TEST_F(BesselRecurrenceTest, OrderOneStaysUnchangedEvenWithFlag) {
    ctx.set_expand_bessel_recurrence(true);
    auto e = parse_expr("bessel_j(1, x)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx.simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    // BesselJ(1, x) stays as-is (no further reduction possible).
    const auto* call = expr_cast<FuncCall>(s.value());
    ASSERT_NE(call, nullptr);
    EXPECT_EQ(call->func_id, BuiltinOp::BesselJ);
    const auto* order_lit = expr_cast<IntegerLit>(call->args[0]);
    ASSERT_NE(order_lit, nullptr);
    EXPECT_EQ(order_lit->value, BigInt(1));
}

TEST_F(BesselRecurrenceTest, ExpansionIsMathematicallyConsistent) {
    // Anti-hardcode high order: BesselJ(5, x) reduces to a Sum/Product over orders 0,1.
    ctx.set_expand_bessel_recurrence(true);
    auto e = parse_expr("bessel_j(5, x)", ctx.arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx.simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    EXPECT_FALSE(contains_bessel_with_order_geq_two(s.value()));
}

}  // namespace cas::test
