#include "cas/calculus.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include <gtest/gtest.h>

namespace cas::calculus {
namespace {

[[nodiscard]] Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

TEST(Calculus, RootSumEmission) {
    symbolic::CASContext ctx;
    // 1/(x^5 + x + 1) -> Irreducible degree 5, should emit RootSum
    auto P = parse_expr("1", ctx.arena()).value();
    auto Q = parse_expr("x^5 + x + 1", ctx.arena()).value();
    Symbol x("x");
    
    auto res = algebra::integrate_rational_lrt(P, Q, x, ctx);
    if (!res.is_ok()) std::cout << "LRT ERROR: " << res.error().message << std::endl;
    ASSERT_TRUE(res.is_ok());
    
    // Check if result contains RootSum
    bool has_rootsum = false;
    std::function<void(ExprPtr)> check_rootsum = [&](ExprPtr e) -> void {
        if (!e) return;
        if (expr_is<FuncCall>(e)) {
            auto call = expr_cast<FuncCall>(e);
            if (call->name == "RootSum" || call->func_id == BuiltinOp::RootSum) {
                has_rootsum = true;
                return;
            }
        }
        if (auto sum = expr_cast<Sum>(e)) {
            for (auto t : sum->terms) check_rootsum(t);
        }
        if (auto prod = expr_cast<Product>(e)) {
            for (auto f : prod->factors) check_rootsum(f);
        }
        if (auto bin = expr_cast<Binary>(e)) {
            check_rootsum(bin->left);
            check_rootsum(bin->right);
        }
        if (auto call = expr_cast<FuncCall>(e)) {
            for (auto arg : call->args) check_rootsum(arg);
        }
    };
    check_rootsum(res.value());
    EXPECT_TRUE(has_rootsum);
}

TEST(Calculus, AsymptoteHorizontal) {
    symbolic::CASContext ctx;
    // f(x) = (2x+1)/(x-1) -> y=2
    auto f = parse_expr("(2*x+1)/(x-1)", ctx.arena()).value();
    Symbol x("x");
    auto res = find_asymptotes(f, x, ctx);
    ASSERT_TRUE(res.is_ok());
    
    bool found_h = false;
    for (const auto& asy : res.value()) {
        std::cout << "ASYMPTOTE FOUND: " << (int)asy.type << " value=" << debug_print(asy.expression) << std::endl;
        if (asy.type == Asymptote::Type::Horizontal) {
            if (expr_is<IntegerLit>(asy.expression) && expr_cast<IntegerLit>(asy.expression)->value == BigInt(2LL)) {
                found_h = true;
            }
        }
    }
    EXPECT_TRUE(found_h);
}

TEST(Calculus, AsymptoteSlant) {
    symbolic::CASContext ctx;
    // f(x) = (x^2+1)/x -> y=x
    auto f = parse_expr("(x^2+1)/x", ctx.arena()).value();
    Symbol x("x");
    auto res = find_asymptotes(f, x, ctx);
    ASSERT_TRUE(res.is_ok());
    
    bool found_s = false;
    for (const auto& asy : res.value()) {
        if (asy.type == Asymptote::Type::Slant) {
            found_s = true;
        }
    }
    EXPECT_TRUE(found_s);
}

TEST(Calculus, AsymptoteVertical) {
    symbolic::CASContext ctx;
    // f(x) = 1/(x-3) -> x=3
    auto f = parse_expr("1/(x-3)", ctx.arena()).value();
    Symbol x("x");
    auto res = find_asymptotes(f, x, ctx);
    ASSERT_TRUE(res.is_ok());
    
    bool found_v = false;
    for (const auto& asy : res.value()) {
        if (asy.type == Asymptote::Type::Vertical) {
            if (expr_is<IntegerLit>(asy.expression) && expr_cast<IntegerLit>(asy.expression)->value == BigInt(3LL)) {
                found_v = true;
            }
        }
    }
    EXPECT_TRUE(found_v);
}

// L1-11: x → -∞ analysis — deduplication of symmetric horizontal asymptote
TEST(Calculus, L1_11_AsymptoteHorizontalNegInfinity) {
    symbolic::CASContext ctx;
    // 1/(x^2+1) → y=0 for both x→+∞ and x→-∞: should appear exactly once
    auto f = parse_expr("1/(x^2+1)", ctx.arena()).value();
    Symbol x("x");
    auto res = find_asymptotes(f, x, ctx);
    ASSERT_TRUE(res.is_ok());
    int h_count = 0;
    for (const auto& a : res.value()) {
        if (a.type == Asymptote::Type::Horizontal) h_count++;
    }
    EXPECT_EQ(h_count, 1) << "Symmetric horizontal must not be duplicated";
}

}
}
