// F5.9-pre — Smoke tests parser+AST per Dirac delta + Heaviside builtins.
//
// Scope: verifica che il parser riconosca "delta(...)", "DiracDelta(...)",
// "theta(...)", "Heaviside(...)" come BuiltinOp e che round-trip via
// builtin_op_name sia stabile.  Le regole semantiche (sifting property,
// diff, simplify) sono scope incrementale follow-up.

#include "cas/ast.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas::calculus {
namespace {

[[nodiscard]] ExprPtr parse_expr(const std::string& s, AstArena& arena) {
    auto t = Lexer(s).tokenize();
    EXPECT_TRUE(t.is_ok()) << s << ": " << t.error().message;
    Parser p(t.value(), arena);
    auto r = p.parse();
    EXPECT_TRUE(r.is_ok()) << s;
    return r.value();
}

class DiracHeavisideSmokeTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

// delta(x) parsa come FuncCall(DiracDelta, [x]).
TEST_F(DiracHeavisideSmokeTest, ParseDeltaLowercase) {
    auto e = parse_expr("delta(x)", ctx.arena());
    const auto* fc = expr_cast<FuncCall>(e);
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->func_id, BuiltinOp::DiracDelta);
    ASSERT_EQ(fc->args.size(), 1U);
}

// DiracDelta(x-1) parsa correttamente.
TEST_F(DiracHeavisideSmokeTest, ParseDiracDeltaMixedCase) {
    auto e = parse_expr("DiracDelta(x - 1)", ctx.arena());
    const auto* fc = expr_cast<FuncCall>(e);
    ASSERT_NE(fc, nullptr);
    EXPECT_EQ(fc->func_id, BuiltinOp::DiracDelta);
}

// theta(x), Heaviside(x) entrambi → BuiltinOp::Heaviside.
TEST_F(DiracHeavisideSmokeTest, ParseHeavisideAndTheta) {
    auto e1 = parse_expr("theta(x)", ctx.arena());
    auto e2 = parse_expr("Heaviside(x)", ctx.arena());
    const auto* fc1 = expr_cast<FuncCall>(e1);
    const auto* fc2 = expr_cast<FuncCall>(e2);
    ASSERT_NE(fc1, nullptr);
    ASSERT_NE(fc2, nullptr);
    EXPECT_EQ(fc1->func_id, BuiltinOp::Heaviside);
    EXPECT_EQ(fc2->func_id, BuiltinOp::Heaviside);
}

// builtin_op_name round-trip.
TEST_F(DiracHeavisideSmokeTest, BuiltinNameRoundTrip) {
    EXPECT_EQ(builtin_op_name(BuiltinOp::DiracDelta), "DiracDelta");
    EXPECT_EQ(builtin_op_name(BuiltinOp::Heaviside), "Heaviside");
    EXPECT_EQ(get_builtin_op("DiracDelta"), BuiltinOp::DiracDelta);
    EXPECT_EQ(get_builtin_op("delta"), BuiltinOp::DiracDelta);
    EXPECT_EQ(get_builtin_op("Heaviside"), BuiltinOp::Heaviside);
    EXPECT_EQ(get_builtin_op("theta"), BuiltinOp::Heaviside);
}

}  // namespace
}  // namespace cas::calculus
