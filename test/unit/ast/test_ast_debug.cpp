#include "cas/ast_debug.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"

#include <gtest/gtest.h>

namespace cas {
namespace {

Result<ExprPtr> parse_debug_input(const std::string& input, AstArena& arena) {
    Lexer lexer(input);
    auto tokens = lexer.tokenize();
    if (tokens.is_error()) {
        return Result<ExprPtr>(tokens.error());
    }

    Parser parser(tokens.value(), arena);
    return parser.parse();
}

TEST(AstDebugTest, PrintsNullExpression) {
    EXPECT_EQ(debug_print(ExprPtr{}), "Null");
}

TEST(AstDebugTest, PrintsDeterministicPrefixTreeForBinaryExpression) {
    AstArena arena;
    const auto parsed = parse_debug_input("sin(x^2 + 1) / (x - 3)", arena);

    ASSERT_TRUE(parsed.is_ok());
    EXPECT_EQ(
        debug_print(parsed.value()),
        "Binary(Div, FuncCall(sin, [Sum([Binary(Pow, Symbol(x), IntegerLit(2)), IntegerLit(1)])]), "
        "Binary(Sub, Symbol(x), IntegerLit(3)))");
}

TEST(AstDebugTest, PrintsCalculusAndRootNodesWithoutEvaluation) {
    AstArena arena;
    const auto x = arena.make<Symbol>(std::string("x"));
    const auto two = arena.make<IntegerLit>(BigInt(2));
    const auto poly = arena.make<Binary>(BinaryOp::Pow, x, two);
    const auto root = arena.make<RootOf>(poly, Symbol{"x"}, std::optional<std::size_t>{1U});
    const auto integral = arena.make<Integral>(poly, Symbol{"x"}, std::nullopt, std::optional<ExprPtr>{two});

    EXPECT_EQ(debug_print(root), "RootOf(Binary(Pow, Symbol(x), IntegerLit(2)), Symbol(x), 1)");
    EXPECT_EQ(debug_print(integral), "Integral(Binary(Pow, Symbol(x), IntegerLit(2)), Symbol(x), None, IntegerLit(2))");
}

}  // namespace
}  // namespace cas
