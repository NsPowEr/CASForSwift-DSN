#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace cas::algebra {
namespace {

[[nodiscard]] Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

TEST(PartialFractions, IrreducibleQuartic) {
    symbolic::CASContext ctx;
    // 1/(x^4 + 1): x^4 + 1 is irreducible over Q.
    // Partial fraction decomposition should be just 1/(x^4 + 1).
    auto expr = parse_expr("1/(x^4 + 1)", ctx.arena());
    ASSERT_TRUE(expr.is_ok());
    Symbol x("x");
    auto terms = partial_fractions(expr.value(), x, ctx);
    ASSERT_TRUE(terms.is_ok());
    
    // It should have 1 term: (0*x^3 + 0*x^2 + 0*x + 1) / (x^4 + 1)
    ASSERT_EQ(terms.value().size(), 1U);
    
    // Verify by substitution
    ExprPtr x_val = ctx.arena().make<IntegerLit>(BigInt(2LL));
    auto orig_val = ctx.simplify(substitute(expr.value(), x, x_val, ctx).value());
    auto term_val = ctx.simplify(substitute(terms.value()[0], x, x_val, ctx).value());
    
    ASSERT_TRUE(orig_val.is_ok() && term_val.is_ok());
    EXPECT_TRUE(structural_equal(orig_val.value(), term_val.value()));
}

} // namespace
} // namespace cas::algebra

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
