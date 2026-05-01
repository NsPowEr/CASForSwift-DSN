#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "cas/formatter.hpp"

#include <iostream>
#include <string>

using namespace cas;
using namespace cas::algebra;
using namespace cas::formatter;

int main() {
    symbolic::CASContext ctx;
    AstArena& arena = ctx.arena();
    TextFormatter formatter;
    
    auto p_tokens = Lexer("2*x^2 - 2").tokenize();
    auto q_tokens = Lexer("4*x - 4").tokenize();
    
    Parser p_parser(p_tokens.value(), arena);
    Parser q_parser(q_tokens.value(), arena);
    
    auto p = p_parser.parse().value();
    auto q = q_parser.parse().value();
    
    Symbol x("x");
    auto gcd_res = polynomial_gcd(p, q, x, ctx);
    
    if (gcd_res.is_error()) {
        std::cerr << "Error: " << gcd_res.error().message << std::endl;
        return 1;
    }
    
    std::cout << "GCD: " << formatter.format(gcd_res.value()) << std::endl;
    
    auto expected_tokens = Lexer("2*x - 2").tokenize();
    auto expected = Parser(expected_tokens.value(), arena).parse().value();
    
    auto diff = arena.make<Binary>(BinaryOp::Sub, gcd_res.value(), expected);
    auto expanded = expand(diff, ctx).value();
    auto simplified = ctx.simplify(expanded).value();
    
    std::cout << "Simplified Diff: " << formatter.format(simplified) << std::endl;
    
    return 0;
}
