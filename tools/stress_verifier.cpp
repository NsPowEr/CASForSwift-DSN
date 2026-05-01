#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include "cas/symbolic.hpp"
#include "cas/parser.hpp"
#include "cas/lexer.hpp"
#include "cas/formatter.hpp"
#include "cas/calculus.hpp"
#include "cas/algebra.hpp"
#include "cas/linalg/Matrix.hpp"

using namespace cas;

struct TestResult {
    std::string output;
    size_t nodes;
    long long duration_ms;
    bool timeout;
    bool error;
    std::string error_msg;
};

TestResult run_test(const std::string& input, symbolic::CASContext& ctx) {
    auto start = std::chrono::steady_clock::now();
    
    Lexer lex(input);
    auto tokens_res = lex.tokenize();
    if (tokens_res.is_error()) {
        return {"", 0, 0, false, true, tokens_res.error().message};
    }
    
    Parser parser(tokens_res.value(), ctx.arena());
    auto expr_res = parser.parse();
    if (expr_res.is_error()) {
        return {"", 0, 0, false, true, expr_res.error().message};
    }
    
    ExprPtr expr = expr_res.value();
    ExprPtr result = nullptr;
    
    try {
        if (const auto* integral = expr_cast<Integral>(expr)) {
            auto res = calculus::integrate(integral->integrand, integral->variable, ctx);
            if (res.is_ok()) result = res.value();
            else return {"", 0, 0, false, true, res.error().message};
        } else if (const auto* limit = expr_cast<Limit>(expr)) {
            auto res = calculus::limit(limit->expression, limit->variable, limit->point, limit->direction, ctx);
            if (res.is_ok()) result = res.value();
            else return {"", 0, 0, false, true, res.error().message};
        } else if (const auto* call = expr_cast<FuncCall>(expr)) {
            if (call->name == "factor") {
                if (call->args.size() == 2) {
                    // Try to see if it's the requested stress test 2
                    // factor(x^4 + 1, sqrt(2))
                    // Actually, we don't have a direct API for extension fields yet in algebra.hpp
                    // but let's see what happens if we try to factor it.
                    // For now, let's just use factor_over_integers on the first arg.
                    // Or if it's the determinant one: factor(det(...))
                    if (const auto* inner_call = expr_cast<FuncCall>(call->args[0])) {
                        if (inner_call->name == "det") {
                             auto simplified = ctx.simplify(call->args[0]);
                             if (simplified.is_ok()) {
                                 auto res = algebra::factor_over_integers(simplified.value(), Symbol("a"), ctx); // Arbitrary var for now
                                 // Handle multivar if possible, but factor_over_integers takes one var.
                                 // Let's just simplify the whole thing if it's a known function.
                                 if (res.is_ok()) {
                                     // Format factorization
                                     std::string out = "";
                                     for (const auto& f : res.value().factors) {
                                         formatter::TextFormatter fmt;
                                         out += "(" + fmt.format(f.factor) + ")^" + std::to_string(f.multiplicity) + " ";
                                     }
                                     auto end = std::chrono::steady_clock::now();
                                     return {out, ctx.arena().size(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(), false, false, ""};
                                 }
                             }
                        }
                    }
                    
                    // Default to simplify
                    auto res = ctx.simplify(expr);
                    if (res.is_ok()) result = res.value();
                    else return {"", 0, 0, false, true, res.error().message};
                } else {
                    auto res = ctx.simplify(expr);
                    if (res.is_ok()) result = res.value();
                    else return {"", 0, 0, false, true, res.error().message};
                }
            } else {
                auto res = ctx.simplify(expr);
                if (res.is_ok()) result = res.value();
                else return {"", 0, 0, false, true, res.error().message};
            }
        } else {
            auto res = ctx.simplify(expr);
            if (res.is_ok()) result = res.value();
            else return {"", 0, 0, false, true, res.error().message};
        }
    } catch (...) {
        return {"", 0, 0, false, true, "Exception during execution"};
    }
    
    auto end = std::chrono::steady_clock::now();
    formatter::TextFormatter fmt;
    return {fmt.format(result), ctx.arena().size(), std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(), false, false, ""};
}

int main() {
    symbolic::CASContext ctx;
    ctx.set_timeout(std::chrono::seconds(10));
    
    std::vector<std::string> tests = {
        "int(1 / (x^8 + 1), x)",
        "factor(x^4 + 1, sqrt(2))",
        "lim((cos(x) - exp(-x^2/2)) / x^4, x, 0)",
        "simplify(sqrt(-1) * sqrt(-1) - sqrt((-1) * (-1)))",
        "factor(det([[1, a, a^2, a^3], [1, b, b^2, b^3], [1, c, c^2, c^3], [1, d, d^2, d^3]]))"
    };
    
    for (size_t i = 0; i < tests.size(); ++i) {
        std::cout << "TEST " << (i + 1) << ": " << tests[i] << std::endl;
        TestResult res = run_test(tests[i], ctx);
        if (res.error) {
            std::cout << "  ERROR: " << res.error_msg << std::endl;
        } else {
            std::cout << "  OUTPUT: " << res.output << std::endl;
            std::cout << "  NODES: " << res.nodes << std::endl;
            std::cout << "  TIME: " << res.duration_ms << " ms" << std::endl;
        }
        std::cout << "----------------------------------------" << std::endl;
    }
    
    return 0;
}
