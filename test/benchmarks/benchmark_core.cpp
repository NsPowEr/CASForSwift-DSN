#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "cas/algebra.hpp"
#include "cas/linalg/Matrix.hpp"

#include <chrono>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

cas::Result<cas::ExprPtr> parse_expr(std::string_view input, cas::AstArena& arena) {
    auto tokens = cas::Lexer(input).tokenize();
    if (tokens.is_error()) {
        return cas::fail<cas::ExprPtr>(tokens.error());
    }

    cas::Parser parser(tokens.value(), arena);
    return parser.parse();
}

double benchmark_parse_implicit_multiplication() {
    constexpr int iterations = 400;
    const std::string input = "2x(y+1)(x+y+z)(a+b+c)";

    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        cas::AstArena arena;
        auto parsed = parse_expr(input, arena);
        if (parsed.is_error()) {
            throw std::runtime_error(parsed.error().message);
        }
    }
    const auto end = Clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count() / static_cast<double>(iterations);
}

double benchmark_round_trip() {
    constexpr int iterations = 200;
    const std::string input = "(x^2 + 2x + 1)/(y - 3) + sin(theta)^2 + cos(theta)^2";

    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        cas::AstArena arena;
        auto parsed = parse_expr(input, arena);
        if (parsed.is_error()) {
            throw std::runtime_error(parsed.error().message);
        }

        auto text = cas::to_round_trip_text(parsed.value());
        if (text.is_error()) {
            throw std::runtime_error(text.error().message);
        }
    }
    const auto end = Clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count() / static_cast<double>(iterations);
}

double benchmark_symbolic_simplify() {
    constexpr int iterations = 250;
    const std::string input = "x + x + 3/4 + 5/6 + sin(t)^2 + cos(t)^2";

    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        cas::AstArena parse_arena;
        auto parsed = parse_expr(input, parse_arena);
        if (parsed.is_error()) {
            throw std::runtime_error(parsed.error().message);
        }

        cas::symbolic::CASContext context;
        auto simplified = context.simplify(parsed.value());
        if (simplified.is_error()) {
            throw std::runtime_error(simplified.error().message);
        }
    }
    const auto end = Clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count() / static_cast<double>(iterations);
}

double benchmark_arena_allocation() {
    constexpr int iterations = 80;
    constexpr int node_count = 2048;

    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        cas::AstArena arena;
        std::vector<cas::ExprPtr> values;
        values.reserve(static_cast<std::size_t>(node_count));

        for (int j = 0; j < node_count; ++j) {
            values.push_back(arena.make<cas::IntegerLit>(cas::BigInt(j + 1)));
        }

        cas::ExprPtr chain = values.front();
        for (int j = 1; j < node_count; ++j) {
            chain = arena.make<cas::Binary>(cas::BinaryOp::Add, chain, values[static_cast<std::size_t>(j)]);
        }

        if (!chain) {
            throw std::runtime_error("arena allocation benchmark produced null expression");
        }
    }
    const auto end = Clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count() / static_cast<double>(iterations);
}

double benchmark_polynomial_expansion_complex() {
    constexpr int iterations = 10;
    const std::string input = "(x + y + z + 1)^4";

    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        cas::AstArena parse_arena;
        auto parsed = parse_expr(input, parse_arena);
        if (parsed.is_error()) throw std::runtime_error(parsed.error().message);
        
        cas::symbolic::CASContext context;
        auto expanded = cas::algebra::expand(parsed.value(), context);
        if (expanded.is_error()) throw std::runtime_error(expanded.error().message);
    }
    const auto end = Clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count() / static_cast<double>(iterations);
}

double benchmark_polynomial_gcd_subresultant() {
    constexpr int iterations = 5;
    const std::string p_str = "(x^3 + 2x^2 + 3x + 4)^3 * (x^2 - x + 1)";
    const std::string q_str = "(x^3 + 2x^2 + 3x + 4)^2 * (x^4 + x^2 + 1)";

    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        cas::AstArena parse_arena;
        auto p = parse_expr(p_str, parse_arena);
        auto q = parse_expr(q_str, parse_arena);
        if (p.is_error()) throw std::runtime_error(p.error().message);
        if (q.is_error()) throw std::runtime_error(q.error().message);
        
        cas::symbolic::CASContext context;
        auto res = cas::algebra::polynomial_gcd(p.value(), q.value(), cas::Symbol{"x"}, context);
        if (res.is_error()) throw std::runtime_error(res.error().message);
    }
    const auto end = Clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count() / static_cast<double>(iterations);
}

double benchmark_matrix_inversion_10x10() {
    constexpr int iterations = 2;
    
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        cas::symbolic::CASContext context;
        auto& arena = context.arena();
        
        cas::linalg::MatrixExpr mat(10, 10);
        for(int r = 0; r < 10; ++r) {
            for(int c = 0; c < 10; ++c) {
                if (r == c) {
                    mat(r, c) = arena.make<cas::Symbol>("x_" + std::to_string(r));
                } else {
                    mat(r, c) = arena.make<cas::IntegerLit>(cas::BigInt(r + c + 1));
                }
            }
        }
        
        auto inv = cas::linalg::inverse(mat, context);
        if (inv.is_error()) throw std::runtime_error(inv.error().message);
    }
    const auto end = Clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count() / static_cast<double>(iterations);
}

void print_row(const char* name, double time_ms) {
    std::cout << name << ' ' << std::fixed << std::setprecision(3) << time_ms << '\n';
}

}  // namespace

int main() {
    try {
        print_row("parse_implicit_mul_ms", benchmark_parse_implicit_multiplication());
        print_row("round_trip_ms", benchmark_round_trip());
        print_row("simplify_basic_ms", benchmark_symbolic_simplify());
        print_row("arena_alloc_ms", benchmark_arena_allocation());
        print_row("poly_expand_complex_ms", benchmark_polynomial_expansion_complex());
        print_row("poly_gcd_subresultant_ms", benchmark_polynomial_gcd_subresultant());
        print_row("matrix_inv_10x10_ms", benchmark_matrix_inversion_10x10());
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "benchmark_failure " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
