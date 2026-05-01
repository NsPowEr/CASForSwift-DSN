#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "../helpers/property_test.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace cas {
namespace {

Result<ExprPtr> parse_input(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) {
        return fail<ExprPtr>(tokens.error());
    }

    Parser parser(tokens.value(), arena);
    return parser.parse();
}

Result<ExprPtr> simplify_input(const std::string& input, AstArena& parse_arena, symbolic::CASContext& context) {
    auto parsed = parse_input(input, parse_arena);
    if (parsed.is_error()) {
        return fail<ExprPtr>(parsed.error());
    }
    return context.simplify(parsed.value());
}

bool simplification_is_idempotent(const std::string& input) {
    AstArena parse_arena;
    symbolic::CASContext context;
    auto first = simplify_input(input, parse_arena, context);
    if (first.is_error()) {
        return false;
    }

    auto second = context.simplify(first.value());
    return second.is_ok() && structural_equal(first.value(), second.value());
}

bool simplified_forms_are_structurally_equal(const std::string& lhs_text, const std::string& rhs_text) {
    AstArena lhs_arena;
    AstArena rhs_arena;
    symbolic::CASContext lhs_context;
    symbolic::CASContext rhs_context;

    auto lhs = simplify_input(lhs_text, lhs_arena, lhs_context);
    auto rhs = simplify_input(rhs_text, rhs_arena, rhs_context);
    return lhs.is_ok() && rhs.is_ok() && structural_equal(lhs.value(), rhs.value());
}

[[nodiscard]] std::string strip_redundant_outer_parens(std::string input) {
    while (input.size() >= 2U && input.front() == '(' && input.back() == ')') {
        int depth = 0;
        bool wraps_all = true;
        for (std::size_t index = 0; index < input.size(); ++index) {
            if (input[index] == '(') {
                ++depth;
            } else if (input[index] == ')') {
                --depth;
                if (depth == 0 && index + 1U < input.size()) {
                    wraps_all = false;
                    break;
                }
            }
            if (depth < 0) {
                wraps_all = false;
                break;
            }
        }
        if (!wraps_all) {
            break;
        }
        input = input.substr(1U, input.size() - 2U);
    }
    return input;
}

[[nodiscard]] std::vector<std::string> shrink_candidates(const std::string& input) {
    std::vector<std::string> candidates;
    candidates.reserve(32U);

    candidates.push_back(strip_redundant_outer_parens(input));
    candidates.push_back("0");
    candidates.push_back("1");
    candidates.push_back("x");

    std::vector<std::size_t> opens;
    for (std::size_t index = 0; index < input.size(); ++index) {
        if (input[index] == '(') {
            opens.push_back(index);
            continue;
        }
        if (input[index] != ')' || opens.empty()) {
            continue;
        }

        const std::size_t open = opens.back();
        opens.pop_back();
        const std::string inner = input.substr(open + 1U, index - open - 1U);
        if (!inner.empty()) {
            candidates.push_back(inner);
        }
        for (const std::string& replacement : {"0", "1", "x"}) {
            candidates.push_back(input.substr(0U, open) + replacement + input.substr(index + 1U));
        }
    }

    return candidates;
}

template <typename StillFails>
[[nodiscard]] std::string shrink_expression(std::string input, StillFails&& still_fails) {
    input = strip_redundant_outer_parens(std::move(input));
    bool changed = true;
    while (changed) {
        changed = false;
        for (std::string candidate : shrink_candidates(input)) {
            candidate = strip_redundant_outer_parens(std::move(candidate));
            if (candidate.empty() || candidate == input) {
                continue;
            }
            if (still_fails(candidate)) {
                input = std::move(candidate);
                changed = true;
                break;
            }
        }
    }
    return input;
}

template <typename StillFails>
[[nodiscard]] std::pair<std::string, std::string> shrink_expression_pair(
    std::string lhs,
    std::string rhs,
    StillFails&& still_fails) {
    bool changed = true;
    while (changed) {
        changed = false;

        const std::string shrunk_lhs = shrink_expression(lhs, [&](const std::string& candidate) {
            return still_fails(candidate, rhs);
        });
        if (shrunk_lhs != lhs) {
            lhs = shrunk_lhs;
            changed = true;
        }

        const std::string shrunk_rhs = shrink_expression(rhs, [&](const std::string& candidate) {
            return still_fails(lhs, candidate);
        });
        if (shrunk_rhs != rhs) {
            rhs = shrunk_rhs;
            changed = true;
        }
    }

    return {lhs, rhs};
}

void expect_structurally_same_after_simplify(
    const std::string& lhs_text,
    const std::string& rhs_text,
    std::string trace_label) {
    SCOPED_TRACE(std::move(trace_label));
    SCOPED_TRACE("lhs=" + lhs_text);
    SCOPED_TRACE("rhs=" + rhs_text);

    if (simplified_forms_are_structurally_equal(lhs_text, rhs_text)) {
        return;
    }

    auto minimal = shrink_expression_pair(lhs_text, rhs_text, [](const std::string& lhs, const std::string& rhs) {
        return !simplified_forms_are_structurally_equal(lhs, rhs);
    });
    ADD_FAILURE()
        << "minimal_lhs=" << minimal.first
        << "\nminimal_rhs=" << minimal.second;
}

void expect_simplification_idempotent(const std::string& input, std::string trace_label) {
    SCOPED_TRACE(std::move(trace_label));
    SCOPED_TRACE("expr=" + input);

    if (simplification_is_idempotent(input)) {
        return;
    }

    const std::string minimal = shrink_expression(input, [](const std::string& candidate) {
        return !simplification_is_idempotent(candidate);
    });
    ADD_FAILURE() << "minimal_expr=" << minimal;
}

std::string generated_symbolic_expr(test::DeterministicRng& rng, int depth) {
    if (depth <= 0) {
        return test::sample_exact_atom(rng);
    }

    switch (rng.next_int(0, 5)) {
        case 0:
            return test::generate_polynomial_expr(rng, depth - 1);
        case 1:
            return "(" + generated_symbolic_expr(rng, depth - 1) + ") + (" +
                generated_symbolic_expr(rng, depth - 1) + ")";
        case 2:
            return "(" + generated_symbolic_expr(rng, depth - 1) + ") * (" +
                generated_symbolic_expr(rng, depth - 1) + ")";
        case 3:
            return "(" + test::sample_symbol(rng) + " + " + test::sample_symbol(rng) + ")^" +
                std::to_string(rng.next_int(0, 3));
        case 4: {
            const std::string symbol = test::sample_symbol(rng);
            return "sin(" + symbol + ")^2 + cos(" + symbol + ")^2";
        }
        case 5:
        default:
            return "sqrt((" + test::sample_symbol(rng) + ")^2)";
    }
}

TEST(SymbolicFuzzTest, GeneratedPolynomialSimplificationIsIdempotent) {
    test::run_seeded_cases(0xF20051DULL, 384U, [](test::DeterministicRng& rng, std::size_t index) {
        const std::string expr = generated_symbolic_expr(rng, 3);
        expect_simplification_idempotent(expr, "seed=0xF20051D iteration=" + std::to_string(index));
    });
}

TEST(SymbolicFuzzTest, GeneratedAdditionAndMultiplicationCanonicalizeCommutativity) {
    test::run_seeded_cases(0xF2C0A11U, 256U, [](test::DeterministicRng& rng, std::size_t index) {
        const std::string lhs = generated_symbolic_expr(rng, 2);
        const std::string rhs = generated_symbolic_expr(rng, 2);
        const std::string trace = "seed=0xF2C0A11 iteration=" + std::to_string(index);

        expect_structurally_same_after_simplify(
            "(" + lhs + ") + (" + rhs + ")",
            "(" + rhs + ") + (" + lhs + ")",
            trace + " additive");
        expect_structurally_same_after_simplify(
            "(" + lhs + ") * (" + rhs + ")",
            "(" + rhs + ") * (" + lhs + ")",
            trace + " multiplicative");
    });
}

TEST(SymbolicFuzzTest, GeneratedDistributiveFormsCanonicalize) {
    test::run_seeded_cases(0xF2D157ULL, 64U, [](test::DeterministicRng& rng, std::size_t index) {
        const std::string a = test::sample_exact_atom(rng);
        const std::string b = test::sample_exact_atom(rng);
        const std::string c = test::sample_exact_atom(rng);

        expect_structurally_same_after_simplify(
            "(" + a + ") * ((" + b + ") + (" + c + "))",
            "((" + a + ") * (" + b + ")) + ((" + a + ") * (" + c + "))",
            "seed=0xF2D157 iteration=" + std::to_string(index));
    });
}

}  // namespace
}  // namespace cas
