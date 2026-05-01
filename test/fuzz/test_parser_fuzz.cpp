#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "../helpers/property_test.hpp"

#include <gtest/gtest.h>

#include <string>

namespace cas {
namespace {

Result<ExprPtr> parse_input(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) {
        return Result<ExprPtr>(tokens.error());
    }

    Parser parser(tokens.value(), arena);
    return parser.parse();
}

std::string generate_parser_stress_input(test::DeterministicRng& rng, int depth) {
    if (depth <= 0) {
        return test::sample_exact_atom(rng);
    }

    switch (rng.next_int(0, 7)) {
        case 0:
            return test::generate_polynomial_expr(rng, depth - 1);
        case 1:
            return "sin(" + generate_parser_stress_input(rng, depth - 1) + ")";
        case 2:
            return "ln(" + generate_parser_stress_input(rng, depth - 1) + ")";
        case 3:
            return "(" + generate_parser_stress_input(rng, depth - 1) + ")/(" +
                generate_parser_stress_input(rng, depth - 1) + ")";
        case 4:
            return "(" + generate_parser_stress_input(rng, depth - 1) + ")^" +
                std::to_string(rng.next_int(0, 5));
        case 5:
            return "[[" + std::to_string(rng.next_int(0, 4)) + "," + std::to_string(rng.next_int(0, 4)) +
                "],[" + std::to_string(rng.next_int(0, 4)) + "," + std::to_string(rng.next_int(0, 4)) + "]]";
        case 6:
            return "RootOf(" + test::sample_symbol(rng) + "^2-" + std::to_string(rng.next_int(1, 5)) + ", " +
                test::sample_symbol(rng) + ")";
        case 7:
        default:
            return "diff(" + generate_parser_stress_input(rng, depth - 1) + ", " + test::sample_symbol(rng) + ")";
    }
}

TEST(ParserFuzzTest, DeterministicGeneratedInputsDoNotCrashLexerOrParser) {
    test::run_seeded_cases(0xCA5F00DULL, 256U, [](test::DeterministicRng& rng, std::size_t) {
        AstArena arena;
        const std::string input = generate_parser_stress_input(rng, 3);
        const auto parsed = parse_input(input, arena);

        if (parsed.is_ok()) {
            EXPECT_TRUE(static_cast<bool>(parsed.value())) << input;
        } else {
            EXPECT_EQ(parsed.error().kind, CASErrorKind::ParseError) << input;
        }
    });
}

}  // namespace
}  // namespace cas
